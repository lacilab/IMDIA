"""Conservative fixed-protonation preparation; not constant-pH MD.
"""
import csv
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

CANONICAL = {
    'ASH': 'ASP', 'GLH': 'GLU', 'HID': 'HIS', 'HIE': 'HIS',
    'HIP': 'HIS', 'LYN': 'LYS', 'CYM': 'CYS', 'CYX': 'CYS',
}
STATES = {
    'ASP': ('ASH', 'ASP'), 'GLU': ('GLH', 'GLU'),
    'LYS': ('LYS', 'LYN'), 'CYS': ('CYS', 'CYM'),
    'HIS': ('HIP', None), 'TYR': ('TYR', None),
    'ARG': ('ARG', None),
}


def normalize_config(config):
    if config is None:
        return {'method': 'none'}
    if not isinstance(config, dict):
        raise ValueError('protonation must be a mapping')
    method = config.get('method', 'none')
    if method not in ('none', 'propka'):
        raise ValueError('protonation.method must be none or propka')
    if method == 'none':
        return {'method': 'none'}
    ph = float(config.get('ph', 7.4))
    window = float(config.get('review_pka_window', 1.0))
    if not math.isfinite(ph) or not 0 <= ph <= 14:
        raise ValueError('pH must be finite and between 0 and 14')
    if not math.isfinite(window) or window < 0:
        raise ValueError('review_pka_window must be finite and nonnegative')
    overrides = config.get('overrides') or []
    if not isinstance(overrides, list):
        raise ValueError('overrides must be a list')
    return dict(method=method, ph=ph, review_pka_window=window,
                overrides=overrides)


def guard_preparation_cache(datadir, input_pdb, config, dna, phos, ffs):
    """Refuse stale/legacy checkpoints; use a fresh run directory to rebuild."""
    config = normalize_config(config)
    datadir, input_pdb = Path(datadir), Path(input_pdb)
    manifest = datadir / 'protonation_preparation.json'
    payload = dict(input_sha256=hashlib.sha256(input_pdb.read_bytes()).hexdigest(),
                   protonation=config, dna=dna, phos=phos, ffs=ffs,
                   implementation=1)
    # Existing legacy workflows may continue unchanged when feature is disabled.
    if config['method'] == 'none' and not manifest.exists():
        return config
    if manifest.exists():
        if json.loads(manifest.read_text()) != payload:
            raise RuntimeError('Preparation settings/input changed. Use a fresh '
                               'run directory containing the input PDB; old '
                               'AMBER/GROMACS checkpoints must not be reused.')
    else:
        existing = [datadir / 'amber', datadir / 'gromacs_common']
        if any(p.exists() and any(p.iterdir()) for p in existing):
            raise RuntimeError('Existing preparation has no pH manifest. Use a '
                               'fresh run directory for pH-aware preparation.')
        manifest.write_text(json.dumps(payload, indent=2) + '\n')
    return config


def residues(lines):
    result = {}
    for line in lines:
        if line.startswith('MODEL'):
            raise ValueError('Select a single PDB model before preparation')
        if line.startswith(('ATOM  ', 'HETATM')):
            key = (line[21].strip(), int(line[22:26]), line[26].strip())
            name = line[17:20].strip()
            if key in result and result[key] != name:
                raise ValueError(f'Ambiguous residue identifier: {key}')
            result[key] = name
    return result


def parse_pkas(text):
    """Parse PROPKA summary; predictor uses unique numeric residue identifiers."""
    if 'SUMMARY OF THIS PREDICTION' not in text:
        raise RuntimeError('PROPKA output has no prediction summary')
    summary = text.split('SUMMARY OF THIS PREDICTION', 1)[1]
    result = {}
    for line in summary.splitlines():
        parts = line.split()
        if len(parts) < 4 or parts[0] not in STATES:
            continue
        try:
            index = int(parts[1])
            # Protein summary: NAME NUMBER [CHAIN] pKa model-pKa.
            pka = float(parts[-2].rstrip('*'))
        except ValueError:
            continue
        if index in result:
            raise RuntimeError(f'Duplicate PROPKA prediction: {index}')
        result[index] = pka
    return result


def assign_protonation(pdbfile, original_pdb, config, logger):
    config = normalize_config(config)
    if config['method'] == 'none':
        return str(pdbfile)
    path = Path(pdbfile).resolve()
    lines = path.read_text().splitlines(keepends=True)
    prepared = residues(lines)
    original = residues(Path(original_pdb).read_text().splitlines(keepends=True))
    # pdb4amber can renumber residues. Map by residue order, with a strict
    # identity check; refuse removals/reordering instead of guessing.
    if len(prepared) != len(original):
        raise ValueError('pdb4amber changed residue count; explicit mapping needed')
    mapping = {}
    for (old_key, old_name), (key, name) in zip(original.items(), prepared.items()):
        if CANONICAL.get(old_name, old_name) != CANONICAL.get(name, name):
            raise ValueError('pdb4amber changed residue order/identity; mapping needed')
        mapping[old_key] = key
    overrides = {}
    for item in config['overrides']:
        old_key = (str(item.get('chain', '')).strip(), int(item['residue']),
                   str(item.get('insertion_code', '')).strip())
        if old_key not in mapping:
            raise ValueError(f'Override residue absent from original PDB: {old_key}')
        key = mapping[old_key]
        if key in overrides:
            raise ValueError(f'Duplicate override: {old_key}')
        base = CANONICAL.get(prepared[key], prepared[key])
        state = str(item['state']).upper()
        allowed = {x for x in STATES.get(base, ()) if x}
        if base == 'HIS':
            allowed.update(('HID', 'HIE'))
        if state not in allowed:
            raise ValueError(f'Unsupported override {state} for {old_key} ({base})')
        if prepared[key] == 'CYX':
            raise ValueError('Do not titrate disulfide CYX residues')
        overrides[key] = state
    # Renumber ONLY the temporary predictor input to eliminate insertion-code
    # collisions. Keep every component and TER record in the assembled complex.
    indices = {key: i for i, key in enumerate(prepared, 1)}
    if len(indices) > 9999:
        raise ValueError('Predictor PDB exceeds 9999 residues; mmCIF support needed')
    predictor = []
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            key = (line[21].strip(), int(line[22:26]), line[26].strip())
            name = CANONICAL.get(prepared[key], prepared[key])
            line = line[:17] + f'{name:>3}' + line[20:22] + f'{indices[key]:4d} ' + line[27:]
        # Original serial/residue references do not apply to predictor input.
        if not line.startswith(('SSBOND', 'CONECT', 'ANISOU', 'LINK')):
            predictor.append(line)
    with tempfile.TemporaryDirectory(prefix='propka_', dir=path.parent) as tmp:
        tmp = Path(tmp)
        (tmp / 'predictor.pdb').write_text(''.join(predictor))
        result = subprocess.run([sys.executable, '-m', 'propka', 'predictor.pdb',
                                 '--display-coupled-residues'], cwd=tmp,
                                text=True, capture_output=True)
        path.with_suffix('.propka.log').write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError('PROPKA failed; inspect the .propka.log file')
        pka_files = sorted(tmp.glob('*.pka'))

        if len(pka_files) != 1:
            generated = ', '.join(p.name for p in tmp.iterdir())
            raise RuntimeError(
                f'Expected one PROPKA .pka file, found {len(pka_files)}. '
                f'Files generated: {generated}. '
                'Inspect the .propka.log file.'
            )

        pka_path = pka_files[0]
        logger.debug(f'Reading PROPKA output: {pka_path.name}')
        text = pka_path.read_text()
        path.with_suffix('.pka').write_text(text)
        pkas = parse_pkas(text)
    decisions, rows, errors = {}, [], []
    old_keys = {v: k for k, v in mapping.items()}
    for key, name in prepared.items():
        base = CANONICAL.get(name, name)
        if base not in STATES or name == 'CYX':
            continue
        pka = pkas.get(indices[key])
        state, reason = name, 'PROPKA'
        if key in overrides:
            state, reason = overrides[key], 'manual override'
        elif pka is None:
            reason = 'missing prediction'
            errors.append(f'{old_keys[key]}: no pKa prediction; supply override')
        else:
            protonated = config['ph'] < pka
            state = STATES[base][0 if protonated else 1]
            if base == 'HIS' and not protonated:
                if name in ('HID', 'HIE'):
                    state, reason = name, 'preserved input neutral tautomer'
                else:
                    errors.append(f'{old_keys[key]}: choose HID or HIE in overrides')
            elif state is None:
                errors.append(f'{old_keys[key]}: unsupported deprotonated {base}')
        review = pka is not None and abs(config['ph'] - pka) <= config['review_pka_window']
        if review:
            logger.warning(f'Review {old_keys[key]}: pKa={pka:.2f}, pH={config["ph"]}')
        decisions[key] = state or name
        rows.append([*old_keys[key], *key, name, pka, state, reason, review])
    report = path.with_suffix('.protonation.csv')
    with report.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['original_chain', 'original_residue', 'original_icode',
                         'prepared_chain', 'prepared_residue', 'prepared_icode',
                         'input_state', 'pka', 'assigned_state', 'reason', 'near_pka'])
        writer.writerows(rows)
    logger.warning('Review PROPKA coupled-group output; assignments assume independent sites. '
                   'Termini, DNA, phosphoresidues, ligands and metals retain input/default chemistry.')
    if errors:
        raise ValueError('Protonation needs review:\n' + '\n'.join(errors) + f'\nReport: {report}')
    # In-place naming preserves existing downstream file naming. Backup first.
    path.with_suffix('.unassigned.pdb').write_text(''.join(lines))
    output = []
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            key = (line[21].strip(), int(line[22:26]), line[26].strip())
            if key in decisions:
                line = line[:17] + f'{decisions[key]:>3}' + line[20:]
        output.append(line)
    path.write_text(''.join(output))
    logger.info(f'Assigned protein side-chain states at pH {config["ph"]}; report: {report}')
    return str(path)
