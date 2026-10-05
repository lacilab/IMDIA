"""Assign fixed protonation states before LEaP adds hydrogens
"""
import csv
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

ATOM_RECORDS = ('ATOM  ', 'HETATM')
STALE_REFERENCES = ('SSBOND', 'CONECT', 'ANISOU', 'LINK')
HISTIDINE_HYDROGENS = {'HID': {'HD1'}, 'HIE': {'HE2'}}

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


def residue_key(line):
    return (line[21].strip(), int(line[22:26]), line[26].strip())


def is_hydrogen(line):
    name = line[12:16].strip()
    element = line[76:78].strip() or name.lstrip('0123456789')[:1]
    return element in ('H', 'D')


def predictor_input(lines, indices, states):
    """Renumber only predictor input; retain components and TER boundaries."""
    output = []
    for line in lines:
        if line.startswith(STALE_REFERENCES):
            continue
        if line.startswith(ATOM_RECORDS):
            key = residue_key(line)
            line = (line[:17] + f'{states[key]:>3}' + line[20:22]
                    + f'{indices[key]:4d} ' + line[27:])
        output.append(line)
    return ''.join(output)


def run_predictor(command, work, logpath, failure_message):
    result = subprocess.run(command, cwd=work, capture_output=True, text=True)
    logpath.write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f'{failure_message}; inspect {logpath}')


def residues(lines):
    result = {}
    for line in lines:
        if line.startswith('MODEL'):
            raise ValueError('Select a single PDB model before preparation')
        if line.startswith(ATOM_RECORDS):
            key = residue_key(line)
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


PROTEIN_NAMES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
    'TYR', 'VAL',
}


def pdb2pqr_histidines(lines, prepared, indices, decisions, targets, path, config, logger):
    """Choose HID/HIE and validate the optimized protein coordinates."""
    unsupported = {name for name in prepared.values()
                   if CANONICAL.get(name, name) not in PROTEIN_NAMES
                   and name not in ('HOH', 'WAT')}
    if unsupported:
        raise ValueError('Automatic histidine assignment currently supports only '
                         'standard proteins and water. Unsupported components: '
                         + ', '.join(sorted(unsupported))
                         + '. Supply explicit histidine overrides or use a '
                         'validated preparation adapter retaining their context.')
    original_atoms = {}
    for line in lines:
        if line.startswith(ATOM_RECORDS):
            key = residue_key(line)
            if line[16].strip():
                raise ValueError('Resolve alternate conformations before automatic assignment')
            name = line[12:16].strip()
            if is_hydrogen(line):
                raise ValueError('Automatic stage expects a hydrogen-free input PDB')
            if CANONICAL.get(prepared[key], prepared[key]) in PROTEIN_NAMES:
                atom_key = (key, name)
                if atom_key in original_atoms:
                    raise ValueError(f'Duplicate input atom: {atom_key}')
                original_atoms[atom_key] = line
    # Only histidines without overrides are free to change tautomer.
    prediction_states = {}
    for key, name in prepared.items():
        state = 'HIS' if key in targets else decisions.get(key, name)
        prediction_states[key] = 'HOH' if state == 'WAT' else state
    # Run PDB2PQR in the pipeline's Python environment.
    with tempfile.TemporaryDirectory(prefix='pdb2pqr_', dir=path.parent) as work:
        work = Path(work)
        (work / 'input.pdb').write_text(predictor_input(lines, indices, prediction_states))
        command = [sys.executable, '-c', 'from pdb2pqr.main import main; main()',
                   '--ff=AMBER', '--ffout=AMBER', '--keep-chain', '--nodebump',
                   '--titration-state-method=propka', '--with-ph', str(config['ph']),
                   '--pdb-output=optimized.pdb', 'input.pdb', 'optimized.pqr']
        run_predictor(command, work, path.with_suffix('.pdb2pqr.log'),
                      'Automatic histidine assignment failed')
        output = work / 'optimized.pdb'
        if not output.is_file():
            raise RuntimeError('PDB2PQR did not write optimized.pdb')
        output_text = output.read_text()
        path.with_suffix('.pdb2pqr.pdb').write_text(output_text)
    by_index = {number: key for key, number in indices.items()}
    output_states, output_hydrogens, coords = {}, {}, {}
    for line in output_text.splitlines():
        if not line.startswith(ATOM_RECORDS):
            continue
        number = int(line[22:26])
        if number not in by_index:
            raise ValueError(f'Unexpected PDB2PQR residue number: {number}')
        key = by_index[number]
        if line[21].strip() != key[0]:
            raise ValueError(f'PDB2PQR changed chain identity at {key}')
        if CANONICAL.get(prepared[key], prepared[key]) not in PROTEIN_NAMES:
            continue  # Original waters remain unchanged in the final PDB.
        state = line[17:20].strip()
        if CANONICAL.get(state, state) != CANONICAL.get(prepared[key], prepared[key]):
            raise ValueError(f'PDB2PQR changed residue identity at {key}')
        if key in output_states and output_states[key] != state:
            raise ValueError(f'Ambiguous PDB2PQR state at {key}')
        output_states[key] = state
        name = line[12:16].strip()
        if is_hydrogen(line):
            output_hydrogens.setdefault(key, set()).add(name)
            continue
        atom_key = (key, name)
        if atom_key in coords:
            raise ValueError(f'Duplicate PDB2PQR atom: {atom_key}')
        xyz = [float(line[i:i+8]) for i in (30, 38, 46)]
        if not all(math.isfinite(x) for x in xyz):
            raise ValueError(f'Invalid PDB2PQR coordinates: {atom_key}')
        coords[atom_key] = line[30:54]
    if coords.keys() != original_atoms.keys():
        raise ValueError('PDB2PQR added/removed or renamed protein heavy atoms. '
                         'Automatic coordinate transfer requires explicit atom mapping.')
    for key, state in decisions.items():
        if key not in targets and key in output_states and output_states[key] != state:
            raise ValueError(f'PDB2PQR changed requested state at {key}: '
                             f'{state} -> {output_states[key]}; assignment needs review')
    assignments = {}
    for key in targets:
        state = output_states.get(key)
        if state not in HISTIDINE_HYDROGENS:
            raise ValueError(f'PDB2PQR did not return neutral HID/HIE for {key}: {state}')
        ring_h = output_hydrogens.get(key, set()) & {'HD1', 'HE2'}
        if ring_h != HISTIDINE_HYDROGENS[state]:
            raise ValueError(f'PDB2PQR state/hydrogens disagree at {key}')
        assignments[key] = state
        logger.info(f'Automatic histidine {key}: {state} (PDB2PQR optimization)')
    changed = sum(original_atoms[k][30:54] != value for k, value in coords.items())
    logger.info(f'Transferring optimized protein heavy-atom coordinates; '
                f'{changed} atom records differ. Original waters/components are retained.')
    return assignments, coords


def run_propka(lines, prepared, indices, path, logger):
    """Run PROPKA and save its log and predictions beside the input PDB."""
    predictor = predictor_input(lines, indices,
                                {key: CANONICAL.get(name, name) for key, name in prepared.items()})
    with tempfile.TemporaryDirectory(prefix='propka_', dir=path.parent) as tmp:
        tmp = Path(tmp)
        (tmp / 'predictor.pdb').write_text(predictor)
        run_predictor([sys.executable, '-m', 'propka', 'predictor.pdb',
                       '--display-coupled-residues'], tmp,
                      path.with_suffix('.propka.log'), 'PROPKA failed')
        pka_files = sorted(tmp.glob('*.pka'))
        if len(pka_files) != 1:
            raise RuntimeError(f'Expected one PROPKA .pka file, found {len(pka_files)}')
        pka_path = pka_files[0]
        logger.debug(f'Reading PROPKA output: {pka_path.name}')
        text = pka_path.read_text()
        path.with_suffix('.pka').write_text(text)
        return parse_pkas(text)


def prepare_overrides(original, prepared, override_items):
    """Map original residue IDs and validate the requested AMBER states."""
    # pdb4amber may renumber residues; overrides use the original IDs.
    if len(prepared) != len(original):
        raise ValueError('pdb4amber changed residue count; explicit mapping needed')
    mapping = {}
    for (old_key, old_name), (key, name) in zip(original.items(), prepared.items()):
        if CANONICAL.get(old_name, old_name) != CANONICAL.get(name, name):
            raise ValueError('pdb4amber changed residue order/identity; mapping needed')
        mapping[old_key] = key
    overrides = {}
    for item in override_items:
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
    return mapping, overrides

def terminal_residues(lines):
    """Find protein endpoints using chain IDs; a terminal GLU would be GLH already and further protonation would make it -2 OXT, not supported yet"""
    terminals = set()
    segment = []

    for line in lines:
        if line.startswith('TER'):
            if segment:
                terminals.update((segment[0], segment[-1]))
                segment = []
            continue

        if not line.startswith(('ATOM  ', 'HETATM')):
            continue

        name = line[17:20].strip()
        if CANONICAL.get(name, name) not in PROTEIN_NAMES:
            continue

        key = residue_key(line)

        if segment and key[0] != segment[-1][0]:
            terminals.update((segment[0], segment[-1]))
            segment = []

        if not segment or key != segment[-1]:
            segment.append(key)

    if segment:
        terminals.update((segment[0], segment[-1]))

    return terminals

def assign_protonation(pdbfile, original_pdb, config, logger):
    """Apply overrides, predict remaining states, and write the PDB and report."""
    config = normalize_config(config)
    if config['method'] == 'none':
        return str(pdbfile)
    path = Path(pdbfile).resolve()
    lines = path.read_text().splitlines(keepends=True)
    prepared = residues(lines)
    original = residues(Path(original_pdb).read_text().splitlines(keepends=True))
    mapping, overrides = prepare_overrides(original, prepared, config['overrides'])

    # Unique predictor IDs avoid chain and insertion-code collisions.
    indices = {key: i for i, key in enumerate(prepared, 1)}
    if len(indices) > 9999:
        raise ValueError('Predictor PDB exceeds 9999 residues; mmCIF support needed')
    pkas = run_propka(lines, prepared, indices, path, logger)

    decisions, rows, errors = {}, [], []
    targets, coords = set(), {}
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
                # Reassess neutral HIS unless YAML explicitly fixes its state.
                targets.add(key)
                reason = 'pending PDB2PQR optimization'
            elif state is None:
                errors.append(f'{old_keys[key]}: unsupported deprotonated {base}')
        review = pka is not None and abs(config['ph'] - pka) <= config['review_pka_window']
        if review:
            logger.warning(f'Review {old_keys[key]}: pKa={pka:.2f}, pH={config["ph"]}')
        decisions[key] = state or name
        rows.append([*old_keys[key], *key, name, pka, state, reason, review, config['ph']])
    for key in sorted(terminal_residues(lines)):
        state = decisions.get(key, prepared[key])

        if state in {'ASH', 'GLH', 'LYN', 'CYM'}:
            errors.append(
                f'{old_keys[key]}: terminal {state} is not validated '
                'in this pipeline. Supply an explicit supported-state '
                'override after reviewing terminal chemistry.'
            )
    if targets and not errors:
        try:
            automatic, coords = pdb2pqr_histidines(
                lines, prepared, indices, decisions, targets, path, config, logger)
            decisions.update(automatic)
            for row in rows:
                key = tuple(row[3:6])
                if key in automatic:
                    row[8] = automatic[key]
                    row[9] = 'PDB2PQR hydrogen-bond optimization'
        except (ValueError, RuntimeError, OSError) as exc:
            errors.append(str(exc))
    report = path.with_suffix('.protonation.csv')
    with report.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['original_chain', 'original_residue', 'original_icode',
                         'prepared_chain', 'prepared_residue', 'prepared_icode',
                         'input_state', 'pka', 'assigned_state', 'reason', 'near_pka', 'target_ph'])
        writer.writerows(rows)
    logger.warning('Review PROPKA coupled-group output; assignments assume independent sites. '
                   'Termini, DNA, phosphoresidues, ligands and metals retain input/default chemistry.')
    if errors:
        raise ValueError('Protonation needs review:\n' + '\n'.join(errors) + f'\nReport: {report}')
    # In-place naming preserves existing downstream file naming. Backup first.
    path.with_suffix('.unassigned.pdb').write_text(''.join(lines))
    output = []
    for line in lines:
        if line.startswith(ATOM_RECORDS):
            key = residue_key(line)
            if key in decisions:
                line = line[:17] + f'{decisions[key]:>3}' + line[20:]
            atom_key = (key, line[12:16].strip())
            if atom_key in coords:
                line = line[:30] + coords[atom_key] + line[54:]
        output.append(line)
    path.write_text(''.join(output))
    logger.info(f'Assigned protein side-chain states at pH {config["ph"]}; report: {report}')
    return str(path)
