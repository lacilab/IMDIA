#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Assign supported protein side-chain states, including residues at either terminus.

PROPKA terminal-group predictions are reported separately; LEaP controls terminal
charges and templates. Modified-residue states retain the existing input policy.
"""
import csv
import contextlib
import io
import logging
import json
import math
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

ATOM_RECORDS = ('ATOM  ', 'HETATM')
STALE_REFERENCES = ('SSBOND', 'CONECT', 'ANISOU', 'LINK')
CANONICAL = {
    'ASH': 'ASP', 'GLH': 'GLU', 'HID': 'HIS', 'HIE': 'HIS',
    'HIP': 'HIS', 'LYN': 'LYS', 'CYM': 'CYS', 'CYX': 'CYS',
}
# Each pair is (protonated, deprotonated). None requires special handling.
STATES = {
    'ASP': ('ASH', 'ASP'), 'GLU': ('GLH', 'GLU'),
    'LYS': ('LYS', 'LYN'), 'CYS': ('CYS', 'CYM'),
    'HIS': ('HIP', None), 'TYR': ('TYR', None), 'ARG': ('ARG', None),
}
PROTEIN_NAMES = set(
    'ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split()
)
# Admission list only; the output reader also checks that these atoms survive.
NUCLEIC_NAMES = {
    base + suffix for base in ('DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'U')
    for suffix in ('', '3', '5', 'N')
}
WATER_NAMES = {'HOH', 'WAT', 'TP3', 'SOL'}
PHOSPHO_RESIDUES = {
    'S1P': {'parent': 'SER', 'phosphate_charge': -1, 'residue_charge': -1, 'attachment': 'OG'},
    'SEP': {'parent': 'SER', 'phosphate_charge': -2, 'residue_charge': -2, 'attachment': 'OG'},
    'T1P': {'parent': 'THR', 'phosphate_charge': -1, 'residue_charge': -1, 'attachment': 'OG'},
    'TPO': {'parent': 'THR', 'phosphate_charge': -2, 'residue_charge': -2, 'attachment': 'OG1'},
    'Y1P': {'parent': 'TYR', 'phosphate_charge': -1, 'residue_charge': -1, 'attachment': 'OH'},
    'PTR': {'parent': 'TYR', 'phosphate_charge': -2, 'residue_charge': -2, 'attachment': 'OH'},
    'H1D': {'parent': None, 'phosphate_charge': -1, 'residue_charge': 0, 'attachment': 'ND1'},
    'H2D': {'parent': None, 'phosphate_charge': -2, 'residue_charge': -1, 'attachment': 'ND1'},
    'H1E': {'parent': None, 'phosphate_charge': -1, 'residue_charge': 0, 'attachment': 'NE2'},
    'H2E': {'parent': None, 'phosphate_charge': -2, 'residue_charge': -1, 'attachment': 'NE2'},
}
PHOSPHO_HISTIDINES = {'H1D', 'H2D', 'H1E', 'H2E'}
SURROGATES = {name: definition['parent'] for name, definition in PHOSPHO_RESIDUES.items()
              if definition['parent'] is not None}
PARENT_ATOMS = {
    'SER': set('N CA C O CB OG'.split()),
    'THR': set('N CA C O CB OG1 CG2'.split()),
    'TYR': set('N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH'.split()),
}
# T1P's library calls its side-chain oxygen OG, whereas THR uses OG1.
# Rename it only in the temporary input; final atoms keep their original names.
SURROGATE_ATOM_ALIASES = {'T1P': {'OG': 'OG1'}, 'TPO': {'OG': 'OG1'}}
# No HIP replacement: an N-phosphate and an N-proton are not equivalent.

PHOSPHATE_ATOMS = {'P', 'O1P', 'O2P', 'O3P', 'OP1', 'OP2', 'OP3'}
PHOSPHATE_ALIASES = {'OP1': 'O1P', 'OP2': 'O2P', 'OP3': 'O3P'}
UNUSUAL_STATES = {'LYN', 'ASH', 'GLH', 'CYM'}
# Expected movements for this PDB2PQR setup. Other transformations require review.
FLIPPABLE_ATOMS = {
    'HIS': {'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
    'ASN': {'CG', 'OD1', 'ND2'}, 'GLN': {'CD', 'OE1', 'NE2'},
}
# Side-chain charges only; do not use their sum as the total system charge.
STATE_CHARGES = {
    'ASH': 0, 'ASP': -1, 'GLH': 0, 'GLU': -1, 'LYS': 1, 'LYN': 0,
    'CYS': 0, 'CYM': -1, 'HIP': 1, 'HID': 0, 'HIE': 0, 'TYR': 0, 'ARG': 1,
}
RESIDUE_NAMES = {
    'ASP': 'Aspartate', 'GLU': 'Glutamate', 'HIS': 'Histidine',
    'LYS': 'Lysine', 'CYS': 'Cysteine', 'TYR': 'Tyrosine', 'ARG': 'Arginine',
}
HISTIDINE_HYDROGENS = {'HID': {'HD1'}, 'HIE': {'HE2'}}
# Keep the original 12 columns in place for existing consumers; append diagnostics.
REPORT_FIELDS = [
    'original_chain', 'original_residue', 'original_icode',
    'prepared_chain', 'prepared_residue', 'prepared_icode',
    'input_state', 'pka', 'assigned_state', 'reason', 'near_pka', 'target_ph',
    'coupled', 'delta_pka', 'flags', 'propka_index', 'pka_source',
    'alternative_pka', 'assigned_charge', 'chain_mismatch', 'tautomer_mismatch',
    'decision_required', 'decision_status', 'decision_reasons',
    'pka_status', 'suppressed_by', 'covalent_partners',
]


def canonical(name):
    return CANONICAL.get(name, name)


def normalize_config(config):
    """Validate protonation options and fill in defaults."""
    if config is None:
        return {'method': 'none'}
    if not isinstance(config, dict):
        raise ValueError('protonation must be a mapping')
    unknown = set(config) - {'method', 'ph', 'review_pka_window', 'overrides',
                             'decision_mode',
                             'acknowledge_phospho_approximation'}
    if unknown:
        raise ValueError('Unknown protonation key(s): ' + ', '.join(sorted(map(str, unknown))))
    method = config.get('method', 'none')
    if method not in ('none', 'propka'):
        raise ValueError('protonation.method must be none or propka')
    decision_mode = config.get('decision_mode', 'automatic')
    if decision_mode not in ('automatic', 'manual'):
        raise ValueError('protonation.decision_mode must be automatic or manual')
    acknowledge = config.get('acknowledge_phospho_approximation', False)
    if not isinstance(acknowledge, bool):
        raise ValueError(
            'protonation.acknowledge_phospho_approximation '
            'must be a YAML boolean'
        )
    if method == 'none':
        return {'method': 'none'}

    ph = float(config.get('ph', 7.4))
    window = float(config.get('review_pka_window', 1.0))
    if not math.isfinite(ph) or not 0 <= ph <= 14:
        raise ValueError('pH must be finite and between 0 and 14')
    if not math.isfinite(window) or window < 0:
        raise ValueError('review_pka_window must be finite and nonnegative')
    overrides = config.get('overrides', [])
    if overrides is None:
        overrides = []
    if not isinstance(overrides, list):
        raise ValueError('overrides must be a list')
    return dict(method=method, ph=ph, review_pka_window=window, overrides=overrides,
                decision_mode=decision_mode,
                acknowledge_phospho_approximation=acknowledge)


def residue_key(line):
    return line[21].strip(), int(line[22:26]), line[26].strip()


def is_hydrogen(line):
    name = line[12:16].strip()
    element = line[76:78].strip() or name.lstrip('0123456789')[:1]
    return element in ('H', 'D')


def residues(lines):
    found = {}
    for line in lines:
        if line.startswith('MODEL'):
            raise ValueError('Select one PDB model and remove its MODEL/ENDMDL records')
        if not line.startswith(ATOM_RECORDS):
            continue
        key, name = residue_key(line), line[17:20].strip()
        if key in found and found[key] != name:
            raise ValueError(f'Ambiguous residue identifier: {key}')
        found[key] = name
    return found


def predictor_input(lines, indices, states, drop=None, atom_aliases=None):
    """Renumber only the temporary predictor input.

    drop lists atoms omitted from the PDB2PQR surrogate. PROPKA receives the
    complete prepared structure, including all original phosphate atoms.
    """
    output = []
    previous_chain = None
    for line in lines:
        if line.startswith(STALE_REFERENCES):
            continue
        if line.startswith('TER'):
            previous_chain = None
        if line.startswith(ATOM_RECORDS):
            key = residue_key(line)
            # PROPKA uses TER, rather than chain-ID changes, to start a new N terminus.
            if previous_chain is not None and key[0] != previous_chain:
                output.append('TER   \n')
            previous_chain = key[0]
            if line[12:16].strip() in (drop or {}).get(key, set()):
                continue
            atom_name = line[12:16].strip()
            alias = (atom_aliases or {}).get(key, {}).get(atom_name)
            if alias is not None:
                line = line[:12] + f'{alias:>4}' + line[16:]
            line = (line[:17] + f'{states[key]:>3}' + line[20:22]
                    + f'{indices[key]:4d} ' + line[27:])
        output.append(line)
    return ''.join(output)


def run_predictor(command, work, logpath, failure_message):
    result = subprocess.run(command, cwd=work, capture_output=True, text=True)
    logpath.write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f'{failure_message}; inspect {logpath}')


def run_propka(lines, prepared, indices, path, logger):
    """Read calculated groups directly, including side chains at protein termini.

    N+/C- groups are recorded separately and never replace side-chain pKas.
    Hidden calculated values are flagged candidates, not independent titrations.
    This stage assigns side chains; LEaP supplies the terminal-group templates.
    """
    import propka.run

    states = {key: canonical(name) for key, name in prepared.items()}
    input_text = predictor_input(lines, indices, states)
    path.with_suffix('.propka.input.pdb').write_text(input_text)
    # Capture PROPKA messages without changing the pipeline logger or cwd.
    propka_logger = logging.getLogger('propka')
    previous = propka_logger.level, propka_logger.propagate
    with path.with_suffix('.propka.log').open('w') as log:
        handler = logging.StreamHandler(log)
        propka_logger.addHandler(handler)
        propka_logger.setLevel(logging.INFO)
        propka_logger.propagate = False
        try:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                molecule = propka.run.single(
                    'predictor.pdb', stream=io.StringIO(input_text), write_pka=False)
                molecule.write_pka(filename=str(path.with_suffix('.pka')))
        finally:
            propka_logger.removeHandler(handler)
            propka_logger.setLevel(previous[0])
            propka_logger.propagate = previous[1]

    by_index = {index: key for key, index in indices.items()}
    if len(by_index) != len(indices):
        raise ValueError('Predictor indices must be globally unique')
    pkas, coupled, details, records = {}, set(), {}, []
    for group in molecule.conformations['AVR'].groups:
        name, index = group.residue_type.strip(), group.atom.res_num
        # Disulfide CYS is non-titratable but carries PROPKA's 99.99 sentinel.
        if not group.titratable and name not in STATES:
            continue
        chain, pka = group.atom.chain_id.strip(), group.pka_value
        key = by_index.get(index)
        if key is None or key[0] != chain or not math.isfinite(pka):
            raise ValueError(f'Invalid PROPKA group: {name} {index} {chain}')
        partner = group.coupled_titrating_group
        record = dict(name=name, index=index, chain=chain, pka=pka,
                      prepared_residue=key[1], prepared_icode=key[2],
                      covalent_partners=[g.label for g in group.covalently_coupled_groups],
                      noncovalent_partners=[g.label for g in group.non_covalently_coupled_groups],
                      suppressed_by=partner.label if partner is not None else '')
        records.append(record)
        if record['covalent_partners'] or record['noncovalent_partners'] or partner is not None:
            coupled.add((name, index, chain))
        if name not in STATES:
            continue
        if canonical(prepared[key]) != name or index in pkas:
            raise ValueError(f'Unexpected or duplicate PROPKA side chain: {name} {index} {chain}')
        pkas[index], details[index] = pka, record
    path.with_suffix('.propka.groups.json').write_text(json.dumps(records, indent=2, allow_nan=False))
    # These files belong to the old, second calculation, which is no longer run.
    for suffix in ('.coupled.pka', '.coupled.propka.log'):
        path.with_suffix(suffix).unlink(missing_ok=True)
    console_message(logger, 'info', 'Side-chain assignments use calculated PROPKA groups; '
                    'coupled sites remain flagged for review. Terminal N+/C- predictions '
                    'are recorded in .propka.groups.json; terminal charges remain with LEaP.')
    return pkas, coupled, details


def prepare_overrides(original, prepared, items):
    """Match original IDs to prepared residues in their existing order.

    Counts and canonical identities must agree. This is not a sequence alignment.
    Disulfide CYX residues cannot be titrated through an override.
    """
    if len(original) != len(prepared):
        raise ValueError('pdb4amber changed residue count; explicit mapping needed')
    mapping = {}
    for (old_key, old_name), (key, name) in zip(original.items(), prepared.items()):
        if canonical(old_name) != canonical(name):
            raise ValueError('pdb4amber changed residue order/identity; mapping needed')
        mapping[old_key] = key

    overrides = {}
    for item in items:
        if not isinstance(item, dict):
            raise ValueError('Each override must be a mapping')
        old_key = (str(item.get('chain', '')).strip(), int(item['residue']),
                   str(item.get('insertion_code', '')).strip())
        if old_key not in mapping:
            raise ValueError(f'Override residue absent from original PDB: {old_key}')
        key = mapping[old_key]
        if key in overrides:
            raise ValueError(f'Duplicate override: {old_key}')
        if prepared[key] == 'CYX':
            raise ValueError('Do not titrate disulfide CYX residues')
        base = canonical(prepared[key])
        allowed = {state for state in STATES.get(base, ()) if state}
        if base == 'HIS':
            allowed.update(HISTIDINE_HYDROGENS)
        state = str(item['state']).upper()
        if state not in allowed:
            raise ValueError(f'Unsupported override {state} for {old_key} ({base})')
        overrides[key] = state
    return mapping, overrides


def terminal_residues(lines):
    """Find protein segment endpoints using TER records and chain changes."""
    terminals, segment = set(), []
    for line in lines:
        if line.startswith('TER'):
            if segment:
                terminals.update((segment[0], segment[-1]))
                segment = []
        elif line.startswith(ATOM_RECORDS) and (
                canonical(line[17:20].strip()) in PROTEIN_NAMES
                or line[17:20].strip() in PHOSPHO_RESIDUES):
            key = residue_key(line)
            if segment and key[0] != segment[-1][0]:
                terminals.update((segment[0], segment[-1]))
                segment = []
            if not segment or key != segment[-1]:
                segment.append(key)
    if segment:
        terminals.update((segment[0], segment[-1]))
    return terminals


def residue_label(key, name):
    chain, number, icode = key
    return f'{name} {chain or "_"}{number}{icode}'


def row_key(row):
    return row['prepared_chain'], row['prepared_residue'], row['prepared_icode']


def choose_states(prepared, indices, pkas, overrides, old_keys, config, logger,
                  coupled=None, details=None):
    """Apply overrides first, then select candidate states from pH and pKa.

    Neutral histidines stay unresolved until PDB2PQR chooses HID or HIE.
    """
    decisions = {}
    targets = set()
    rows = []
    errors = []
    ph = config['ph']
    for key, name in prepared.items():
        base = canonical(name)
        if base not in STATES or name == 'CYX':
            continue
        index = indices[key]
        pka = pkas.get(index)
        detail = (details or {}).get(index, {})
        suppressed_by = detail.get('suppressed_by', '')
        state, reason = None, 'PROPKA calculated prediction'
        label = residue_label(old_keys[key], name)
        if key in overrides:
            state, reason = overrides[key], 'manual override'
        elif pka is None:
            reason = 'missing prediction'
            errors.append(f'{label}: no pKa prediction; supply override')
        else:
            protonated = ph < pka
            state = STATES[base][0 if protonated else 1]
            if base == 'HIS' and not protonated:
                targets.add(key)
                reason = 'pending PDB2PQR optimization'
            elif state is None:
                reason = f'unsupported deprotonated {base}'
                errors.append(f'{label}: {reason}')
        if state is not None:
            decisions[key] = state
        rows.append({
            'original_chain': old_keys[key][0], 'original_residue': old_keys[key][1],
            'original_icode': old_keys[key][2], 'prepared_chain': key[0],
            'prepared_residue': key[1], 'prepared_icode': key[2],
            'input_state': name, 'pka': pka, 'assigned_state': state,
            'reason': reason, 'near_pka': pka is not None and
            abs(ph - pka) <= config['review_pka_window'], 'target_ph': ph,
            'coupled': key in (coupled or set()),
            'delta_pka': abs(ph - pka) if pka is not None else None,
            'flags': '', 'propka_index': index,
            'pka_source': ('PROPKA calculated (suppressed in text)' if suppressed_by else
                           'PROPKA calculated'),
            'alternative_pka': None,  # Kept for existing CSV readers; no second calculation.
            'assigned_charge': STATE_CHARGES.get(state),
            'chain_mismatch': False, 'tautomer_mismatch': False,
            'decision_required': False, 'decision_status': '', 'decision_reasons': '',
            'pka_status': ('calculated-suppressed' if suppressed_by else
                           'calculated' if pka is not None else 'missing'),
            'suppressed_by': suppressed_by,
            'covalent_partners': '; '.join(detail.get('covalent_partners', [])),
        })
    return decisions, targets, rows, errors


def flag_rows(rows, surrogate_targets=(), surrogate_names=None):
    """Flag unusual states and number-matched differences between chains.

    Matching residue numbers does not establish equivalent sites. These flags
    never force symmetry or reject a chemically valid state.
    """
    groups = {}
    flags = {}
    for row in rows:
        key = row_key(row)
        tags = []
        if row['near_pka']:
            tags.append('near-pka')
        if row['assigned_state'] in UNUSUAL_STATES:
            tags.append('nonstandard-state')
        if row['coupled']:
            tags.append('coupled')
        if row.get('pka_status') == 'calculated-suppressed':
            tags.append('suppressed-covalent-coupling')
        if row['assigned_state'] is None:
            tags.append('unresolved')
        if key in surrogate_targets:
            tags.append('phospho-surrogate-approximation')
            if 'PTR' in (surrogate_names if surrogate_names is not None else {'PTR'}):
                tags.append('PTR-surrogate-approximation')  # compatibility with older CSV consumers
        row['assigned_charge'] = STATE_CHARGES.get(row['assigned_state'])
        row['chain_mismatch'] = row['tautomer_mismatch'] = False
        flags[key] = tags
        groups.setdefault((row['original_residue'], row['original_icode']), []).append(row)
    for group in groups.values():
        if (len({r['original_chain'] for r in group}) < 2 or
                len({canonical(r['input_state']) for r in group}) != 1):
            continue
        assigned = [r for r in group if r['assigned_state'] is not None]
        charge_difference = len({r['assigned_charge'] for r in assigned}) > 1
        tautomer_difference = (not charge_difference and
                               {r['assigned_state'] for r in assigned} == {'HID', 'HIE'})
        for row in assigned:
            row['chain_mismatch'] = charge_difference
            row['tautomer_mismatch'] = tautomer_difference
            if charge_difference or tautomer_difference:
                peers = ','.join(f'{r["original_chain"] or "_"}:{r["assigned_state"]}'
                                 for r in assigned if r is not row)
                kind = 'charge-difference' if charge_difference else 'tautomer-difference'
                flags[row_key(row)].append(f'heuristic-{kind}({peers})')
    for row in rows:
        row['flags'] = ' '.join(flags[row_key(row)])


def readable_residue(row):
    """Format the original residue ID for the console."""
    base = canonical(row['input_state'])
    key = row['original_chain'], row['original_residue'], row['original_icode']
    return residue_label(key, RESIDUE_NAMES.get(base, base))


def console_message(logger, level, message):
    """Wrap console text to a readable width."""
    for line in textwrap.wrap(message, width=86, break_long_words=True,
                              break_on_hyphens=False, subsequent_indent='    ') or ['']:
        getattr(logger, level)(line)


class ProtonationReviewRequired(ValueError):
    """A residue choice or experimental-method opt-in is pending."""


def review_decisions(rows, config):
    """Mark choices that need an override in manual mode.

    Keeping the predicted state counts as a choice. Original uncertainty flags
    remain, even after an explicit override.
    """
    pending = []
    for row in rows:
        reasons = []
        if row['near_pka']:
            reasons.append('predicted pKa is within the review window')
        if row['coupled']:
            reasons.append('coupled protonation site: review interacting states together')
        if row.get('suppressed_by'):
            reasons.append('calculated side-chain pKa was suppressed in the PROPKA text '
                           f'report due to coupling with {row["suppressed_by"]}')
        if row['pka'] is None:
            reasons.append('no pKa prediction')
        if row['assigned_state'] is None:
            reasons.append('candidate assignment is unresolved')
        explicit = row['reason'] == 'manual override'
        row['decision_reasons'] = '; '.join(reasons)
        row['decision_required'] = (config.get('decision_mode', 'automatic') == 'manual'
                                    and bool(reasons) and not explicit)
        if explicit:
            row['decision_status'] = 'explicit override'
        elif row['decision_required']:
            row['decision_status'] = 'awaiting explicit choice'
            pending.append(row)
        elif row['assigned_state'] is None:
            row['decision_status'] = 'unresolved'
        else:
            row['decision_status'] = 'automatic selection'
    return pending


def write_review_status(path, rows, config, pending, errors, surrogate_used,
                        acknowledgment_required, execution_status, phospho_residues=None):
    """Record pending choices and completion status without claiming chemical validation."""
    phospho_residues = phospho_residues or []
    uses_ptr = surrogate_used and any(row['name'] == 'PTR' for row in phospho_residues)
    record = {
        'scope': 'protein side-chain protonation assignment only',
        'decision_mode': config['decision_mode'],
        'target_ph': config['ph'],
        'review_pka_window': config['review_pka_window'],
        'execution_status': execution_status,
        'chemical_accuracy_validated': False,
        'scientific_status': ('review_required' if any(r['flags'] for r in rows)
                              or phospho_residues else 'not_validated'),
        'production_readiness': 'not_assessed',
        'assignment_policy': 'calculated PROPKA group candidates, including values suppressed '
                             'in the text report; not proof of a preferred joint state',
        'pka_geometry_source': 'prepared structure before PDB2PQR; no self-consistent iteration',
        'phosphate_state_source': 'input residue names; not predicted or validated here',
        'structural_check_scope': ('residue mapping; optimizer atom identity and permitted '
                                   'movements only when optimization runs'),
        'pending_decisions': [
            {'chain': r['original_chain'], 'residue': r['original_residue'],
             'insertion_code': r['original_icode'], 'candidate_state': r['assigned_state'],
             'pka': r['pka'], 'alternative_pka': r['alternative_pka'],
             'reason': r['decision_reasons']}
            for r in pending
        ],
        'explicit_overrides': config['overrides'],
        'phospho_surrogate_used_for_automatic_histidines': surrogate_used,
        'phospho_approximation_acknowledged': config['acknowledge_phospho_approximation'],
        'phospho_acknowledgment_required': acknowledgment_required,
        'phosphorylated_residues': phospho_residues,
        'acknowledgment_is_scientific_validation': False,
        'residues_with_review_flags': sum(bool(r['flags']) for r in rows),
        'errors': list(errors),
    }
    output = path.with_suffix('.protonation.review.json')
    output.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    return output


def log_reviews(rows, config, report, template, logger):
    """Print counts and blockers; keep residue-level explanations in the CSV."""
    unresolved = [r for r in rows if r['assigned_state'] is None]
    flagged = sum(bool(r['flags']) for r in rows)
    neutral = sum(r['assigned_state'] in {'HID', 'HIE'} for r in rows)
    positive = sum(r['assigned_state'] == 'HIP' for r in rows)
    console_message(logger, 'info',
                    f'Protonation pH {config["ph"]:.2f} ({config["decision_mode"]}): '
                    f'{len(rows)} sites; histidines {neutral} neutral, {positive} positive.')
    if flagged:
        categories = {
            'near-pKa': sum(bool(r['near_pka']) for r in rows),
            'coupled': sum(bool(r['coupled']) for r in rows),
            'suppressed-in-text': sum(bool(r.get('suppressed_by')) for r in rows),
            'unusual-state': sum(r['assigned_state'] in UNUSUAL_STATES for r in rows),
            'chain-difference': sum(r['chain_mismatch'] or r['tautomer_mismatch'] for r in rows),
            'phosphate-approximation': sum('phospho-surrogate-approximation' in r['flags']
                                           for r in rows),
        }
        counts = ', '.join(f'{name}={count}' for name, count in categories.items() if count)
        console_message(logger, 'warning', f'{flagged} sites have review flags'
                        + (f' ({counts})' if counts else '') + '; categories may overlap. See CSV.')
    if unresolved:
        labels = ', '.join(readable_residue(r) for r in sorted(
            unresolved, key=lambda r: (r['original_chain'], r['original_residue'], r['original_icode'])))
        console_message(logger, 'error', 'Unresolved: ' + labels + '.')
    console_message(logger, 'info', f'Assignments: {report.resolve()}')
    console_message(logger, 'info', f'Overrides template (not applied): {template.name}')


def write_overrides_template(path, rows, original_name, config):
    """Write commented entries using original residue IDs.

    Users must select and uncomment entries themselves. CHOOSE_STATE deliberately
    cannot be used unchanged for an unresolved assignment.
    """
    output = path.with_suffix('.overrides.yml')
    text = [
        f'# Original input: {original_name}; pH {config["ph"]}; mode {config["decision_mode"]}.',
        '# Review candidates in the CSV; choosing a state is not chemical validation.',
        '# IDs refer to the ORIGINAL PDB. Select states and uncomment entries below.',
        '# Copy selected entries into protonation.overrides in params.yml.',
        '# In manual mode, decision_required=True needs an explicit choice, even if unchanged.',
        '# CHOOSE_STATE must be replaced. All entries are initially commented out.',
        'protonation:', '  overrides: []',
        '    # Replace [] with an indented list of selected entries below.',
    ]
    for row in rows:
        if not row['flags'] and 'PDB2PQR' not in row['reason']:
            continue
        base = canonical(row['input_state'])
        allowed = {s for s in STATES[base] if s}
        if base == 'HIS':
            allowed.update(HISTIDINE_HYDROGENS)
        old = row['original_chain'], row['original_residue'], row['original_icode']
        state = row['assigned_state']
        text.append(f'    # {residue_label(old, base)} assigned={state or "UNRESOLVED"} '
                    f'pKa={row["pka"]} [{row["flags"]}] allowed={", ".join(sorted(allowed))} '
                    f'decision_required={row.get("decision_required", False)} '
                    f'status={row.get("decision_status", "")}')
        if row.get('pka_status') == 'calculated-suppressed':
            text.append(f'    # Calculated pKa={row["pka"]:.2f} was suppressed in text; '
                        f'coupled to {row["suppressed_by"]}. Not an independent prediction.')
        state_text = json.dumps(state) if state else 'CHOOSE_STATE'
        text.append('    # - {chain: ' + json.dumps(old[0]) +
                    f', residue: {old[1]}, insertion_code: ' + json.dumps(old[2]) +
                    ', state: ' + state_text + '}')
    output.write_text('\n'.join(text) + '\n')
    return output


def xyz(line):
    values = tuple(float(line[i:i+8]) for i in (30, 38, 46))
    if not all(math.isfinite(v) for v in values):
        raise ValueError(f'Invalid coordinates at {residue_key(line)}')
    return values


def pdb2pqr_plan(lines, prepared, terminals, old_keys, allow_phospho_approximation=False):
    """Prepare the histidine input without silently dropping components.

    Return temporary names, withheld phosphate atoms, protein atoms to transfer,
    and the complete input inventory. Incomplete parents stop here.
    """
    atoms = {}
    submitted = {}
    withheld = {}
    transferable = {}
    unsupported = []
    for key, name in prepared.items():
        if canonical(name) in PROTEIN_NAMES:
            submitted[key] = name
        elif name in WATER_NAMES:
            submitted[key] = 'HOH'
        elif name in NUCLEIC_NAMES:
            submitted[key] = name
        elif name in PHOSPHO_HISTIDINES:
            raise ValueError(
                f'{residue_label(old_keys[key], name)}: automatic histidine optimization '
                'with phosphohistidine is unsupported; no HIP surrogate is permitted. '
                'Use reviewed neutral-HIS overrides or a method retaining the full chemistry.')
        elif name in SURROGATES:
            submitted[key] = SURROGATES[name]
        else:
            unsupported.append(residue_label(old_keys[key], name))
    if unsupported:
        raise ValueError('PDB2PQR cannot type: ' + ', '.join(unsupported) +
                         '. Use a compatible preparation method or reviewed neutral-HIS '
                         'overrides from .overrides.yml. Do not remove required components '
                         'merely to bypass this check.')
    for line in lines:
        if not line.startswith(ATOM_RECORDS):
            continue
        key = residue_key(line)
        label = residue_label(old_keys[key], prepared[key])
        if line[16].strip():
            raise ValueError(f'{label}: resolve alternate conformations')
        if is_hydrogen(line):
            raise ValueError(f'{label}: automatic stage expects hydrogen-free input')
        atom_key = key, line[12:16].strip()
        if atom_key in atoms:
            raise ValueError(f'Duplicate input atom: {label} {atom_key[1]}')
        xyz(line)
        atoms[atom_key] = line
        if canonical(prepared[key]) in PROTEIN_NAMES:
            transferable[atom_key] = line
    for key, name in prepared.items():
        if name not in SURROGATES:
            continue
        present = {atom for (residue, atom) in atoms if residue == key}
        aliases = SURROGATE_ATOM_ALIASES.get(name, {})
        parent_name = SURROGATES[name]
        parent = PARENT_ATOMS[parent_name] | ({'OXT'} if key in terminals and 'OXT' in present else set())
        normalized = [PHOSPHATE_ALIASES.get(aliases.get(atom, atom), aliases.get(atom, atom))
                      for atom in present]
        if len(set(normalized)) != len(normalized):
            raise ValueError(f'{residue_label(old_keys[key], name)}: duplicate parent atom aliases')
        normalized = set(normalized)
        extra = normalized - parent
        if not parent <= normalized or extra != {'P', 'O1P', 'O2P', 'O3P'}:
            raise ValueError(f'{residue_label(old_keys[key], name)}: expected complete '
                             f'{parent_name} parent plus P and three phosphate oxygens; '
                             f'missing={sorted(parent-normalized)}, extra={sorted(extra)}')
        withheld[key] = {atom for atom in present if atom in PHOSPHATE_ATOMS}
    if withheld and not allow_phospho_approximation:
        raise ProtonationReviewRequired(
            'Phosphate removal is disabled. Experimental Ser/Thr/Tyr parents require '
            'acknowledge_phospho_approximation: true, or use reviewed histidine overrides.')
    return submitted, withheld, transferable, atoms


def near_withheld(lines, targets, withheld, cutoff=6.0):
    """Find close contacts with omitted phosphate atoms; this is not a safety test."""
    target_atoms, removed = {}, []
    for line in lines:
        if not line.startswith(ATOM_RECORDS):
            continue
        key = residue_key(line)
        if key in targets:
            target_atoms.setdefault(key, []).append(xyz(line))
        if line[12:16].strip() in withheld.get(key, set()):
            removed.append((key, xyz(line)))
    hits = {}
    for key, points in target_atoms.items():
        for other, point in removed:
            distance = min(math.dist(p, point) for p in points)
            if distance <= cutoff:
                hits[key, other] = min(distance, hits.get((key, other), math.inf))
    return hits


def nucleic_atom_name(name):
    """Normalize only known phosphate and sugar aliases for the retention check."""
    return {'O1P': 'OP1', 'O2P': 'OP2', 'O3P': 'OP3'}.get(name, name).replace('*', "'")


def nucleic_base(name):
    if name in NUCLEIC_NAMES:
        return name[:-1] if name[-1] in '35N' else name
    return name


def read_pdb2pqr(text, prepared, indices, transferable, submitted, original_atoms):
    """Check protein and environment atoms before accepting the optimizer output.

    DNA retention is checked directly, not inferred from a returned histidine state.
    The 0.002-angstrom tolerance allows coordinate rounding, not rearrangement.
    """
    def environment_atom_name(key, name):
        if prepared[key] in NUCLEIC_NAMES:
            return nucleic_atom_name(name)
        if prepared[key] in WATER_NAMES and name in {"O", "OW"}:
            return "O"
        return name
    by_index = {index: key for key, index in indices.items()}
    if len(by_index) != len(indices):
        raise ValueError('Predictor indices must be globally unique')
    protein_keys = {key for key, _ in transferable}
    states = {}
    hydrogens = {}
    coords = {}
    environment = {}
    seen = set()
    for line in text.splitlines():
        if not line.startswith(ATOM_RECORDS):
            continue
        index = int(line[22:26])
        if index not in by_index:
            raise ValueError(f'Unexpected PDB2PQR residue number: {index}')
        key = by_index[index]
        if line[21].strip() != key[0] or line[26].strip():
            raise ValueError(f'PDB2PQR changed residue identity at {key}')
        if line[16].strip():
            raise ValueError(f'PDB2PQR returned ambiguous alternate location at {key}')
        state = line[17:20].strip()
        requested = submitted[key]
        if key in protein_keys:
            valid = canonical(state) == canonical(prepared[key])
        elif prepared[key] in NUCLEIC_NAMES:
            valid = nucleic_base(state) == nucleic_base(requested)
        elif prepared[key] in WATER_NAMES:
            valid = state in WATER_NAMES
        else:
            valid = state == requested
        if not valid:
            raise ValueError(f'PDB2PQR changed residue identity at {key}: {requested} -> {state}')
        if key in states and states[key] != state:
            raise ValueError(f'Ambiguous PDB2PQR state at {key}')
        states[key] = state
        name = line[12:16].strip()
        if prepared[key] in WATER_NAMES and is_hydrogen(line):
            continue
        atom_key = key, name
        if atom_key in seen:
            raise ValueError(f'Duplicate PDB2PQR atom: {atom_key}')
        seen.add(atom_key)
        if is_hydrogen(line):
            hydrogens.setdefault(key, set()).add(name)
            continue
        xyz(line)
        if key in protein_keys:
            coords[atom_key] = line[30:54]
        else:
            norm = environment_atom_name(key, name)
            env_key = key, norm
            if env_key in environment:
                raise ValueError(f'Duplicate normalized environment atom: {env_key}')
            environment[env_key] = line
    if coords.keys() != transferable.keys():
        missing = sorted(transferable.keys() - coords.keys())
        extra = sorted(coords.keys() - transferable.keys())
        raise ValueError(f'PDB2PQR changed protein heavy atoms; missing={missing}, extra={extra}')
    expected = {}
    for (key, name), line in original_atoms.items():
        if key in protein_keys:
            continue
        if prepared[key] in SURROGATES:
            norm = SURROGATE_ATOM_ALIASES.get(prepared[key], {}).get(name, name)
            if norm not in PARENT_ATOMS[submitted[key]] | {'OXT'}:
                continue
        else:
            norm = environment_atom_name(key, name)
        if (key, norm) in expected:
            raise ValueError(f'Duplicate normalized input environment atom: {(key, norm)}')
        expected[key, norm] = line
    if environment.keys() != expected.keys():
        raise ValueError('PDB2PQR dropped/added/renamed environment heavy atoms; '
                         f'missing={sorted(expected.keys()-environment.keys())}, '
                         f'extra={sorted(environment.keys()-expected.keys())}')
    # Environment coordinates are retained from the original, so require consistency.
    for atom_key, line in expected.items():
        if math.dist(xyz(line), xyz(environment[atom_key])) > 0.002:
            raise ValueError(f'PDB2PQR moved an environment atom: {atom_key}')
    return states, hydrogens, coords


def _run_pdb2pqr(input_text, path):
    """Run hydrogen optimization only; keep the submitted input and tool output."""
    # Keep the actual predictor input so the temporary chemistry can be inspected.
    path.with_suffix('.pdb2pqr.input.pdb').write_text(input_text)
    with tempfile.TemporaryDirectory(prefix='pdb2pqr_', dir=path.parent) as folder:
        work = Path(folder)
        (work / 'input.pdb').write_text(input_text)
        # Do not rerun PROPKA on the surrogate or overwrite explicit charge decisions.
        # Unpatched HIS is neutral after PDB2PQR H-bond optimization; validate below.
        command = [sys.executable, '-c', 'from pdb2pqr.main import main; main()',
                   '--ff=AMBER', '--ffout=AMBER', '--keep-chain', '--nodebump',
                   '--pdb-output=optimized.pdb', 'input.pdb', 'optimized.pqr']
        run_predictor(command, work, path.with_suffix('.pdb2pqr.log'),
                      'Automatic histidine assignment failed')
        output = work / 'optimized.pdb'
        if not output.is_file():
            raise RuntimeError('PDB2PQR did not write optimized.pdb')
        text = output.read_text()
        path.with_suffix('.pdb2pqr.pdb').write_text(text)
    return text


def _audit_coordinate_changes(path, original_atoms, coords, prepared, old_keys, logger):
    """Check movements supported by this implementation; this is not an energy test."""
    changes = []
    for atom_key, coordinate in coords.items():
        before = xyz(original_atoms[atom_key])
        after = tuple(float(coordinate[i:i+8]) for i in (0, 8, 16))
        displacement = math.dist(before, after)
        if displacement > 0.002:
            key, atom = atom_key
            if atom not in FLIPPABLE_ATOMS.get(canonical(prepared[key]), set()):
                raise ValueError(f'PDB2PQR movement is outside the supported flip checks: {atom_key}')
            changes.append((*old_keys[key], prepared[key], atom, displacement))
    audit = path.with_suffix('.pdb2pqr.changes.csv')
    with audit.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['original_chain', 'original_residue', 'original_icode',
                         'residue_name', 'atom', 'displacement_angstrom'])
        writer.writerows(changes)
    moved = {entry[:3] for entry in changes}
    console_message(logger, 'warning' if changes else 'info',
                    f'PDB2PQR moved {len(changes)} heavy atoms in {len(moved)} residues '
                    f'by more than 0.002 angstroms; audit: {audit.name}.')


def pdb2pqr_histidines(lines, prepared, indices, decisions, targets, path, config,
                      logger, old_keys):
    """Choose neutral histidine tautomers and audit any side-chain flips.

    PDB2PQR may use a default tautomer. A returned state does not prove a unique
    hydrogen-bond arrangement; validate the ring hydrogens and preserve its geometry.
    """
    submitted, withheld, original_atoms, all_atoms = pdb2pqr_plan(
        lines, prepared, terminal_residues(lines), old_keys,
        allow_phospho_approximation=config.get('acknowledge_phospho_approximation', False))
    prediction_states = dict(submitted)
    for key in {key for key, _ in original_atoms}:
        prediction_states[key] = 'HIS' if key in targets else decisions.get(key, prepared[key])
    if withheld:
        labels = ', '.join(
            f'{residue_label(old_keys[key], prepared[key])} -> {SURROGATES[prepared[key]]}'
            for key in sorted(withheld, key=lambda k: old_keys[k]))
        console_message(logger, 'warning', f'Histidine step uses temporary parents: {labels}.')
        console_message(logger, 'warning', 'Phosphate atoms are omitted for this step only; '
                        'the effect is unvalidated. Distance alone does not establish safety.')
        console_message(logger, 'info', 'Final phosphorylated residues will retain their original atoms and names.')
    for (key, other), distance in near_withheld(lines, targets, withheld).items():
        console_message(logger, 'warning', f'{residue_label(old_keys[key], "Histidine")}: '
                        f'{distance:.2f} angstroms from omitted phosphate at '
                        f'{residue_label(old_keys[other], prepared[other])}. Check proton position.')
    atom_aliases = {key: SURROGATE_ATOM_ALIASES.get(prepared[key], {}) for key in withheld}
    input_text = predictor_input(lines, indices, prediction_states, drop=withheld,
                                 atom_aliases=atom_aliases)
    text = _run_pdb2pqr(input_text, path)
    states, hydrogens, coords = read_pdb2pqr(
        text, prepared, indices, original_atoms, prediction_states, all_atoms)
    for key in {key for key, _ in original_atoms}:
        requested = prediction_states[key]
        if key not in targets and states.get(key) != requested:
            raise ValueError(f'PDB2PQR changed requested state at {key}: '
                             f'{requested} -> {states.get(key)}')
    assignments = {}
    for key in targets:
        state = states.get(key)
        if state not in HISTIDINE_HYDROGENS:
            raise ValueError(f'PDB2PQR did not return neutral HID/HIE for {key}: {state}')
        ring_h = hydrogens.get(key, set()) & {'HD1', 'HE2'}
        if ring_h != HISTIDINE_HYDROGENS[state]:
            raise ValueError(f'PDB2PQR state/hydrogens disagree at {key}')
        assignments[key] = state
    console_message(logger, 'info', f'PDB2PQR assigned {len(assignments)} neutral histidines '
                    '(HID: proton on ND1; HIE: proton on NE2).')
    console_message(logger, 'info', 'Default tautomer choices remain possible; see CSV.')
    _audit_coordinate_changes(
        path, original_atoms, coords, prepared, old_keys, logger)
    return assignments, coords


def write_report(path, rows):
    report = path.with_suffix('.protonation.csv')
    with report.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return report


def write_assigned_pdb(path, lines, decisions, coords):
    """Write selected states and audited protein coordinates; back up the input first."""
    output = []
    for line in lines:
        if line.startswith(ATOM_RECORDS):
            key = residue_key(line)
            if key in decisions:
                line = line[:17] + f'{decisions[key]:>3}' + line[20:]
            atom_key = key, line[12:16].strip()
            if atom_key in coords:
                line = line[:30] + coords[atom_key] + line[54:]
        output.append(line)
    path.with_suffix('.unassigned.pdb').write_text(''.join(lines))
    path.write_text(''.join(output))


def _require_manual_choices(pending, acknowledgment_required, report, template,
                            status_path, logger):
    """Stop before writing when choices or phosphate-method opt-in are pending."""
    if not pending and not acknowledgment_required:
        return
    if pending:
        console_message(logger, 'warning',
                        f'Manual mode: {len(pending)} sites need explicit choices: '
                        + ', '.join(readable_residue(r) for r in pending) + '.')
        console_message(logger, 'warning',
                        f'Copy reviewed entries from {template.name} into '
                        'protonation.overrides in params.yml; keeping a candidate counts as a choice.')
    if acknowledgment_required:
        console_message(
            logger, 'warning',
            'Phosphate-surrogate histidine optimization is disabled; '
            'the required acknowledgment has not been given. '
            'To enable it, set protonation.acknowledge_phospho_approximation '
            'to true in params.yml. This accepts an unvalidated approximation. '
            'Alternatively, supply reviewed neutral-HIS overrides or use a '
            'preparation method retaining phosphate chemistry.'
        )
    raise ProtonationReviewRequired(
        'Protonation preparation needs explicit choices or method opt-in; '
        'the assigned structure was not written. '
        f'Report: {report}; overrides template: {template}; status: {status_path}')


def assign_protonation(pdbfile, original_pdb, config, logger):
    """Prepare candidate states, write review reports, then commit the structure.

    The public entry point and return path are unchanged. Manual review stops
    before writing the assigned PDB. Early predictor/mapping failures may have
    only a tool log rather than a complete report.
    """
    config = normalize_config(config)
    if config['method'] == 'none':
        return str(pdbfile)
    path = Path(pdbfile).resolve()
    lines = path.read_text().splitlines(keepends=True)
    prepared_residues = residues(lines)
    original_residues = residues(Path(original_pdb).read_text().splitlines(keepends=True))
    mapping, overrides = prepare_overrides(original_residues, prepared_residues, config['overrides'])
    original_keys = {key: old_key for old_key, key in mapping.items()}
    indices = {key: index for index, key in enumerate(prepared_residues, 1)}
    if len(indices) > 9999:
        raise ValueError('Predictor PDB exceeds 9999 residues; mmCIF support needed')
    # Remove stale files so they cannot be mistaken for this run's results.
    for suffix in ('.pdb2pqr.input.pdb', '.pdb2pqr.pdb', '.pdb2pqr.log',
                   '.pdb2pqr.changes.csv', '.protonation.review.json'):
        path.with_suffix(suffix).unlink(missing_ok=True)
    pka_values, coupled_groups, pka_details = run_propka(
        lines, prepared_residues, indices, path, logger)
    coupled = {key for key, name in prepared_residues.items()
               if (canonical(name), indices[key], key[0]) in coupled_groups and name != 'CYX'}
    decisions, histidine_targets, rows, errors = choose_states(
        prepared_residues, indices, pka_values, overrides, original_keys, config,
        logger, coupled, pka_details)
    terminal_states = [residue_label(original_keys[key], decisions[key])
                       for key in terminal_residues(lines)
                       if decisions.get(key) in UNUSUAL_STATES]
    if terminal_states:
        console_message(logger, 'warning', 'Selected terminal side-chain states: '
                        + ', '.join(sorted(terminal_states)) + '. LEaP must support the '
                        'corresponding terminal templates; this step does not validate them.')
    phospho_names = set(prepared_residues.values()) & PHOSPHO_RESIDUES.keys()
    unsupported_phosphohistidines = phospho_names & PHOSPHO_HISTIDINES
    needs_phospho_surrogate = bool(histidine_targets and phospho_names & SURROGATES.keys())
    acknowledgment_required = False
    protein_coordinates = {}
    if histidine_targets and not errors and unsupported_phosphohistidines:
        labels = ', '.join(residue_label(original_keys[key], name)
                           for key, name in prepared_residues.items()
                           if name in unsupported_phosphohistidines)
        errors.append('Automatic histidine optimization cannot retain the chemistry of '
                      f'{labels}. No phosphohistidine-to-HIP replacement is allowed. '
                      'Use reviewed neutral-HIS overrides or a method that keeps the full chemistry.')
    if histidine_targets and not errors and needs_phospho_surrogate:
        acknowledgment_required = not config['acknowledge_phospho_approximation']
    if histidine_targets and not errors and not acknowledgment_required:
        try:
            histidine_states, protein_coordinates = pdb2pqr_histidines(
                lines, prepared_residues, indices, decisions, histidine_targets,
                path, config, logger, original_keys)
            decisions.update(histidine_states)
            for row in rows:
                if row_key(row) in histidine_states:
                    row['assigned_state'] = histidine_states[row_key(row)]
                    row['reason'] = 'PDB2PQR hydrogen-bond optimization'
        except (ValueError, RuntimeError, OSError) as exc:
            errors.append(str(exc))
            protein_coordinates = {}
            for row in rows:
                if row_key(row) in histidine_targets:
                    row['assigned_state'] = None
                    row['reason'] = 'failed PDB2PQR optimization: ' + str(exc)
    elif histidine_targets:
        reason = ('experimental phosphate simplification is not enabled'
                  if acknowledgment_required else 'earlier preparation errors')
        for row in rows:
            if row_key(row) in histidine_targets:
                row['reason'] = 'skipped PDB2PQR optimization: ' + reason
    surrogate_used = needs_phospho_surrogate and any(
        row_key(row) in histidine_targets
        and row['reason'] == 'PDB2PQR hydrogen-bond optimization' for row in rows)
    # Flag all successful tautomer choices made under the altered environment.
    # A distance cutoff cannot determine which choices are unaffected.
    surrogate_targets = histidine_targets if surrogate_used else set()
    flag_rows(rows, surrogate_targets, phospho_names)
    phospho_records = [
        dict(name=name, original_chain=original_keys[key][0],
             original_residue=original_keys[key][1], original_icode=original_keys[key][2],
             **PHOSPHO_RESIDUES[name])
        for key, name in prepared_residues.items() if name in PHOSPHO_RESIDUES
    ]
    pending_choices = review_decisions(rows, config)
    report = write_report(path, rows)
    template = write_overrides_template(path, rows, Path(original_pdb).name, config)
    log_reviews(rows, config, report, template, logger)
    if errors:
        execution_status = 'failed'
    elif acknowledgment_required:
        execution_status = 'awaiting_approximation_opt_in'
    elif pending_choices:
        execution_status = 'awaiting_manual_decisions'
    else:
        execution_status = 'ready_to_write'
    status_path = write_review_status(path, rows, config, pending_choices, errors, surrogate_used,
                                      acknowledgment_required, execution_status, phospho_records)
    console_message(logger, 'info', f'Decision status: {status_path.name}')
    if errors:
        for error in errors:
            logger.error("%s", error)
        raise ValueError('Protonation needs review:\n' + '\n'.join(errors) +
                         f'\nReport: {report}\nOverrides template: {template}')
    _require_manual_choices(
        pending_choices, acknowledgment_required, report, template, status_path, logger)
    unresolved = [row for row in rows if row['assigned_state'] is None]
    if unresolved:
        raise ValueError(f'Unresolved protonation assignments; inspect {report}')
    # Record success only after the assigned PDB has been written.
    try:
        write_assigned_pdb(path, lines, decisions, protein_coordinates)
    except OSError as exc:
        write_review_status(path, rows, config, pending_choices, [str(exc)], surrogate_used,
                            acknowledgment_required, 'write_failed', phospho_records)
        raise
    write_review_status(path, rows, config, pending_choices, [], surrogate_used,
                        acknowledgment_required, 'assigned_structure_written', phospho_records)
    console_message(logger, 'info', 'Fixed side-chain protonation written; '
                    'check subsequent LEaP/topology results for system preparation status.')
    if any(row['flags'] for row in rows):
        console_message(logger, 'warning', 'Review flagged choices before production use.')
    return str(path)
