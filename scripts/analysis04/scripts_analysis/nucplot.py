#!/usr/bin/env python3
"""DNA-focused hydrogen bond plot: full-duplex NUCPLOT-style schematic only.

Adapted from https://github.com/RomanLas/NUCPLOT

It reuses data already produced by the hbonds pipeline and never recomputes bonds:

  - ``hbonds_{method}_ext.csv`` : one row per unique (donor, hydrogen, acceptor)
                                   bond with occupancy (``frame%``/``count``),
                                   residue info, ``*_type`` (protein/DNA) and
                                   ``*_section`` (for DNA atoms: ``base`` /
                                   ``backbone``).
  - a structure file             : any format readable by MDAnalysis, used to
                                   get the full DNA sequence (including
                                   uncontacted nucleotides) and true
                                   Watson-Crick pairing.

Output: one figure -> outfiles['hbonds_dna_nucplot_full']
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Ellipse
from matplotlib.lines import Line2D
from scripts.utils.helpers import setup_logger
import MDAnalysis as mda
import warnings


# ----- small helpers -----

def format_residue(chain, res_num):
    return f"{chain}:{res_num}"


def _dna_side_frame(df):
    """Return the protein-DNA subset of ``df`` with normalised DNA-side columns.

    Adds: dna_chain, dna_res_num, dna_res_name, dna_section, dna_label,
    prot_chain, prot_res_num, prot_res_name, prot_label, dna_full_label.
    """
    # Compare case-insensitively: the producer (hbonds_extend.py) writes
    # 'Protein'/'DNA' (capitalised), so a literal 'protein' match would find nothing.
    donor_type = df['donor_type'].astype(str).str.lower()
    acceptor_type = df['acceptor_type'].astype(str).str.lower()
    pd_mask = (
        ((donor_type == 'protein') & (acceptor_type == 'dna')) |
        ((donor_type == 'dna') & (acceptor_type == 'protein'))
    )
    sub = df[pd_mask].copy()
    if sub.empty:
        return sub

    donor_is_dna = sub['donor_type'].astype(str).str.lower() == 'dna'
    sub['dna_chain'] = np.where(donor_is_dna, sub['donor_chain'], sub['acceptor_chain'])
    sub['dna_res_num'] = np.where(
        donor_is_dna, sub['donor_res_num_adjusted'], sub['acceptor_res_num_adjusted']
    ).astype(int)
    sub['dna_res_name'] = np.where(donor_is_dna, sub['donor_res_name'], sub['acceptor_res_name'])
    sub['dna_section'] = np.where(donor_is_dna, sub['donor_section'], sub['acceptor_section'])
    sub['dna_atom'] = np.where(donor_is_dna, sub['donor_atom_name'], sub['acceptor_atom_name'])
    sub['dna_label'] = [format_residue(c, r) for c, r in zip(sub['dna_chain'], sub['dna_res_num'])]

    # ...and the protein side (the opposite end of every protein-DNA bond).
    sub['prot_chain'] = np.where(donor_is_dna, sub['acceptor_chain'], sub['donor_chain'])
    sub['prot_res_num'] = np.where(
        donor_is_dna, sub['acceptor_res_num_adjusted'], sub['donor_res_num_adjusted']
    ).astype(int)
    sub['prot_res_name'] = np.where(donor_is_dna, sub['acceptor_res_name'], sub['donor_res_name'])
    sub['prot_label'] = [
        f"{format_residue(c, r)} {n}" for c, r, n in
        zip(sub['prot_chain'], sub['prot_res_num'], sub['prot_res_name'])
    ]
    sub['dna_full_label'] = [
        f"{format_residue(c, r)} {n}" for c, r, n in
        zip(sub['dna_chain'], sub['dna_res_num'], sub['dna_res_name'])
    ]
    return sub


# Pastel fills per base identity.
BASE_COLORS = {'A': '#a6cee3', 'T': '#fb9a99', 'G': '#b2df8a', 'C': '#cab2d6', 'U': '#fdbf6f'}


def _base_letter(resname):
    r = str(resname).upper().lstrip('D')
    return r[0] if r else '?'


_COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
# H-bonds per Watson-Crick pair (for the rung drawing).
_WC_BONDS = {frozenset(('A', 'T')): 2, frozenset(('G', 'C')): 3}


def _read_dna_strands(pdb):
    """Return {chain: [(resid, base_letter), ...] sorted by resid} for all DNA."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        u = mda.Universe(pdb)
    strands = {}
    for seg in u.select_atoms('nucleic').segments:
        res = seg.atoms.select_atoms('nucleic').residues
        if len(res) == 0:
            continue
        strands[seg.segid] = sorted(
            [(int(r.resid), _base_letter(r.resname)) for r in res], key=lambda t: t[0]
        )
    return strands


def _pair_strands(s1, s2):
    """Watson-Crick pairing of two antiparallel strands by best complementary
    registration. Returns {resid_in_s1: resid_in_s2}. Robust to distorted frames
    because it uses sequence complementarity, not 3D distance."""
    r2 = list(reversed(s2))
    best = None
    for off in range(-len(r2) + 1, len(s1)):
        matches = 0
        for i in range(len(s1)):
            j = i - off
            if 0 <= j < len(r2) and _COMPLEMENT.get(s1[i][1]) == r2[j][1]:
                matches += 1
        if best is None or matches > best[0]:
            best = (matches, off)
    pairs = {}
    if best:
        off = best[1]
        for i in range(len(s1)):
            j = i - off
            if 0 <= j < len(r2) and _COMPLEMENT.get(s1[i][1]) == r2[j][1]:
                pairs[s1[i][0]] = r2[j][0]
    return pairs


def plot_nucplot_full(dna, pdb, outfile, logger,
                      min_occupancy=5.0):
    """Full-sequence NUCPLOT-style schematic with real Watson-Crick pairing.

    Draws every DNA nucleotide of both strands as an antiparallel ladder (true WC
    rungs: 2 lines A-T, 3 lines G-C), highlights the nucleotides that are actually
    contacted, and connects protein residues - coloured by protein chain, one chain
    above the duplex and the next below - to the exact base (solid line) or backbone
    (dashed line) atom they hydrogen-bond to. Width scales with occupancy; pairs
    below ``min_occupancy`` (%) are hidden.
    """


    strands = _read_dna_strands(pdb)
    csv_chains = set(dna['dna_chain'].unique())
    # Keep DNA chains that appear in the contacts (fall back to all if no overlap).
    chains = [c for c in sorted(strands) if c in csv_chains] or sorted(strands)
    if not chains:
        raise ValueError("No DNA strands found in structure for the schematic.")
    top_chain = chains[0]
    bottom_chain = chains[1] if len(chains) > 1 else None

    df = dna.copy()
    df['base_letter'] = df['dna_res_name'].map(_base_letter)
    agg = (df.groupby(['prot_label', 'prot_chain', 'dna_chain', 'dna_res_num',
                       'dna_section'])['frame%'].sum().reset_index())
    agg = agg[agg['frame%'] >= min_occupancy]
    contacted = set(zip(agg['dna_chain'], agg['dna_res_num']))

    # --- geometry ---
    dx = 1.0
    y_top, y_bot = 2.0, 0.0
    base_off, phos_off = 0.55, 0.55
    y_res_top, y_res_bot = y_top + phos_off + 1.5, y_bot - phos_off - 1.5

    pairs = _pair_strands(strands[top_chain], strands[bottom_chain]) if bottom_chain else {}

    # Column layout: top strand 5'->3' across columns; paired bottom nt sits beneath.
    nt_pos = {}
    col = 0
    bottom_placed = set()
    for resid, letter in strands[top_chain]:
        x = col * dx
        nt_pos[(top_chain, resid)] = {'x': x, 'base': (x, y_top - base_off),
                                      'phos': (x, y_top + phos_off), 'letter': letter,
                                      'strand': 'top', 'resid': resid}
        if resid in pairs:
            br = pairs[resid]
            bl = dict(strands[bottom_chain])[br]
            nt_pos[(bottom_chain, br)] = {'x': x, 'base': (x, y_bot + base_off),
                                          'phos': (x, y_bot - phos_off), 'letter': bl,
                                          'strand': 'bottom', 'resid': br}
            bottom_placed.add(br)
        col += 1
    # Any unpaired bottom nucleotides get appended columns.
    if bottom_chain:
        for resid, letter in strands[bottom_chain]:
            if resid in bottom_placed:
                continue
            x = col * dx
            nt_pos[(bottom_chain, resid)] = {'x': x, 'base': (x, y_bot + base_off),
                                             'phos': (x, y_bot - phos_off), 'letter': letter,
                                             'strand': 'bottom', 'resid': resid}
            col += 1
    width = max(1, col - 1) * dx

    fig, ax = plt.subplots(figsize=(max(12, col * 0.62), 9))

    # --- WC rungs ---
    if bottom_chain:
        for tr, br in pairs.items():
            t = nt_pos[(top_chain, tr)]['base']
            b = nt_pos[(bottom_chain, br)]['base']
            n = _WC_BONDS.get(frozenset((nt_pos[(top_chain, tr)]['letter'],
                                         nt_pos[(bottom_chain, br)]['letter'])), 2)
            for k in range(n):
                xo = (k - (n - 1) / 2) * 0.07
                ax.plot([t[0] + xo, b[0] + xo], [t[1], b[1]], color='0.6', lw=1.0, zorder=0)

    # --- strands: backbone, bases, phosphates ---
    for chain, ylevel in [(top_chain, y_top), (bottom_chain, y_bot)]:
        if chain is None:
            continue
        seq = strands[chain]
        phos = [nt_pos[(chain, r)]['phos'] for r, _ in seq if (chain, r) in nt_pos]
        ax.plot([p[0] for p in phos], [p[1] for p in phos], color='0.4', lw=2.5,
                zorder=1, solid_capstyle='round')
        for resid, _ in seq:
            nt = nt_pos[(chain, resid)]
            bx, by = nt['base']; px, py = nt['phos']
            ax.plot([bx, px], [by, py], color='0.6', lw=1.2, zorder=1)
            ax.add_patch(Circle((px, py), 0.16, facecolor='white', edgecolor='0.4',
                                lw=1.2, zorder=2))
            ax.text(px, py, 'P', ha='center', va='center', fontsize=6.5, color='0.3', zorder=3)
            is_hit = (chain, resid) in contacted
            fc = BASE_COLORS.get(nt['letter'], '#dddddd') if is_hit else '#f0f0f0'
            ec = '0.2' if is_hit else '0.7'
            box = FancyBboxPatch((bx - 0.34, by - 0.26), 0.68, 0.52,
                                 boxstyle='round,pad=0.02,rounding_size=0.08',
                                 facecolor=fc, edgecolor=ec,
                                 lw=1.4 if is_hit else 0.8, zorder=2)
            ax.add_patch(box)
            ax.text(bx, by, f"{nt['letter']}{resid}", ha='center', va='center',
                    fontsize=7.5, color='black' if is_hit else '0.5', zorder=3,
                    fontweight='bold' if is_hit else 'normal')
        # 5'/3' end labels (by residue number)
        resids = [r for r, _ in seq if (chain, r) in nt_pos]
        for end_resid, tag in [(min(resids), "5'"), (max(resids), "3'")]:
            ex = nt_pos[(chain, end_resid)]['phos'][0]
            ax.text(ex, ylevel + (0.95 if ylevel == y_top else -0.95),
                    f"{chain} {tag}", ha='center', va='center', fontsize=9, fontweight='bold')

    # --- protein residues: colour per chain, chain A above / chain B below ---
    prot_chains = sorted(agg['prot_chain'].unique())
    palette = plt.cm.tab10(np.linspace(0, 1, 10))
    prot_color = {c: palette[i % 10] for i, c in enumerate(prot_chains)}
    chain_side = {c: ('top' if i % 2 == 0 else 'bottom') for i, c in enumerate(prot_chains)}

    res_info = {}
    for r in agg['prot_label'].unique():
        sub = agg[agg['prot_label'] == r]
        xs = [nt_pos[(c, rn)]['x'] for c, rn in zip(sub['dna_chain'], sub['dna_res_num'])
              if (c, rn) in nt_pos]
        if xs:
            res_info[r] = {'x': float(np.mean(xs)), 'chain': sub['prot_chain'].iloc[0]}

    res_xy = {}
    for side, yv in [('top', y_res_top), ('bottom', y_res_bot)]:
        lst = sorted([r for r in res_info if chain_side[res_info[r]['chain']] == side],
                     key=lambda r: res_info[r]['x'])
        if not lst:
            continue
        xs = np.linspace(0, width, len(lst)) if len(lst) > 1 else [width / 2]
        for r, xx in zip(lst, xs):
            res_xy[r] = (xx, yv)

    for r, (xx, yy) in res_xy.items():
        chain, rest = r.split(':', 1)
        num, name = rest.split()
        col_fill = prot_color[chain]
        ax.add_patch(Ellipse((xx, yy), 0.95, 0.62, facecolor=col_fill, alpha=0.45,
                             edgecolor=col_fill, lw=1.6, zorder=2))
        ax.text(xx, yy, f"{name}\n{chain}{num}", ha='center', va='center',
                fontsize=7.5, zorder=3)

    # --- H-bond connections (colour = protein chain, style = base/backbone) ---
    occ_max = agg['frame%'].max() if not agg.empty else 1.0
    for _, row in agg.sort_values('frame%').iterrows():
        if row['prot_label'] not in res_xy:
            continue
        p = res_xy[row['prot_label']]
        nt = nt_pos.get((row['dna_chain'], row['dna_res_num']))
        if nt is None:
            continue
        target = nt['base'] if row['dna_section'] == 'base' else nt['phos']
        lw = 0.6 + 3.5 * (row['frame%'] / occ_max)
        ls = '-' if row['dna_section'] == 'base' else (0, (4, 2))
        ax.plot([p[0], target[0]], [p[1], target[1]], color=prot_color[row['prot_chain']],
                lw=lw, alpha=0.6, linestyle=ls, zorder=1)

    # --- legend ---
    legend = [Line2D([0], [0], color=prot_color[c], lw=3, label=f'Protein chain {c}')
              for c in prot_chains]
    legend += [
        Line2D([0], [0], color='0.3', lw=2, linestyle='-', label='-> base'),
        Line2D([0], [0], color='0.3', lw=2, linestyle=(0, (4, 2)), label='-> backbone'),
        Line2D([0], [0], color='gray', lw=0.6, label=f'Occupancy thin->thick = Min {min_occupancy:.0f}% - Max {occ_max:.0f}%'),
    ]
    ax.legend(handles=legend, loc='upper center', bbox_to_anchor=(0.5, -0.01),
              ncol=min(len(legend), 5), frameon=False, fontsize=9)
    ax.set_xlim(-1.5, width + 1.0)
    ax.set_ylim(y_res_bot - 0.6, y_res_top + 0.6)
    ax.set_aspect('equal')
    ax.set_axis_off()
    ax.set_title(f'Protein-DNA hydrogen-bond schematic '
                 f'(>= {min_occupancy:.0f}% occupancy)', fontsize=13)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

    if logger and Path(outfile).is_file():
        logger.info("DNA full-duplex NUCPLOT scheme saved as %s (%d pairs shown)", outfile, len(agg))
    else:
        logger.warning(f"DNA NUCPLOT scheme could not be correctly generated.")
        logger.warning(f"Check log {logfile}")

# ----------------------------------- main ------------------------------------

def main(logfile: str = None, hbonds_ext: str = None, 
        outfile: str = None, min_occupancy: float = 5.0,
        pdb: str = None):

    logger = setup_logger("analysis", logfile)

    if not pdb:
        raise ValueError("plot_nucplot_full requires a PDB structure_file.")

    hbonds_df = pd.read_csv(hbonds_ext)
    dna = _dna_side_frame(hbonds_df)

    if dna.empty:
        logger.warning(
            "DNA nucplot-full plot - no protein-DNA hydrogen bonds found in %s; "
            "skipping plot.", hbonds_ext
        )
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No protein-DNA hydrogen bonds", ha='center', va='center')
        ax.set_axis_off()
        plt.savefig(outfile)
        plt.close()
        return

    plot_nucplot_full(dna, pdb, outfile, logger, min_occupancy)


if __name__ == "__main__":
    main()