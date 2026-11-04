"""Text diagrams of nucleic-acid strands and their Watson-Crick pairs.

Each register (duplex or hairpin stem) is drawn as two aligned rows:

    position     1   2   3   4
    B 5'->3'     A   C   A   G   ...
                 -   |   :   |
    C 3'->5'         G   T   C   ...
    position         18  17  16

``|`` observed pair, ``:`` inferred pair. For a hairpin the top row is the 5' arm and
the bottom row the 3' arm (read right to left), so the loop sits at the right end.
"""

from .infer import nucleic_residues

MARK = {"observed": "|", "inferred": ":"}
COL = 4      # characters per column
LABEL = 14   # width of the row labels


def _position_map(residues):
    """One line mapping 1-based positions to residue numbers."""
    ids = [f"{r['res_id']}{r['ins_code']}" for r in residues]
    contiguous = all(r["ins_code"] == "" for r in residues) and all(
        b["res_id"] - a["res_id"] == 1 for a, b in zip(residues, residues[1:]))
    if contiguous:
        return f"positions 1-{len(ids)} = residues {ids[0]}-{ids[-1]}"
    return "position:residue " + " ".join(f"{k}:{i}" for k, i in enumerate(ids, start=1))


def _row(label, cells):
    return (label.ljust(LABEL) + "".join(c.ljust(COL) for c in cells)).rstrip()


def _draw_register(chains, marks, c1, c2, s, width):
    res1, res2 = chains[c1]["residues"], chains[c2]["residues"]
    if c1 != c2:
        top = list(range(len(res1)))
        bottom = list(range(len(res2)))
        kind = f"Duplex {c1}/{c2}"
    else:  # hairpin: split the chain at the middle of the register
        top = [k for k in range(len(res1)) if 2 * k <= s]
        bottom = [k for k in range(len(res1)) if 2 * k > s]
        kind = f"Hairpin in chain {c1}"

    cols = [k for k in top] + [s - k for k in bottom]
    first, last = min(cols), max(cols)
    top_at = {k: k for k in top}
    bottom_at = {s - k: k for k in bottom}

    lines = [f"### {kind} (register: position sum {s + 2})", ""]
    for start in range(first, last + 1, width):
        block = range(start, min(start + width, last + 1))
        pos_t, let_t, mark, let_b, pos_b = [], [], [], [], []
        for col in block:
            t, b = top_at.get(col), bottom_at.get(col)
            pos_t.append("" if t is None else str(t + 1))
            let_t.append("" if t is None else (res1[t]["letter"] or "?"))
            pos_b.append("" if b is None else str(b + 1))
            let_b.append("" if b is None else (res2[b]["letter"] or "?"))
            if t is None or b is None:
                mark.append("")
            else:
                mark.append(MARK.get(marks.get(((c1, t), (c2, b))), ""))
        lines += ["```",
                  _row("position", pos_t),
                  _row(f"{c1} 5'->3'", let_t),
                  _row("", mark),
                  _row(f"{c2} 3'->5'", let_b),
                  _row("position", pos_b),
                  "```", ""]
    return lines


def pairs_diagram(atoms, nucleic_chains, pairs, info, title="", width=25):
    """Markdown text with one aligned diagram per register.

    Parameters
    ----------
    atoms : AtomArray
        The model the pairs were computed on.
    nucleic_chains : list of str
        Nucleic-acid chain IDs (DNA, RNA or hybrid).
    pairs, info
        As returned by :func:`find_wc_pairs`.
    title : str
        Heading for the document (e.g. the PDB ID).
    width : int
        Columns per block; long duplexes are wrapped.
    """
    lines = [f"# {title} base pairs".strip(), ""]
    if not nucleic_chains:
        return "\n".join(lines + ["No nucleic-acid chains: no base-pair diagram.", ""])

    chains = nucleic_residues(atoms, nucleic_chains)
    index_of = {(r["chain"], r["res_id"], r["ins_code"]): k
                for c in chains.values() for k, r in enumerate(c["residues"])}
    marks = {}
    for row in pairs.itertuples():
        a = (row.chain1, index_of[(row.chain1, row.res_id1, row.ins_code1)])
        b = (row.chain2, index_of[(row.chain2, row.res_id2, row.ins_code2)])
        marks[(a, b)] = marks[(b, a)] = row.evidence

    n_obs = (pairs["evidence"] == "observed").sum() if len(pairs) else 0
    n_inf = (pairs["evidence"] == "inferred").sum() if len(pairs) else 0
    lines += [f"Method: branch {info['method']} | pairs: {len(pairs)} "
              f"({n_obs} observed `|`, {n_inf} inferred `:`)", ""]
    lines += ["Residue numbering:", ""]
    lines += [f"- chain {c}: {_position_map(v['residues'])}" for c, v in chains.items()]
    lines.append("")

    drawn = set()
    for reg in info["registers"]:
        if reg["n_observed"] < 2:
            continue
        c1, c2, s = reg["chain1"], reg["chain2"], reg["sum"]
        lines += _draw_register(chains, marks, c1, c2, s, width)
        drawn.add((c1, c2, s))

    # Pairs that belong to no drawn register (e.g. a single isolated pair).
    loose = []
    for (a, b), ev in marks.items():
        (c1, i), (c2, j) = sorted((a, b))
        if (a, b) == ((c1, i), (c2, j)) and (c1, c2, i + j) not in drawn:
            r1, r2 = chains[c1]["residues"][i], chains[c2]["residues"][j]
            loose.append(f"- {c1}{r1['res_id']}{r1['ins_code']} ({r1['letter']}) - "
                         f"{c2}{r2['res_id']}{r2['ins_code']} ({r2['letter']}): {ev}")
    if loose:
        lines += ["### Pairs outside any drawn register", ""] + loose + [""]

    unpaired = [c for c in chains if not any(a[0] == c for a, _ in marks)]
    for c in unpaired:
        seq = "".join(r["letter"] or "?" for r in chains[c]["residues"])
        lines += [f"Chain {c} has no Watson-Crick pairs: 5'-{seq}-3'", ""]
    return "\n".join(lines)
