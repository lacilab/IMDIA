Use 1PGB (1 domain of streptococcal protein G) at pH 3.0, 4.0 and 7.4 to test acidic side-chain assignment and the terminal guard.

| Test | Expected result | Ionic changes
| **pH 3.0, no overrides** | Stops before LEaP: terminal GLU56 would become unvalidated GLH | # Fail |
| **pH 3.0, GLU56 override** | Continues; internal acidic residues follow their predicted pKa | # Na 14, Cl 17 |
| **pH 4.0, no overrides** | GLU56 stays GLU because its observed pKa is 3.13; internal residues show mixed states | # Na 14, Cl 14 |
| **pH 7.4, no overrides** | Based on your observed pKa values, all ASP/GLU side chains remain negatively charged | # Na 18, Cl 14 |


## Test protonation after --prep run:

rg ',(ASP|GLU),' usecase_1PGB_pH/amber/*.protonation.csv

rg '; residue\s+[0-9]+\s+(ASH|GLH|ASP|GLU)\b' \
  usecase_1PGB_pH/gromacs_common/toppar/molecule_*.itp