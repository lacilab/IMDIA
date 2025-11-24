import os
import subprocess

def main(datadir: str, filename: str,
    with_dna = False, with_phosphate = False,
    range_a = None, range_b = None,
    ffs = None):
    here = os.path.dirname(os.path.abspath(__file__))
    sh_path = os.path.join(here, "prep_00_all.sh")

    dna_flag = "--with-dna" if with_dna else ""
    phosphate_flag = "--with-phosphate" if with_phosphate else ""
    range_a_flag = "--range_a" if range_a is not None else ""
    rang_a_val = range_a if range_a is not None else ""
    range_b_flag = "--range_b" if range_b is not None else ""
    range_b_val = range_b if range_b is not None else ""
    ffs_flag = "--ffs" if ffs is not None else ""
    ffs_val = ffs if ffs is not None else ""

    cmd = [
        sh_path,
        datadir,
        filename,
        dna_flag,
        range_a_flag, rang_a_val,
        range_b_flag, range_b_val,
        phosphate_flag,
        ffs_flag, ffs_val
    ]

    print("Running command for prep script:")
    # remove empty strings from cmd before joining
    cmd = [c for c in cmd if c != ""]
    print(" ".join(cmd))
    subprocess.call(cmd)