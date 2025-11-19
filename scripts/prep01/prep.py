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

    print ("Running prep script...")
    print( f"sh_path: {sh_path}",
           f"datadir: {datadir}",
           f"filename: {filename}",
           f"dna_flag: {dna_flag}",
           f"range_a_flag: {range_a_flag}", f"rang_a_val: {rang_a_val}",
           f"range_b_flag: {range_b_flag}", f"range_b_val: {range_b_val}",
           f"phosphate_flag: {phosphate_flag}",
           f"ffs_flag: {ffs_flag}", f"ffs_val: {ffs_val}",
         )

    subprocess.run([
        sh_path,
        datadir,
        filename,
        dna_flag,
        range_a_flag, rang_a_val,
        range_b_flag, range_b_val,
        phosphate_flag,
        ffs_flag, ffs_val
    ])