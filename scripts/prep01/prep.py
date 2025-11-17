import os
import subprocess

def main():
    # Directory where runscript.py is located
    here = os.path.dirname(os.path.abspath(__file__))
    sh_path = os.path.join(here, "prep_00_all.sh")
    subprocess.call([
        sh_path,
        "/data/leuven/359/vsc35916/thesis_main/data/8RM6/",
        "8RM6_mod.pdb",
        "--with-dna",
        "--range_a", "2,-1",
        "--range_b", "2,-1"
    ])