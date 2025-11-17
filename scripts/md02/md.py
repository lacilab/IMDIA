import os
import subprocess

def main():
    # Directory where runscript.py is located
    here = os.path.dirname(os.path.abspath(__file__))
    sh_path = os.path.join(here, "md.sh")
    subprocess.call([
        sh_path,
        "/data/leuven/359/vsc35916/thesis_main/data/8RM6/gromacs_common/",
        "1",
        "--with-dna",
        "--ntmpi", "2",
        "--ntomp", "16",
        "--ngpu", "1",
        "--gpulist", "01",
    ])