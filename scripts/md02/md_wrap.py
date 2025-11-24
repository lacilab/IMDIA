import cmd
import os
import subprocess

def main(gmxdir: str, repn: int, 
    with_dna = False, with_phosphate = False,
    ntmpi = None, ntomp = None, ngpu = None, gpulist = None, 
    em = None, nvt = None, npt = None, md = None, 
    ):
    here = os.path.dirname(os.path.abspath(__file__))
    sh_path = os.path.join(here, "md.py")

    repn_flag = "--repn" if repn is not None else ""
    repn_val = repn if repn is not None else ""
    dna_flag = "--with-dna" if with_dna else ""
    phosphate_flag = "--with-phosphate" if with_phosphate else ""
    ntmpi_flag = "--ntmpi" if ntmpi is not None else ""
    ntmpi_val = ntmpi if ntmpi is not None else ""
    ntomp_flag = "--ntomp" if ntomp is not None else ""
    ntomp_val = ntomp if ntomp is not None else ""
    ngpu_flag = "--ngpu" if ngpu is not None else ""
    ngpu_val = ngpu if ngpu is not None else ""
    gpulist_flag = "--gpulist" if gpulist is not None else ""
    gpulist_val = gpulist if gpulist is not None else ""
    em_flag = "--em" if em is not None else ""
    em_val = em if em is not None else ""       
    nvt_flag = "--nvt" if nvt is not None else ""
    nvt_val = nvt if nvt is not None else ""       
    npt_flag = "--npt" if npt is not None else ""
    npt_val = npt if npt is not None else ""       
    md_flag = "--md" if md is not None else ""
    md_val = md if md is not None else ""
    
    cmd = [
        sh_path,
        "--gmxdir", gmxdir,
        repn_flag, repn_val,
        dna_flag,
        phosphate_flag,
        ntmpi_flag, ntmpi_val,
        ntomp_flag, ntomp_val,
        ngpu_flag, ngpu_val,
        gpulist_flag, gpulist_val,
        em_flag, em_val,
        nvt_flag, nvt_val,
        npt_flag, npt_val,
        md_flag, md_val
    ]

    print("Running command for md script:")
    # remove empty strings from cmd before joining
    cmd = [c for c in cmd if c != ""]
    print(" ".join(cmd))
    subprocess.call(cmd)