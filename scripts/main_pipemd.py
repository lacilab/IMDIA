import prep01.prep
import md02.md_wrap

def main(datadir, filename, gmxdir, repn,
        with_dna = False, with_phosphate = False, range_a = None, range_b = None, ffs = None,
        ntmpi = None, ntomp = None, ngpu = None, gpulist = None, 
        em = None, nvt = None, npt = None, md = None):

        prep01.prep.main(datadir, filename, with_dna, with_phosphate, range_a, range_b, ffs)
        md02.md_wrap.main(gmxdir, repn, with_dna, with_phosphate, ntmpi, ntomp, ngpu, gpulist, em, nvt, npt, md)

main(datadir = '../data/8RM6/',
        filename = '8RM6_mod.pdb',
        with_dna = True,
        with_phosphate = False,
        range_a = '2,-1',
        range_b = '2,-1',
        gmxdir = '../data/8RM6/gromacs_common',
        repn = "1",
        ntmpi = "2",
        ntomp = "16",
        ngpu = "1",
        gpulist = '01',
        md = 'md02/md_10ns.mdp'
    )