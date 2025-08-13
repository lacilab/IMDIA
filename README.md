This repository contains the code developed for the thesis *Implementation of a molecular dynamics pipeline and integration with omics data to predict the impact of missense mutations for STAT1 gain-of-function* presented in fulfillment of the requirements for the degree of Master of Science: Bioinformatics at KU Leuven. 

**mut_choice** contains the results of the preliminary FoldX analysis that was used to support the choice of mutations for this project.

**motif_analysis** contains the motif analysis performed with HOMER based on the differentially expressed genes identified through RNAseq:

* **4h** - 4h IFNg estimulation on U3A cells (Giovannozi *et al.*, 2021, https://doi.org/10.3389/fimmu.2021.632997)
* **1h** - 1h IFNg estimulation on THP-1 cells

**scripts** contains the semi-automated pipeline to perform molecular dynamics simulations on STAT1 bound to DNA, and analyze the results. The pipeline consists of four main steps:

1. System preparation
2. Molecular dynamics simulation
3. Post-processing
4. Analysis

The scripts are prepared to work on the following directory system:

```text
/
├── data/
│   ├── {structure}
│   │   ├── {mutation}
│   │   │   └── {prep}.pdb
├── scripts/
│   ├── 01_prep
│   ├── 02_md
│   ├── 03_post
│   └── 04_analysis
```

From this initial directory system, the scripts generate the necessary folders to achieve this structure:

```text
/
├── data/
│   ├── {structure}
│   │   ├── {mutation}
│   │   │   ├── amber
│   │   │   ├── gromacs_common
│   │   │   ├── rep{i}
│   │   │   │   ├── analysis
│   │   │   │   │   ├── 000_raw/
│   │   │   │   │   ├── 001_qc/
│   │   │   │   │   ├── 002_hbonds/
│   │   │   │   │   ├── 003_pca/
│   │   │   │   │   ├── 004_contacts/
│   │   │   │   │   ├── 005_clustering/
│   │   │   │   │   └── 006_gbsa/
│   │   │   │   └── gromacs
├── scripts
│   ├── 01_prep
│   ├── 02_md
│   ├── 03_post
│   └── 04_analysis
```

In `scripts\01_prep`, run `bash prep_00_all.sh data/{structure}/{mutation} {prep}` to execute these steps:

1. Manual structure revision and preparation
2. Force field selection and solvation
3. Conversion to GROMACS
4. Topology extraction
5. Neutralization
6. Distance restrains

In `scripts\02_md`, run `bash md.sh  data/{structure}/{mutation}/gromacs_common replicate_n` to run the MD simulation.

In `scripts\03_post`, run `bash wrap.sh  data/{structure}/{mutation}/rep{i}/gromacs` to process the MD trajectory.

In `scripts\04_analysis`, run `snakemake` (or `snakemake --cores=12`, or the available number of cores) to execute the snakemake analysis. Once the analysis is completed, run `snakemake --report` to generate the report.
