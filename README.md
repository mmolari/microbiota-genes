# microbiota-associated gene distribution in *E. coli*

A simple analysis to survery the distribution of a panel of genes of interest across the *E. coli* species using the Horesh et al. 2021 pangenome dataset ([10.1099/mgen.0.000499](https://doi.org/10.1099/mgen.0.000499)). These focal genes were shown to increase *E. coli* fitness specifically **in the presence of gut microbiota**.

You can find in:

- [`notes/methods.md`](notes/methods.md): an overview of the pipeline and the methods used.
- [`notes/results.md`](notes/results.md): resulting figures and their description.

## Setup

To reproduce the analysis, the only host-side requirements are [conda](https://docs.conda.io/) and [snakemake](https://snakemake.readthedocs.io/) (tested on snakemake 9.20.0). All other dependencies (BLAST, pandas, biopython, etc.) are pulled in automatically by `--use-conda`. The pipeline downloads the Horesh dataset from figshare on first run.

## Run

```bash
snakemake -c4 --use-conda
```

The main results produced are the figures (in `results/figs/`) and the focal-gene presence/absence matrix (in `results/presence_absence.csv`)

## Citation

The dataset analyzed here is from Horesh et al. 2021:

> Horesh, Gal, et al. "A comprehensive and high-quality collection of *Escherichia coli* genomes and their genes." *Microbial Genomics* 7.2 (2021): 000499. [doi:10.1099/mgen.0.000499](https://doi.org/10.1099/mgen.0.000499)