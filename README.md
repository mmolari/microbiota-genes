# microbiota-genes

Snakemake pipeline that surveys a panel of *focal genes* across the *E. coli* species using the Horesh et al. 2021 pangenome dataset ([10.1099/mgen.0.000499](https://doi.org/10.1099/mgen.0.000499)). Focal genes are a curated set of genes that were shown to confer a fitness advantage to *E. coli* isolates colonizing the mouse gut, **only in the presence of a complex microbiota**.

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

The dataset used here is from Horesh et al. 2021:

> Horesh, Gal, et al. "A comprehensive and high-quality collection of *Escherichia coli* genomes and their genes." *Microbial Genomics* 7.2 (2021): 000499. [doi:10.1099/mgen.0.000499](https://doi.org/10.1099/mgen.0.000499)