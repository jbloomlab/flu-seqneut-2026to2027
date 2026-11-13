# Full library design

Snakemake pipeline for selecting influenza HA strains for the 2026-2027 neutralization library.

**Status: scaffold.** The directory layout, environment, and cluster profile are in place. The
`Snakefile`, the `scripts/` implementations, and the input data under `data/` still need to be
brought over and adapted from the `full_design` subdirectory of `flu-seqneut-2026`, which is the
reference implementation for every step described below.

## Layout

* [config.yaml](./config.yaml) — all input paths and options (templated, with TODOs to fill in)
* [environment.yml](./environment.yml) — `conda` environment for the pipeline
* [Snakefile](./Snakefile) — workflow definition (to be added)
* [data/](./data/) — input haplotype TSVs, selection files, outgroups, frozen GenBank matches
* [scripts/](./scripts/) — the scripts each rule runs (to be added)
* [results/](./results/) — all pipeline output
* [profiles/config.yaml](./profiles/config.yaml) — SLURM profile for the Hutch cluster

## Expected workflow

1. **Validate** recent haplotypes, flagging any with ambiguous or missing amino acids
   → `results/recent_haplotype_validation/`
2. **Curate** the library, reporting Hamming and epitope-site distances between selected and
   non-selected haplotypes → `results/curated_library/`
3. **Annotate sites** in both sequential and HA1/HA2 numbering → `results/site_annotations/`
4. **Build trees** as Auspice JSONs for interactive inspection at
   [auspice.us](https://auspice.us/) → `results/tree/`
5. **Aggregate** the final strain set, with GenBank accessions and nucleotide sequences
   → `results/aggregated_library_strains/library_strains.tsv`

Step 5 produces the input consumed by [../construct_order/](../construct_order/).

Selection is iterative: inspect the curated distances and the trees, adjust which haplotypes are
selected in `config.yaml`, and re-run until the library composition looks right.

## Submodules

Two submodules are needed by this pipeline and are added at the top level of the repository, not
here:

* [match_prot_to_genbank_nt](https://github.com/jbloom/match_prot_to_genbank_nt/) — matches protein
  sequences to the closest GenBank nucleotide sequence
* [nextstrain-prot-titers-tree](https://github.com/jbloomlab/nextstrain-prot-titers-tree) — builds
  Nextstrain trees from protein alignments

## Running

Build and activate the environment, then run:

```
conda env create -f environment.yml
conda activate library_design
snakemake --use-conda -j <ncpus>
```

Or on the Hutch cluster:

```
sbatch run_Hutch_cluster.bash
```
