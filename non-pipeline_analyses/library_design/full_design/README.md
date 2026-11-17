# Full library design

Snakemake pipeline for selecting influenza HA strains for the 2026-2027 neutralization library.

Adapted from the `full_design` subdirectory of `flu-seqneut-2026`, which remains the reference
implementation for every step described below.

## Layout

* [config.yaml](./config.yaml) — all input paths and options
* [environment.yml](./environment.yml) — `conda` environment for the pipeline
* [Snakefile](./Snakefile) — workflow definition
* [data/](./data/) — input haplotype TSVs, selection files, outgroups, frozen GenBank matches
* [scripts/](./scripts/) — the scripts each rule runs
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

## Choosing haplotypes

Haplotypes enter the library three ways, all configured in [config.yaml](./config.yaml):

* `select_recent_haplotype_files` — the normal route. Plain lists (`data/{person}_{subtype}.txt`) of
  `derived_haplotype` names to pick from the recent-haplotypes TSV. A name not found there is
  **warned about and ignored**, not an error, so check the validation report after editing.
* `additional_haplotypes` — haplotypes needing a hand-corrected sequence, currently those whose
  representative strain has ambiguous (`X`) residues. Paired with
  `override_select_recent_haplotypes` to drop the uncorrected version.
* `provisional_haplotypes` — haplotypes with **no representative sequence at all**, because no
  isolate carrying them has been deposited. See below.

### Provisional haplotypes

Listed in `data/provisional_{subtype}.txt`, one `derived_haplotype` per line. The
`build_provisional_haplotypes` rule reconstructs a sequence for each by applying its HA1 mutations
to the subclade backbone (the bare-clade row of the recent-haplotypes TSV), writing a TSV that
joins `additional_haplotypes` as input to curation. Edit only the `.txt`; the TSV is pipeline
output.

HA2 is taken from the backbone unmodified, because `derived_haplotype` encodes HA1 mutations only
and so cannot constrain HA2. It is the modal HA2 within its clade by a wide margin (90.1% of K,
92.9% of D.3.1.1 by sequence count), but it is a **design decision, not a derivation**: 583 H3N2 and
284 H1N1 observed haplotypes carry an HA2 differing from their backbone's. A haplotype expected on a
lineage with a different HA2 cannot be expressed through the `.txt` list — give it a hand-built row
in an `additional_haplotypes` TSV instead.

Each is marked three ways so the status survives into `construct_order`'s `final_library.csv`: a
`provisional` column, a `_provisional` strain-name suffix, and `selection_file`.

To check the reconstruction against every haplotype that *does* have an observed sequence:

```
python3 scripts/build_provisional_haplotypes.py \
    --vcm-tsv data/2027-Feb-VCM-seqneut-library-strain-selection-H3N2.tsv --self-test
```

As of the 2027-Feb input: 2465/2592 exact for H3N2 and 1344/1380 for H1N1, with zero genuine
disagreements — every mismatch is where the observed sequence has an ambiguous `X`, so the
reconstruction is more complete than the deposit rather than wrong.

Adding a provisional haplotype requires re-freezing the GenBank match, since the reconstructed
sequences are not in GenBank: set `freeze_genbank_match: false`, run, copy each regenerated
`results/genbank_match/{subtype}/match_prot_to_genbank_nt.csv` to
`data/frozen_genbank_match/{subtype}.csv`, then set it back to `true`. Without this,
`aggregate_library_strains` fails on the missing strains.

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
