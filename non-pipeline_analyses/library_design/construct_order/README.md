# Construct order

Designs the HA insert sequences submitted to [GenScript](https://www.genscript.com) for gene
fragment synthesis and cloning, for the 2026-2027 library.

**Status: scaffold.** The directory layout, environment, and cluster profile are in place. The
`workflow/Snakefile`, the scripts under `workflow/scripts/`, and the reference data under `data/`
still need to be brought over and adapted from the `construct_order` subdirectory of
`flu-seqneut-2026`, which is the reference implementation.

Inserts are designed to sit between the BsmBI-v2 cut sites in Bloom lab vector 5764, beginning after
the 19th codon (H3N2) or 20th codon (H1N1) of the WSN upstream signal peptide and continuing through
the end of the HA coding region, followed by a double stop codon and a 16-nucleotide barcode.

**Note:** the scripts are hardcoded for vector 5764 and make assumptions that will not generalize to
other backbones. A GenBank file of the vector belongs in `data/`.

## Layout

* [config.yaml](./config.yaml) — insert-design parameters and the list of orders (templated)
* [environment.yml](./environment.yml) — `conda` environment for the pipeline
* `workflow/Snakefile` — workflow definition (to be added)
* `workflow/scripts/` — the scripts each rule runs (to be added)
* [data/](./data/) — vector GenBank file, construct log, vaccine annotations, prior barcodes
* [results/](./results/) — all pipeline output
* [logs/](./logs/) — per-rule logs
* [profiles/config.yaml](./profiles/config.yaml) — SLURM profile for the Hutch cluster

## Expected outputs

* `results/orders/{date}_{order}_order.csv` — minimal CSV of shortname and insert sequence for
  upload to the ordering portal; includes a platemap for constructs ordered in plate format
* `results/plasmid_log/{date}_{order}.csv` — everything needed to record the new plasmids in the
  seqneut plasmid log
* `results/genbank/` — one GenBank file per plasmid
* `results/excluded_from_order/{date}_{order}_excluded.csv` — entries dropped because the same HA
  ectodomain sequence already exists in a prior plasmid or elsewhere in the order set
* `results/final_library/final_library.csv` — the final sequence set and metadata for the library

## Configuration

Orders are added under the `orders` key in [config.yaml](./config.yaml); the commented-out block
there lists the required fields. Orders are processed in order, each aware of the barcodes and
sequences consumed by the ones before it, so append new orders at the end.

The input to each order is the `library_strains.tsv` produced by the corresponding design
subdirectory, e.g. [../full_design/](../full_design/).

## Insert design

H3 constructs: WSN upstream signal peptide, HA ectodomain from a circulating strain, endodomain from
a recent H3 consensus sequence, and a WSN C-terminal domain.

H1 constructs: WSN upstream signal peptide, HA ectodomain from a circulating strain, with both the
endodomain and C-terminal domain from WSN.

## Barcode design

Barcodes avoid those used by prior libraries (identified via `past_barcodes_to_avoid` in
[config.yaml](./config.yaml)) and avoid barcodes starting with `GG`, which sequence poorly.

## Running

```
conda env create -f environment.yml
conda activate constructs
```

Then on the Hutch cluster:

```
sbatch run_snakemake.bash
```
