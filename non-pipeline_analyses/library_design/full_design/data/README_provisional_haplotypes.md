# Provisional haplotypes

`provisional_h3n2.txt` and `provisional_h1n1pdm.txt` list haplotypes added to the library that have
**no representative sequence** in the VCM strain-selection TSVs — either because no isolate carrying
the haplotype has been deposited, or because the haplotype is defined from outside that input. Their
HA1 sequences are reconstructed from the subclade backbone rather than observed, which is why they
are marked provisional.

These are distinct from `manually_fixed_recent_*_haplotypes.tsv`, which re-adds haplotypes that
*do* have a representative strain whose sequence contained ambiguous (`X`) residues.

## Adding a haplotype

**These two `.txt` files are the only thing you edit.** One `derived_haplotype` per line; blank
lines and `#` comments are ignored. They are named by `provisional_haplotypes` in
[../config.yaml](../config.yaml).

The `build_provisional_haplotypes` rule expands each list into
`results/provisional_haplotypes/{subtype}_provisional_haplotypes.tsv`, which joins the files listed
under `additional_haplotypes` as input to `curate_library_sequences`. Snakemake runs it
automatically, so adding a line here and re-running the pipeline is the whole workflow — there is no
separate generation step to remember, and no committed TSV that can fall out of sync with the list.

Every row of an `additional_haplotypes` file is forced to `selected_haplotype = True` by
`scripts/curate_library_sequences.py`, so only list haplotypes you actually want in the library.

The generated TSV is pipeline output, not input: don't hand-edit it. The generation step exists
because `curate_library_sequences.py` hard-requires `representative_strain`,
`representative_strain_ha1_sequence`, and `representative_ha2_sequence` in the file it reads, so a
bare list cannot be fed to `additional_haplotypes` directly.

To run the build by hand, e.g. to inspect the output before a full pipeline run:

```
python3 scripts/build_provisional_haplotypes.py \
    --vcm-tsv data/2027-Feb-VCM-seqneut-library-strain-selection-H3N2.tsv \
    --haplotype-list data/provisional_h3n2.txt \
    --output results/provisional_haplotypes/H3N2_provisional_haplotypes.tsv
```

### Haplotypes that already have an observed sequence

The build **fails** on a haplotype that already has a representative sequence in the VCM TSV, and
tells you to use a selection file instead. Such a haplotype is not provisional: listing it here
would discard its real strain name, count, and collection date, give it a GenBank accession that is
only a near match, and select it twice over — with `construct_order`'s sequence dedup dropping one
copy by row order rather than by which record is better.

### Checking the reconstruction

`--self-test` applies the reconstruction to every haplotype that *does* have an observed sequence
and reports how often it reproduces the observed HA1:

```
python3 scripts/build_provisional_haplotypes.py \
    --vcm-tsv data/2027-Feb-VCM-seqneut-library-strain-selection-H3N2.tsv --self-test
```

As of the 2027-Feb VCM input: 2465/2592 exact for H3N2 and 1344/1380 for H1N1, with **zero** genuine
disagreements — every mismatch is at a position where the observed sequence has an ambiguous `X`,
meaning the reconstruction is more complete than the deposit rather than wrong.

## How "provisional" is recorded

Three independent markers, so the status survives every downstream join:

1. **`provisional` column** — set to `TRUE` on every row here. Unlisted elsewhere, so it arrives as
   empty for all non-provisional haplotypes after the concat in `curate_library_sequences.py`.
2. **`representative_strain` suffix** — named `{derived_haplotype}_provisional`, e.g.
   `K:N6D,K135E_provisional`. This is the marker that reaches furthest: `representative_strain`
   becomes `strain` in the plasmid log, the GenBank record name, and the GenBank filename, none of
   which carry the `provisional` column.
3. **`selection_file`** — set automatically to this file's basename by
   `curate_library_sequences.py`, and carried through `library_strains.tsv` into
   `construct_order`'s `final_library.csv`. No action needed.

## Columns

`derived_haplotype`, `representative_strain`, and `latest_sequence` are hard requirements of
`scripts/aggregate_library_strains.py`. The rest exist because
`nextstrain_tree.*.color_by_metadata` in [../config.yaml](../config.yaml) requires the columns to be
present (`scripts/prepare_nextstrain_tree.py` checks existence, not values), so they may be left
empty.

| Column | Fill | Notes |
|---|---|---|
| `derived_haplotype` | required | `{subclade}:{mut},{mut}` — e.g. `K:N6D,K135E`. 1-indexed HA1 numbering. Mutation **order is significant** if you also reference this haplotype in `construct_order`'s `extra_barcode_haplotypes`. |
| `count` | `0` | No observed sequences. Note `count` uses a log color scale in the tree config; 0 may render oddly. |
| `emerging_haplotype` | optional | Leave empty unless meaningful. |
| `n_countries`, `n_regions` | `0` | |
| `latest_sequence` | required | No observed sequence, so there is no real date. Must be non-null for `aggregate_library_strains.py`; it is normalized to `""` there and used as `date` in the tree metadata. |
| `representative_strain` | required | `{derived_haplotype}_provisional`. Must be unique across both subtypes. |
| `representative_strain_ha1_sequence` | required | Reconstructed: subclade backbone HA1 with the haplotype's mutations applied. 329 aa (H3N2) / 327 aa (H1N1). |
| `representative_ha2_sequence` | required | Subclade backbone HA2, unmodified — see below. 221 aa (H3N2) / 222 aa (H1N1). |
| `days_since_latest_sequence` | empty | |
| `relative_fitness_derived_haplotype` | empty | |
| `provisional` | `TRUE` | |

Sequences must contain only the 20 standard amino acids and match the length of every other
haplotype, or `curate_library_sequences.py` fails.

## Subclade backbones

Taken from the bare-clade rows (no `:` in `derived_haplotype`) of the VCM TSVs:

| Subtype | Subclade | Backbone strain | Sequence count |
|---|---|---|---|
| H3N2 | `K` | A/Romania/439/2026 | 4156 |
| H3N2 | `J.2.4.2` | A/Zacapa/FLU-227/2026 | 480 |
| H3N2 | `J.2.3` | A/MatoGrossoDoSul/LEIAL8184/2026 | 16 |
| H3N2 | `J.2.4` | A/Madagascar/1149/2026 | 2 |
| H3N2 | `J.2.5` | A/Cayman/26InfM314/2026 | 1 |
| H1N1 | `D.3.1.1` | A/Victoria/199/2026 | 1641 |
| H1N1 | `D.3.1` | A/SouthAfrica/PATH-CERI-C076784/2026 | 268 |

A haplotype in any other subclade has no backbone row and cannot be reconstructed this way.

### Why HA1 reconstruction is trustworthy

Applying backbone + mutations to every haplotype in the VCM TSVs that *does* have a representative
sequence reproduces the observed HA1 exactly for 2465/2592 (H3N2) and 1344/1380 (H1N1). Every
mismatch is an `X` in the observed sequence — the reconstruction is more complete than the deposited
sequence at those sites. There were **zero** cases where the backbone residue disagreed with a
mutation's stated wild-type, and zero mutations at sites beyond HA1.

### Why HA2 is the backbone's, unmodified

`derived_haplotype` encodes HA1 mutations only, so it does not determine HA2 — 583 H3N2 and 284 H1N1
haplotypes carry an HA2 differing from their subclade backbone's. There is no way to infer the
"correct" HA2 for a haplotype that has never been observed, so the backbone HA2 is used as the
deliberate default. It is the modal HA2 within its clade by a wide margin — 90.1% of K-clade and
92.9% of D.3.1.1 sequences by count, with the next most common variants only 1–2 residues away.

**This is a design decision, not a derivation.** If a provisional haplotype is expected to sit on a
lineage carrying a non-backbone HA2, override `representative_ha2_sequence` on that row and note why
here.

## Downstream consequences of editing these files

* **Re-freeze the GenBank match.** Provisional sequences are not in GenBank, so they are absent from
  `data/frozen_genbank_match/{subtype}.csv` and `aggregate_library_strains.py` hard-errors on them.
  Set `freeze_genbank_match: false`, run, copy the regenerated
  `results/genbank_match/{subtype}/match_prot_to_genbank_nt.csv` back in as `{subtype}.csv`, then set
  it true again. `match_prot_to_genbank_nt` resolves a never-observed protein to the nearest
  full-coverage accession and records the added codon changes in `aa_muts_added_to_accession`, which
  surfaces as `accession_w_aa_muts_added` — so the library's `genbank_accession` for a provisional
  strain refers to a *related* accession, not an exact match.
* **Construct-order assignments shift.** `construct_order` assigns shortnames, barcodes, and plate
  wells positionally from curated-TSV row order, so adding rows renumbers constructs. Only
  regenerate an order that has not yet been submitted for synthesis.
* **Check `extra_barcode_haplotypes`.** That list in `../../construct_order/config.yaml` is a hard
  error if an entry matches no strain in the order.
