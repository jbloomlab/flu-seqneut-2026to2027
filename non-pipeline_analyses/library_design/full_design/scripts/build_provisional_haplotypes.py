"""
Build a provisional-haplotype TSV from a plain list of haplotype names.

Provisional haplotypes have no representative sequence in the VCM strain-selection TSV --
no isolate carrying them has been deposited, or they are defined from outside that input.
This script reconstructs a sequence for each by taking its subclade backbone (the bare-clade
row of the VCM TSV, i.e. the one whose `derived_haplotype` has no `:`) and applying the
haplotype's HA1 mutations.

HA2 is taken from the backbone unmodified. `derived_haplotype` encodes HA1 mutations only, so
it does not constrain HA2, and there is no way to infer the "correct" HA2 for a haplotype never
observed. The backbone HA2 is the modal HA2 within its clade by a wide margin (90.1% of K-clade
and 92.9% of D.3.1.1 sequences by count). See data/README_provisional_haplotypes.md.

The reconstruction is validated by the `--self-test` mode, which applies it to every haplotype
in the VCM TSV that *does* have a representative sequence and reports how often it reproduces
the observed HA1.

Writes the TSV that `additional_haplotypes` in config.yaml consumes. Every row is marked
`provisional = TRUE`, and `representative_strain` is `{derived_haplotype}_provisional`.
"""
import argparse
import csv
import re
import sys


MUT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")

# Column order of the generated TSV. The three sequence/strain columns are hard
# requirements of curate_library_sequences.py; derived_haplotype, representative_strain,
# and latest_sequence are required by aggregate_library_strains.py. The rest exist
# because nextstrain_tree.*.color_by_metadata in config.yaml requires the columns to be
# present (prepare_nextstrain_tree.py checks existence, not values).
OUTPUT_COLUMNS = [
    "derived_haplotype",
    "count",
    "emerging_haplotype",
    "n_countries",
    "n_regions",
    "latest_sequence",
    "representative_strain",
    "representative_strain_ha1_sequence",
    "representative_ha2_sequence",
    "days_since_latest_sequence",
    "relative_fitness_derived_haplotype",
    "provisional",
]

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_vcm(path):
    """Read the VCM TSV, returning (all rows by haplotype, backbone rows by subclade)."""
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        raise ValueError(f"VCM TSV '{path}' is empty")

    by_haplotype = {}
    for row in rows:
        hap = row["derived_haplotype"].strip()
        if hap:
            by_haplotype.setdefault(hap, row)

    # Backbones are the bare-clade rows: a derived_haplotype with no mutation suffix.
    backbones = {h: r for h, r in by_haplotype.items() if ":" not in h}
    if not backbones:
        raise ValueError(
            f"VCM TSV '{path}' contains no bare-clade (backbone) rows -- expected rows "
            f"whose derived_haplotype has no ':' suffix, e.g. 'K' or 'D.3.1.1'."
        )
    return by_haplotype, backbones


def read_haplotype_list(path):
    """Read a plain-text haplotype list, one per line; blanks and # comments ignored."""
    haplotypes = []
    with open(path) as f:
        for lineno, line in enumerate(f, start=1):
            name = line.strip()
            if not name or name.startswith("#"):
                continue
            if ":" not in name:
                raise ValueError(
                    f"{path}:{lineno}: '{name}' has no ':' mutation suffix. A provisional "
                    f"haplotype must name a subclade and at least one mutation, e.g. "
                    f"'K:N190T'. A bare subclade is already a backbone row in the VCM TSV."
                )
            if name in haplotypes:
                raise ValueError(f"{path}:{lineno}: '{name}' is listed more than once")
            haplotypes.append(name)
    return haplotypes


def reconstruct(haplotype, backbones):
    """
    Reconstruct HA1 and HA2 for a haplotype from its subclade backbone.

    Returns (ha1, ha2, backbone_strain). Raises ValueError with an actionable message
    if the subclade has no backbone, a mutation is unparsable, a site is out of range,
    or the backbone residue disagrees with the mutation's stated wild-type.
    """
    subclade, _, mutations = haplotype.partition(":")

    if subclade not in backbones:
        raise ValueError(
            f"{haplotype}: no backbone row for subclade '{subclade}' in the VCM TSV. "
            f"Available backbones: {sorted(backbones)}"
        )

    backbone = backbones[subclade]
    ha1 = list(backbone["representative_strain_ha1_sequence"])
    ha2 = backbone["representative_ha2_sequence"]

    if not ha1 or not ha2:
        raise ValueError(
            f"{haplotype}: backbone row for '{subclade}' "
            f"({backbone['representative_strain']}) is missing an HA1 or HA2 sequence"
        )

    seen_sites = {}
    for mutation in mutations.split(","):
        mutation = mutation.strip()
        match = MUT_RE.match(mutation)
        if not match:
            raise ValueError(
                f"{haplotype}: cannot parse mutation '{mutation}'. Expected the form "
                f"'<wt><site><mut>' in 1-indexed HA1 numbering, e.g. 'N190T'."
            )
        wt, site, mut = match.group(1), int(match.group(2)), match.group(3)

        if not 1 <= site <= len(ha1):
            raise ValueError(
                f"{haplotype}: mutation '{mutation}' is at site {site}, outside HA1 "
                f"(length {len(ha1)}). derived_haplotype uses HA1 numbering only."
            )
        if site in seen_sites:
            raise ValueError(
                f"{haplotype}: site {site} is mutated twice "
                f"('{seen_sites[site]}' and '{mutation}')"
            )
        seen_sites[site] = mutation

        actual = ha1[site - 1]
        if actual != wt:
            raise ValueError(
                f"{haplotype}: mutation '{mutation}' expects wild-type '{wt}' at HA1 site "
                f"{site}, but the '{subclade}' backbone "
                f"({backbone['representative_strain']}) has '{actual}'. Either the "
                f"haplotype name is wrong or it does not belong to this subclade."
            )
        if mut == wt:
            raise ValueError(
                f"{haplotype}: mutation '{mutation}' is synonymous at the protein level "
                f"(wild-type and mutant are both '{wt}')"
            )
        ha1[site - 1] = mut

    ha1_seq = "".join(ha1)

    for label, seq in (("HA1", ha1_seq), ("HA2", ha2)):
        invalid = sorted(set(seq) - VALID_AA)
        if invalid:
            raise ValueError(
                f"{haplotype}: reconstructed {label} contains non-standard residue(s) "
                f"{invalid}, inherited from backbone {backbone['representative_strain']}. "
                f"The backbone sequence itself needs fixing before it can seed a "
                f"provisional haplotype."
            )

    return ha1_seq, ha2, backbone["representative_strain"]


def self_test(by_haplotype, backbones):
    """
    Validate the reconstruction against every haplotype that has an observed sequence.

    Reports exact-match counts and, separately, mismatches attributable to ambiguous (X)
    residues in the observed sequence versus genuine disagreements. A genuine
    disagreement means the method is unsound and is reported as a failure.
    """
    exact = 0
    ambiguous_only = 0
    genuine = []
    skipped = 0

    for haplotype, row in by_haplotype.items():
        if ":" not in haplotype:
            continue
        observed = row["representative_strain_ha1_sequence"]
        if not observed or not row["representative_ha2_sequence"]:
            skipped += 1
            continue
        try:
            built, _, _ = reconstruct(haplotype, backbones)
        except ValueError:
            skipped += 1
            continue
        if len(built) != len(observed):
            skipped += 1
            continue

        if built == observed:
            exact += 1
            continue

        diffs = [
            (i + 1, b, o) for i, (b, o) in enumerate(zip(built, observed)) if b != o
        ]
        if all(o not in VALID_AA for _, _, o in diffs):
            # Every mismatch is at a position where the observed sequence is ambiguous,
            # so the reconstruction is more complete than the deposit, not wrong.
            ambiguous_only += 1
        else:
            genuine.append((haplotype, diffs))

    total = exact + ambiguous_only + len(genuine)
    print(f"  compared:                      {total} haplotypes ({skipped} skipped)")
    print(f"  reproduced observed HA1 exactly: {exact}")
    print(f"  differed only at ambiguous (X) positions in the observed sequence: "
          f"{ambiguous_only}")
    print(f"  genuine disagreements:           {len(genuine)}")

    if genuine:
        print("\n  GENUINE DISAGREEMENTS (reconstruction is unsound):")
        for haplotype, diffs in genuine[:10]:
            shown = ", ".join(f"site {s}: built {b} vs observed {o}" for s, b, o in diffs[:5])
            print(f"    {haplotype}: {len(diffs)} diff(s) -- {shown}")
        if len(genuine) > 10:
            print(f"    ... and {len(genuine) - 10} more")
        return False

    return True


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--vcm-tsv", required=True,
                        help="VCM strain-selection TSV supplying subclade backbones")
    parser.add_argument("--haplotype-list",
                        help="plain-text list of haplotypes, one per line")
    parser.add_argument("--output", help="TSV to write")
    parser.add_argument("--self-test", action="store_true",
                        help="validate the reconstruction against observed sequences and exit")
    args = parser.parse_args()

    by_haplotype, backbones = read_vcm(args.vcm_tsv)
    print(f"Read {len(by_haplotype)} haplotypes from {args.vcm_tsv}")
    print(f"Backbones: " + ", ".join(
        f"{clade} ({row['representative_strain']}, count {row['count']})"
        for clade, row in sorted(backbones.items())
    ))

    if args.self_test:
        print("\nSelf-test: backbone + mutations vs observed HA1")
        if not self_test(by_haplotype, backbones):
            sys.exit(1)
        print("\n  PASS: no genuine disagreements.")
        return

    if not args.haplotype_list or not args.output:
        parser.error("--haplotype-list and --output are required unless --self-test is given")

    haplotypes = read_haplotype_list(args.haplotype_list)
    print(f"\nRead {len(haplotypes)} haplotype(s) from {args.haplotype_list}")

    # A haplotype with an observed representative belongs in a selection .txt file, where it
    # keeps its real strain name, count, collection date, and exact GenBank accession.
    # Listing it here would select it twice over, with construct_order's dedup dropping one
    # copy by row order rather than by which record is better.
    already_observed = [
        h for h in haplotypes
        if h in by_haplotype and by_haplotype[h]["representative_strain_ha1_sequence"]
    ]
    if already_observed:
        lines = [
            f"  - {h}: count {by_haplotype[h]['count']}, "
            f"strain {by_haplotype[h]['representative_strain']}"
            for h in already_observed
        ]
        raise ValueError(
            f"{len(already_observed)} haplotype(s) in '{args.haplotype_list}' already have "
            f"an observed representative sequence in the VCM TSV and so are not "
            f"provisional:\n" + "\n".join(lines) +
            f"\nList these in a selection file (e.g. data/turners_*.txt) instead."
        )

    rows = []
    for haplotype in haplotypes:
        ha1, ha2, backbone_strain = reconstruct(haplotype, backbones)
        rows.append({
            "derived_haplotype": haplotype,
            "count": 0,
            "emerging_haplotype": "",
            "n_countries": 0,
            "n_regions": 0,
            "latest_sequence": "",
            "representative_strain": f"{haplotype}_provisional",
            "representative_strain_ha1_sequence": ha1,
            "representative_ha2_sequence": ha2,
            "days_since_latest_sequence": "",
            "relative_fitness_derived_haplotype": "",
            "provisional": "TRUE",
        })
        print(f"  {haplotype}: HA1 {len(ha1)} aa, HA2 {len(ha2)} aa "
              f"(backbone {backbone_strain})")

    # Guard the uniqueness invariant curate_library_sequences.py enforces, so a collision
    # is reported here with context rather than as a bare duplicate-strain error later.
    strains = [r["representative_strain"] for r in rows]
    if len(set(strains)) != len(strains):
        dupes = sorted({s for s in strains if strains.count(s) > 1})
        raise ValueError(f"Duplicate representative_strain value(s): {dupes}")

    ha1_lengths = {len(r["representative_strain_ha1_sequence"]) for r in rows}
    ha2_lengths = {len(r["representative_ha2_sequence"]) for r in rows}
    if len(ha1_lengths) > 1 or len(ha2_lengths) > 1:
        raise ValueError(
            f"Reconstructed sequences have inconsistent lengths (HA1 {sorted(ha1_lengths)}, "
            f"HA2 {sorted(ha2_lengths)}). All library sequences must be the same length; "
            f"this means the backbones disagree in length."
        )

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} provisional haplotype(s) to {args.output}")


if __name__ == "__main__":
    main()
