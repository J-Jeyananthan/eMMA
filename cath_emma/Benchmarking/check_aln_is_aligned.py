import os
import click
from Bio import SeqIO


@click.command()
@click.option(
    "--aln-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
    help="Directory containing .aln files to inspect",
)
@click.option(
    "--aln-suffix",
    default="aln",
    help="Suffix of alignment files to check (default: aln)",
)
def check_aln_is_aligned(aln_dir, aln_suffix):
    """For each .aln file, report whether sequences are equal-length and gapped
    (i.e. a real alignment) or ragged/ungapped (i.e. raw unaligned FASTA)."""
    aln_files = sorted(f for f in os.listdir(aln_dir) if f.endswith(f".{aln_suffix}"))

    n_aligned = 0
    n_unaligned = 0
    n_single = 0

    for aln_file in aln_files:
        path = os.path.join(aln_dir, aln_file)
        records = list(SeqIO.parse(path, "fasta"))

        if len(records) <= 1:
            n_single += 1
            continue

        lengths = {len(r.seq) for r in records}
        has_gaps = any("-" in str(r.seq) for r in records)
        is_aligned = len(lengths) == 1 and has_gaps

        if is_aligned:
            n_aligned += 1
        else:
            n_unaligned += 1
            print(
                f"{aln_file}: NOT aligned "
                f"(n_seqs={len(records)}, n_distinct_lengths={len(lengths)}, has_gaps={has_gaps})"
            )

    print("---")
    print(f"Total files: {len(aln_files)}")
    print(f"Single-sequence (skipped): {n_single}")
    print(f"Aligned (equal-length, gapped): {n_aligned}")
    print(f"NOT aligned (ragged/ungapped): {n_unaligned}")


if __name__ == "__main__":
    check_aln_is_aligned()
