import logging
import os

import click
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s"
)
LOG = logging.getLogger(__name__)

GAP_CHARS = "-."


def write_degapped_fasta(records, out_path):
    with open(out_path, "w") as out_fh:
        for record in records:
            seq = str(record.seq).translate({ord(c): None for c in GAP_CHARS})
            out_fh.write(f">{record.description}\n")
            out_fh.write(f"{seq}\n")


@click.command()
@click.option(
    "--input-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
    help="Directory containing input alignment files (one per FunFam)",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False, resolve_path=True),
    help="Directory to write de-gapped .fa files",
)
@click.option(
    "--input-suffix",
    default="aln",
    help="Suffix of input alignment files to process (default: aln)",
)
def degap_aln_dir(input_dir, output_dir, input_suffix):
    """Strip gap characters ('-' and '.') from each sequence in every alignment
    file in INPUT_DIR, writing raw FASTA (.fa) files to OUTPUT_DIR with the same
    basename. Ragged/unaligned input files pass through unchanged (de-gapping
    is a no-op when there are no gap characters)."""
    os.makedirs(output_dir, exist_ok=True)

    input_files = sorted(f for f in os.listdir(input_dir) if f.endswith(f".{input_suffix}"))
    if not input_files:
        LOG.warning(f"No *.{input_suffix} files found in {input_dir}")
        return

    for input_file in input_files:
        input_path = os.path.join(input_dir, input_file)
        basename = os.path.splitext(input_file)[0]
        output_path = os.path.join(output_dir, f"{basename}.fa")

        records = list(SeqIO.parse(input_path, "fasta"))
        if not records:
            LOG.warning(f"Skipping empty file: {input_file}")
            continue

        write_degapped_fasta(records, output_path)

    LOG.info("DONE.")


if __name__ == "__main__":
    degap_aln_dir()
