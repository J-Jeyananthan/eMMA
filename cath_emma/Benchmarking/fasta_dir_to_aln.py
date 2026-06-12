import logging
import os
import shutil
import subprocess
from io import StringIO

import click
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s"
)
LOG = logging.getLogger(__name__)

# Mirrors Cath::Gemma::Tool::Aligner::make_alignment_file
MAFFT_PARAMS_SLOW_HIGH_QUAL = ["--amino", "--anysymbol", "--localpair", "--maxiterate", "1000", "--quiet"]
MAFFT_PARAMS_FAST_LOW_QUAL = ["--amino", "--anysymbol", "--parttree", "--retree", "1", "--quiet"]
MANY_SEQUENCES_THRESHOLD = 200


def write_flattened_fasta(records, out_path):
    with open(out_path, "w") as out_fh:
        for record in records:
            out_fh.write(f">{record.description}\n")
            out_fh.write(f"{str(record.seq)}\n")


def align_fasta_file(fasta_path, aln_path, mafft_exe):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    num_sequences = len(records)

    if num_sequences == 0:
        LOG.warning(f"Skipping empty FASTA file: {fasta_path}")
        return

    if num_sequences == 1:
        LOG.info(f"Copying single-sequence file {fasta_path} -> {aln_path}")
        shutil.copy(fasta_path, aln_path)
        return

    mafft_params = (
        MAFFT_PARAMS_SLOW_HIGH_QUAL
        if num_sequences <= MANY_SEQUENCES_THRESHOLD
        else MAFFT_PARAMS_FAST_LOW_QUAL
    )

    command = [mafft_exe, *mafft_params, str(fasta_path)]
    LOG.info(f"About to mafft-align {num_sequences} sequences for {fasta_path}")
    LOG.debug("Running: " + " ".join(command))

    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"mafft command {' '.join(command)} failed with:\n"
            f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
        )
    if result.stderr:
        LOG.warning(f"mafft produced stderr output for {fasta_path}:\n{result.stderr}")

    # Flatten the alignment to single-line-per-sequence FASTA, as done with
    # Bio::SeqIO (-width => 32000) in Aligner.pm
    aligned_records = list(SeqIO.parse(StringIO(result.stdout), "fasta"))
    write_flattened_fasta(aligned_records, aln_path)

    total_residues = sum(len(r.seq) for r in aligned_records)
    total_gaps = sum(str(r.seq).count("-") for r in aligned_records)
    gap_percentage = (100 * total_gaps / total_residues) if total_residues else 0
    LOG.info(
        f"Finished mafft-aligning {num_sequences} sequences for {fasta_path} "
        f"(gap%={gap_percentage:.1f})"
    )


@click.command()
@click.option(
    "--input-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
    help="Directory containing input FASTA files (one per FunFam)",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False, resolve_path=True),
    help="Directory to write .aln files",
)
@click.option(
    "--fasta-suffix",
    default="fa",
    help="Suffix of input FASTA files to process (default: fa)",
)
@click.option(
    "--mafft-exe",
    default="mafft",
    help="Path to the mafft executable (default: mafft on PATH)",
)
def fasta_dir_to_aln(input_dir, output_dir, fasta_suffix, mafft_exe):
    """Align each FASTA file in INPUT_DIR with mafft, using the same protocol
    as Cath::Gemma::Tool::Aligner::make_alignment_file, writing .aln files to OUTPUT_DIR"""
    os.makedirs(output_dir, exist_ok=True)

    fasta_files = sorted(f for f in os.listdir(input_dir) if f.endswith(f".{fasta_suffix}"))
    if not fasta_files:
        LOG.warning(f"No *.{fasta_suffix} files found in {input_dir}")
        return

    for fasta_file in fasta_files:
        fasta_path = os.path.join(input_dir, fasta_file)
        basename = os.path.splitext(fasta_file)[0]
        aln_path = os.path.join(output_dir, f"{basename}.aln")

        if os.path.exists(aln_path) and os.path.getsize(aln_path) > 0:
            LOG.info(f"Skipping {fasta_file}: {aln_path} already exists")
            continue

        align_fasta_file(fasta_path, aln_path, mafft_exe)

    LOG.info("DONE.")


if __name__ == "__main__":
    fasta_dir_to_aln()
