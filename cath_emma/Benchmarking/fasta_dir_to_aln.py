import logging
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO

import click
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s"
)
LOG = logging.getLogger(__name__)

# Mirrors Funfhmmer::Align::align (cath-funfhmmer), pinned to a single thread
# per mafft call so that --jobs can run many of these concurrently without
# every process trying to grab every core.
MAFFT_PARAMS_HIGH_QUAL = ["--anysymbol", "--amino", "--quiet", "--localpair", "--maxiterate", "1000", "--thread", "1"]
MAFFT_PARAMS_MID_QUAL = ["--anysymbol", "--amino", "--quiet", "--maxiterate", "2", "--thread", "1"]
MAFFT_PARAMS_LOW_QUAL = ["--anysymbol", "--amino", "--quiet", "--retree", "1", "--thread", "1"]
HIGH_QUAL_MAX_SEQUENCES = 500
MID_QUAL_MAX_SEQUENCES = 2000


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

    if num_sequences <= HIGH_QUAL_MAX_SEQUENCES:
        mafft_params = MAFFT_PARAMS_HIGH_QUAL
    elif num_sequences <= MID_QUAL_MAX_SEQUENCES:
        mafft_params = MAFFT_PARAMS_MID_QUAL
    else:
        mafft_params = MAFFT_PARAMS_LOW_QUAL

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

    # Flatten the alignment to single-line-per-sequence FASTA
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
@click.option(
    "--jobs",
    "-j",
    default=1,
    type=int,
    help="Number of mafft alignments to run concurrently (default: 1)",
)
def fasta_dir_to_aln(input_dir, output_dir, fasta_suffix, mafft_exe, jobs):
    """Align each FASTA file in INPUT_DIR with mafft, using the same protocol
    as Funfhmmer::Align::align, writing .aln files to OUTPUT_DIR"""
    os.makedirs(output_dir, exist_ok=True)

    fasta_files = sorted(f for f in os.listdir(input_dir) if f.endswith(f".{fasta_suffix}"))
    if not fasta_files:
        LOG.warning(f"No *.{fasta_suffix} files found in {input_dir}")
        return

    todo = []
    for fasta_file in fasta_files:
        fasta_path = os.path.join(input_dir, fasta_file)
        basename = os.path.splitext(fasta_file)[0]
        aln_path = os.path.join(output_dir, f"{basename}.aln")

        if os.path.exists(aln_path) and os.path.getsize(aln_path) > 0:
            LOG.info(f"Skipping {fasta_file}: {aln_path} already exists")
            continue

        todo.append((fasta_path, aln_path))

    if jobs <= 1:
        for fasta_path, aln_path in todo:
            align_fasta_file(fasta_path, aln_path, mafft_exe)
    else:
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            futures = {
                executor.submit(align_fasta_file, fasta_path, aln_path, mafft_exe): fasta_path
                for fasta_path, aln_path in todo
            }
            for future in as_completed(futures):
                future.result()

    LOG.info("DONE.")


if __name__ == "__main__":
    fasta_dir_to_aln()
