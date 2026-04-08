import os
import argparse

def filter_ec_file(funfams_dir, ec_file, output_ec_file):
    fasta_files = [f for f in os.listdir(funfams_dir) if f.endswith(".faa") or f.endswith(".aln")]

    uniprot_ids = set()
    for filename in fasta_files:
        with open(os.path.join(funfams_dir, filename), "r") as f:
            for line in f:
                if line.startswith('>'):
                    uniprot_ids.add(line.split('/')[0][1:].strip())

    print(f"Collected {len(uniprot_ids)} unique UniProt IDs from {len(fasta_files)} alignment files.")

    matched = 0
    with open(ec_file, "r") as ec_in, open(output_ec_file, "w") as ec_out:
        next(ec_in)  # skip header
        for line in ec_in:
            fields = line.split()
            uniprot_id = fields[0]
            if uniprot_id in uniprot_ids:
                ec_out.write(f"{fields[0]},{fields[1]}\n")
                matched += 1

    print(f"Wrote {matched} EC annotations to {output_ec_file}.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter a large EC annotation file to only IDs present in FunFam alignments.')
    parser.add_argument('--funfams-dir', required=True, help='Path to the folder containing FASTA alignments.')
    parser.add_argument('--ec-file', required=True, help='Path to the full EC annotation file.')
    parser.add_argument('--output-ec-file', required=True, help='Path to write the filtered EC annotation file.')
    args = parser.parse_args()
    filter_ec_file(args.funfams_dir, args.ec_file, args.output_ec_file)
