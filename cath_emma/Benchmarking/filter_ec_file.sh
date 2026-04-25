#!/bin/bash
# Usage: bash filter_ec_file.sh <funfams_dir> <ec_file> <output_ec_file>

IDS=$(grep -h "^>" "$1"/*.faa "$1"/*.aln 2>/dev/null | sed 's/>//; s/\/.*//')

DUPES=$(echo "$IDS" | sort | uniq -d | wc -l)
echo "Duplicate IDs: $DUPES"

echo "$IDS" \
    | sort -u \
    | awk 'FNR==NR{ids[$1]=1; next} NR>1 && ids[$1]{print $1","$2}' - "$2" \
    > "$3"
