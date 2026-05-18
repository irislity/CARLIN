#!/bin/bash
# Analyzes spatial barcode TSV files:
# 1. Counts unique spatial barcodes (excluding NA)
# 2. For each unique spatial barcode, counts how many unique corrected sequences it has

RESULTS_DIR="$(dirname "$0")/results"

for tsv in "$RESULTS_DIR"/*.tsv; do
    filename=$(basename "$tsv")
    echo "=== $filename ==="

    echo "Total unique spatial_barcodes (excl NA):"
    awk -F'\t' 'NR>1 && $6 != "NA" {print $6}' "$tsv" | sort -u | wc -l

    echo ""
    echo "Unique corrected sequences per spatial_barcode (top 20 by count):"
    awk -F'\t' 'NR>1 && $6 != "NA" && $4 != "NA" {print $6, $4}' "$tsv" \
        | sort -u \
        | awk '{count[$1]++} END {for (sb in count) print sb, count[sb]}' \
        | sort -t' ' -k2 -rn \
        | head -20

    echo ""
    echo "Distribution (# of spatial_barcodes with N corrected sequences):"
    printf "%-30s %s\n" "corrected_sequences_count" "num_spatial_barcodes"
    awk -F'\t' 'NR>1 && $6 != "NA" && $4 != "NA" {print $6, $4}' "$tsv" \
        | sort -u \
        | awk '{count[$1]++} END {for (sb in count) print count[sb]}' \
        | sort -n \
        | uniq -c \
        | awk '{printf "%-30s %s\n", $2, $1}'

    echo ""
    echo "-------------------------------------------"
    echo ""
done
