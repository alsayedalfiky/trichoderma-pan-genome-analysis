#!/bin/bash
set -euo pipefail

# =============================================================================
# eggNOG Functional Annotation Pipeline
# Exact code used in the research analysis (adapted for reusability)
#
# PURPOSE: Comprehensive functional annotation using eggNOG-mapper
# INPUT: Predicted protein sequences (FASTA)
# OUTPUT: eggNOG annotations, COG/GO/KEGG summaries, and a distribution table
# TOOLS: eggNOG-mapper v2.+, DIAMOND
# PARAMETERS: -m diamond, --seed_ortholog_evalue 1e-05, 12 CPUs
# =============================================================================

# --- Default configuration (edit these or override with arguments) ---
INPUT_PROTEINS="input/proteins.faa"               # Protein FASTA file
OUTPUT_DIR="output/eggnog_results"                 # Base output directory
OUTPUT_PREFIX="eggnog_annotation"                  # Prefix for output files
CPU=12                                              # Number of CPUs
EVALUE="1e-05"                                      # Seed ortholog e-value

# --- Optionally allow command-line overrides ---
# Usage: ./script.sh [input.faa] [output_dir] [prefix] [cpus] [evalue]
if [[ $# -ge 1 ]]; then INPUT_PROTEINS="$1"; fi
if [[ $# -ge 2 ]]; then OUTPUT_DIR="$2"; fi
if [[ $# -ge 3 ]]; then OUTPUT_PREFIX="$3"; fi
if [[ $# -ge 4 ]]; then CPU="$4"; fi
if [[ $# -ge 5 ]]; then EVALUE="$5"; fi

# --- Check prerequisites ---
if ! command -v emapper.py &> /dev/null; then
    echo "ERROR: emapper.py (eggNOG-mapper) not found in PATH." >&2
    echo "       Please install eggNOG-mapper v2.+ (https://github.com/eggnogdb/eggnog-mapper)" >&2
    exit 1
fi

# Check that input file exists
if [[ ! -f "$INPUT_PROTEINS" ]]; then
    echo "ERROR: Input protein file not found: $INPUT_PROTEINS" >&2
    exit 1
fi

# --- Create output directory ---
mkdir -p "$OUTPUT_DIR"

# --- Log start ---
echo "Starting eggNOG-mapper at $(date)"
echo "Input: $INPUT_PROTEINS"
echo "Output directory: $OUTPUT_DIR"
echo "Prefix: $OUTPUT_PREFIX"
echo "CPUs: $CPU, e-value: $EVALUE"
echo "----------------------------------------"

# --- Run eggNOG-mapper ---
emapper.py \
    -i "$INPUT_PROTEINS" \
    -o "$OUTPUT_PREFIX" \
    --output_dir "$OUTPUT_DIR" \
    -m diamond \
    --cpu "$CPU" \
    --seed_ortholog_evalue "$EVALUE" \
    2>&1 | tee "$OUTPUT_DIR/${OUTPUT_PREFIX}.log"

echo "----------------------------------------"
echo "eggNOG-mapper finished at $(date)"
echo "Generating summary statistics..."

# --- Summary statistics (COG, GO, KEGG) ---
ANNOT_FILE="$OUTPUT_DIR/${OUTPUT_PREFIX}.emapper.annotations"
SUMMARY_FILE="$OUTPUT_DIR/${OUTPUT_PREFIX}_summary.txt"

if [[ ! -f "$ANNOT_FILE" ]]; then
    echo "ERROR: Annotation file not found: $ANNOT_FILE" >&2
    exit 1
fi

# Extract data rows (skip comments and header)
DATA_ROWS=$(grep -v '^##' "$ANNOT_FILE" | grep -v '^#query' || true)
TOTAL_PROTEINS=$(echo "$DATA_ROWS" | wc -l)

# Function to calculate percentage
calc_pct() {
    echo "scale=1; $1 * 100 / $TOTAL_PROTEINS" | bc
}

# Start summary
{
    echo "=== EGGNOG ANNOTATION SUMMARY ==="
    echo "Total proteins analyzed: $TOTAL_PROTEINS"
    echo ""

    # COG annotations (column 7 is COG functional categories)
    PROTEINS_WITH_COG=$(echo "$DATA_ROWS" | cut -f7 | grep -vE '^(-|--)$' | grep -E '[A-Z]' | wc -l)
    PCT_COG=$(calc_pct "$PROTEINS_WITH_COG")
    echo "Proteins with COG annotation: $PROTEINS_WITH_COG ($PCT_COG%)"

    # COG category distribution (split multi-letter codes)
    echo ""
    echo "COG CATEGORY DISTRIBUTION (multi-letter codes split):"
    echo "$DATA_ROWS" | cut -f7 | grep -E '[A-Z]' | grep -vE '^(-|--)$' | \
        grep -o '[A-Z]' | sort | uniq -c | sort -nr | \
        awk '{printf "  %s: %d\n", $2, $1}'

    echo ""
    echo "DETAILED COG BREAKDOWN:"
    SINGLE_LETTER=$(echo "$DATA_ROWS" | cut -f7 | grep -E '^[A-Z]$' | wc -l)
    MULTI_LETTER=$(echo "$DATA_ROWS" | cut -f7 | grep -E '[A-Z][A-Z]+' | wc -l)
    echo "  Single-letter COG proteins: $SINGLE_LETTER"
    echo "  Multi-letter COG proteins: $MULTI_LETTER"
    echo "  Total COG-annotated proteins: $PROTEINS_WITH_COG"
    echo ""

    # GO annotations (column 10)
    PROTEINS_WITH_GO=$(echo "$DATA_ROWS" | cut -f10 | grep 'GO:' | wc -l)
    PCT_GO=$(calc_pct "$PROTEINS_WITH_GO")
    echo "Proteins with GO terms: $PROTEINS_WITH_GO ($PCT_GO%)"

    # KEGG annotations (column 12)
    PROTEINS_WITH_KEGG=$(echo "$DATA_ROWS" | cut -f12 | grep 'ko:' | wc -l)
    PCT_KEGG=$(calc_pct "$PROTEINS_WITH_KEGG")
    echo "Proteins with KEGG pathways: $PROTEINS_WITH_KEGG ($PCT_KEGG%)"

} | tee "$SUMMARY_FILE"

echo "----------------------------------------"
echo "Summary saved to: $SUMMARY_FILE"
echo "eggNOG pipeline completed successfully at $(date)"
