#!/bin/bash
set -euo pipefail

# =============================================================================
# Funannotate Gene Prediction Pipeline
# Exact code used in the research analysis (adapted for reusability)
# 
# PURPOSE: Structural gene annotation of Trichoderma genomes using evidence-based prediction
# INPUT: Masked genome assembly, transcript evidence, protein evidence
# OUTPUT: Predicted genes, proteins, and functional annotations
# TOOLS: Funannotate v1.8.+
# PARAMETERS: --augustus_species fusarium_graminearum, 10-12 CPUs
# =============================================================================

# --- Default configuration (edit these or override with command-line arguments) ---
BASE_DIR="."                         # Base directory for input/output (change if needed)
GENOME="${BASE_DIR}/input/genome.fasta"
OUTDIR="${BASE_DIR}/output"
SPECIES="Trichoderma species"
STRAIN="Strain-1"
TRANSCRIPTS="${BASE_DIR}/input/transcripts.fa"
PROTEINS="${BASE_DIR}/input/proteins.faa"
AUGUSTUS_SPECIES="fusarium_graminearum"
CPUS=12

# --- Command-line overrides (up to 6 arguments: genome outdir species strain transcripts proteins) ---
if [[ $# -ge 1 ]]; then GENOME="$1"; fi
if [[ $# -ge 2 ]]; then OUTDIR="$2"; fi
if [[ $# -ge 3 ]]; then SPECIES="$3"; fi
if [[ $# -ge 4 ]]; then STRAIN="$4"; fi
if [[ $# -ge 5 ]]; then TRANSCRIPTS="$5"; fi
if [[ $# -ge 6 ]]; then PROTEINS="$6"; fi

# --- Check prerequisites ---
if ! command -v funannotate &> /dev/null; then
    echo "ERROR: funannotate not found in PATH. Please install Funannotate v1.8.+." >&2
    exit 1
fi

# --- Create output directory if it doesn't exist ---
mkdir -p "$OUTDIR"

# --- Run funannotate predict ---
echo "Starting funannotate predict for $SPECIES $STRAIN"
funannotate predict \
    -i "$GENOME" \
    -o "$OUTDIR" \
    --species "$SPECIES" \
    --strain "$STRAIN" \
    --transcript_evidence "$TRANSCRIPTS" \
    --protein_evidence "$PROTEINS" \
    --augustus_species "$AUGUSTUS_SPECIES" \
    --cpus "$CPUS" \
    2>&1 | tee "$OUTDIR/funannotate_predict.log"

echo "Funannotate predict completed. Log saved to $OUTDIR/funannotate_predict.log"
