#!/bin/bash
set -euo pipefail

# =============================================================================
# NCBI Genome Download and Processing
# Exact code used in the research analysis (adapted for reusability)
#
# PURPOSE: Download reference genomes from NCBI for multiple Trichoderma species
# INPUT: NCBI accession numbers (hardcoded below or passed as arguments)
# OUTPUT: .fna genomic files in the 'genomes/' directory
# TOOLS: ncbi-datasets-cli
# NOTES: Downloads 6 Trichoderma genomes for comparative analysis
# =============================================================================

# --- Configuration (edit these variables or override with arguments) ---
# List of accessions (space-separated)
ACCESSIONS="GCA_036169685.1 GCA_002006585.1 GCA_026401105.1 GCA_026259275.1 GCA_003025115.1 GCA_010015525.1"
# Add more accessions as needed, e.g.,:
# ACCESSIONS="$ACCESSIONS GCA_XXXXXXXXX.X"

# Output directory for genome files
OUTDIR="genomes"

# --- Optionally allow accessions as command-line arguments ---
# If arguments are provided, they replace the default list.
if [[ $# -gt 0 ]]; then
    ACCESSIONS="$*"
fi

# --- Check prerequisites ---
if ! command -v datasets &> /dev/null; then
    echo "ERROR: 'datasets' CLI not found. Please install ncbi-datasets-cli first." >&2
    echo "       See: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/" >&2
    exit 1
fi

if ! command -v unzip &> /dev/null; then
    echo "ERROR: 'unzip' not found. Please install it (e.g., sudo apt install unzip)." >&2
    exit 1
fi

# --- Create output directory ---
mkdir -p "$OUTDIR"

# --- Loop over accessions ---
for ACC in $ACCESSIONS; do
    echo "Processing $ACC..."

    # Download genome data (creates a zip file named after the accession)
    datasets download genome accession "$ACC" --include genome --filename "${ACC}.zip"

    # Unzip into a temporary directory to avoid clutter
    unzip -q -o "${ACC}.zip" -d "tmp_${ACC}"

    # Locate the genomic FASTA file (usually named 'genomic.fna' inside the accession folder)
    GENOME_FILE=$(find "tmp_${ACC}" -type f -name "*.fna" | head -n 1)

    if [[ -z "$GENOME_FILE" ]]; then
        echo "Warning: No .fna file found for $ACC. Skipping."
        rm -rf "tmp_${ACC}" "${ACC}.zip"
        continue
    fi

    # Move and rename the genome file to the output directory
    mv "$GENOME_FILE" "${OUTDIR}/${ACC}_genomic.fna"

    # Clean up temporary files
    rm -rf "tmp_${ACC}" "${ACC}.zip"

    echo "Finished $ACC → ${OUTDIR}/${ACC}_genomic.fna"
done

echo "All done. Genome files are in the '${OUTDIR}' directory."
