#!/bin/bash
set -euo pipefail

# =============================================================================
# Repeat Analysis Pipeline: RepeatMasker and RepeatModeler
# Exact code used in the research analysis (adapted for reusability)
#
# PURPOSE: Identify and mask repetitive elements in Trichoderma genomes
# INPUT: Genome assemblies in FASTA format (list or directory)
# OUTPUT: Repeat libraries, masked genomes, repeat annotations (GFF, HTML)
# TOOLS: RepeatModeler, RepeatMasker, RMBlast, TRF, CD-HIT
# PARAMETERS: -LTRStruct, 12 threads, cd-hit-est 0.95 identity
# =============================================================================

# --- Default configuration (edit these or override with arguments) ---
BASE_DIR="."                                         # Root directory containing genomes/
GENOME_DIR="${BASE_DIR}/genomes/fasta"                # Where individual genome FASTA files are stored
OUTPUT_BASE="${BASE_DIR}/repeat_analysis"            # Base output directory
SPECIES_NAMES=()                                      # List of species (auto-detected from files if empty)

# Step control (set to "yes" or "no")
RUN_REPEATMODELER="yes"
RUN_COMBINE_LIBRARIES="yes"
RUN_REPEATMASKER="yes"

# RepeatModeler parameters
RMODELER_THREADS=12
RMODELER_LTRSTRUCT="--LTRStruct"                     # Include LTR structural analysis

# CD-HIT parameters (for deduplication)
CDHIT_IDENTITY=0.95
CDHIT_WORDSIZE=8                                      # -n 8 for 0.95 identity
CDHIT_THREADS=4
CDHIT_MEMORY=16000                                    # Memory in MB

# RepeatMasker parameters
RMASKER_THREADS=10
RMASKER_ENGINE="ncbi"                                 # Use RMBlast

# --- Parse command-line arguments (simple overrides) ---
# Usage: ./script.sh [genome_dir] [output_base] [species1 species2 ...]
if [[ $# -ge 1 && -d "$1" ]]; then
    GENOME_DIR="$1"
    shift
fi
if [[ $# -ge 1 && ! "$1" =~ ^- ]]; then
    OUTPUT_BASE="$1"
    shift
fi
if [[ $# -gt 0 ]]; then
    SPECIES_NAMES=("$@")   # Remaining arguments are species names
fi

# --- Create output directory structure ---
mkdir -p "$OUTPUT_BASE"/{repeatmodeler_output,combined_repeat_library,repeatmasker_output}

# --- Log start ---
LOG_FILE="$OUTPUT_BASE/repeat_pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "=== Repeat analysis pipeline started at $(date) ==="
echo "Genome directory: $GENOME_DIR"
echo "Output base: $OUTPUT_BASE"
echo "----------------------------------------"

# -----------------------------------------------------------------------------
# Helper function: check required tools
# -----------------------------------------------------------------------------
check_tool() {
    if ! command -v "$1" &> /dev/null; then
        echo "ERROR: $1 not found in PATH. Please install it." >&2
        exit 1
    fi
}

check_tool "BuildDatabase"
check_tool "RepeatModeler"
check_tool "RepeatMasker"
check_tool "cd-hit-est"

# Optional: check for RMBlast/TRF configuration? Not easily done, but warn.
if [[ "$RMASKER_ENGINE" == "ncbi" ]] && ! RepeatMasker -e ncbi -version &>/dev/null; then
    echo "WARNING: RepeatMasker may not be configured to use RMBlast. Check your installation." >&2
fi

# -----------------------------------------------------------------------------
# Determine list of species
# -----------------------------------------------------------------------------
if [[ ${#SPECIES_NAMES[@]} -eq 0 ]]; then
    # Auto-detect from genome files (assuming .fna, .fasta, .fa)
    cd "$GENOME_DIR"
    SPECIES_NAMES=($(ls *.fna *.fasta *.fa 2>/dev/null | sed -E 's/\.[^.]+$//' | sort -u))
    cd - > /dev/null
    if [[ ${#SPECIES_NAMES[@]} -eq 0 ]]; then
        echo "ERROR: No genome files found in $GENOME_DIR" >&2
        exit 1
    fi
    echo "Auto-detected species: ${SPECIES_NAMES[*]}"
else
    echo "Using provided species list: ${SPECIES_NAMES[*]}"
fi

# -----------------------------------------------------------------------------
# Step 1: Run RepeatModeler for each species
# -----------------------------------------------------------------------------
if [[ "$RUN_REPEATMODELER" == "yes" ]]; then
    echo ">>> Step 1: Running RepeatModeler on each genome..."
    for species in "${SPECIES_NAMES[@]}"; do
        # Locate genome file (try common extensions)
        genome_file=""
        for ext in fna fasta fa; do
            if [[ -f "$GENOME_DIR/${species}.${ext}" ]]; then
                genome_file="$GENOME_DIR/${species}.${ext}"
                break
            fi
        done
        if [[ -z "$genome_file" ]]; then
            echo "WARNING: Genome file for $species not found. Skipping."
            continue
        fi

        outdir="$OUTPUT_BASE/repeatmodeler_output/$species"
        mkdir -p "$outdir"
        db_name="$outdir/${species}_db"

        echo "  Processing $species ..."
        BuildDatabase -name "$db_name" "$genome_file"
        RepeatModeler -database "$db_name" -threads "$RMODELER_THREADS" "$RMODELER_LTRSTRUCT" > "$outdir/run.log" 2>&1

        # Check that the families file was created
        families_file="$outdir/${species}_db-families.fa"
        if [[ ! -f "$families_file" ]]; then
            echo "WARNING: RepeatModeler did not produce $families_file for $species." >&2
        else
            echo "  RepeatModeler finished for $species. Families: $(grep -c '^>' "$families_file")"
        fi
    done
else
    echo ">>> Step 1: Skipping RepeatModeler."
fi

# -----------------------------------------------------------------------------
# Step 2: Combine and deduplicate repeat families
# -----------------------------------------------------------------------------
COMBINED_RAW="$OUTPUT_BASE/combined_repeat_library/Trichoderma_combined_repeats.fa"
COMBINED_NR="$OUTPUT_BASE/combined_repeat_library/Trichoderma_combined_repeats_nr.fa"

if [[ "$RUN_COMBINE_LIBRARIES" == "yes" ]]; then
    echo ">>> Step 2: Combining repeat families..."

    # Gather all families files
    families_files=()
    for species in "${SPECIES_NAMES[@]}"; do
        ff="$OUTPUT_BASE/repeatmodeler_output/$species/${species}_db-families.fa"
        if [[ -f "$ff" ]]; then
            families_files+=("$ff")
        fi
    done

    if [[ ${#families_files[@]} -eq 0 ]]; then
        echo "ERROR: No families files found. Cannot combine." >&2
        exit 1
    fi

    # Concatenate
    cat "${families_files[@]}" > "$COMBINED_RAW"
    echo "Combined raw library: $(grep -c '^>' "$COMBINED_RAW") sequences."

    # Deduplicate with cd-hit-est
    echo "Deduplicating with cd-hit-est (identity $CDHIT_IDENTITY)..."
    cd-hit-est \
        -i "$COMBINED_RAW" \
        -o "$COMBINED_NR" \
        -c "$CDHIT_IDENTITY" \
        -n "$CDHIT_WORDSIZE" \
        -T "$CDHIT_THREADS" \
        -M "$CDHIT_MEMORY" \
        > "$OUTPUT_BASE/combined_repeat_library/cd-hit.log" 2>&1

    echo "Non‑redundant library: $(grep -c '^>' "$COMBINED_NR") sequences."
else
    echo ">>> Step 2: Skipping library combination."
fi

# -----------------------------------------------------------------------------
# Step 3: Run RepeatMasker on each genome using the combined library
# -----------------------------------------------------------------------------
if [[ "$RUN_REPEATMASKER" == "yes" ]]; then
    echo ">>> Step 3: Masking genomes with RepeatMasker..."

    # Ensure combined library exists
    if [[ ! -f "$COMBINED_NR" ]]; then
        echo "ERROR: Non‑redundant library not found: $COMBINED_NR" >&2
        exit 1
    fi

    for species in "${SPECIES_NAMES[@]}"; do
        # Find the genome file again (use same logic as step 1)
        genome_file=""
        for ext in fna fasta fa; do
            if [[ -f "$GENOME_DIR/${species}.${ext}" ]]; then
                genome_file="$GENOME_DIR/${species}.${ext}"
                break
            fi
        done
        if [[ -z "$genome_file" ]]; then
            echo "WARNING: Genome file for $species not found. Skipping RepeatMasker."
            continue
        fi

        # Create output directory for this species
        mask_outdir="$OUTPUT_BASE/repeatmasker_output/$species"
        mkdir -p "$mask_outdir"

        echo "  Masking $species ..."
        RepeatMasker \
            -pa "$RMASKER_THREADS" \
            -lib "$COMBINED_NR" \
            -dir "$mask_outdir" \
            -gff -html -nolow -xsmall \
            -e "$RMASKER_ENGINE" \
            "$genome_file" \
            > "$mask_outdir/RepeatMasker.log" 2>&1

        echo "  Finished $species. Results in $mask_outdir"
    done
else
    echo ">>> Step 3: Skipping RepeatMasker."
fi

echo "----------------------------------------"
echo "=== Repeat analysis pipeline completed at $(date) ==="
echo "Log saved to: $LOG_FILE"
