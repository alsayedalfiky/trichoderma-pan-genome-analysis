#!/bin/bash
set -euo pipefail

# =============================================================================
# dbCAN CAZyme Annotation Pipeline
# Exact code used in the research analysis (adapted for reusability)
#
# PURPOSE: Identify carbohydrate-active enzymes in Trichoderma proteomes
# INPUT: Protein sequences (FASTA)
# OUTPUT: CAZyme families, validated lists, class distribution, Excel summary
# TOOLS: DIAMOND, HMMER, seqtk, awk, ssconvert (optional)
# PARAMETERS: DIAMOND E-value < 1e-11, identity > 30%, coverage > 70%;
#             HMMER domain E-value < 1e-15
# =============================================================================

# --- Default configuration (edit these or override with arguments) ---
PROTEOME="input/proteins.faa"                       # Protein FASTA file
OUTPUT_BASE="output/dbcan_results"                   # Base output directory
SPECIES_DIR=""                                        # Species subfolder (auto-derived if empty)
DIAMOND_DB="databases/dbCAN2/diamond_db.dmnd"        # DIAMOND dbCAN2 database
HMM_DB="databases/dbCAN2/dbCAN-HMMdb-V13.txt"        # HMMER dbCAN HMM database
CPU=12                                                # CPUs for DIAMOND and HMMER

# Step control (set to "yes" or "no")
RUN_DIAMOND="yes"
RUN_HMMER="yes"
RUN_SUMMARY="yes"
RUN_EXCEL="yes"                                       # Requires ssconvert (Gnumeric)

# --- Parse command-line arguments (simple override of input and output) ---
# Usage: ./script.sh [proteome.faa] [output_base] [species_dir]
if [[ $# -ge 1 ]]; then PROTEOME="$1"; fi
if [[ $# -ge 2 ]]; then OUTPUT_BASE="$2"; fi
if [[ $# -ge 3 ]]; then SPECIES_DIR="$3"; fi

# --- Derive species directory from input filename if not provided ---
if [[ -z "$SPECIES_DIR" ]]; then
    BASENAME=$(basename "$PROTEOME" .fa)
    BASENAME=${BASENAME%.faa}   # handle .faa as well
    SPECIES_DIR="$BASENAME"
fi

# Full output path for this species
OUTPUT_DIR="$OUTPUT_BASE/$SPECIES_DIR"

# --- Check prerequisites ---
MISSING_TOOLS=()
for tool in diamond hmmsearch seqtk awk; do
    if ! command -v "$tool" &> /dev/null; then
        MISSING_TOOLS+=("$tool")
    fi
done
if [[ "$RUN_EXCEL" == "yes" ]] && ! command -v ssconvert &> /dev/null; then
    echo "WARNING: ssconvert not found. Excel conversion will be skipped." >&2
    RUN_EXCEL="no"
fi

if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
    echo "ERROR: Missing required tools: ${MISSING_TOOLS[*]}" >&2
    echo "Please install them and ensure they are in PATH." >&2
    exit 1
fi

# Check that input file exists
if [[ ! -f "$PROTEOME" ]]; then
    echo "ERROR: Input protein file not found: $PROTEOME" >&2
    exit 1
fi

# Check database files if those steps are enabled
if [[ "$RUN_DIAMOND" == "yes" && ! -f "$DIAMOND_DB" ]]; then
    echo "ERROR: DIAMOND database not found: $DIAMOND_DB" >&2
    exit 1
fi
if [[ "$RUN_HMMER" == "yes" && ! -f "$HMM_DB" ]]; then
    echo "ERROR: HMMER database not found: $HMM_DB" >&2
    exit 1
fi

# --- Create output directory ---
mkdir -p "$OUTPUT_DIR"

# --- Log start ---
LOG_FILE="$OUTPUT_DIR/dbcan_pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "=== dbCAN pipeline started at $(date) ==="
echo "Proteome: $PROTEOME"
echo "Output directory: $OUTPUT_DIR"
echo "Species: $SPECIES_DIR"
echo "CPU: $CPU"
echo "----------------------------------------"

cd "$OUTPUT_DIR"

# -----------------------------------------------------------------------------
# Step 1: DIAMOND search against CAZy database
# -----------------------------------------------------------------------------
if [[ "$RUN_DIAMOND" == "yes" ]]; then
    echo ">>> Step 1: Running DIAMOND BLASTp against CAZy database..."
    
    diamond blastp \
        --query "$PROTEOME" \
        --db "$DIAMOND_DB" \
        --out diamond.full.tsv \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen \
        --threads "$CPU"
    
    echo "DIAMOND search complete. Total hits: $(wc -l < diamond.full.tsv)"
    
    # Apply strict filter: E-value < 1e-11, identity > 30%, coverage > 70%
    echo "Filtering DIAMOND hits..."
    awk -F'\t' '$11 < 1e-11 && $3 > 30 && (($8 - $7 + 1) / $13 * 100) > 70' diamond.full.tsv \
        > diamond.filtered.strict.tsv
    
    echo "Filtered hits: $(wc -l < diamond.filtered.strict.tsv)"
    
    # Select best hit per protein (highest bit score)
    sort -k1,1 -k12,12nr diamond.filtered.strict.tsv | awk '!seen[$1]++' > diamond.final_cazymes.tsv
    
    echo "Best hits per protein: $(wc -l < diamond.final_cazymes.tsv)"
    
    # Extract potential CAZyme IDs and sequences for HMMER validation
    cut -f1 diamond.final_cazymes.tsv > potential_cazyme_ids.txt
    seqtk subseq "$PROTEOME" potential_cazyme_ids.txt > potential_cazymes.fa
else
    echo ">>> Step 1: Skipping DIAMOND (using existing files if present)."
fi

# -----------------------------------------------------------------------------
# Step 2: HMMER validation
# -----------------------------------------------------------------------------
if [[ "$RUN_HMMER" == "yes" && -f potential_cazymes.fa ]]; then
    echo ">>> Step 2: Running HMMER validation..."
    
    hmmsearch --domtblout hmmer_validation.tsv --cpu "$CPU" "$HMM_DB" potential_cazymes.fa > hmmer.log
    
    # Extract protein IDs with domain E-value < 1e-15
    awk '$12 < 1e-15 && $12 != "-"' hmmer_validation.tsv | cut -f1 -d ' ' | sort -u > validated_cazyme_ids.txt
    
    echo "HMMER validation complete. Validated proteins: $(wc -l < validated_cazyme_ids.txt)"
    
    # Create final validated CAZyme table (intersect with diamond final hits)
    awk 'NR==FNR {valid[$1]; next} $1 in valid' validated_cazyme_ids.txt diamond.final_cazymes.tsv \
        > final_validated_cazymes.tsv
    
    echo "Final validated CAZymes: $(wc -l < final_validated_cazymes.tsv)"
    
    # Extract validated protein sequences
    seqtk subseq "$PROTEOME" validated_cazyme_ids.txt > "${SPECIES_DIR}_CAZymes.final.fa"
else
    echo ">>> Step 2: Skipping HMMER validation."
fi

# -----------------------------------------------------------------------------
# Step 3: Summary statistics
# -----------------------------------------------------------------------------
if [[ "$RUN_SUMMARY" == "yes" && -f final_validated_cazymes.tsv ]]; then
    echo ">>> Step 3: Generating CAZyme class distribution..."
    
    # Count families by class (GH, GT, PL, CE, AA, CBM)
    awk -F'\t' '{
        split($2, a, "|");
        fam = a[length(a)];
        if (fam ~ /^GH/) gh++;
        else if (fam ~ /^GT/) gt++;
        else if (fam ~ /^PL/) pl++;
        else if (fam ~ /^CE/) ce++;
        else if (fam ~ /^AA/) aa++;
        else if (fam ~ /^CBM/) cbm++;
    } END {
        print "GH:", gh+0;
        print "GT:", gt+0;
        print "PL:", pl+0;
        print "CE:", ce+0;
        print "AA:", aa+0;
        print "CBM:", cbm+0;
        print "Total:", gh+gt+pl+ce+aa+cbm;
    }' final_validated_cazymes.tsv > "${SPECIES_DIR}_CAZyme_class_distribution.txt"
    
    echo "Distribution saved."
    cat "${SPECIES_DIR}_CAZyme_class_distribution.txt"
    
    # Also create a detailed TSV with best hits for each validated protein
    # Use original diamond.full.tsv to get all info, then filter by validated IDs
    if [[ -f diamond.full.tsv ]]; then
        awk 'NR==FNR {valid[$1]; next} $1 in valid' validated_cazyme_ids.txt diamond.full.tsv \
            | sort -k1,1 -k12,12nr \
            | awk '!seen[$1]++' \
            | awk -F'\t' 'BEGIN {OFS="\t"; print "ID", "SseqID", "Qstart", "Qend", "Sstart", "Send", "evalue", "score", "pident", "qcovs", "CAZy_Family"} {
                split($2, a, "|");
                family = a[length(a)];
                qcovs = ($8 - $7 + 1) / $13 * 100;
                print $1, $2, $7, $8, $9, $10, $11, $12, $3, qcovs, family
              }' > "${SPECIES_DIR}_CAZymes_detailed_BEST.tsv"
        
        echo "Detailed best-hit TSV created."
    fi
else
    echo ">>> Step 3: Skipping summary generation."
fi

# -----------------------------------------------------------------------------
# Step 4: Convert TSV to Excel (if ssconvert available)
# -----------------------------------------------------------------------------
if [[ "$RUN_EXCEL" == "yes" && -f "${SPECIES_DIR}_CAZymes_detailed_BEST.tsv" ]]; then
    echo ">>> Step 4: Converting TSV to Excel format..."
    
    # ssconvert sometimes struggles with large TSV; ensure clean input
    # Use a temporary CSV if needed, but direct conversion usually works
    ssconvert "${SPECIES_DIR}_CAZymes_detailed_BEST.tsv" "${SPECIES_DIR}_CAZymes_detailed_BEST.xlsx" 2>/dev/null \
        || echo "Warning: ssconvert failed; Excel file not created."
    
    if [[ -f "${SPECIES_DIR}_CAZymes_detailed_BEST.xlsx" ]]; then
        echo "Excel file created: ${SPECIES_DIR}_CAZymes_detailed_BEST.xlsx"
    fi
fi

echo "----------------------------------------"
echo "=== dbCAN pipeline completed at $(date) ==="
echo "Log saved to: $LOG_FILE"
