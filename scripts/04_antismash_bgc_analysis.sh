#!/bin/bash
set -euo pipefail

# =============================================================================
# antiSMASH Biosynthetic Gene Cluster Analysis Pipeline
# Exact code used in the research analysis (adapted for reusability)
#
# PURPOSE: Identify secondary metabolite biosynthetic gene clusters in Trichoderma genomes
# INPUT: GenBank file (raw or pre-filtered), optional MIBiG and Pfam databases
# OUTPUT: BGC predictions, extracted proteins, MIBiG BLAST results, Pfam domains
# TOOLS: antiSMASH v6.0, Python 3, jq, BLAST+, HMMER
# =============================================================================

# --- Default configuration (edit these or override with arguments) ---
INPUT_GBK="input/genome.gbk"                    # Raw GenBank file
OUTPUT_BASE="output/antismash_results"           # Base output directory
SPECIES_DIR=""                                   # Species subfolder (auto-derived if empty)
TAXON="fungi"                                     # antiSMASH taxon
CPUS=12                                           # CPUs for antiSMASH and BLAST

# Database paths (set these if you want to run BLAST and Pfam steps)
MIBIG_DB=""                                       # Path to MIBiG BLAST database (without .phr etc.)
PFAM_DB=""                                        # Path to Pfam HMM database (Pfam-A.hmm)

# Step control (set to "yes" or "no")
RUN_FILTER="yes"                                  # Filter GenBank for CDS-containing contigs
RUN_ANTISMASH="yes"                               # Run antiSMASH
RUN_PROTEIN_EXTRACT="yes"                          # Extract proteins from antiSMASH GBKs
RUN_MIBIG_BLAST="yes"                              # BLAST against MIBiG (requires MIBIG_DB)
RUN_PFAM="yes"                                     # Pfam domain search (requires PFAM_DB)

# --- Parse command-line arguments (simple override of input and output) ---
# Usage: ./script.sh [input.gbk] [output_dir] [species_dir]
if [[ $# -ge 1 ]]; then INPUT_GBK="$1"; fi
if [[ $# -ge 2 ]]; then OUTPUT_BASE="$2"; fi
if [[ $# -ge 3 ]]; then SPECIES_DIR="$3"; fi

# --- Derive species directory from input filename if not provided ---
if [[ -z "$SPECIES_DIR" ]]; then
    BASENAME=$(basename "$INPUT_GBK" .gbk)
    # Remove common prefixes like "Trichoderma_" or just use the base name
    SPECIES_DIR="$BASENAME"
fi

# Full output path for this species
OUTPUT_DIR="$OUTPUT_BASE/$SPECIES_DIR"

# --- Check prerequisites ---
MISSING_TOOLS=()
for tool in python3 jq antismash blastp hmmscan; do
    if ! command -v "$tool" &> /dev/null; then
        MISSING_TOOLS+=("$tool")
    fi
done

if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
    echo "ERROR: Missing required tools: ${MISSING_TOOLS[*]}" >&2
    echo "Please install them and ensure they are in PATH." >&2
    exit 1
fi

# Check that input file exists
if [[ ! -f "$INPUT_GBK" ]]; then
    echo "ERROR: Input GenBank file not found: $INPUT_GBK" >&2
    exit 1
fi

# Check database paths if those steps are enabled
if [[ "$RUN_MIBIG_BLAST" == "yes" && ( -z "$MIBIG_DB" || ! -f "${MIBIG_DB}.phr" ) ]]; then
    echo "ERROR: MIBiG BLAST database not found or incomplete: $MIBIG_DB" >&2
    echo "       Please set MIBIG_DB correctly or disable RUN_MIBIG_BLAST." >&2
    exit 1
fi

if [[ "$RUN_PFAM" == "yes" && ( -z "$PFAM_DB" || ! -f "$PFAM_DB" ) ]]; then
    echo "ERROR: Pfam HMM database not found: $PFAM_DB" >&2
    echo "       Please set PFAM_DB correctly or disable RUN_PFAM." >&2
    exit 1
fi

# --- Create output directory ---
mkdir -p "$OUTPUT_DIR"

# --- Log start ---
LOG_FILE="$OUTPUT_DIR/antismash_pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "=== antiSMASH pipeline started at $(date) ==="
echo "Input GBK: $INPUT_GBK"
echo "Output directory: $OUTPUT_DIR"
echo "Species: $SPECIES_DIR"
echo "CPUs: $CPUS"
echo "----------------------------------------"

# -----------------------------------------------------------------------------
# Step 1: Filter GenBank for contigs with CDS features (optional)
# -----------------------------------------------------------------------------
FILTERED_GBK="$OUTPUT_DIR/$(basename "${INPUT_GBK%.gbk}").filtered.gbk"

if [[ "$RUN_FILTER" == "yes" ]]; then
    echo ">>> Step 1: Filtering GenBank file to keep only contigs with CDS features..."
    
    # Create a temporary Python script
    cat > "$OUTPUT_DIR/filter_genbank.py" << 'EOF'
#!/usr/bin/env python3
from Bio import SeqIO
import sys

def filter_genbank(input_file, output_file):
    print(f"Reading records from {input_file}...")
    all_records = list(SeqIO.parse(input_file, "genbank"))
    print(f"Found {len(all_records)} total contigs.")

    filtered_records = []
    for record in all_records:
        cds_count = sum(1 for feature in record.features if feature.type == "CDS")
        if cds_count > 0:
            filtered_records.append(record)
            print(f"  Keeping contig '{record.id}' with {cds_count} CDS features.")
        else:
            print(f"  Discarding contig '{record.id}' with 0 CDS features.")

    print(f"Keeping {len(filtered_records)} contigs with CDS features.")
    if filtered_records:
        SeqIO.write(filtered_records, output_file, "genbank")
        print(f"Successfully wrote filtered records to {output_file}")
    else:
        print("Error: No records with CDS features were found!", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filter_genbank.py <input.gbk> <output.gbk>")
        sys.exit(1)
    filter_genbank(sys.argv[1], sys.argv[2])
EOF

    python3 "$OUTPUT_DIR/filter_genbank.py" "$INPUT_GBK" "$FILTERED_GBK"
    echo "Filtering complete."
else
    echo ">>> Step 1: Skipping filter (using raw input as filtered)."
    FILTERED_GBK="$INPUT_GBK"
fi

# -----------------------------------------------------------------------------
# Step 2: Run antiSMASH
# -----------------------------------------------------------------------------
ANTISMASH_DIR="$OUTPUT_DIR/antismash"
mkdir -p "$ANTISMASH_DIR"

if [[ "$RUN_ANTISMASH" == "yes" ]]; then
    echo ">>> Step 2: Running antiSMASH..."
    antismash \
        "$FILTERED_GBK" \
        --taxon "$TAXON" \
        --cb-general \
        --cc-mibig \
        --cpus "$CPUS" \
        --output-dir "$ANTISMASH_DIR" \
        --logfile "$ANTISMASH_DIR/antismash.log"
    echo "antiSMASH finished."
else
    echo ">>> Step 2: Skipping antiSMASH (assume results exist in $ANTISMASH_DIR)."
fi

# -----------------------------------------------------------------------------
# Step 3: Extract proteins from antiSMASH region GBK files
# -----------------------------------------------------------------------------
PROTEIN_FASTA="$OUTPUT_DIR/BGC_proteins.faa"
if [[ "$RUN_PROTEIN_EXTRACT" == "yes" ]]; then
    echo ">>> Step 3: Extracting protein sequences from antiSMASH region GBK files..."
    
    # Find all .gbk files in the antiSMASH directory (excluding the filtered input)
    mapfile -t GBK_FILES < <(find "$ANTISMASH_DIR" -maxdepth 1 -name "*.gbk" ! -name "*.filtered.gbk")
    
    if [[ ${#GBK_FILES[@]} -eq 0 ]]; then
        echo "WARNING: No region GBK files found in $ANTISMASH_DIR. Skipping protein extraction."
    else
        # Create extraction script
        cat > "$OUTPUT_DIR/extract_proteins.py" << 'EOF'
#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO

gbk_files = sys.argv[1:-1]
output_faa = sys.argv[-1]

extracted = set()
with open(output_faa, "w") as out_f:
    for gbk in gbk_files:
        for record in SeqIO.parse(gbk, "genbank"):
            for feature in record.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    prot_id = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                    if prot_id not in extracted:
                        translation = feature.qualifiers["translation"][0]
                        out_f.write(f">{prot_id}\n{translation}\n")
                        extracted.add(prot_id)
print(f"Extracted {len(extracted)} proteins to {output_faa}")
EOF

        python3 "$OUTPUT_DIR/extract_proteins.py" "${GBK_FILES[@]}" "$PROTEIN_FASTA"
        echo "Extraction complete. Found $(grep -c '^>' "$PROTEIN_FASTA") proteins."
    fi
else
    echo ">>> Step 3: Skipping protein extraction."
fi

# -----------------------------------------------------------------------------
# Step 4: Analyze antiSMASH results (BGC types, hybrids) using jq
# -----------------------------------------------------------------------------
JSON_FILE=$(find "$ANTISMASH_DIR" -name "*.json" | head -n 1)
if [[ -f "$JSON_FILE" ]]; then
    echo ">>> Analyzing antiSMASH JSON results..."
    ANALYSIS_FILE="$OUTPUT_DIR/antismash_analysis.txt"
    {
        echo "=== antiSMASH Region Analysis ==="
        echo "JSON file: $JSON_FILE"
        echo ""
        
        # All BGC types count
        echo "BGC type counts:"
        jq -r '.records[].features[] | select(.type == "region") | .qualifiers.product[]' "$JSON_FILE" \
            | sort | uniq -c | sort -nr
        echo ""
        
        # Hybrid regions (multiple unique products)
        echo "Hybrid regions (multiple unique products):"
        jq -r '.records[].features[] | select(.type == "region") | {region: .qualifiers.region_number[0], product: .qualifiers.product[]} | "\(.region)\t\(.product)"' "$JSON_FILE" \
            | sort -u \
            | awk '
                {
                    products[$1] = (products[$1] ? products[$1] "; " : "") $2
                }
                END {
                    for (region in products) {
                        split(products[region], arr, "; ");
                        delete seen;
                        for (i in arr) seen[arr[i]]++;
                        if (length(seen) > 1) {
                            printf "Region %s: %s\n", region, products[region];
                        }
                    }
                }'
        echo ""
        
        # Count of hybrid regions
        hybrid_count=$(jq -r '.records[].features[] | select(.type == "region") | {region: .qualifiers.region_number[0], product: .qualifiers.product[]} | "\(.region)\t\(.product)"' "$JSON_FILE" \
            | sort -u \
            | awk '
                {
                    products[$1] = (products[$1] ? products[$1] "; " : "") $2
                }
                END {
                    count=0;
                    for (region in products) {
                        split(products[region], arr, "; ");
                        delete seen;
                        for (i in arr) seen[arr[i]]++;
                        if (length(seen) > 1) count++;
                    }
                    print count;
                }')
        echo "Total hybrid regions: $hybrid_count"
    } | tee "$ANALYSIS_FILE"
else
    echo "WARNING: No JSON file found in $ANTISMASH_DIR; skipping antiSMASH analysis."
fi

# -----------------------------------------------------------------------------
# Step 5: BLAST against MIBiG database (optional)
# -----------------------------------------------------------------------------
if [[ "$RUN_MIBIG_BLAST" == "yes" && -f "$PROTEIN_FASTA" ]]; then
    echo ">>> Step 5: BLASTing proteins against MIBiG database..."
    cd "$OUTPUT_DIR"
    
    BLAST_OUT="mibig_blast.tsv"
    blastp \
        -query "$PROTEIN_FASTA" \
        -db "$MIBIG_DB" \
        -out "$BLAST_OUT" \
        -evalue 1e-5 \
        -num_threads "$CPUS" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    
    # Get top hit per query
    awk '!seen[$1]++' "$BLAST_OUT" > mibig_blast_top.tsv
    
    # Filter by significance (≥30% identity, ≤1e-10, length ≥100)
    awk '$3 >= 30 && $11 <= 1e-10 && $4 >= 100' mibig_blast_top.tsv > mibig_blast_filtered.tsv
    
    echo "BLAST summary:"
    echo "  Total hits: $(wc -l < "$BLAST_OUT")"
    echo "  Unique query hits: $(wc -l < mibig_blast_top.tsv)"
    echo "  High-confidence hits: $(wc -l < mibig_blast_filtered.tsv)"
    
    # Convert to Excel-friendly format
    if [[ -s mibig_blast_filtered.tsv ]]; then
        cat > "$OUTPUT_DIR/blast_to_excel.py" << 'EOF'
#!/usr/bin/env python3
import pandas as pd

df = pd.read_csv('mibig_blast_filtered.tsv', sep='\t', header=None)
blast_columns = [
    'Query ID', 'Subject ID Full', 'Percentage Identity', 'Alignment Length',
    'Mismatches', 'Gap Opens', 'Query Start', 'Query End',
    'Subject Start', 'Subject End', 'E-value', 'Bit Score'
]
df.columns = blast_columns

def parse_mibig_id(full_id):
    parts = str(full_id).split('|')
    clean_id = parts[0] if parts else full_id
    annotation = parts[-2] if len(parts) > 2 else 'Unknown'
    return clean_id, annotation

df[['MIBiG Subject ID', 'Annotation']] = df['Subject ID Full'].apply(
    lambda x: pd.Series(parse_mibig_id(x))
)

final_df = df[[
    'Query ID', 'MIBiG Subject ID', 'Percentage Identity',
    'Alignment Length', 'E-value', 'Bit Score', 'Annotation'
]].rename(columns={'Percentage Identity': '% Identity'})

final_df.to_excel('mibig_blast_filtered.xlsx', index=False, engine='openpyxl')
print("Excel file created: mibig_blast_filtered.xlsx")
EOF
        python3 "$OUTPUT_DIR/blast_to_excel.py"
    else
        echo "No high-confidence hits to convert."
    fi
    cd - > /dev/null
elif [[ "$RUN_MIBIG_BLAST" == "yes" ]]; then
    echo "WARNING: Protein FASTA not found; skipping MIBiG BLAST."
fi

# -----------------------------------------------------------------------------
# Step 6: Pfam domain analysis with HMMER (optional)
# -----------------------------------------------------------------------------
if [[ "$RUN_PFAM" == "yes" && -f "$PROTEIN_FASTA" ]]; then
    echo ">>> Step 6: Running Pfam domain search with HMMER..."
    cd "$OUTPUT_DIR"
    
    hmmscan --cpu "$CPUS" --domtblout pfam.domtblout "$PFAM_DB" "$PROTEIN_FASTA" > hmmscan.log
    
    # Count significant domains (E-value < 1e-5)
    total_domains=$(grep -v '^#' pfam.domtblout | awk '$7 < 1e-5' | wc -l)
    unique_proteins=$(grep -v '^#' pfam.domtblout | awk '$7 < 1e-5 {print $4}' | sort -u | wc -l)
    
    echo "Pfam results:"
    echo "  Significant domain hits: $total_domains"
    echo "  Proteins with at least one domain: $unique_proteins"
    
    cd - > /dev/null
elif [[ "$RUN_PFAM" == "yes" ]]; then
    echo "WARNING: Protein FASTA not found; skipping Pfam."
fi

# -----------------------------------------------------------------------------
# Cleanup temporary Python scripts (optional)
# -----------------------------------------------------------------------------
rm -f "$OUTPUT_DIR"/filter_genbank.py "$OUTPUT_DIR"/extract_proteins.py "$OUTPUT_DIR"/blast_to_excel.py

echo "----------------------------------------"
echo "=== antiSMASH pipeline completed at $(date) ==="
echo "Log saved to: $LOG_FILE"
