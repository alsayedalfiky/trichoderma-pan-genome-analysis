#!/bin/bash
set -euo pipefail

# =============================================================================
# RNA-Seq Analysis Pipeline for Trichoderma Species
# Exact code used in the research analysis (adapted for reusability)
#
# PURPOSE: Transcriptome assembly and validation using RNA-seq data
# INPUT: RNA-seq FASTQ files, masked genome assemblies
# OUTPUT: Transcript assemblies, alignment files, quality reports
# TOOLS: HISAT2, StringTie, samtools, FastQC, gffcompare, gffread
# =============================================================================

# --- Default configuration (edit these or override with arguments) ---
BASE_DIR="."                                         # Repository root
RNA_DIR="${BASE_DIR}/RNA_Seq"                         # RNA-Seq working directory
GENOME_DIR="${BASE_DIR}/genomes/repeatmasker_output"  # Masked genomes from Step 07

# Species to process (space‑separated)
SPECIES_NAMES=("T_harzianum" "T_asperellum" "T_longibrachiatum")  # adjust as needed

# Step control (set to "yes" or "no")
RUN_FASTQC="yes"
RUN_HISAT2_INDEX="yes"
RUN_ALIGNMENT="yes"
RUN_SAMTOOLS="yes"
RUN_STRINGTIE="yes"
RUN_MERGE="yes"
RUN_VALIDATION="yes"

# Threads
THREADS=10

# --- Parse command-line arguments (override base dir, RNA dir, species list) ---
# Usage: ./script.sh [base_dir] [rna_dir] [species1 species2 ...]
if [[ $# -ge 1 && -d "$1" ]]; then
    BASE_DIR="$1"
    shift
fi
if [[ $# -ge 1 && -d "$1" ]]; then
    RNA_DIR="$1"
    shift
fi
if [[ $# -gt 0 ]]; then
    SPECIES_NAMES=("$@")
fi

# --- Create output directories ---
mkdir -p "$RNA_DIR"/{fastqc_output,hisat2_index,alignments,transcripts,merged,validation}

# --- Log start ---
LOG_FILE="$RNA_DIR/rnaseq_pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "=== RNA-Seq pipeline started at $(date) ==="
echo "Base directory: $BASE_DIR"
echo "RNA directory: $RNA_DIR"
echo "Species: ${SPECIES_NAMES[*]}"
echo "Threads: $THREADS"
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

REQUIRED_TOOLS=(fastqc hisat2 samtools stringtie gffcompare gffread)
for tool in "${REQUIRED_TOOLS[@]}"; do
    check_tool "$tool"
done

# -----------------------------------------------------------------------------
# Step 1: FastQC quality control
# -----------------------------------------------------------------------------
if [[ "$RUN_FASTQC" == "yes" ]]; then
    echo ">>> Step 1: Running FastQC on all FASTQ files..."
    for species in "${SPECIES_NAMES[@]}"; do
        sample_dir="$RNA_DIR/rnaseq_raw/$species"
        if [[ ! -d "$sample_dir" ]]; then
            echo "WARNING: Sample directory $sample_dir not found. Skipping $species FastQC."
            continue
        fi
        outdir="$RNA_DIR/fastqc_output/$species"
        mkdir -p "$outdir"
        # Find all FASTQ files (assuming _1.fastq and _2.fastq, but we run on both)
        find "$sample_dir" -name "*.fastq" -o -name "*.fq" | while read -r fq; do
            fastqc "$fq" -o "$outdir" -t "$THREADS"
        done
    done
    echo "FastQC completed."
fi

# -----------------------------------------------------------------------------
# Step 2: HISAT2 genome indexing (per species)
# -----------------------------------------------------------------------------
if [[ "$RUN_HISAT2_INDEX" == "yes" ]]; then
    echo ">>> Step 2: Building HISAT2 indices..."
    for species in "${SPECIES_NAMES[@]}"; do
        # Find the masked genome file (from RepeatMasker output)
        genome_file=$(find "$GENOME_DIR/$species" -name "*.fna.masked" -o -name "*.fasta.masked" | head -n1)
        if [[ -z "$genome_file" ]]; then
            echo "WARNING: Masked genome not found for $species. Skipping index."
            continue
        fi
        index_dir="$RNA_DIR/hisat2_index/$species"
        mkdir -p "$index_dir"
        echo "  Building index for $species..."
        hisat2-build "$genome_file" "$index_dir/genome_index" -p "$THREADS"
    done
    echo "HISAT2 indexing completed."
fi

# -----------------------------------------------------------------------------
# Step 3: Read alignment with HISAT2
# -----------------------------------------------------------------------------
if [[ "$RUN_ALIGNMENT" == "yes" ]]; then
    echo ">>> Step 3: Aligning reads with HISAT2..."
    for species in "${SPECIES_NAMES[@]}"; do
        index_dir="$RNA_DIR/hisat2_index/$species"
        if [[ ! -f "$index_dir/genome_index.1.ht2" ]]; then
            echo "WARNING: HISAT2 index missing for $species. Skipping alignment."
            continue
        fi
        # Assume each species has a subfolder in rnaseq_raw with sample folders
        # For simplicity, we look for sample folders directly under rnaseq_raw/$species/
        sample_base="$RNA_DIR/rnaseq_raw/$species"
        if [[ ! -d "$sample_base" ]]; then
            echo "WARNING: Sample base $sample_base not found. Skipping $species alignment."
            continue
        fi
        # Iterate over subdirectories (each sample/condition)
        for sample_dir in "$sample_base"/*/ ; do
            [[ -d "$sample_dir" ]] || continue
            sample=$(basename "$sample_dir")
            # Find paired FASTQ files (assuming naming convention: sample_1.fastq, sample_2.fastq)
            r1=$(find "$sample_dir" -name "*_1.fastq" -o -name "*_1.fq" | head -n1)
            r2=$(find "$sample_dir" -name "*_2.fastq" -o -name "*_2.fq" | head -n1)
            if [[ -z "$r1" || -z "$r2" ]]; then
                echo "WARNING: Paired FASTQ files not found for $sample_dir. Skipping."
                continue
            fi
            outdir="$RNA_DIR/alignments/$species"
            mkdir -p "$outdir"
            sam="$outdir/${sample}.sam"
            echo "  Aligning $sample (species $species)..."
            hisat2 -x "$index_dir/genome_index" -1 "$r1" -2 "$r2" -S "$sam" -p "$THREADS"
        done
    done
    echo "Alignment completed."
fi

# -----------------------------------------------------------------------------
# Step 4: SAM/BAM processing (samtools view, sort, index)
# -----------------------------------------------------------------------------
if [[ "$RUN_SAMTOOLS" == "yes" ]]; then
    echo ">>> Step 4: Converting and sorting BAM files..."
    for species in "${SPECIES_NAMES[@]}"; do
        align_dir="$RNA_DIR/alignments/$species"
        if [[ ! -d "$align_dir" ]]; then
            echo "WARNING: Alignment directory $align_dir not found. Skipping $species."
            continue
        fi
        for sam in "$align_dir"/*.sam; do
            [[ -f "$sam" ]] || continue
            base=$(basename "$sam" .sam)
            bam="$align_dir/$base.bam"
            sorted="$align_dir/$base.sorted.bam"
            echo "  Processing $base ..."
            samtools view -bS "$sam" -o "$bam"
            samtools sort -@ "$THREADS" "$bam" -o "$sorted"
            samtools index "$sorted"
            # Optionally remove the intermediate BAM and SAM to save space
            # rm -f "$bam" "$sam"
        done
    done
    echo "SAM/BAM processing completed."
fi

# -----------------------------------------------------------------------------
# Step 5: Transcript assembly with StringTie (per sample)
# -----------------------------------------------------------------------------
if [[ "$RUN_STRINGTIE" == "yes" ]]; then
    echo ">>> Step 5: Assembling transcripts with StringTie..."
    for species in "${SPECIES_NAMES[@]}"; do
        align_dir="$RNA_DIR/alignments/$species"
        if [[ ! -d "$align_dir" ]]; then
            echo "WARNING: Alignment directory $align_dir not found. Skipping $species."
            continue
        fi
        outdir="$RNA_DIR/transcripts/$species"
        mkdir -p "$outdir"
        # Process each sorted BAM
        for bam in "$align_dir"/*.sorted.bam; do
            [[ -f "$bam" ]] || continue
            base=$(basename "$bam" .sorted.bam)
            gtf="$outdir/$base.gtf"
            echo "  Assembling $base ..."
            stringtie "$bam" -o "$gtf" -p "$THREADS" -l "$base"
        done
    done
    echo "StringTie assembly completed."
fi

# -----------------------------------------------------------------------------
# Step 6: Merge GTFs for each species
# -----------------------------------------------------------------------------
if [[ "$RUN_MERGE" == "yes" ]]; then
    echo ">>> Step 6: Merging GTFs per species..."
    for species in "${SPECIES_NAMES[@]}"; do
        transcript_dir="$RNA_DIR/transcripts/$species"
        if [[ ! -d "$transcript_dir" ]]; then
            echo "WARNING: Transcript directory $transcript_dir not found. Skipping $species."
            continue
        fi
        # Create list of GTF files
        gtf_list="$RNA_DIR/merged/${species}_gtf_list.txt"
        ls "$transcript_dir"/*.gtf > "$gtf_list" 2>/dev/null || true
        if [[ ! -s "$gtf_list" ]]; then
            echo "WARNING: No GTF files found for $species. Skipping merge."
            continue
        fi
        merged_gtf="$RNA_DIR/merged/${species}_merged.gtf"
        echo "  Merging $species GTFs..."
        stringtie --merge -p "$THREADS" -o "$merged_gtf" "$gtf_list"
    done
    echo "Merging completed."
fi

# -----------------------------------------------------------------------------
# Step 7: Transcriptome validation (gffcompare and gffread)
# -----------------------------------------------------------------------------
if [[ "$RUN_VALIDATION" == "yes" ]]; then
    echo ">>> Step 7: Validating merged transcriptomes..."
    for species in "${SPECIES_NAMES[@]}"; do
        merged_gtf="$RNA_DIR/merged/${species}_merged.gtf"
        if [[ ! -f "$merged_gtf" ]]; then
            echo "WARNING: Merged GTF not found for $species. Skipping validation."
            continue
        fi
        # Find masked genome for this species
        genome_file=$(find "$GENOME_DIR/$species" -name "*.fna.masked" -o -name "*.fasta.masked" | head -n1)
        if [[ -z "$genome_file" ]]; then
            echo "WARNING: Masked genome not found for $species; gffread will skip."
            genome_opt=""
        else
            genome_opt="-g $genome_file"
        fi

        val_dir="$RNA_DIR/validation/$species"
        mkdir -p "$val_dir"

        # gffread structural validation
        echo "  Running gffread for $species..."
        gffread "$merged_gtf" $genome_opt -E > "$val_dir/gffread.errors.txt" 2>&1

        # gffcompare statistics (self‑comparison)
        echo "  Running gffcompare for $species..."
        gffcompare -o "$val_dir/compare" "$merged_gtf"
        if [[ -f "$val_dir/compare.stats" ]]; then
            echo "    gffcompare stats:"
            head -n 20 "$val_dir/compare.stats"
        fi
    done
    echo "Validation completed."
fi

echo "----------------------------------------"
echo "=== RNA-Seq pipeline completed at $(date) ==="
echo "Log saved to: $LOG_FILE"
