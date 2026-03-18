#!/bin/bash
set -euo pipefail

# =============================================================================
# FastANI Analysis for Trichoderma Species - Complete Workflow
# Exact code used in the research analysis (adapted for reusability)
#
# PURPOSE: Calculate Average Nucleotide Identity between Trichoderma genomes
# INPUT: List of genome assemblies (FASTA paths) or directory containing genomes
# OUTPUT: ANI matrix, heatmap visualizations (PDF/PNG), clustering plots
# TOOLS: FastANI, R (with reshape2, viridis, ape, pheatmap, ggplot2)
# PARAMETERS: All-vs-all comparison, 8 threads
# =============================================================================

# --- Default configuration (edit these or override with arguments) ---
GENOME_LIST="input/genome_list.txt"          # File with one genome path per line, or
GENOME_DIR=""                                 # Directory containing .fna/.fasta genomes (alternative)
OUTPUT_BASE="output/fastani"                  # Base output directory
SPECIES_PREFIX="Trichoderma"                  # Prefix for output files
THREADS=8                                      # Number of threads for FastANI
MIN_ANI=0                                      # Minimum ANI to report (0 = all)
FRAGMENT_LEN=0                                 # Fragment length (0 = auto)

# Step control
RUN_FASTANI="yes"
RUN_HEATMAP="yes"                              # Run R heatmap scripts
RUN_TREE_HEATMAP="yes"                         # Run tree-based heatmap

# --- Parse command-line arguments (simple overrides) ---
# Usage: ./script.sh [genome_list|genome_dir] [output_base] [threads]
if [[ $# -ge 1 ]]; then
    if [[ -f "$1" ]]; then
        GENOME_LIST="$1"
    elif [[ -d "$1" ]]; then
        GENOME_DIR="$1"
    else
        echo "ERROR: First argument must be a file (genome list) or directory (genomes)." >&2
        exit 1
    fi
fi
if [[ $# -ge 2 ]]; then OUTPUT_BASE="$2"; fi
if [[ $# -ge 3 ]]; then THREADS="$3"; fi

# --- Create output directory ---
OUTPUT_DIR="$OUTPUT_BASE"
mkdir -p "$OUTPUT_DIR"

# --- Log start ---
LOG_FILE="$OUTPUT_DIR/fastani_pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "=== FastANI pipeline started at $(date) ==="
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "----------------------------------------"

# -----------------------------------------------------------------------------
# Prepare genome list
# -----------------------------------------------------------------------------
FINAL_GENOME_LIST="$OUTPUT_DIR/genome_list.txt"

if [[ -n "$GENOME_DIR" ]]; then
    echo "Using genome directory: $GENOME_DIR"
    # Find all .fna, .fasta, .fa files (adjust extensions as needed)
    find "$GENOME_DIR" -maxdepth 1 -type f \( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \) > "$FINAL_GENOME_LIST"
    if [[ ! -s "$FINAL_GENOME_LIST" ]]; then
        echo "ERROR: No genome files found in $GENOME_DIR" >&2
        exit 1
    fi
    echo "Found $(wc -l < "$FINAL_GENOME_LIST") genome files."
elif [[ -f "$GENOME_LIST" ]]; then
    cp "$GENOME_LIST" "$FINAL_GENOME_LIST"
else
    echo "ERROR: Neither GENOME_LIST file nor GENOME_DIR provided/exists." >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Check prerequisites
# -----------------------------------------------------------------------------
MISSING_TOOLS=()
for tool in fastani; do
    if ! command -v "$tool" &> /dev/null; then
        MISSING_TOOLS+=("$tool")
    fi
done

if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
    echo "ERROR: Missing required tools: ${MISSING_TOOLS[*]}" >&2
    echo "Please install FastANI (e.g., conda install -c bioconda fastani)." >&2
    exit 1
fi

# Check R if heatmaps are enabled
if [[ "$RUN_HEATMAP" == "yes" || "$RUN_TREE_HEATMAP" == "yes" ]]; then
    if ! command -v Rscript &> /dev/null; then
        echo "WARNING: Rscript not found. Disabling heatmap generation." >&2
        RUN_HEATMAP="no"
        RUN_TREE_HEATMAP="no"
    else
        # Check required R packages
        R_PACKAGES=(reshape2 viridis ape pheatmap ggplot2)
        MISSING_PKGS=()
        for pkg in "${R_PACKAGES[@]}"; do
            if ! Rscript -e "suppressPackageStartupMessages(require('$pkg', quietly=TRUE))" 2>/dev/null; then
                MISSING_PKGS+=("$pkg")
            fi
        done
        if [[ ${#MISSING_PKGS[@]} -gt 0 ]]; then
            echo "WARNING: Missing R packages: ${MISSING_PKGS[*]}. Install them or disable heatmaps." >&2
            echo "Disabling heatmap generation." >&2
            RUN_HEATMAP="no"
            RUN_TREE_HEATMAP="no"
        fi
    fi
fi

# -----------------------------------------------------------------------------
# Step 1: Run FastANI
# -----------------------------------------------------------------------------
FASTANI_OUT="$OUTPUT_DIR/fastani_results.txt"

if [[ "$RUN_FASTANI" == "yes" ]]; then
    echo ">>> Step 1: Running FastANI all-vs-all..."
    fastani --ql "$FINAL_GENOME_LIST" --rl "$FINAL_GENOME_LIST" \
            -o "$FASTANI_OUT" -t "$THREADS" \
            ${MIN_ANI} --minFraction $MIN_ANI \
            ${FRAGMENT_LEN:+--fragLen $FRAGMENT_LEN}
    
    echo "FastANI completed. Results saved to $FASTANI_OUT"
    echo "Total comparisons: $(wc -l < "$FASTANI_OUT")"
    
    # Quick summary of ANI values
    echo "ANI summary (lowest 10):"
    awk '{print $3}' "$FASTANI_OUT" | sort -n | head -10
    echo "ANI summary (highest 10):"
    awk '{print $3}' "$FASTANI_OUT" | sort -nr | head -10
else
    echo ">>> Step 1: Skipping FastANI (using existing results if present)."
    if [[ ! -f "$FASTANI_OUT" ]]; then
        echo "ERROR: FastANI results not found at $FASTANI_OUT" >&2
        exit 1
    fi
fi

# -----------------------------------------------------------------------------
# Step 2: Generate basic heatmap (make_heatmap.R)
# -----------------------------------------------------------------------------
if [[ "$RUN_HEATMAP" == "yes" ]]; then
    echo ">>> Step 2: Generating basic ANI heatmap..."
    
    # Create R script as here-document
    cat > "$OUTPUT_DIR/make_heatmap.R" << 'EOF'
# Load required libraries
library(reshape2)
library(ggplot2)
library(viridis)

# 1. Read the FastANI results
cat("Reading FastANI results...\n")
data <- read.table("fastani_results.txt", header = FALSE, sep = "\t")
colnames(data) <- c("Query", "Reference", "ANI", "Frags", "Total")

# 2. Clean up the species names (remove paths and extensions)
cat("Cleaning up species names...\n")
clean_name <- function(x) {
    x <- gsub(".*/", "", x)           # Remove path
    x <- gsub("\\..*", "", x)          # Remove file extension
    x <- gsub("_genomic", "", x)       # Remove "_genomic"
    x <- gsub("GCA_.*?_", "", x)       # Remove GCA numbers
    return(x)
}
data$Query <- clean_name(data$Query)
data$Reference <- clean_name(data$Reference)

# 3. Convert to a matrix format for the heatmap
cat("Creating matrix for heatmap...\n)
ani_matrix <- acast(data, Query ~ Reference, value.var = "ANI")

# 4. Create the heatmap
cat("Creating heatmap...\n")
heatmap_plot <- ggplot(data, aes(x=Reference, y=Query, fill=ANI)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(ANI, 1)), color = "black", size = 2.5) +
  scale_fill_viridis(option = "plasma", name = "ANI (%)", limits = c(60, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Reference Genome", 
       y = "Query Genome", 
       title = "Average Nucleotide Identity (ANI) between Trichoderma species")

# 5. Save the plot
cat("Saving heatmap...\n")
ggsave("fastani_heatmap.pdf", heatmap_plot, width = 12, height = 10, dpi = 300)
ggsave("fastani_heatmap.png", heatmap_plot, width = 12, height = 10, dpi = 300)

cat("Done! Check fastani_heatmap.pdf and fastani_heatmap.png\n")
EOF

    cd "$OUTPUT_DIR"
    Rscript make_heatmap.R
    cd - > /dev/null
else
    echo ">>> Step 2: Skipping basic heatmap."
fi

# -----------------------------------------------------------------------------
# Step 3: Generate tree-based heatmap (tree_heatmap.R)
# -----------------------------------------------------------------------------
if [[ "$RUN_TREE_HEATMAP" == "yes" ]]; then
    echo ">>> Step 3: Generating tree-based ANI heatmap..."

    cat > "$OUTPUT_DIR/tree_heatmap.R" << 'EOF'
# Load required libraries
library(reshape2)
library(viridis)
library(ape)
library(pheatmap)

# 1. Read the FastANI results
cat("Reading FastANI results...\n")
data <- read.table("fastani_results.txt", header = FALSE, sep = "\t")
colnames(data) <- c("Query", "Reference", "ANI", "Frags", "Total")

# 2. Clean up the species names
cat("Cleaning species names...\n")
clean_name <- function(x) {
    x <- gsub(".*/", "", x)
    x <- gsub("\\..*", "", x)
    x <- gsub("_genomic", "", x)
    x <- gsub("GCA_.*?_", "", x)
    return(x)
}
data$Query <- clean_name(data$Query)
data$Reference <- clean_name(data$Reference)

# 3. Convert to a matrix format
cat("Creating ANI matrix...\n)
ani_matrix <- acast(data, Query ~ Reference, value.var = "ANI")

# 4. FILL MISSING VALUES with a low ANI value (e.g., 50) instead of 0
cat("Filling missing values...\n")
ani_matrix[is.na(ani_matrix)] <- 50  # Biologically reasonable

# 5. Create the heatmap WITH automatic clustering
cat("Creating heatmap with automatic clustering...\n")
png("ani_heatmap_auto_cluster.png", width = 3000, height = 2500, res = 300)
pheatmap(ani_matrix,
         color = viridis(100),
         fontsize_row = 8,
         fontsize_col = 8,
         main = "ANI Heatmap with Automatic Clustering",
         clustering_method = "average")
dev.off()

# 6. Alternative: Simple heatmap without dendrograms
cat("Creating simple heatmap...\n")
png("ani_heatmap_simple.png", width = 3000, height = 2500, res = 300)
pheatmap(ani_matrix,
         color = viridis(100),
         fontsize_row = 8,
         fontsize_col = 8,
         main = "ANI Heatmap (Simple)",
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

cat("Done! Check ani_heatmap_auto_cluster.png and ani_heatmap_simple.png\n")
EOF

    cd "$OUTPUT_DIR"
    Rscript tree_heatmap.R
    cd - > /dev/null
else
    echo ">>> Step 3: Skipping tree-based heatmap."
fi

echo "----------------------------------------"
echo "=== FastANI pipeline completed at $(date) ==="
echo "Log saved to: $LOG_FILE"
echo "Results are in: $OUTPUT_DIR"
