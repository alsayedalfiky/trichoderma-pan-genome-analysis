# TrichoBase: A Unified Genomic Resource for Trichoderma

TrichoBase is a standardized computational framework for high-resolution comparative genomics of the genus Trichoderma. This repository contains the versioned pipeline used to establish a consistent functional and structural annotation framework across six key species.

##  Quick Start: Reproducibility
To ensure absolute reproducibility, we provide versioned Conda environments. To replicate the main functional annotation environment:

```bash
# Clone the repository
git clone [https://github.com/alsayedalfiky/trichoderma-pan-genome-analysis.git](https://github.com/alsayedalfiky/trichoderma-pan-genome-analysis.git)
cd trichoderma-pan-genome-analysis

# Create the environment from the master blueprint
conda env create -f environment.yml
conda activate func_funannotate
```

## Repository Structure
The project is organized into modular directories to separate logic from environment requirements:
- /scripts: Numbered shell (.sh) and Python (.py) scripts (01–10) covering the end-to-end pipeline.

- /envs: Specialized Conda .yml blueprints for RNA-Seq, RepeatMasking, and full versioned lockfiles.

- environment.yml: The master environment file for core functional annotation.

### Pipeline Workflow
The analysis follows a multi-step modular progression. While parameters are optimized for Trichoderma, the scripts are adaptable for other fungal pangenomes. The scripts are available in the scripts folder in this repository. 

## Key Software Versions
- Python v3.11.13

- Funannotate v1.8.17

- antiSMASH v8.0.2

- OrthoFinder v3.1.0

- CAFE v5.1.0

##  How to Run the Pipeline

Each script in the `/scripts` directory is modular and can be run independently. You can use default relative paths or provide custom arguments.

### Step 01: Genome Acquisition
This script automates the download and organization of the 6 *Trichoderma* reference genomes from NCBI.

**Usage:**
```bash
# Download the 6 study accessions into /genomes
bash scripts/01_ncbi_genome_download.sh

# Or download custom accessions by passing them as arguments
bash scripts/01_ncbi_genome_download.sh GCA_000001234.1 GCA_000005678.1

### Step 02: Gene Prediction (Script 02)
This pipeline runs `funannotate predict` using evidence-based parameters optimized for *Trichoderma*.

**Default input files (place them in the `input/` folder at the repository root):**
* `genome.fasta` – Masked assembly.
* `transcripts.fa` – Transcript evidence (FASTA).
* `proteins.faa` – Protein evidence (FASTA).

**Default output:** Results will be written to the `output/` directory, including a log file `funannotate_predict.log`.

**Prerequisites:** [Funannotate](https://github.com/nextgenusfs/funannotate) v1.8+ must be installed and available in your PATH.

**Usage:**
```bash
# Option A: Use default paths (after placing files in input/)
bash scripts/02_funannotate_annotation.sh

# Option B: Specify custom paths and metadata
bash scripts/02_funannotate_annotation.sh <genome> <outdir> <species> <strain> <transcripts> <proteins>
