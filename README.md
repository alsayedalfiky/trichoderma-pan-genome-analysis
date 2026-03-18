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
