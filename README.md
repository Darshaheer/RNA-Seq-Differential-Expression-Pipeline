# RNA-Seq-Differential-Expression-Pipeline
This project demonstrates preprocessing, alignment, quantification, and differential expression analysis from Cancer and Normal RNA-Seq data.

## Repository Structure
RNA-Seq-Differential-Expression-Pipeline/
│── RNA-Seq.sh                # Bash pipeline (QC → Trimming → Alignment → Counting)
│── DESeq-Analysis.R          # R pipeline (Differential Expression + Plotting)
│── Results/                  # Analysis outputs
│   ├── FINAL FILE.xlsx       # Differential expression results
│   ├── featurecounts.txt     # Raw featureCounts output
│   ├── featurecounts_mod.txt # Processed counts file
│   ├── Heatmap.png           # Sample distance heatmap
│   ├── PCA Plot.png          # Principal Component Analysis
│   ├── MA Plot.png           # MA plot of DE genes
│   ├── Volcano Plot.png      # Volcano plot of DE genes
│── README.md                 # Project documentation
│── LICENSE                   # License file (if included)
│── .gitignore                # Ignore unnecessary files

## Features
  Quality control with FastQC
  Trimming adapters and low-quality reads (Fastp)
  Genome indexing & alignment (BWA + Samtools)
  Quantification of read counts (featureCounts)
  Differential expression analysis (DESeq2 in R)
  Data visualization (PCA, Heatmap, MA Plot, Volcano Plot)
  Export results to Excel

## Requirements
  Linux Tools
  FastQC
  Fastp
  BWA
  Samtools
  Subread (for featureCounts)
  Install on Ubuntu/Debian:
    sudo apt-get install fastqc fastp bwa samtools subread
  R Packages
  DESeq2
  pheatmap
  apeglm
  RColorBrewer
  EnhancedVolcano
  openxlsx
  Install in R:
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(c("DESeq2", "apeglm", "EnhancedVolcano", "pheatmap"))
    install.packages(c("RColorBrewer", "openxlsx"))

## Usage
  1. Run the preprocessing, indexing, alignment & quantification:
  bash rnaseq_pipeline.sh
  2. Perform DESeq2 analysis in R:
  source("rnaseq_deseq.R")

## Example Outputs
  QC Reports → FastQC_Results/
  Trimmed FASTQ files → Preprocessed_files/
  Sorted BAMs → Aligned_Results/
  Feature Counts → Results/Counts/featurecounts_0_mod.txt
  Plots → PCA, Heatmap, MA Plot, Volcano Plot
  Excel File → FINAL FILE.xlsx

## Notes
  Replace Input_fastq/ with your FASTQ data.
  Ensure your reference genome (genome.fa) and annotation file (annotation.gtf) are in Reference_data/.
  Metadata (metadata.tsv) and counts file (counts.tsv) are required for DESeq2.

## License
  This project is open-source under the MIT License. You are free to use and modify it for research or educational purposes.
