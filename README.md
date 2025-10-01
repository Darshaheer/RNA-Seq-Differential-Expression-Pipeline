# RNA-Seq-Differential-Expression-Pipeline
This project demonstrates preprocessing, alignment, quantification, and differential expression analysis from Cancer and Normal RNA-Seq data.

## Repository Structure
RNA-Seq-Pipeline/
│── Input_fastq/              # Raw FASTQ files
│── Reference_data/           # Reference genome (.fa) and annotation (.gtf)
│── FastQC_Results/           # QC reports
│── Preprocessed_files/       # Trimmed FASTQ files
│── GenomeIndex/              # Genome index for alignment
│── Aligned_Results/          # BAM/SAM alignment results
│── Results/Counts/           # featureCounts output
│── rnaseq_pipeline.sh        # Main Bash pipeline script
│── rnaseq_deseq.R            # R script for DESeq2 analysis
│── README.md                 # Documentation
│── LICENSE                   # License file
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
