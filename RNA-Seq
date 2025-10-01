#RNA-Seq Analysis
#!/bin/bash
# -----------------------1. Preprocessing Script -------------------------------------
# Quality Check
mkdir -p FastQC_Results Preprocessed_files

for R1_File in ./Input_fastq/*_1.fastq; do
	filename=$(basename $R1_File _1.fastq)
	R2_File=./Input_fastq/${filename}_2.fastq

	fastqc $R1_File $R2_File --outdir FastQC_Results

# Trimming 		
	fastp \
	-i $R1_File \
	-I $R2_File \
	-o Preprocessed_files/${filename}_1_trimmed.fastq \
	-O Preprocessed_files/${filename}_2_trimmed.fastq \
	-h Preprocessed_files/${filename}.html \
	-j Preprocessed_files/${filename}.json
done
echo "The PreProcessing is Done"

# -----------------------2. Indexing Script -------------------------------------
Reference_genome=../Reference_data/genome.fa
Genome_index=./GenomeIndex

mkdir -p $Genome_index

cd $Genome_index
bwa index $Reference_genome
echo "The Indexing is Done"

# -----------------------3. Alignment Script -------------------------------------
Input_files=./Preprocessed_files
Genomic_index_file=./GenomeIndex/genome.fa
Output_DIR=./Aligned_Results

mkdir -p $Output_DIR

for R1_File in $Input_files/*_1_trimmed.fastq; do
	filename=$(basename $R1_File _1_trimmed.fastq)
	R2_File=$Input_files/${filename}_2_trimmed.fastq

	if [[ -f $R1_File && -f $R2_File ]]; then
		bwa mem -t 24 $Genomic_index_file $R1_File $R2_File > $Output_DIR/${filename}.sam
		
		samtools view -b -S $Output_DIR/${filename}.sam > $Output_DIR/${filename}.bam
		samtools sort -o $Output_DIR/${filename}_sorted.bam $Output_DIR/${filename}.bam
		samtools index $Output_DIR/${filename}_sorted.bam
		
		rm $Output_DIR/${filename}.sam
	else
		echo "The paired file was not found"
	fi
done
echo "The Alignment is Complete"

# -----------------------4. Quantification Script -------------------------------------
BAM_DIR=./Aligned_Results
Ann_File=./Reference_data/annotation.gtf
OUTPUT_DIR=./Results/Counts

mkdir -p $OUTPUT_DIR

featureCounts -T 4 -s 0 -p -a $Ann_File -o $OUTPUT_DIR/featurecounts_0.txt $BAM_DIR/*_sorted.bam

cut -f1,7- $OUTPUT_DIR/featurecounts_0.txt > $OUTPUT_DIR/featurecounts_0_mod.txt

echo "Reading the counts is complete"

