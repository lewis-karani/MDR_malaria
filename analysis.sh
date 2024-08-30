#!/bin/bash

# Define paths and variables
REF_FASTA="/mnt/d/Design_beds/MDR-MALARIA/Plasmodium_falciparum_reference.fasta"
PRIMER_BED="/mnt/d/Design_beds/MDR-MALARIA/Plasmodium_falciparum_bed.bed"
READS_DIR="/mnt/d/Design_beds/MDR-MALARIA/"
OUTPUT_DIR="/mnt/d/Design_beds/MDR-MALARIA/results/"  

# Make sure tools are installed and in PATH

# Index the reference genome
bwa index $REF_FASTA
samtools faidx $REF_FASTA

# Step 1: Loop through each sample and process the FASTQ files
for SAMPLE in KKGR03 KKGR06 KKGR12 KKGR15 KKGR25 KKGR26 KKGR27 Reference; do

  # Define paths to FASTQ files
  READS_FASTQ_R1="$READS_DIR/${SAMPLE}.fastq"
  TRIMMED_R1="$OUTPUT_DIR/${SAMPLE}_trimmed.fastq"
  OUTPUT_PREFIX="$OUTPUT_DIR/$SAMPLE"

 # Step 2: Trimming using Trimmomatic SE
  trimmomatic SE -threads 8 \
      $READS_FASTQ_R1 $TRIMMED_R1 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

  # Step 3: Quality Control 
  fastqc $TRIMMED_R1 -o $OUTPUT_DIR

  # Step 4: Align the reads to the reference genome
  bwa mem -t 8 $REF_FASTA $TRIMMED_R1 > ${OUTPUT_PREFIX}.sam

  # Step 5: Convert SAM to BAM, sort, and index the BAM file
  samtools view -Sb ${OUTPUT_PREFIX}.sam > ${OUTPUT_PREFIX}.bam
  samtools sort ${OUTPUT_PREFIX}.bam -o ${OUTPUT_PREFIX}_sorted.bam
  samtools index ${OUTPUT_PREFIX}_sorted.bam

  # Step 6: Call variants focusing on the specified regions using GATK HaplotypeCaller
  gatk HaplotypeCaller \
      -R $REF_FASTA \
      -I ${OUTPUT_PREFIX}_sorted.bam \
      --intervals $PRIMER_BED \
      -O ${OUTPUT_PREFIX}_raw_variants.vcf

  # Step 7: Filter the variants
  gatk VariantFiltration \
      -R $REF_FASTA \
      -V ${OUTPUT_PREFIX}_raw_variants.vcf \
      --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
      --filter-name "FILTER" \
      -O ${OUTPUT_PREFIX}_filtered_variants.vcf

   # Step 8: Generate a BAM file focused on the regions of interest for visualization
  samtools view -b -L $PRIMER_BED ${OUTPUT_PREFIX}_sorted.bam > ${OUTPUT_PREFIX}_targeted.bam
  samtools sort ${OUTPUT_PREFIX}_targeted.bam -o ${OUTPUT_PREFIX}_targeted_sorted.bam
  samtools index ${OUTPUT_PREFIX}_targeted_sorted.bam

  # Step 9: Generate a VCF file with only the variants within the targeted regions
  gatk SelectVariants \
      -R $REF_FASTA \
      -V ${OUTPUT_PREFIX}_filtered_variants.vcf \
      --intervals $PRIMER_BED \
      -O ${OUTPUT_PREFIX}_targeted_variants.vcf

  # Cleanup intermediate files if needed
  rm ${OUTPUT_PREFIX}.sam ${OUTPUT_PREFIX}.bam ${OUTPUT_PREFIX}_targeted.bam

done

# Step 10: IGV Preparation - Copy files to the IGV data directory
IGV_DIR="/mnt/d/Design_beds/MDR-MALARIA/igv/"
mkdir -p $IGV_DIR
cp $OUTPUT_DIR/*_targeted_sorted.bam $IGV_DIR/
cp $OUTPUT_DIR/*_targeted_sorted.bam.bai $IGV_DIR/
cp $OUTPUT_DIR/*_targeted_variants.vcf $IGV_DIR/

echo "Alignment and variant files prepared for IGV visualization."
echo "You can now open IGV and load the BAM files and VCF files for visualization."
