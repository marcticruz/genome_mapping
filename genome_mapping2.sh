#!/bin/sh
#SBATCH --time=45:00:00
#SBATCH --mem=20G
#SBATCH --job-name=genome_mapping
#SBATCH --output=/scratch/mdacruz/bank_vole_raw_genome/1GB_100GB_bank_vole_files/genome_mapping_reduced.out
#SBATCH --error=/scratch/mdacruz/bank_vole_raw_genome/1GB_100GB_bank_vole_files/genome_mapping_reduced.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marcos.cruz@ou.edu
#SBATCH --ntasks=1 
#SBATCH --partition=normal
 
## genome mapping workflow

module load BBMap/38.36-intel-2016a # Load BBMap program

module load BWA/0.7.13-intel-2016a  # Load BWA program

module load SAMtools/1.9-foss-2018b # Load SAMtools program

module load GATK/3.8-0-Java-1.8.0_141  # Load GATK program 

module load VCFtools/0.1.14-intel-2016a # Load VCFtools program

Java/1.7.0_21 # Load java program  

samtools faidx /scratch/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_11Jun2018_KbcOz2.fasta # Execute samtools faidx to index the reference genome

bwa index /scratch/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_11Jun2018_KbcOz2.fasta # Execute bwa index to index the reference genome   

java -jar picard.jar CreateSequenceDictionary R=/scratch/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_11Jun2018_KbcOz2.fasta O=/scratch/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_11Jun2018_KbcOz2.dict # Execute picard.jar CreateSequenceDictionary to create a sequence dictionary for the reference sequences

cd /scratch/mdacruz/bank_vole_raw_genome/1GB_100GB_bank_vole_files # Define directory

for name in $(cat sample_ids_reduced.txt);

do
	bbduk.sh in=${name}_R1_001.fastq in2=${name}_R2_001.fastq out=${name}_trimmed_F.fastq.gz out2=${name}_trimmed_R.fastq.gz ref=adapter_reference.fasta k=13 ktrim=r qtrim=t trimq=10 minlength=100 # Execute bbduk.sh to remove adapters from the samples 
 
	bbduk.sh in1=${name}_trimmed_F.fastq.gz  in2=${name}_trimmed_R.fastq.gz ref=/scratch/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_mitochondrial_genome.fasta out1=${name}_trimmed_nomtdna_F.fastq.gz out2=${name}_trimmed_nomtdna_R.fastq.gz # Execute bbduk.sh to remove mitochondrial DNA 

	bwa mem -t 10 /scratch/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_11Jun2018_KbcOz2.fasta ./${name}_trimmed_nomtdna_F.fastq.gz ./${name}_trimmed_nomtdna_R.fastq.gz  > ./${name}.sam 2> ./${name}_genome_assembly_err.txt # Execute bwa mem to map the Bank Vole sequences to a reference genome. Store the results in a sam file

	samtools view -bS ./${name}.sam > ./${name}.bam 2> ${name}_bam_err.txt # Convert sam files to bam files

	samtools flagstat ${name}.bam > ${name}_flagstat_out.txt 2> ${name}_flagstat_err.txt # Execute samtools flagstat to count the number of reads that were mapped to the reference genome.
 
	java -Xmx2g -jar /scratch/mdacruz/bank_vole_raw_genome/picard.jar CleanSam INPUT=./${name}.bam OUTPUT=./${name}_cleaned.bam # Execute picard.jar CleaSam to filter all the reads that were not mapped to the reference
 
	java -Xmx2g -jar /scratch/mdacruz/bank_vole_raw_genome/picard.jar SortSam INPUT=${name}_cleaned.bam OUTPUT=${name}_cleaned_sorted.bam SORT_ORDER=coordinate # Execute SortSam to sort bam files by coordinates of the bam file

	java -Xmx2g -jar /scratch/mdacruz/bank_vole_raw_genome/picard.jar AddOrReplaceReadGroups INPUT=${name}_cleaned_sorted.bam OUTPUT=${name}_readgroups_sorted_cleaned.bam RGID=1 RGLB=library1 RGPL=illumina RGPU=unit1 RGSM=${name} # Execute picard.jar AddOrReplaceReadGroups to mark heterozygous sites with a specific tag

	java -Xmx2g -jar /scratch/mdacruz/bank_vole_raw_genome/picard.jar MarkDuplicates REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 M=${name}_markdups_metric_file.txt INPUT=${name}_readgroups_sorted_cleaned.bam OUTPUT=${name}_final.bam # Execute picard.jar MarkDuplicates to mark and remove duplicate reads

	java -Xmx2g -jar /scratch/mdacruz/bank_vole_raw_genome/picard.jar BuildBamIndex INPUT=./${name}_final.bam # Create indexes for the bam file
     
    java -Xmx2g -jar /scratch/mdacruz/bank_vole_raw_genome/picard.jar CollectMultipleMetrics INPUT=${name}_final.bam OUTPUT=${name}_final.bam_metrics.txt PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution # Execute picard.jar CollectMultipleMetrics to output information about the reading process 

done

cat *_final.bam > bank_vole_bam_files.bam  # Redirect the content of all the bam files to a single bam file

java -Xmx45g -jar GenomeAnalysisTK.jar  -T HaplotypeCaller -nct 6  -R /scratch/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_11Jun2018_KbcOz2.fasta -I ./bank_vole_bam_files.bam -ERC GVCF -o ./bank_vole_raw.g.vcf # Execute GenomeAnalysisTK.jar -T HaplotypeCaller to identify SNP
 
java -Xmx45g -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R /scratch/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_11Jun2018_KbcOz2.fasta -V ./bank_vole_raw.g.vcf -o ./bank_vole.vcf # Execute GenomeAnalysisTK.jar -T GenotypeGVCFs to genotype the samples

vcftools --vcf ./bank_vole.vcf --min-alleles 2 --max-alleles 2 --remove-indels --max-missing 1 --maf 0.15 --max-maf 0.83 --min-meanDP 7 --max-meanDP 25 --thin 2000 --recode -c  > ./bank_vole_snps.vcf # Execute VCFtools to filter out low-quality and uninformative SNP
