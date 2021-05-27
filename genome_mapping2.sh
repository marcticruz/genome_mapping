samtools faidx /scratch/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_11Jun2018_KbcOz2.fasta # Execute samtools faidx to index the reference genome

bwa index /scratch/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_11Jun2018_KbcOz2.fasta # Execute bwa index to index the reference genome   

java -jar /scratch/mdacruz/bank_vole_raw_genome/picard.jar CreateSequenceDictionary R=/scratch/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_11Jun2018_KbcOz2.fasta O=/scratch/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_11Jun2018_KbcOz2.dict # Execute picard.jar CreateSequenceDictionary to create a sequence dictionary for the reference sequences

cd /work/mdacruz/bank_vole_raw_genome/1GB_100GB_bank_vole_files/ # Define directory

for name in $(cat sample_ids_unique_sorted_reduced.txt);

do

bbduk.sh in1=${name}_R1_001.fastq in2=${name}_R2_001.fastq out1=${name}_trimmed_F.fastq.gz out2=${name}_trimmed_R.fastq.gz ref=adapter_reference.fasta k=13 ktrim=r qtrim=t trimq=10 minlength=100 usejni=t -Xmx28g # Execute bbduk.sh to remove adapters from the samples 

bbduk.sh in1=${name}_trimmed_F.fastq.gz  in2=${name}_trimmed_R.fastq.gz ref=/work/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_mitochondrial_genome.fasta out1=${name}_trimmed_nomtdna_F.fastq.gz out2=${name}_trimmed_nomtdna_R.fastq.gz usejni=t -Xmx28g # Execute bbduk.sh to remove mitochondrial DNA   

bwa mem /work/mdacruz/bank_vole_raw_genome/bank_vole_reference/bank_vole_11Jun2018_KbcOz2.fasta ./${name}_trimmed_nomtdna_F.fastq.gz ./${name}_trimmed_nomtdna_R.fastq.gz -t 8  > ./${name}.sam 2> ./${name}_genome_assembly_err.txt # Execute bwa mem to map the Bank Vole sequences to a reference genome. Store the results in a sam file

rm *trimmed*

samtools view -bS ./${name}.sam  --threads 8  > ./${name}.bam 2> ${name}_bam_err.txt # Convert sam files to bam files
#
samtools flagstat ${name}.bam --threads 8 > ${name}_flagstat_out.txt 2> ${name}_flagstat_err.txt # Execute samtools flagstat to count the number of reads that were mapped to the reference genome.
# 
java -Xmx28g  -jar /work/mdacruz/bank_vole_raw_genome/picard.jar CleanSam INPUT=./${name}.bam OUTPUT=./${name}_cleaned.bam # Execute picard.jar CleaSam to mark all the reads that were not mapped to the reference
# 
java -Xmx28g -jar /work/mdacruz/bank_vole_raw_genome/picard.jar SortSam INPUT=${name}_cleaned.bam OUTPUT=${name}_cleaned_sorted.bam SORT_ORDER=coordinate # Execute SortSam to sort bam files by coordinates of the bam file
#
java -Xmx28g -jar /work/mdacruz/bank_vole_raw_genome/picard.jar AddOrReplaceReadGroups INPUT=${name}_cleaned_sorted.bam OUTPUT=${name}_readgroups_sorted_cleaned.bam RGID=1 RGLB=library1 RGPL=illumina RGPU=unit1 RGSM=${name} # Execute picard.jar AddOrReplaceReadGroups to mark heterozygous sites with a specific tag
#
java -Xmx28g -jar /work/mdacruz/bank_vole_raw_genome/picard.jar MarkDuplicates REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 M=${name}_markdups_metric_file.txt INPUT=${name}_readgroups_sorted_cleaned.bam OUTPUT=${name}_final.bam # Execute picard.jar MarkDuplicates to mark and remove duplicate reads
#
java -Xmx28g -jar /work/mdacruz/bank_vole_raw_genome/picard.jar BuildBamIndex INPUT=./${name}_final.bam # Create indexes for the bam file
#     
java -Xmx28g -jar /work/mdacruz/bank_vole_raw_genome/picard.jar CollectMultipleMetrics INPUT=${name}_final.bam OUTPUT=${name}_final.bam_metrics.txt PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution # Execute picard.jar CollectMultipleMetrics to output information about the mapping process 

mv *_final.bam ./final_bam_files/ # move
rm *.sam # delete files
rm *.bam # delete files

done
~                                                                                                                                               
