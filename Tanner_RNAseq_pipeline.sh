#!/bin/sh

#  Tanner_RNAseq_pipeline.sh
#  
#
#  Created by Tanner Anderson on 4/26/22.
#
# The following is the pipeline used for the differential expression analysis of RNAseq data from 96 banked rhesus macaque hippocampus samples

### STEP 1: Downloading the data to talapas
# The following code allows for dowlnoad of the data from PRIMeSeq to talapas directly
wget -e robots=off --wait=1 -r -nH --cut-dirs=1 --no-parent --http-user=tander10@uoregon.edu --http-password=3YnNgdxbxhVnSET  --reject "index.html*"  'https://prime-seq.ohsu.edu/_webdav/Internal/BioinformaticsCore/SteveKohama/1/%40files/raw/RNA220216SK/220311_A01058_0220_BHMGWTDSX3/RNA220216SK/?listing=html'
#The download is also a little messy and requires some reorganizing/renaming of files. This can be done with basic bash commands but the following code is useful for getting rid of the 'listing=html' extension on the files
for f in *\?listing=html; do mv -v "${f}" "${f%%?listing=html}"; done

### STEP 2: Performing quality checks and trimming poor quality bases
# Next we need to perform some quality checks to assess the integrity of our samples using FastQC
#The following wesbite is useful when learning to interpret FastQC and trimmomatic outputs:https://datacarpentry.org/wrangling-genomics/02-quality-control/index.html
# First need to make a list of sample names without the fastq.gz extension using the cut command. Importantly, the sample list file needs to be in the same directory where you submit the array job. 
find *_001.fastq.gz | cut -d _ -f 1-4 > sampleList
# Then we can begin the array job
sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p sampleList` ## This is how slurm knows to iterate over your sample list; each "array task ID" corresponds to the processing of one file
module load fastqc/0.11.5

fastqc -o /home/tander10/sternerlab/RNAseq/RNA_fastq_files/FastQC /home/tander10/sternerlab/RNAseq/RNA_fastq_files/${sampleID}*.fastq*
# This will output a directory filled with FastQC reports for each sample. We need to take the summary reports for each sample and concatenate them so we know what to trim out later on. Note: The code above refers to an array job being performed so all of our samples run through fastqc at once (saves time!!) [FastQCbyArray.sh]
# Now we need to unzip all of the fastqc outputs to have access to the summary.txt files [Unzip_FastQC.sh]
for filename in /home/tander10/sternerlab/RNAseq/RNA_fastq_files/FastQC/*.zip
do
    unzip $filename
done
# Now we need to contatenate all of the summary.txt files into one massive summary file (This command can be done in the same directory where fastqc folders are located)
cat */summary.txt > fastqc_summaries.txt
# We can also make an aggregate qc report using multiqc. First, naviagte to the directory where all the fastqc reports are located. Then run the following script:
module load easybuild
module load ifort/2017.1.132-GCC-6.3.0-2.27
module load impi/2017.1.132
module load MultiQC/1.3-Python-3.6.1
multiqc . --verbose
# Now we will need to start trimming some poor quality bases to ensure the inegrity of the samples we are using for mapping
# First, we will make a txt file that contains example sample names like how we did earlier but excluding R1 or R2 so we can specify that better in trimmomatic.
find *R1_001.fastq.gz | cut -d _ -f 1-3 > Trimm_input
# Next, we can submit our array job. The following code was included in the job we submitted:
module load racs-eb
module load Trimmomatic/0.36-Java-1.8.0_162

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p Trimm_input`

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE \
-threads 8 \
-phred33 \
${sampleID}_R1_001.fastq.gz \
${sampleID}_R2_001.fastq.gz \
${sampleID}_F_paired.fq.gz \
${sampleID}_F_unpaired.fq.gz \
${sampleID}_R_paired.fq.gz \
${sampleID}_R_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36 \

# Once our trimming is complete, we need to do some reorganizing/decluttering
mkdir trimmed_fastqc

mv *_F_paired.fq.gz trimmed_fastq
mv *_F_unpaired.fq.gz trimmed_fastq
mv *_R_paired.fq.gz trimmed_fastq
mv *_R_unpaired.fq.gz trimmed_fastq

# Now that we have trimmed our fastq files, we will perform another round of qc to make sure that the trimming did what it was supposed to do.
# First we will make a sample list of our trimmed files like before.
find *_paired.fq.gz | cut -d _ -f 1-4 > sampleList
# Then we will submit our fastqc array job of our trimmed reads
sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p sampleList` ## This is how slurm knows to iterate over your sample list; each "array task ID" corresponds to the processing of one file
mkdir fastqc_output
module load fastqc/0.11.5

fastqc -o /home/tander10/sternerlab/RNAseq/RNA_fastq_files/trimmed_fastq/fastqc_output /home/tander10/sternerlab/RNAseq/RNA_fastq_files/trimmed_fastqc/${sampleID}*_paired.fq.gz

# Once we have our fastqc output, we will again make an aggregate report by performing multiqc
module load easybuild
module load ifort/2017.1.132-GCC-6.3.0-2.27
module load impi/2017.1.132
module load MultiQC/1.3-Python-3.6.1
multiqc . --verbose
# If the multiqc report looks good, then we are free to move on to the mapping step

### STEP 3: Mapping our trimmed reads to the reference genome
# Mapping involves two separate steps: generating a genome index & mapping reads to the genome
#First we will generate the genome index using the rhesus macaque reference genome Mmul_10
module load easybuild  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load easybuild  ifort/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load STAR/2.5.3a
STAR --runMode genomeGenerate --genomeDir /home/tander10/sternerlab/RNAseq/Mmul_10_Directory --genomeFastaFiles /home/tander10/sternerlab/RNAseq/Mmul_10_genomic.fna --sjdbGTFfile /home/tander10/sternerlab/RNAseq/Mmul_10_genomic.gtf --runThreadN 12 --limitGenomeGenerateRAM=51900000000
# Once the index has been generated, we can then begin our array job to map our trimmed reads to the reference genome
# First, we are going to make a list of our sample names that will be included in the array job
find *F_paired.fq.gz | cut -d _ -f 1-3 > STAR_input
# You will need to make sure that this list is in the same directory where you submit the job
# Next, we can submit our array job [Array_mapping.sh]
sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p STAR_input`

module load easybuild  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load easybuild  ifort/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load STAR/2.5.3a

STAR --runMode alignReads \
--runThreadN 12 \
--genomeDir /home/tander10/sternerlab/RNAseq/Mmul_10_Directory \
--readFilesCommand gunzip -c \
--readFilesIn /home/tander10/sternerlab/RNAseq/RNA_fastq_files/trimmed_fastq/${sampleID}_F_paired.fq.gz /home/tander10/sternerlab/RNAseq/RNA_fastq_files/trimmed_fastq/${sampleID}_R_paired.fq.gz \
--twopassMode Basic \
--outFilterMismatchNoverLmax 0.05 \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /home/tander10/sternerlab/RNAseq/mapped_reads/${sampleID}_mapped
#This will output bam files that are aligned and sorted by coordinate

### STEP 4: Making a count matrix
# We will use featurecounts to quantify read counts and produce a count matrix
# First, we need to clean up our GTF file and remove any blank gene_id identifiers
grep 'gene_id ""' Mmul_10_genomic.gtf
grep -v 'gene_id ""' Mmul_10_genomic.gtf > Mmul_10_genomic_fixed.gtf
# Next, we can submit our featurecounts job. featurecounts is an R program so we will need to submit it a special way in talapas.
# The following is the Rscript that we will need to have prewritten so our job can run it:

#!/usr/bin/env Rscript

library(Rsubread)

sortedBAM.files <- list.files(path = "/home/tander10/sternerlab/RNAseq/mapped_reads", pattern = "mappedAligned.sortedByCoord.out.bam$", full.names = TRUE)

featureCounts(files=sortedBAM.files, annot.ext="/home/tander10/sternerlab/RNAseq/Mmul_10_genomic_fixed.gtf", isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id", isPairedEnd=TRUE)

# Once that script is in the proper directory, then we can submit our job:
module load gcc
module load R/4.0.2

Rscript runfeaturecounts.R

# Running an R script on talapas for a high number of samples take a very long time (even with threading), therefore we can run an array job in bash for featurecounts and then concatenate all of the outputs at the end into one large count matrix. We start with making a list of our sample names
find *_mappedAligned.sortedByCoord.out.bam | cut -d _ -f 1-3 > feature_input
# Once we have a list of our sample names, we can submit our array job
sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p feature_input`

module load subread/1.6.0

featureCounts -p -t exon -g gene_id -a /home/tander10/sternerlab/RNAseq/Mmul_10_genomic_fixed.gtf -o /home/tander10/sternerlab/RNAseq/mapped_reads/featurecounts_output/${sampleID}_featureCount.txt /home/tander10/sternerlab/RNAseq/mapped_reads/${sampleID}_mappedAligned.sortedByCoord.out.bam
# Once this has completed, we can concatenate all of the ouput txt files into one large count matrix
module load easybuild  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load easybuild  ifort/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load parallel/20160622

ls -1  *featureCount.txt | parallel 'cat {} | sed '1d' | cut -f7 {} > {/.}_clean.txt'
ls -1  *featureCount.txt | head -1 | xargs cut -f1 > genes.txt
paste genes.txt *featureCount_clean.txt > output.txt

# Now that we have our final count matrix, we can then move the matrix into an R environment to perform the differential expression analysis
