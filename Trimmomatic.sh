#!/bin/sh
#SBATCH --partition=short
#SBATCH --job-name=Trimmomatic
#SBATCH --mail-user=tander10@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/home/tander10/sternerlab/RNAseq/Trimmomatic_%A_%a.out ## %A_%a means add "slurmJobID_ArrayID" to the end of the file name
#SBATCH --error=/home/tander10/sternerlab/RNAseq/Trimmomatic_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --account=sternerlab
#SBATCH --array=1-96%50

cd /home/tander10/sternerlab/RNAseq/RNA_fastq_files

mkdir trimmed_fastqc

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
