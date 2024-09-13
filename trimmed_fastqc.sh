#!/bin/sh
#SBATCH --partition=short
#SBATCH --job-name=trimmed_FastQC
#SBATCH --mail-user=tander10@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/home/tander10/sternerlab/RNAseq/FastQC_%A_%a.out ## %A_%a means add "slurmJobID_ArrayID" to the end of the file name
#SBATCH --error=/home/tander10/sternerlab/RNAseq/FastQC_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --account=sternerlab
#SBATCH --array=1-192%100

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p sampleList` ## This is how slurm knows to iterate over your sample list; each "array task ID" corresponds to the processing of one file

module load fastqc/0.11.5

fastqc -o /home/tander10/sternerlab/RNAseq/RNA_fastq_files/trimmed_fastqc/fastqc_output /home/tander10/sternerlab/RNAseq/RNA_fastq_files/trimmed_fastqc/${sampleID}*_paired.fq.gz
