#!/bin/sh
#SBATCH --partition=long
#SBATCH --job-name=STARmap
#SBATCH --mail-user=tander10@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/home/tander10/sternerlab/RNAseq/STARmap_%A_%a.out ## %A_%a means add "slurmJobID_ArrayID" to the end of the file name
#SBATCH --error=/home/tander10/sternerlab/RNAseq/STARmap_%A_%a.err
#SBATCH --time=10-00:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --account=sternerlab
#SBATCH --array=1-96%50

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
