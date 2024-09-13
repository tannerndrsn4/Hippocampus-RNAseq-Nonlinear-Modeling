#!/bin/sh
#SBATCH --partition=compute
#SBATCH --job-name=featurecounts
#SBATCH --mail-user=tander10@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/home/tander10/sternerlab/Wisconsin_Liver_RNAseq/featurecounts_%A_%a.out ## %A_%a means add "slurmJobID_ArrayID" to the end of the file name
#SBATCH --error=/home/tander10/sternerlab/Wisconsin_Liver_RNAseq/featurecounts_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=25G
#SBATCH --nodes=1
#SBATCH --account=sternerlab
#SBATCH --array=1-63%50

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p feature_input`

module load subread/1.6.0

featureCounts -p -t exon -g gene_id -a /home/tander10/sternerlab/Wisconsin_Liver_RNAseq/mmul10/Macaca_mulatta.Mmul_10.104.gtf -o /home/tander10/sternerlab/Wisconsin_Liver_RNAseq/STAR_output/featurecounts_output/${sampleID}_featureCount.txt /home/tander10/sternerlab/RNAseq/mapped_reads/${sampleID}_mappedAligned.sortedByCoord.out.bam
