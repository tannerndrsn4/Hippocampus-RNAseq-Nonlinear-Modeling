#!/bin/sh
#SBATCH --partition=short
#SBATCH --job-name=Genome_index
#SBATCH --mail-user=tander10@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/home/tander10/sternerlab/RNAseq/Genome_index.output
#SBATCH --error=/home/tander10/sternerlab/RNAseq/Genome_index.err
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=10
#SBATCH --account=sternerlab

module load easybuild  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load easybuild  ifort/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load STAR/2.5.3a
STAR --runMode genomeGenerate --genomeDir /home/tander10/sternerlab/RNAseq/Mmul_10_Directory --genomeFastaFiles /home/tander10/sternerlab/RNAseq/Mmul_10_genomic.fna --sjdbGTFfile /home/tander10/sternerlab/RNAseq/Mmul_10_genomic.gtf --runThreadN 12 --limitGenomeGenerateRAM=51900000000
