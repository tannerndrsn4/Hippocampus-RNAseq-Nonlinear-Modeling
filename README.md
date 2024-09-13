# Hippocampus-RNAseq-Nonlinear-Modeling
Differential expression analysis pipeline of RNAseq data derived from banked rhesus macaque hippocampus samples. A nonlinear modeling approach allowed us to identify important midlife transitions. 

1. _Quality control_ w/ FastQC -- [FastQCbyArray.sh](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/FastQCbyArray.sh)
2. _Trimming_ w/ Trimmomatic -- [Trimmomatic.sh](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/Trimmomatic.sh)
3. _Additional quality control_ w/ FastQC -- [trimmed_fastqc.sh](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/trimmed_fastqc.sh)
4. _Mapping to Mmul_1_0 w/ STAR -- [Array_mapping.sh](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/Array_mapping.sh)
5. _Generate count matrix_ w/ featureCounts -- [featurecounts_array.sh](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/featurecounts_array.sh)


