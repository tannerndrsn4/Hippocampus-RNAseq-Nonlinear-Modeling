# Hippocampus-RNAseq-Nonlinear-Modeling
Differential expression analysis pipeline of RNAseq data derived from banked rhesus macaque hippocampus samples. A nonlinear modeling approach allowed us to identify important midlife transitions. 

Steps **1-5** were completed using shell scripts.

1. **Quality control** w/ FastQC -- [FastQCbyArray.sh](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/FastQCbyArray.sh)
2. **Trimming** w/ Trimmomatic -- [Trimmomatic.sh](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/Trimmomatic.sh)
3. **Additional quality control** w/ FastQC -- [trimmed_fastqc.sh](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/trimmed_fastqc.sh)
4. **Mapping to Mmul_10** w/ STAR -- [Array_mapping.sh](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/Array_mapping.sh)
5. **Generate count matrix** w/ featureCounts -- [featurecounts_array.sh](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/featurecounts_array.sh)

Steps **6-10** were completed using RStudio.

6. **Read count normalization** w/ edgeR & limma -- [read_count_normalization](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/read_count_normalization) 
7. **PCA analysis & Sample information** figure generation -- [PCA_sampleinfo](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/PCA_sampleinfo)
8. **Nonlinear modeling** differential expression analyses using ARIMA approach -- [nonlinear_modeling](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/nonlinear_modeling)
9. **Gene dysregulation analysis** -- [dysregulation](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/dysregulation)
10. **Gene ontology enrichment analysis** w/ enrichR -- [gene_ontology](https://github.com/tannerndrsn4/Hippocampus-RNAseq-Nonlinear-Modeling/blob/main/gene_ontology)

All analyses were completed on UO's computer cluster Talapas.
