# Final Project Outline
## Title
Differential gene expression in TCGA within stage 1 lung cancers occurring in lower and upper lobe using DeSEQ2, controlling for race and ethnicity.
## Author
Bohan Zhang
## Overview of project
I will identify differentially expressed genes between lung cancers occuring in the upper and lower lobe of the lung. This analysis will utilize the package DeSEQ2 and follow the specific [vignette](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc DESeq2.html).
For this analysis, I'll use the TCGA cohort and have identified 215 ht-seq counts files for tumors that fit within my cohort with 81 patients with lung cancer in the lower lobe and 134 patients with lung cancer in the upper lober. Within the analysis, I will control for race and ethnicity.
## Data
I will use the data from [National Cancer Institute GDC Data Portal](https://portal.gdc.cancer.gov/repository). Examining clinical data,there are 215 tumor samples. 81 people's lung cancers occurred in lower lobes, and 134 people's lung cancers occurred in upper lobes. The specific files areavailable are [here](https://github.com/DZBohan/final_project/blob/main/clinical.tsv).
## Milestone 1
I will filter and download the HT-SEQ count files from [GDC Data Protal](https://portal.gdc.cancer.gov/repository). Since the number of files is greater than the number of cases under the same filtering condition, I need to do a second filtering of HT-SEQ files based on the case ID to get the 215 HT-SEQ count files that have cases corresponding to them. Then, I will load 215 data into vignette and create the file representing the set of samples through HT-SEQ steps.
## Milestone 2
An initial completion of vignette. I will complete an entire first draft of analysis analyzed through the vignette. After that, I will send the results of the first draft to the professor for feedback.
# Deliverable
A complete repository with clear documentation and description of your analysis and results.
