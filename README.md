# Prostate cancer detection through unbiased capture of methylated cell-free DNA

## Graphical Abstract

![image (2)](https://github.com/MassieLab/cfMBD-seq-for-Prostate-cancer-detection/assets/82373498/8413d3db-2ebd-4e13-92a7-a0ce35f83a1f)

## Abstract 

Prostate cancer screening using prostate-specific antigen (PSA) has been shown to reduce mortality but with substantial overdiagnosis and uncertainty in risk prediction, leading to unnecessary biopsies, and a greater health burden. The identification of a highly specific biomarker using liquid biopsies, represents a critical unmet need in the diagnostic pathway for prostate cancer. In this study, we employed a method that enriches for methylated cell-free DNA fragments coupled with a machine learning algorithm which enabled the detection of metastatic and localised cancers with an AUC of 0.96 and 0.74, respectively. Furthermore, we show that the methylation patterns at certain regions reflect epigenetic and transcriptomic changes at the tissue level. Notably, these differentially methylated regions are significantly enriched for biologically relevant pathways associated with the regulation of cellular proliferation and TGF-beta signalling. This demonstrates the real-world efficacy of circulating tumour DNA methylation for prostate cancer detection and prognostication.


## Data

- DMRS_ALL.csv : list of all 30k dmrs identified during bootstapping
- annotated_DMRS_900.csv : annotated list of the 900 top DMRs used for Machine Learning
- genes_all_900.csv : Beta value table for all patients at the 900 DMRs

## ML 

Python scripts for machine learning models

## R

Plot script
