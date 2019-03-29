# Ribolog
A suite of regression-based tools for RNA-seq and Ribo-seq data analysis

The core function of this package is a test for tranlational efficiency ratio (TER) based on 
ribosome profiling data (module 4: TER_BASE). Module 1: CELP identifies positions of translational 
pause (stalling) and corrects RPF counts to eliminate the stalling bias. Module 2: PREP normalizes 
and combines RNA and RPF datasets and shapes them into a format ready for QC and anlaysis. Module 3:
QC_PCA performs princiapl component analysis on RNA, RPF and TE (translational efficiency) to visualize
reproducibility among replicates and help with hypothesis generation with respect to biological effects.
Modules 5 and 6 are intended for more advanced users. Module 5: TER_ALT provides alternative options 
for the TER test including empirical significance testing and meta-analysis. Module 6: MODEL_SEL 
provides tools for model selection and feature selection to identify factors that influence 
translation rates of individual genes.
