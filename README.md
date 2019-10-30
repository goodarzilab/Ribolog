![Logo-r](https://github.com/Goodarzilab/Ribolog/blob/master/vignettes/Logo-r.jpg)

# Ribolog
A suite of regression-based tools for Ribosome profiling data analysis

## Module 1: CELP (Consistent Excess of Loess Preds) 
Identifies positions of translational pause (stalling) 
and corrects RPF counts to eliminate the impact of stalling bias. The output of CELP can be used to model the
factors that influence translational dynamics. 

## Module 2: PREP 
Normalizes and combines RNA and RPF datasets and shapes them into a format ready for quality control (QC) and translational
efficiency ratio (TER) anlaysis. 

## Module 3: QC 
Includes three powerful tools to quantify and visualize reproducibility among replicates and inform hypothesis generation with respect to biological effects: 
princiapl component analysis (PCA) of TEs, proportion of null features (non-differentially translated transcripts)
and correlation of equivalent TER tests. 

## Module 4: TER 
Tests the size and significance of differential translation
rates among biological samples. Although better results are always obtained with sufficient replicates, Ribolog is able to peform the TER test with only one replicate per sample. The TER test is not restricted to pairwise comparisons; any number of samples described by several indenpedent variables can be compared in a single model. 

The Ribolog workflow is described in great detail in the package vignettes (RIBOLOG.pdf in the vignettes folder). 

## Installing Ribolog
Run the following command in R:

`install_github("Goodarzilab/Ribolog", dependencies = TRUE, build_vignettes = FALSE, build_manual = FALSE)`

Rendering the vignettes during installation requires bam files that are not uploaded onto this repository. The knitted .pdf file should be downloaded directly from the vignettes folder instead.

Ribolog is still a work in progress. Four other modules are being prepared and will be released in near future.
Feel free to contact us (using the issues tab on this page) about problems with current functions and modules or if you you have special requests or ideas about additional analyses you would like to do with ribosome profiling data. 
