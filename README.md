![Logo-r](https://github.com/Goodarzilab/Ribolog/blob/master/vignettes/Logo3.png)

# Ribolog
A suite of regression-based tools for Ribosome profiling data analysis

## Installing Ribolog

### Install directory in R

```R
install.packages('BiocManager')
BiocManager::install("Goodarzilab/Ribolog")
```

### Install Ribolog with a Conda Environment in Bash
Run the code below in your terminal to build and activate a conda environment with all the dependencies of Ribolog inside it.

```sh
conda env create -f 'https://raw.githubusercontent.com/goodarzilab/Ribolog/master/environment.yml'
conda activate ribolog # Now you're inside the conda environment
R -e "BiocManager::install('Goodarzilab/Ribolog', dependencies = FALSE)"
```

## Tutorial 

### Please refer to the ["Ribolog_Full_Tutorial.ipynb"](Ribolog_Full_Tutorial.ipynb) file for a full tutorial of the package.

## Module details

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
rates among biological samples. Although better results are always obtained with sufficient replicates, __Ribolog__ is able to peform the TER test with only one replicate per sample. The TER test is not restricted to pairwise comparisons; any number of samples described by a list of attributes (covariates) can be compared in a single model.

## Module 5: Empirical significance testing and Meta-analysis
Offers tools for two important slightly advanced statistical tasks: 1) Empirical null hypothesis testing to reduce false positives in replicated datasets.
2) Meta-analysis to integrate
results of biologically related and statistically correlated experiments.

The Ribolog workflow is described in great detail in the package vignettes (RIBOLOG.pdf in the vignettes folder).

![Logo-r](https://github.com/Goodarzilab/Ribolog/blob/master/vignettes/Ribolog_workflow.v5.png)

Rendering the vignettes during installation requires bam files that are not uploaded onto this repository. The knitted .pdf file should be downloaded directly from the vignettes folder instead.

__Ribolog__ was developed by Hossein Asgharian and Sohit Miglani at UCSF supervised by Hani Goodarzi and Adam Olshen. More modules are being prepared and will be released in near future.

For questions and comments, email us at:  
- hossein.asgharian@gmail.com
- sohitmiglani@gmail.com
- hani.goodarzi@ucsf.edu
