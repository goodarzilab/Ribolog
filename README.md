![Logo-r](https://github.com/Goodarzilab/Ribolog/blob/master/vignettes/Logo3.png)

# Ribolog
A suite of regression-based tools for Ribosome profiling data analysis

## Installing Ribolog

### Install directly in R

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

## Module 1: Pre-processing
Converts aligned, sorted, and indexed BAM files into per-read P-site assignments: read-length distributions,
P-site offset detection, and reading-frame/periodicity QC. Produces the `reads_psite_list` object used by
later modules.

## Module 2: PREP and CELP
Corrects RPF counts for codon-level stalling bias (CELP: Consistent Excess of Loess Preds) and analyzes
ribosome dwell times per codon/amino acid, then normalizes and filters RNA/RPF read counts into the
transcript-by-sample count matrices used by Modules 3-6.

## Module 3: QC
Three tools for assessing replicate reproducibility and biological signal before testing: PCA of
translational efficiency, the proportion of null (non-differentially-translated) features, and
correlograms of equivalent replicate-vs-replicate TER tests.

## Module 4: TER
Tests the size and significance of differential translational efficiency between biological samples via
logistic regression. Although better results are always obtained with sufficient replicates, __Ribolog__
is able to perform the TER test with only one replicate per sample. The TER test is not restricted to
pairwise comparisons; any number of samples described by several attributes (covariates) can be compared
in a single model.

## Module 5: Empirical Null Testing and Meta-analysis
Two advanced statistical tools for TER testing: 1) empirical null hypothesis testing, which derives the
null distribution directly from replicate-vs-replicate comparisons in your own data instead of assuming a
theoretical one, reducing false positives from batch effects, mapping noise, or overdispersion. 2)
meta-analysis, which combines correlated test results (e.g. from separate datasets, or rep-by-rep
sub-tests) into a single consensus effect size and p-value.

## Module 6: ORF Usage and Stop Codon Readthrough
Tests differential usage of upstream ORFs and stop-codon readthrough by comparing the distribution of
P-site reads across the 5'UTR, CDS, and 3'UTR regions of each transcript between biological samples.

The Ribolog workflow is described in great detail in the package vignettes (RIBOLOG.pdf in the vignettes folder).

![Ribolog workflow diagram](https://github.com/Goodarzilab/Ribolog/blob/master/vignettes/Ribolog_workflow.svg)

Rendering the vignettes during installation requires bam files that are not uploaded onto this repository. The knitted .pdf file should be downloaded directly from the vignettes folder instead.

__Ribolog__ was developed by Hossein Asgharian and Sohit Miglani at UCSF supervised by Hani Goodarzi and Adam Olshen. More modules are being prepared and will be released in near future.

For questions and comments, email us at:  
- hossein.asgharian@gmail.com
- sohitmiglani@gmail.com
- hani.goodarzi@ucsf.edu
