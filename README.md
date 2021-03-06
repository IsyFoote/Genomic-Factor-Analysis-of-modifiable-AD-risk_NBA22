Code walkthrough for genomic factor analysis of modifiable risk for Alzheimer’s disease
================
Isabelle Foote
2022-06-13

- [Overview](#overview)
- [GWAS data identification and formatting](#gwas-data-identification-and-formatting)
  * [Criteria for inclusion](#criteria-for-inclusion)
    + [Links to the GWAS data used in our current study](#links-to-the-gwas-data-used-in-our-current-study)
  * [Example code for formatting UK Biobank GWAS data](#example-code-for-formatting-uk-biobank-gwas-data)
    + [Set your working directory and load the required packages](#set-your-working-directory-and-load-the-required-packages)
    + [Make an edited version of the variant data to merge with the GWAS datasets](#make-an-edited-version-of-the-variant-data-to-merge-with-the-gwas-datasets)
    + [Format UK Biobank GWAS data ready for munging](#format-uk-biobank-gwas-data-ready-for-munging)
- [LD score regression](#ld-score-regression)
    + [Install LDSC software](#install-ldsc-software)
  * [Munge data](#munge-data)
  * [Calculating univariate SNP-based heritability](#calculating-univariate-snp-based-heritability)
  * [Calculating pairwise genetic correlations](#calculating-pairwise-genetic-correlations)
  * [Plotting the data](#plotting-the-data)
    + [Genetic correlation matrix](#genetic-correlation-matrix)
    + [Undirected weighted graph](#undirected-weighted-graph)
- [Genomic SEM](#genomic-sem)
    + [Install *GenomicSEM*](#install--genomicsem-)
    + [Set working directory and load the *GenomicSEM* package](#set-working-directory-and-load-the--genomicsem--package)
  * [Munge data](#munge-data-1)
  * [Exploratory factor analysis performed in the odd autosomes](#exploratory-factor-analysis-performed-in-the-odd-autosomes)
    + [Load the additional required packages](#load-the-additional-required-packages)
    + [Conduct multivariable LD score regression (odd autosomes)](#conduct-multivariable-ld-score-regression--odd-autosomes-)
    + [Exploratory factor analysis](#exploratory-factor-analysis)
    + [EFA loadings plot](#efa-loadings-plot)
  * [Confirmatory factor analysis performed in the even autosomes](#confirmatory-factor-analysis-performed-in-the-even-autosomes)
    + [Conduct multivariable LD score regression (even autosomes)](#conduct-multivariable-ld-score-regression--even-autosomes-)
    + [Confirmatory factor analysis](#confirmatory-factor-analysis)
    + [Plotting CFA models](#plotting-cfa-models)
  * [Follow-up genome-wide analyses](#follow-up-genome-wide-analyses)
    + [Genome-wide multivariable LD score regression](#genome-wide-multivariable-ld-score-regression)
    + [Genome-wide CFA](#genome-wide-cfa)
    + [Additional SEM modelling](#additional-sem-modelling)
    + [Common factor model](#common-factor-model)
    + [Hierarchal model](#hierarchal-model)
    + [Bifactor model](#bifactor-model)
  * [Models without the APOE region](#models-without-the-apoe-region)
- [Sensitivity analyses](#sensitivity-analyses)
    + [LD score regression in different chromosomal groupings](#ld-score-regression-in-different-chromosomal-groupings)
  * [EFA comparison across chromosomal groupings](#efa-comparison-across-chromosomal-groupings)
  * [CFA comparison across chromosomal groupings](#cfa-comparison-across-chromosomal-groupings)
- [Final remarks](#final-remarks)

# Overview

The following document outlines the code that was used to carry out the
work described in the study titled **‘The shared genetic architecture of
modifiable risk for Alzheimer’s disease: a genomic structural equation
modelling study’** by Isabelle Foote, Ben Jacobs, Georgina Mathlin,
Cameron Watson, Phazha Bothongo, Sheena Waters, Ruth Dobson, Alastair
Noyce, Kam Bhui, Ania Korszun and Charles Marshall.

The published version of the paper can be found
[here](https://www.sciencedirect.com/science/article/pii/S0197458022001312).

The preprint for the study can be found on
[MedRxiv](https://www.medrxiv.org/content/10.1101/2021.02.23.21252211v1).

This markdown report and code was written by Isabelle Foote.

The current document provides step-by-step instructions to replicate our
work and utilises open-source software and code that was created by the
original developers of the methods used.

Further details relating to software/code can be found on the
corresponding original GitHub pages:

-   **[LD score regression](https://github.com/bulik/ldsc)** (Brendan
    Bulik-Sullivan and colleagues)
-   **[Genomic SEM](https://github.com/GenomicSEM/GenomicSEM)** (Andrew
    Grotzinger and colleagues)

Further details on the methods used can be found in the following
publications, upon which the current study protocol was based:

**LD score regression**

-   Bulik-Sullivan B., et al. (2015) *LD Score regression distinguishes
    confounding from polygenicity in genome-wide association studies*.
    Nature Genetics. 47: p. 291-295 <https://doi.org/10.1038/ng.3211>.
-   Bulik-Sullivan B., et al. (2015) *An atlas of genetic correlations
    across human diseases and traits*. Nature Genetics. 47(11):
    p. 1236-1241 <https://doi.org/10.1038/ng.3406>.
-   Zheng J., et al. (2016) *LD Hub: a centralized database and web
    interface to perform LD score regression that maximizes the
    potential of summary level GWAS data for SNP heritability and
    genetic correlation analysis*. Bioinformatics. 33(2): p. 272-279
    <https://doi.org/10.1093/bioinformatics/btw613>.

**Genomic SEM**

-   Grotzinger AD., et al. (2019) *Genomic structural equation modelling
    provides insights into the multivariate genetic architecture of
    complex traits*. Nat Hum Behav. 3(5): p. 513-525
    <https://doi.org/10.1038/s41562-019-0566-x>.

# GWAS data identification and formatting

The LD score regression and genomic SEM methods we use in this study are
based off GWAS summary statistics data, so the first step is to identify
the appropriate data for your analysis. For a detailed rationale
relating to the choice of traits we included within our study please
refer to our
[paper](https://www.sciencedirect.com/science/article/pii/S0197458022001312).

In short, we included GWAS data for Alzheimer’s disease and it’s major
potentially modifiable risk factors (less education, hearing loss,
hypertension, high alcohol intake, obesity, smoking, depression, social
isolation, physical inactivity, type 2 diabetes, sleep disturbance and
socioeconomic deprivation).

Some useful resources to find datasets are:

-   [GWAS Catalog](https://www.ebi.ac.uk/gwas/)
-   [GWAS Atlas](https://atlas.ctglab.nl/)
-   [Neale lab](http://www.nealelab.is/uk-biobank) (for GWAS of UK
    Biobank traits)
-   The [UK Biobank Data
    Showcase](https://biobank.ndph.ox.ac.uk/showcase/) is an easy way to
    look up what variables were measured within the UK Biobank study
-   Relevant consortia for your topic of interest e.g. the [Psychiatric
    Genomics Consortium](https://www.med.unc.edu/pgc/download-results/)
    for GWASs of psychiatric disorders

## Criteria for inclusion

In order for your analysis to be reliable you need to ensure that the
GWAS data you include is appropriate and meets the advised criteria
explained in the original methods papers cited above and on the relevant
wikis.

For the current study we ensured that the data met the following
criteria:

-   The samples were of unrelated individuals of European ancestry
-   Each dataset had a sample size \>5000 (ideally over 10,000)
-   The summary statistics had not been adjusted for heritable
    covariates during association testing
-   The association test had not used linear mixed modelling (LMM)
    methods
-   Each trait had a SNP-based heritability Z score \>4

Often the SNP-based heritability is reported in the orginal GWAS paper
so to get an idea about whether your trait is appropriate this can be a
good starting point. [LD Hub](http://ldsc.broadinstitute.org/) also
contains heritability and genetic correlation results for hundreds of
traits. For the Neale lab UK Biobank GWASs you can check details of
their heritability in their [heritability
browser](https://nealelab.github.io/UKBB_ldsc/index.html) and pairwise
genetic correlations in their [correlation
browser](https://ukbb-rg.hail.is/).

This is useful as an initial checking step, however we recommend
performing your own LD score regression as an initial part of your own
analysis pipeline to ensure that the results you get are in line with
previously run results and to facilitate a deeper understanding of your
own study results.

### Links to the GWAS data used in our current study

The GWAS data we used in our current study is publicly available and
available to download via the following links:

-   **Alzheimer’s disease**
    -   Paper: Kunkle BW., et al. (2019) *Genetic meta-analysis of
        diagnosed Alzheimer’s disease identifies new risk loci and
        implicates Aβ, tau, immunity and lipid processing*. Nat
        Genetics. 51: p. 414-430
        <https://doi.org/10.1038/s41588-019-0358-2>.
    -   Data: [IGAP](https://www.niagads.org/datasets/ng00075)
-   **Major depressive disorder**
    -   Paper: Wray NR., et al. (2018) *Genome-wide association analyses
        identify 44 risk variants and refine the genetic architecture of
        major depression*. Nat Genetics. 50(5): p. 668-681
        <https://doi.org/10.1038/s41588-018-0090-3>.
    -   Data: [PGC](https://www.med.unc.edu/pgc/download-results/mdd/)
-   **Insomnia**
    -   Paper: Jansen PR., et al. (2019) *Genome-wide analysis of
        insomnia in 1,331,010 individuals identifies new risk loci and
        functional pathways*. Nat Genetics. 51(3): p. 394-403
        <https://doi.org/10.1038/s41588-018-0333-3>.
    -   Data: [CTG Lab](https://ctg.cncr.nl/software/summary_statistics)
-   **Type 2 diabetes**
    -   Paper: Scott RA., et al. (2017) *An Expanded Genome-Wide
        Association Study of Type 2 Diabetes in Europeans*. Diabetes.
        66(11): p. 2888 <https://doi.org/10.2337/db16-1253>.
    -   Data: [DIAGRAM](http://diagram-consortium.org/downloads.html)
-   **UK Biobank traits**
    -   Traits: Loneliness, less social activity, hearing difficulty,
        less education, physical inactivity, smoking, alcohol intake,
        BMI, deprivation status, systolic blood pressure
    -   Website: <http://www.nealelab.is/uk-biobank>
    -   Data: [Neale
        lab](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=178908679)

## Example code for formatting UK Biobank GWAS data

There are certain columns and information that are required in the GWAS
datasets so that the quality control filtering step (known as munging)
can be completed before conducting LD score regression and genomic SEM.
These requirements are outlined in detail below in the *data munging*
sections.

The LDSC and GenomicSEM `munge` functions can recognise a range of
different column names so often you can just input the downloaded
dataset directly into this function and it will run.

However, there may be occasions where the GWAS data cannot be read by
the `munge` function because the software cannot recognise a certain
column name or if there is data missing that it needs. One of the main
examples of this happening is if you use the GWAS summary statistics of
UK Biobank traits that were conducted by the Neale lab. Therefore, we
provide a detailed example of how to appropriately format these datasets
for use with the LDSC and GenomicSEM software. This code was run in R
Studio (Version 1.2.5033) using R version 3.6.3 on a personal desktop.

We downloaded our UK Biobank GWAS traits from the [Neale lab
manifest](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=227859291)
of their round 2 European GWAS data. This is best done using the `wget`
command in the terminal.

### Set your working directory and load the required packages

``` r
setwd("/path/to/directory")
library(dplyr)
library(data.table)
```

### Make an edited version of the variant data to merge with the GWAS datasets

The individual GWAS files from the Neale lab only include some of the
information that you need to carry out the `munge` step (sample size,
MAF and effect statistics). Therefore, the Neale lab also provide a
general variant file that contains annotation information for all the
variants measured within a single file in the same order that they
appear in the GWAS files. By merging this file with the **variant**
column in the trait GWAS you can then get the rsIDs, A1 (i.e. the alt
allele), A2 (i.e. the ref allele) and INFO information that you need for
your analysis. See
[here](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=227859291)
for further information on the variant file and to download it.

**Note:** Since the UK Biobank files are very large we found it best to
initially make an edited variant file containing only the columns that
you need for this formatting using the code below as this significantly
reduces runtime in the later steps.

In addition, we use the `fread` function from the *data.table* R package
to read in the datasets in our below example. This requires the
downloaded files to be unzipped so that they are .tsv files rather than
.bgz ones (i.e. the downloaded format from the Neale lab website). Using
the base R `read.table` function also works but is *much* slower.

``` r
variant_data <- fread("variants.tsv")
variant_data <- as.data.table(variant_data)
head(variant_data)

variant_data_edited <- variant_data %>% select(variant, ref, alt, rsid, info)
head(variant_data_edited)

write.table(variant_data_edited, file = "variant_data_edited_UKB.txt", row.names = FALSE,
            quote = FALSE)
```

### Format UK Biobank GWAS data ready for munging

Now we can use this edited variant file to merge with the GWASs of our
traits of interest. Please see an example script for the trait *feelings
of loneliness/isolation* below. In addition to merging the data we use
the *dplyr* R package to `rename` the columns to names that are
recognised by the `munge` functions and we `select` only the required
columns to reduce the size of the edited file.

**Note:** The `merge` step takes a while to run.

``` r
### Read in original GWAS sumstats file using fread 

loneliness_UKB_data <- fread("2020.gwas.imputed_v3.both_sexes.tsv")
loneliness_UKB_data <- as.data.table(loneliness_UKB_data)
head(loneliness_UKB_data)

## read in our edited Neale variants file

variant_data <- fread("variant_data_edited_UKB.txt")
variant_data <- as.data.table(variant_data)
head(variant_data)

## join the phenotype file with the variant file by the variant column

joined_df <- merge(loneliness_UKB_data, variant_data , by.x = "variant", 
                   by.y = "variant", all.x = TRUE, all.y = FALSE)

head(joined_df)

## rename the columns & keep only those required for munging

loneliness_UKB_data_edited <- joined_df %>% rename(A1 = "alt") %>% 
  rename(A2 = "ref") %>% 
  rename(MAF = "minor_AF") %>% 
  rename(N = "n_complete_samples") %>% 
  rename(BETA = "beta") %>% 
  rename(SE = "se") %>% 
  rename(P = "pval") %>% 
  rename(SNP = "rsid") %>% 
  rename(INFO = "info") %>%
  select(A1, A2, MAF, N, BETA, SE, P, SNP, INFO)

head(loneliness_UKB_data_edited)

## Write new table as txt file munging step

write.table(loneliness_UKB_data_edited, file = "loneliness_UKB_GWAS_edit.txt", 
            row.names = FALSE, quote = FALSE)
```

# LD score regression

To ensure that our traits were sufficiently heritable to be meaningfully
modelled using genomic SEM we initially performed LD score regression to
calculate the SNP-based heritability of each trait. We also used
cross-trait LD score regression to assess the level of inter-correlation
present between each pair of traits.

This stage of the analysis utilised the command line *LDSC* software
developed by Brendan Bulik-Sullivan and colleagues (for in depth details
please refer to their [original wiki](https://github.com/bulik/ldsc)).
We performed the next few steps in the terminal on an Ubuntu personal
desktop (Linux OS).

### Install LDSC software

The LDSC software is run in the command line and relies on a collection
of Python dependencies to run. In order to install LDSC you first need
to download [Anaconda](https://www.anaconda.com/products/individual) so
that you can create an environment that includes the package
dependencies specifically required for the software. More instructions
on this can be found on the [original LDSC
wiki](https://github.com/bulik/ldsc).

However, we found that if you are using a more recent version of
Anaconda (we used Anaconda3 version 2020.07), you may run into some
functionality issues with the LDSC functions if you use the environment
file specified in the original software. We therefore used an adapted
environment file (environment_edit.yml) within our analysis.

**Adapted environment file**

In a text editor we wrote the following text and saved it as
**LDSC_env.yml**:

``` bash
name: environment_edit
channels:
- bioconda
dependencies:
- python=2.7
- bitarray=0.8
- nose=1.3
- pybedtools=0.7
- pip
- pip:
  - scipy==0.18
  - pandas==0.19
  - numpy==1.16
```

We then clone the repository for LDSC and create the environment with
the adapted package dependencies.

``` bash
git clone https://github.com/bulik/ldsc.git
cd ldsc
conda env create --file LDSC_env.yml
source activate LDSC_env
```

## Munge data

First we must munge the GWAS data so that it can be recognised by the
ldsc software and the data is filtered to only include well-imputed SNPs
with a MAF \>0.01. See the **Reformatting Summary Statistics** section
[here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)
for further details.

In short, for the `munge` function to work we need the following columns
in the downloaded GWAS summary statistics as a minimum:

-   SNP rsIDs
-   Allele 1 (i.e. the effect allele)
-   Allele 2 (i.e. the non-effect allele)
-   Sample size
-   The regression effect (e.g. beta, OR, logOR, Z-score)
-   P values of the effect

It is also useful to have the MAF and INFO columns for filtering, but
these are not crucial since we filter to a reference genome. Here we
used Hapmap 3 SNPs with the MHC region removed as our reference genome.
This data can be downloaded from:
<https://utexas.app.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v>.

``` bash
cd /media/user/Elements/R_markdown_test/

# AD

./ldsc/munge_sumstats.py \
--sumstats Kunkle_etal_Stage1_results.txt \
--N 63926 \
--out ad_ldsc \
--merge-alleles w_hm3.noMHC.snplist \

# MDD

./ldsc/munge_sumstats.py \
--sumstats MDD2018_ex23andMe.gz \
--N 173005 \
--out mdd_ldsc \
--merge-alleles w_hm3.noMHC.snplist \

# Insomnia

./ldsc/munge_sumstats.py \
--sumstats Insomnia_sumstats_Jansenetal.txt.gz \
--N 386533 \
--out insomnia_ldsc \
--merge-alleles w_hm3.noMHC.snplist \

# Loneliness

./ldsc/munge_sumstats.py \
--sumstats loneliness_UKB_GWAS_edit.txt \
--N 355583 \
--out loneliness_ldsc \
--merge-alleles w_hm3.noMHC.snplist \


# Low social activity

./ldsc/munge_sumstats.py \
--sumstats no_social_activity_UKB_GWAS_edit.txt \
--N 360063 \
--out no_social_activity_ldsc \
--merge-alleles w_hm3.noMHC.snplist \

# Hearing difficulty

./ldsc/munge_sumstats.py \
--sumstats hearing_diff_UKB_GWAS_edit.txt \
--N 353983 \
--out hearing_difficulty_ldsc \
--merge-alleles w_hm3.noMHC.snplist \

# Less education

./ldsc/munge_sumstats.py \
--sumstats low_education_UKB_GWAS_edit.txt \
--N 357549 \
--out low_education_ldsc \
--merge-alleles w_hm3.noMHC.snplist \

# Physical inactivity

./ldsc/munge_sumstats.py \
--sumstats physical_inactivity_UKB_GWAS_edit.txt \
--N 359263 \
--out physical_inactivity_ldsc \
--merge-alleles w_hm3.noMHC.snplist \

# Smoking

./ldsc/munge_sumstats.py \
--sumstats current_smoker_UKB_GWAS_edit.txt \
--N 359706 \
--out current_smoker_ldsc \
--merge-alleles w_hm3.noMHC.snplist \

# Alcohol intake frequency

./ldsc/munge_sumstats.py \
--sumstats alcohol_intake_freq_UKB_GWAS_edit.txt \
--N 360726 \
--out alcohol_intake_freq_ldsc \
--merge-alleles w_hm3.noMHC.snplist \

# BMI

./ldsc/munge_sumstats.py \
--sumstats bmi_UKB_GWAS_edit.txt \
--N 354831 \
--out bmi_ldsc \
--merge-alleles w_hm3.noMHC.snplist \

# Deprivation status

./ldsc/munge_sumstats.py \
--sumstats deprivation_UKB_GWAS_edit.txt \
--N 360763 \
--out deprivation_ldsc \
--merge-alleles w_hm3.noMHC.snplist \

# Systolic blood pressure

./ldsc/munge_sumstats.py \
--sumstats systolic_bp_UKB_GWAS_edit.txt \
--N 340159 \
--out systolic_bp_ldsc \
--merge-alleles w_hm3.noMHC.snplist \

# Type 2 diabetes

./ldsc/munge_sumstats.py \
--sumstats diabetes_with_rsids.txt \
--N 159208 \
--out type_2_diabetes_ldsc \
--merge-alleles w_hm3.noMHC.snplist
```

## Calculating univariate SNP-based heritability

Now we want to calculate the univariate SNP-based heritability of each
trait so that we can assess whether it is sufficiently heritable to be
included in our main analysis (i.e. whether it has a h2 Z score \>4). We
use the `ldsc.py` function and the `--h2` argument within the LDSC
software. See
[here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)
for further details.

You need to use LD scores and weights that match your data’s ancestral
background. Here we used pre-calculated European LD scores and weights
for our analysis from the Broad Institute computed from 1000 Genomes
data. These can be downloaded from:
<https://alkesgroup.broadinstitute.org/LDSCORE/>.

This group has made a range of different LD scores and weights that are
publically available but if you need to make some of your own see the
[LDSC
wiki](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial).

For **binary** traits you need to specify the sample and population
prevalences (see `sample-prev` and `population-prev` in the code below).

Sample prevalence can be calculated by dividing the *number of cases* by
the *total sample size*.

Population prevalences are based on a relevant estimate that is widely
cited and reliable for your population of interest (see the
**Supplementary Methods** of our
[paper](https://www.sciencedirect.com/science/article/pii/S0197458022001312)
for details on how we identified our population prevalence estimates).

You do not need to provide the sample and population prevalence for
**continuous** or **categorical** traits.

``` bash
# AD 

./ldsc/ldsc.py \
--h2 ad_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ad_h2 \
--samp-prev 0.34 \
--pop-prev 0.05 \

# MDD

./ldsc/ldsc.py \
--h2 mdd_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out mdd_h2 \
--samp-prev 0.35 \
--pop-prev 0.13 

# Insomnia
./ldsc/ldsc.py \
--h2 insomnia_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out insomnia_h2 \
--samp-prev 0.28 \
--pop-prev 0.32 

# Loneliness
./ldsc/ldsc.py \
--h2 loneliness_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out loneliness_h2 \
--samp-prev 0.18 \
--pop-prev 0.08 

# Low social activity
./ldsc/ldsc.py \
--h2 no_social_activity_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out no_social_activity_h2 \
--samp-prev 0.30 \
--pop-prev 0.30 

# Hearing difficulty
./ldsc/ldsc.py \
--h2 hearing_difficulty_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out hearing_difficulty_h2 \
--samp-prev 0.38 \
--pop-prev 0.39 

# Less education
./ldsc/ldsc.py \
--h2 low_education_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out low_education_h2 \
--samp-prev 0.17 \
--pop-prev 0.27 

# Physical inactivity
./ldsc/ldsc.py \
--h2 physical_inactivity_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out physical_inactivity_h2 \
--samp-prev 0.06 \
--pop-prev 0.18

# Smoking
./ldsc/ldsc.py \
--h2 current_smoker_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out current_smoker_h2 \
--samp-prev 0.10 \
--pop-prev 0.27

# Alcohol intake frequency
./ldsc/ldsc.py \
--h2 alcohol_intake_freq_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out alcohol_intake_freq_h2 

# BMI
./ldsc/ldsc.py \
--h2 bmi_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out bmi_h2 

# Deprivation status
./ldsc/ldsc.py \
--h2 deprivation_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out deprivation_h2 

# Systolic bp
./ldsc/ldsc.py \
--h2 systolic_bp_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out systolic_bp_h2 

# Type 2 diabetes
./ldsc/ldsc.py \
--h2 type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out type_2_diabetes_h2 \
--samp-prev 0.17 \
--pop-prev 0.06
```

## Calculating pairwise genetic correlations

Next, we want to calculate the pairwise genetic correlations between
each pair of traits to ascertain how much genetic correlation is present
between our traits. We use the `ldsc.py` function with the `--rg`
argument. See
[here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)
for further details.

As with the previous step we use pre-calculated European LD scores and
weights for our analysis from the Broad Institute computed from 1000
Genomes data downloaded from:
<https://alkesgroup.broadinstitute.org/LDSCORE/>.

``` bash
# AD
./ldsc/ldsc.py \
--rg ad_ldsc.sumstats.gz,mdd_ldsc.sumstats.gz,insomnia_ldsc.sumstats.gz,loneliness_ldsc.sumstats.gz,no_social_activity_ldsc.sumstats.gz,\
hearing_difficulty_ldsc.sumstats.gz,low_education_ldsc.sumstats.gz,physical_inactivity_ldsc.sumstats.gz,\
current_smoker_ldsc.sumstats.gz,alcohol_intake_freq_ldsc.sumstats.gz,bmi_ldsc.sumstats.gz,deprivation_ldsc.sumstats.gz,\
systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ad_rg \
--samp-prev 0.34,0.35,0.28,0.18,0.30,0.38,0.17,0.06,0.10,nan,nan,nan,nan,0.17 \
--pop-prev 0.05,0.13,0.32,0.08,0.30,0.39,0.27,0.18,0.27,nan,nan,nan,nan,0.06 

# MDD
./ldsc/ldsc.py \
--rg mdd_ldsc.sumstats.gz,insomnia_ldsc.sumstats.gz,loneliness_ldsc.sumstats.gz,no_social_activity_ldsc.sumstats.gz,\
hearing_difficulty_ldsc.sumstats.gz,low_education_ldsc.sumstats.gz,physical_inactivity_ldsc.sumstats.gz,\
current_smoker_ldsc.sumstats.gz,alcohol_intake_freq_ldsc.sumstats.gz,bmi_ldsc.sumstats.gz,deprivation_ldsc.sumstats.gz,\
systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out mdd_rg \
--samp-prev 0.35,0.28,0.18,0.30,0.38,0.17,0.06,0.10,nan,nan,nan,nan,0.17 \
--pop-prev 0.13,0.32,0.08,0.30,0.39,0.27,0.18,0.27,nan,nan,nan,nan,0.06

# insomnia
./ldsc/ldsc.py \
--rg insomnia_ldsc.sumstats.gz,loneliness_ldsc.sumstats.gz,no_social_activity_ldsc.sumstats.gz,\
hearing_difficulty_ldsc.sumstats.gz,low_education_ldsc.sumstats.gz,physical_inactivity_ldsc.sumstats.gz,\
current_smoker_ldsc.sumstats.gz,alcohol_intake_freq_ldsc.sumstats.gz,bmi_ldsc.sumstats.gz,deprivation_ldsc.sumstats.gz,\
systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out insomnia_rg \
--samp-prev 0.28,0.18,0.30,0.38,0.17,0.06,0.10,nan,nan,nan,nan,0.17 \
--pop-prev 0.32,0.08,0.30,0.39,0.27,0.18,0.27,nan,nan,nan,nan,0.06

# loneliness
./ldsc/ldsc.py \
--rg loneliness_ldsc.sumstats.gz,no_social_activity_ldsc.sumstats.gz,\
hearing_difficulty_ldsc.sumstats.gz,low_education_ldsc.sumstats.gz,physical_inactivity_ldsc.sumstats.gz,\
current_smoker_ldsc.sumstats.gz,alcohol_intake_freq_ldsc.sumstats.gz,bmi_ldsc.sumstats.gz,deprivation_ldsc.sumstats.gz,\
systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out loneliness_rg \
--samp-prev 0.18,0.30,0.38,0.17,0.06,0.10,nan,nan,nan,nan,0.17 \
--pop-prev 0.08,0.30,0.39,0.27,0.18,0.27,nan,nan,nan,nan,0.06

# low social activity
./ldsc/ldsc.py \
--rg no_social_activity_ldsc.sumstats.gz,\
hearing_difficulty_ldsc.sumstats.gz,low_education_ldsc.sumstats.gz,physical_inactivity_ldsc.sumstats.gz,\
current_smoker_ldsc.sumstats.gz,alcohol_intake_freq_ldsc.sumstats.gz,bmi_ldsc.sumstats.gz,deprivation_ldsc.sumstats.gz,\
systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out no_social_activity_rg \
--samp-prev 0.30,0.38,0.17,0.06,0.10,nan,nan,nan,nan,0.17 \
--pop-prev 0.30,0.39,0.27,0.18,0.27,nan,nan,nan,nan,0.06

# hearing difficulty
./ldsc/ldsc.py \
--rg hearing_difficulty_ldsc.sumstats.gz,low_education_ldsc.sumstats.gz,physical_inactivity_ldsc.sumstats.gz,\
current_smoker_ldsc.sumstats.gz,alcohol_intake_freq_ldsc.sumstats.gz,bmi_ldsc.sumstats.gz,deprivation_ldsc.sumstats.gz,\
systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out hearing_difficulty_rg \
--samp-prev 0.38,0.17,0.06,0.10,nan,nan,nan,nan,0.17 \
--pop-prev 0.39,0.27,0.18,0.27,nan,nan,nan,nan,0.06

# less education
./ldsc/ldsc.py \
--rg low_education_ldsc.sumstats.gz,physical_inactivity_ldsc.sumstats.gz,\
current_smoker_ldsc.sumstats.gz,alcohol_intake_freq_ldsc.sumstats.gz,bmi_ldsc.sumstats.gz,deprivation_ldsc.sumstats.gz,\
systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out low_education_rg \
--samp-prev 0.17,0.06,0.10,nan,nan,nan,nan,0.17 \
--pop-prev 0.27,0.18,0.27,nan,nan,nan,nan,0.06

# physical inactivity
./ldsc/ldsc.py \
--rg physical_inactivity_ldsc.sumstats.gz,\
current_smoker_ldsc.sumstats.gz,alcohol_intake_freq_ldsc.sumstats.gz,bmi_ldsc.sumstats.gz,deprivation_ldsc.sumstats.gz,\
systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out physical_inactivity_rg \
--samp-prev 0.06,0.10,nan,nan,nan,nan,0.17 \
--pop-prev 0.18,0.27,nan,nan,nan,nan,0.06

# smoking
./ldsc/ldsc.py \
--rg current_smoker_ldsc.sumstats.gz,alcohol_intake_freq_ldsc.sumstats.gz,bmi_ldsc.sumstats.gz,deprivation_ldsc.sumstats.gz,\
systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out current_smoker_rg \
--samp-prev 0.10,nan,nan,nan,nan,0.17 \
--pop-prev 0.27,nan,nan,nan,nan,0.06

# alcohol intake freq
./ldsc/ldsc.py \
--rg alcohol_intake_freq_ldsc.sumstats.gz,bmi_ldsc.sumstats.gz,deprivation_ldsc.sumstats.gz,\
systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out alcohol_intake_freq_rg \
--samp-prev nan,nan,nan,nan,0.17 \
--pop-prev nan,nan,nan,nan,0.06

# BMI
./ldsc/ldsc.py \
--rg bmi_ldsc.sumstats.gz,deprivation_ldsc.sumstats.gz,\
systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out bmi_rg \
--samp-prev nan,nan,nan,0.17 \
--pop-prev nan,nan,nan,0.06

# deprivation
./ldsc/ldsc.py \
--rg deprivation_ldsc.sumstats.gz,\
systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out deprivation_rg \
--samp-prev nan,nan,0.17 \
--pop-prev nan,nan,0.06

# systolic bp
./ldsc/ldsc.py \
--rg systolic_bp_ldsc.sumstats.gz,type_2_diabetes_ldsc.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out systolic_bp_rg \
--samp-prev nan,0.17 \
--pop-prev nan,0.06
```

The LDSC output gives us a table of the genetic correlations at the end
of the .log file that we can easily copy into Excel.

## Plotting the data

The remaining analytic steps were all performed in R Studio (Version
1.2.5033) using R version 3.6.3 on a Windows OS desktop.

### Genetic correlation matrix

Visualising pairwise genetic correlations between lots of traits in a
table is not an easy way to study the results. Plotting your genetic
correlations in a correlation matrix plot provides far better visuals of
the correlations allowing for better interpretation of your results. A
great R package for this is *corrplot* and we provide the code script
below for the figure we used in our
[paper](https://www.sciencedirect.com/science/article/pii/S0197458022001312).

For a detailed description of the different aesthetic options available
in this package see
[here](https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html)
for a wide range of vignettes.

First we load the required R packages and set our working directory. We
use *tidyverse* to turn our .csv file data into appropriate matrices for
the `corrplot` function.

``` r
setwd("E:/GWAS_sumstats")
library(tidyverse)
library(corrplot)
```

In our matrix plot we wanted to include both the correlation
coefficients between the traits but also information about the
statistical significance of the correlation. Therefore we initially made
2 separate .csv files based on the output from the LDSC analysis.

**.csv file 1** This file consisted of 3 columns named *trait_1*,
*trait_2* and *rg*. We wrote abbreviated trait names so that the labels
on the matrix plot will be easier to read (see image below). **Rg**
refers to the correlation coefficient estimate. As the central line on
the matrix plot is the trait correlated with itself we give these
variables a value of 1.

**.csv file 2** This file consisted of 3 columns named *trait_1*,
*trait_2* and *pval* with data in the same order as the 1st file. P
values of traits with themselves are specified as 0.

We then use these two files to construct a rg data matrix and p value
data matrix in R using the *tidyverse* package to use as the input for
our correlation plot (`corrplot` only accepts a data matrix as its
input).

``` r
# Load in the rg csv and make a duplicate for other side of the matrix
rg <- read.csv(file = "LDSC_all_gen_corr_no_pval.csv", header = T, sep = ",")
rg

rg2<-rg

# Swap the column names for trait1 and trait 2 
# (because the other half of the triangle is just a transposed version of the extra half we need)
names(rg2)<-c("trait_2", "trait_1", "rg")
rg2
  
# Drop the columns with an rg of 1 that we don't need as these are already in the first dataset
# and as they represent the matrix intersection we don't need the values twice

rg3<-rg2 %>% mutate(na_if(trait_1, trait_2)) %>% 
  drop_na() %>% 
  select(trait_2, trait_1, rg) #removing the extra column we used to complete this step
rg3

# Combine to make a matrix 
all<-rbind(rg, rg3) %>% #combine the two sides of the matrix
  pivot_wider(names_from=trait_1, values_from=rg) %>% # turn data into wide format
  column_to_rownames(var="trait_2") #transfer the trait_2 column to be the row names

all

# Save matrix
write.table(all, file = "ldsc_rg_matrix.txt")

# Repeat to make p value matrix
options(scipen=10)
pval <- read.csv(file = "LDSC_all_pval.csv", header = T, sep = ",")
pval

pval2<-pval

names(pval2)<-c("trait_2", "trait_1", "pval")
pval2

pval3<-pval2 %>% mutate(na_if(trait_1, trait_2)) %>% 
  drop_na() %>% 
  select(trait_2, trait_1, pval)
pval3

all_pval<-rbind(pval, pval3) %>% 
  pivot_wider(names_from=trait_1, values_from=pval) %>% 
  column_to_rownames(var="trait_2")

all_pval

write.table(all_pval, file = "ldsc_pval_matrix.txt")

# Convert to matrices for corrplot
all_matrix <- as.matrix(all)
all_matrix

pval_matrix <- as.matrix(all_pval)
pval_matrix
```

Now we can plot a matrix that depicts the genetic correlation
coefficients on the lower triangle of the matrix and the p values on the
upper triangle (marked by asterisks if they are statistically
significant). We use a Bonferroni corrected p value threshold to adjust
for multiple testing.

**Bonferroni calculation:**

**Step 1:** Calculate how many pairwise tests are being performed by:
![k=NumberOfTraits\*(NumberOfTraits+1)/2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%3DNumberOfTraits%2A%28NumberOfTraits%2B1%29%2F2 "k=NumberOfTraits*(NumberOfTraits+1)/2")

**Step 2:** Then calculate p value adjusting for the number of tests:
![P=0.05/k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P%3D0.05%2Fk "P=0.05/k")

In our case there were 14 traits giving us 105 pairwise tests with a
corrected p value of \<4.76E-4.

``` r
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# choose colour pallete

# save pdf of corr matrix with asterisks for significant values and no numeric values
pdf(file = "LDSC_all_pval.pdf", height = 5, width = 5.2)
corr_plot <- corrplot(all_matrix, method = "square", col = col(200),
                          type = "full",  order = "original",
                          tl.pos = "lt", tl.col = "black", tl.srt = 90, tl.cex = .7, # Text label color & rotation
                          cl.pos = "r", cl.cex = .7, cl.ratio = 0.1,
                          # Combine with significance
                          p.mat = pval_matrix, sig.level = 4.76E-4, insig = "label_sig",
                          pch.col = "black", pch.cex = 1.2, 
                          # hide correlation coefficient on the principal diagonal
                          diag = TRUE)
dev.off()

### PDF corr matrix where the upper matrix has p values asterisks 
# and the lower displays correlation coefficients
pdf(file = "LDSC_all_pval_with_values.pdf", height = 5, width = 5.2)
corrplot(all_matrix, method = "color", col = col(200),
         addCoef.col = "black", number.cex = .6, number.digits = 2, # Add coefficient of correlation
         type = "lower",  order = "original",
         tl.pos = "lt", tl.col = "black", tl.srt = 90, tl.cex = .7, # Text label color and rotation
         cl.pos = "r", cl.cex = .7, cl.ratio = 0.1, 
         # Combine with significance
         diag = TRUE)
corrplot(all_matrix, add = TRUE, type = "upper", method = "square", order = "original", col = col(200),
         diag = FALSE, tl.pos = "n", cl.pos = "n",
         p.mat = pval_matrix, sig.level = 4.76E-4, insig = "label_sig", 
         pch.col = "black", pch.cex = 1.2)
dev.off()
```

**Output image:**

<p align="center">
  <img width="500" src="/images/ldscmatrix_bivariate.jpg">
</p>

### Undirected weighted graph

The matrix plot gives us a nice image of pairwise correlations but it
can still be hard to distinguish how these inter-relate between multiple
traits. Therefore, we decided to plot a weighted undirected graph of the
statistically significant genetic correlations to visualise how the
traits inter-correlate with one another. For a nice example of this type
of graph where they include the correlation coefficient on the edges
please see **Figure 1B** in the [genomicSEM
paper](https://doi.org/10.1016/j.cell.2019.11.020) by Lee et al (2019).

We use the `qplot` function in the *qgraph* R package to make this plot
with data functions from the *tidyverse* package.

``` r
library(qgraph)
```

As with the LDSC plot, we require an input .csv file with 3 columns
(*trait_1*, *trait_2* and *rg*). However, for rg values that are not
statistically significant we write these values as 0 so that the
non-significant relationships are not included within our graph.

We then convert the data into a matrix form and save the plot as a pdf.

``` r
# Load in the .csv file
rg <- read.csv(file = "ldsc_sig_rg_network.csv", header = T, sep = ",")
rg

# Make duplicate to be the other side of the matrix

rg2<-rg

# Swap the column names for trait 1 and trait 2 
names(rg2)<-c("trait_2", "trait_1", "rg")
rg2

# Drop the columns with rg of 1 that we don't need as they are already in the first dataset
# and as they represent the intersection, we don't need the values twice
rg3<-rg2 %>% mutate(na_if(trait_1, trait_2)) %>% 
  drop_na() %>% 
  select(trait_2, trait_1, rg)
rg3

# Combine as a matrix 
all<-rbind(rg, rg3) %>% 
  pivot_wider(names_from=trait_1, values_from=rg) %>% # Convert data to wide format
  column_to_rownames(var="trait_2") # Transfer the trait_2 column to be the row names

all

#### make matrix
all_matrix <- as.matrix(all)
all_matrix


# Make a weighted undirected graph of significant correlations

pdf(file = "undirected_weighted_graph.pdf", height = 4, width = 6)
qplot <- qgraph(all_matrix, layout = "spring", minimum = 0, posCol = "darkblue", 
                negCol = "red", maximum = 1, labels = colnames(all_matrix), 
                label.scale.equal = FALSE, cut = NULL, edge.labels = FALSE, 
                edge.label.bg = FALSE)
dev.off()
```

<p align="center">
  <img width="600" src="/images/undirected_weighted_graph_kunkle.jpg">
</p>

# Genomic SEM

Now we conduct genomic factor analysis using the *GenomicSEM* R package,
but there are many other potential analyses one can perform (see the
**[original wiki](https://github.com/GenomicSEM/GenomicSEM/wiki)** for
further details).

### Install *GenomicSEM*

``` r
install.packages("devtools")
library(devtools)
install_github("GenomicSEM/GenomicSEM")
```

### Set working directory and load the *GenomicSEM* package

``` r
setwd("E:/GWAS_sumstats")
library(GenomicSEM)
```

## Munge data

We munged the data again using the `munge` function in the *GenomicSEM*
package since this version produces more conservative filtering results
than the original LDSC software (it checks whether both the A1 **and**
A2 alleles match the reference file instead of just checking A1, and
removes SNPs with missing MAF values).

As with the original LDSC software, certain columns are required for
this function to work.

The required columns are:

-   SNP rsID
-   A1 allele (i.e. the effect allele)
-   A2 allele (i.e. the non-effect allele)
-   The regression effect (can be logistic or linear)
-   The p value of the effect

As with the LDSC software, it is preferential to also have MAF and INFO
columns but is not essential since we filter the SNPs to a reference
list of well imputed SNP with MAF information. In contrast, the files do
not require a specific N column since we provide the overall sample size
as an argument in the `munge` function (but if there is a column present
then the function automatically reads this instead of your provided
values).

In our analysis we used Hapmap 3 SNPs with the MHC region removed as our
reference genome. This data can be downloaded from:
<https://utexas.app.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v>.

``` r
munge(files = c("Kunkle_etal_Stage1_results.txt", "MDD2018_ex23andMe.gz", "Insomnia_sumstats_Jansenetal.txt.gz", 
                "loneliness_UKB_GWAS_edit.txt", "no_social_activity_UKB_GWAS_edit.txt", 
                "hearing_diff_UKB_GWAS_edit.txt", "low_education_UKB_GWAS_edit.txt",
                "physical_inactivity_UKB_GWAS_edit.txt", "current_smoker_UKB_GWAS_edit.txt", 
                "alcohol_intake_freq_UKB_GWAS_edit.txt", "bmi_UKB_GWAS_edit.txt", 
                "deprivation_UKB_GWAS_edit.txt", "systolic_bp_UKB_GWAS_edit.txt", 
                "diabetes_with_rsids.txt"), # a list of the dataset files
      hm3 = "w_hm3.noMHC.snplist", ## file of the Hapmap3 SNPs with MHC region removed
      trait.names = c("AD", "MDD", "insomnia", "loneliness_UKB", "no_social_activity_UKB", 
                      "hearing_difficulty_UKB", "low_education_UKB", "physical_inactivity_UKB", 
                      "current_smoker_UKB", "alcohol_intake_UKB", "BMI_UKB", "deprivation_UKB", 
                      "systolic_BP_UKB", "Diabetes"), # list of names of your traits
      N=c(63926, 173005, 386533, 355583, 360063, 353983, 357549, 359263, 359706, 360726, 354831, 
          360763, 340159, 159208), # list of total sample sizes of your traits
      info.filter = 0.9, # imputation filter (default used)
      maf.filter = 0.01) # minor allele frequency filter (default used)
```

## Exploratory factor analysis performed in the odd autosomes

We then performed exploratory factor analysis in the odd autosomes using
the munged summary statistics files as our input.

### Load the additional required packages

``` r
library(Matrix)
library(stats)
```

### Conduct multivariable LD score regression (odd autosomes)

We ran multivariable LD-Score regression to calculate the genetic
covariance matrix (*S*) and associated sampling covariance matrix (*V*)
using the odd autosomes only. When you run an analysis where you are
performing both exploratory and confirmatory factor analysis splitting
the data in this way can guard against model overfitting. If you have
non-overlapping replication datasets for each of your traits then the
best practice would be to run EFA in the discovery data and the CFA in
the replication data, but this was not possible for all of the current
study traits.

As with the initial LD score regression performed using the LDSC
software, we use pre-calculated European LD scores and weights for our
analysis from the Broad Institute computed from 1000 Genomes data
downloaded from: <https://alkesgroup.broadinstitute.org/LDSCORE/>.

We use the same population and sample prevalence details as above.

**Note:** This step can take a while to run if you have lots of traits,
but for 14 traits it took only *2 minutes and 41 seconds* on a personal
desktop computer.

``` r
ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/" #folder of ld scores & weights

traits <- c("AD.sumstats.gz", "MDD.sumstats.gz", "insomnia.sumstats.gz",
            "loneliness_UKB.sumstats.gz", "no_social_activity_UKB.sumstats.gz",
            "hearing_difficulty_UKB.sumstats.gz", "low_education_UKB.sumstats.gz",
            "physical_inactivity_UKB.sumstats.gz", "current_smoker_UKB.sumstats.gz", 
            "alcohol_intake_UKB.sumstats.gz", "BMI_UKB.sumstats.gz", 
            "deprivation_UKB.sumstats.gz", "systolic_BP_UKB.sumstats.gz",
            "Diabetes.sumstats.gz") #munged sumstat files

sample.prev <- c(0.34,0.35,0.28,0.18,0.30,0.38,0.17,0.06,0.10,NA,NA,NA,NA,0.17)

population.prev <- c(0.05,0.13,0.32,0.08,0.30,0.39,0.27,0.18,0.27,NA,NA,NA,NA,0.06) 

trait.names <- c("AD", "MDD", "INS", "LON", "LSA", "HD", "LED",
                 "LPA", "SMK", "ALC", "BMI", "DEP", "SBP", "T2DM")

LDSCoutput_odd <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, 
                   ldsc.log="E:/GSEM_study/GSEM_AD_RISK_FACTORS_STUDY_ODD_LDSC", 
                   select = "ODD") #specify odd autosomes


##optional command to save the ldsc output in case you want to use it in a later R session. 

save(LDSCoutput_odd, file= "E:/GSEM_study/GSEM_AD_RISK_FACTORS_STUDY_ODD_LDSC.RData")
```

**Note:** If you are running an analysis with a lot of traits you may
have a k \>200. In such cases, you will need to specify this using the
`n.blocks` argument in the `ldsc` function because otherwise you will
not have enough blocks for the jackknifing procedure of ldsc (number of
blocks has to be higher than k by at least 1).

Example calculation for 14 traits (same as the calculation to calculate
the number of pairwise tests for Bonferroni correction):

![k=14\*(14+1)/2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%3D14%2A%2814%2B1%29%2F2 "k=14*(14+1)/2")

### Exploratory factor analysis

The covariance matrix from the multivariable LDSC output is used as the
input for factor analysis. We initially smooth the covariance matrix
using the `nearPD` function in the *Matrix* R package so that if the
covariance matrix is non-positive definite it is smoothed prior to
performing EFA. In our example the matrix is positive definite but we
run the function anyway to make sure and to gain an appropriately
labelled matrix so that the trait names are labelled in our factor
analysis output.

Then we use the smoothed covariance matrix to test different numbers of
factors using the `factanal` function in the R *stats* package to see
which model will be best to take forward to confirmatory factor
analysis. We tested a 2-, 3- and 4-factor model in our analysis and used
promax rotation (this allows for intercorrelation between the latent
factors).

``` r
Ssmooth<-as.matrix((nearPD(LDSCoutput_odd$S, corr = FALSE))$mat)

# 2-factor model

EFA<-factanal(covmat = Ssmooth, factors = 2, rotation = "promax")

print(EFA, sort=TRUE) # We sort the results by the magnitude of their loadings
print(EFA,digits=3,cutoff=.20,sort=TRUE) #specify the cut-off you want to use

# 3-factor model

EFA2<-factanal(covmat = Ssmooth, factors = 3, rotation = "promax")

print(EFA2, sort=TRUE)
print(EFA2,digits=3,cutoff=.20,sort=TRUE)

# 4-factor model

EFA3<-factanal(covmat = Ssmooth, factors = 4, rotation = "promax")

print(EFA3, sort=TRUE)
print(EFA3,digits=3,cutoff=.20,sort=TRUE)
```

The output prints details about the factor loadings, the factor scores,
the factor correlations and the uniqueness scores for each individual
trait. From this analysis we chose to base our subsequent analyses on a
3-factor model. Please see our
[paper](https://www.sciencedirect.com/science/article/pii/S0197458022001312)
for further details about the methods and how we chose how many factors
to extract based on the results from the above analysis.

### EFA loadings plot

The factor loadings table can be a little hard to visualise and
interpret, so we decided to make a plot of the loadings for our 3-factor
model based on an adapted version of code script written by Dan Mirman
(see [here](https://rpubs.com/danmirman/plotting_factor_analysis)). This
step uses the *ggplot2* R package.

**Load in the required additional packages**

``` r
library(ggplot2)
library(tidyr)
```

We made a .csv file in Excel of the factor loadings from the EFA step to
use as input for this code. However, we flipped the order of the traits
in the .csv file so that they appear in the correct order in the graph
(otherwise they will appear in reverse order).

``` r
# Load in csv file

EFA_loadings_sorted <- read.csv(file = "EFA_loadings_odd_ordered.csv")
EFA_loadings_sorted

# This step ensures that the traits stay in the order you want them in the final plot 
# (note this is still in reverse order here), otherwise they are plotted in a random order

EFA_loadings_sorted$Trait<- factor(EFA_loadings_sorted$Trait, 
                                   levels = EFA_loadings_sorted$Trait)

# Makes the plot long so column of factors 
# (use the gather function of the tidyr package)

EFA_loading_long <- gather(EFA_loadings_sorted, key="Factor", 
                           value="Loading", c("Factor_1", "Factor_2", "Factor_3"))

EFA_loading_long

######## Ensure the factor facets come out in the order you want

EFA_loading_long$Factor <- factor(EFA_loading_long$Factor, 
                                  levels = c("Factor_1", "Factor_2", "Factor_3"))

# Add a column to categorise the loadings by loading strength

EFA_loading_long$breaks <- cut(EFA_loading_long$Loading, 
                               breaks = c(-Inf, 0, .2, .4, Inf)) 

EFA_loading_long

# Create a name vector for labelling facets

facet_names <- c(
  `Factor_1` = "Factor 1",
  `Factor_2` = "Factor 2",
  `Factor_3` = "Factor 3")

# Make your plot with a chosen theme & save as a pdf

theme_set(theme_bw()) # set a chosen theme
pdf(file = "EFA_loadings_plot_3factors_odd.pdf", height = 5, width = 9)
EFA_loadings_plot <- ggplot(EFA_loading_long, aes(Trait, abs(Loading), fill=breaks)) + 
  facet_wrap(~ Factor, nrow=1, 
             labeller = as_labeller(facet_names)) + #place the factors in separate facets
  geom_bar(stat = "identity") + #make the bars
  coord_flip() + #flip the axes so the labels can be horizontal  
  #define the fill color gradient: blue=positive, red=negative
  scale_fill_manual(name = "Loading Strength", 
                    labels = c("Negative", "0 to .20", ".20 to .40", "Greater than .40"), 
                    values = c("#870E0EAA", "#256C33AA", "#316DB2AA", "#041A58AA")) + 
  ### tell what colours to assign to each break
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0, 1, by = .2)) + # change the scale of the y axis
  ylab("Standardised Factor Loading") + #name the y-axis 
  theme(panel.background = element_rect(fill = "snow3", colour = "snow3", size = 2, 
                                        linetype = "solid"), 
    panel.grid.major = element_line(size = 0.4, linetype = 'solid', colour = "white"), 
    panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "white"),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text.x = element_text(colour = "black"),
    axis.text = element_text(colour = "gray20"),
    panel.spacing.x=unit(1, "lines")) +
  geom_hline(yintercept=c(.20), linetype="dashed", color = "black") # add cut-off line
EFA_loadings_plot
dev.off()
```

The above code gives us our 3-factor EFA graph:

<p align="center">
  <img width="700" src="/images/EFA_loadings_plot_Kunkle.jpg">
</p>

## Confirmatory factor analysis performed in the even autosomes

Next, we want to test the fit of our model using confirmatory factor
analysis.

### Conduct multivariable LD score regression (even autosomes)

Because we are now performing the CFA in even autosome data we must
re-run the multivariable LD score regression step to estimate the
genetic covariance matrix for the even autosome.

``` r
LDSCoutput_even <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, 
                   ldsc.log="E:/GSEM_study/GSEM_AD_RISK_FACTORS_STUDY_EVEN_LDSC", 
                   select = "EVEN") #specify even autosomes

save(LDSCoutput_even, file= "E:/GSEM_study/GSEM_AD_RISK_FACTORS_STUDY_EVEN_LDSC.RData")
```

### Confirmatory factor analysis

Using the output from the even autosome LD score regression we can then
test different models using the `usermodel` function in the *GenomicSEM*
R package. We specify the parameters of our model based on the output of
our chosen EFA model (in this case the 3-factor model), including only
the factor loadings that pass our pre-specified cut-off level (≥0.20).
Model parameters are written in *lavaan* syntax.

We test a number of different CFA models to see how including negative
loadings, cross-loadings and higher loading cut-off points influence the
model fit statistics.

The output gives you the model fit statistics and the estimated
parameters for your models which can be drawn as path diagrams to better
visualise your model. The **unstandardised** estimates are equivalent to
covariance values whereas the **standardised** estimates are equivalent
to correlation coefficients.

We use the default `DWLS` estimation method but the models can also be
estimated using maximum likelhood `ML` estimation if you prefer.

``` r
# Model 1: CFA of 3-factor model: including traits with >.20 positive loadings and cross-loadings (DWLS method),
# no negative loadings; factor correlations present 

CFAofEFA <- 'F1 =~ NA*LPA + SMK + DEP + LED + MDD + INS + HD + BMI
F2 =~ NA*LPA + LSA + LED + ALC + LON + AD + INS + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput <- usermodel(LDSCoutput_even, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput

# Model 2: CFA of 3-factor model: including traits with >.20 loadings and cross-loadings (DWLS method),
# negative loadings included; factor correlations present 

CFAofEFA2 <- 'F1 =~ NA*LPA + SMK + DEP + LED + MDD + INS + HD + BMI + AD
F2 =~ NA*LPA + LSA + LED + ALC + LON + AD + INS + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS + SBP
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput2<- usermodel(LDSCoutput_even, estimation = "DWLS", model = CFAofEFA2, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput2

# Model 3: CFA of 3-factor model: including traits with >.20 loadings and no cross-loadings (DWLS method),
# negative loadings included; factor correlations present 

CFAofEFA3 <- 'F1 =~ NA*LPA + SMK + DEP
F2 =~ NA*LSA + LED + ALC + AD + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS + SBP
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput3 <- usermodel(LDSCoutput_even, estimation = "DWLS", model = CFAofEFA3, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput3

# Model 4: CFA of 3-factor model: including traits with >.20 loadings and no cross-loadings (DWLS method),
# negative loadings excluded; factor correlations present 

CFAofEFA4 <- 'F1 =~ NA*LPA + SMK + DEP
F2 =~ NA*LSA + LED + ALC + AD + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput4 <- usermodel(LDSCoutput_even, estimation = "DWLS", model = CFAofEFA4, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput4

# Model 5: CFA of 3-factor model: including traits with >.25 positive loadings and cross-loadings (DWLS method),
# no negative loadings; factor correlations present (KEPT INSOMNIA IN AS .24)

CFAofEFA5 <- 'F1 =~ NA*LPA + SMK + DEP + LED
F2 =~ NA*LPA + LSA + LED + ALC + LON + AD + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput5 <- usermodel(LDSCoutput_even, estimation = "DWLS", model = CFAofEFA5, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput5

# Model 6: CFA of 3-factor model: including traits with >.25 positive loadings and cross-loadings (DWLS method),
# negative loadings incl; factor correlations present (KEPT INSOMNIA IN AS .24)

CFAofEFA6 <- 'F1 =~ NA*LPA + SMK + DEP + LED
F2 =~ NA*LPA + LSA + LED + ALC + LON + AD + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS + SBP
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'


CFAoutput6 <- usermodel(LDSCoutput_even, estimation = "DWLS", model = CFAofEFA6, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput6

## .30 CUTOFF WAS THE SAME AS .25 SO DIDN'T RUN THAT
```

### Plotting CFA models

Although there are specific R packages available that plot path diagrams
for SEM automatically e.g. *semPlot*, we found that the clearest and
best way is to draw out your path diagrams in PowerPoint.

Below is an example of a path diagram drawn using this method based on
the standardised results of the CFA for *model 5* in the even autosomes
(see code above). This was the best fitting model that we tested so we
took it forward for subsequent analysis in data from all autosomes.

<p align="center">
  <img width="700" src="/images/CFA_even_best_model.png">
</p>

## Follow-up genome-wide analyses

Now we want to measure the fit of our best-fitting model and other more
complex SEM models using data from across the genome, so we run a final
multivariable LD score regression in all chromosomes.

### Genome-wide multivariable LD score regression

``` r
LDSCoutput_all <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, 
                   ldsc.log="E:/GSEM_study/GSEM_AD_RISK_FACTORS_STUDY_ALL_LDSC")

save(LDSCoutput_all, file= "E:/GSEM_study/GSEM_AD_RISK_FACTORS_STUDY_ALL_LDSC.RData")
```

### Genome-wide CFA

We then run the 3-factor model (based off Model 5 in the code script of
the even CFA) using the LDSC output estimated from all the chromosomes.

``` r
CFAofEFAall <- 'F1 =~ NA*LPA + SMK + DEP + LED
F2 =~ NA*LPA + LSA + LED + ALC + LON + AD + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput_all <- usermodel(LDSCoutput_all, estimation = "DWLS", model = CFAofEFAall, 
                           CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

CFAoutput_all
```

*The path diagram of the best-fitting CFA model (all autosomes):*

<p align="center">
  <img width="700" src="/images/CFA_all_best_model.png">
</p>

### Additional SEM modelling

The genome-wide 3-factor CFA model fits the data moderately well.
However, there is high inter-factor correlation present, which makes it
hard to meaningfully interpret distinct patterns of covariance.
Therefore, we tested a number of additional SEM models in the
*GenomicSEM* R package to try and disentangle the covariance
relationships in more detail. We provide examples of these models below
with brief explanations and example path diagrams.

### Common factor model

A common factor model is a model where there is a single latent
construct representing common variance between all of the traits in your
data. In our example although the traits load into a common factor, the
model fit statistics are poor, which suggests that the data is better
represented by a model with multiple latent constructs.

We run this model using the `commonfactor` function in *GenomicSEM*.

``` r
CommonFactor<- commonfactor(covstruc = LDSCoutput_all, estimation="DWLS")

CommonFactor
```

*The path diagram of the common factor model:*

<p align="center">
  <img width="700" src="/images/CommonFactor_all_SBP.png">
</p>

**Note:** Here in this walkthrough we use the LDSC output that includes
systolic blood pressure in the common factor model. However, in our
[paper](https://www.sciencedirect.com/science/article/pii/S0197458022001312)
we excluded systolic blood pressure in our genome-wide LD score
regression to estimate the parameters without it included.

### Hierarchal model

A hierarchal model refers to a structural equation model that has two
levels of latent constructs (i.e. a latent construct of shared
covariance between the specified lower level latent constructs).

Here we run a hierarchal model of our 3-factor CFA model.

``` r
hierarchal <- 'F1 =~ NA*LPA + SMK + DEP + LED
F2 =~ NA*LPA + LSA + LED + ALC + LON + AD + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS
F4 =~ F1 + F2 + F3'

hierarchal <- usermodel(LDSCoutput_all, estimation = "DWLS", model = hierarchal, 
                         CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
hierarchal
```

*The path diagram of the hierarchal model:*

<p align="center">
  <img width="700" src="/images/hierarchal_all.png">
</p>

**Note:** As you can see in our path diagram of the hierarchal model we
do not include the standard errors of our parameters. This is because
the model includes an endogenous variable so the STD_All column of the
results output is used to specify the standardised values. However, the
*GenomicSEM* package is unable to calculate standard errors for these
estimates (see the section titled **‘A Note on Standardized Output in
Genomic SEM’** on the original [wiki
page](https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects)
for further explanation).

### Bifactor model

Finally, we test the fit of a bifactor model. This model produces
orthogonal factors (i.e. uncorrelated factors) by including an
overarching general latent factor that all the traits load into, as well
as additional latent factors that only distinct sub-clusters of traits
load into. The benefit of these models is that they often provide better
model fit and allow for an enhanced interpretation of what the
sub-clusters might be capturing.

In our analysis, the bifactor model provided a good fit to our data and
the sub-clusters were specified based on the 3 factors in our CFA.

**Note:** For the bifactor model we had to constrain smoking to ensure
it didn’t output as a negative residual variance (see [original
wiki](https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects)
for further details on this).

``` r
bifactor <- 'F1 =~ NA*LPA + SMK + DEP + LED
F2 =~ NA*LPA + LSA + LED + ALC + LON + AD + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS
F4 =~ NA*AD + MDD + INS + LON + LSA + HD + LED + LPA + SMK + ALC + BMI + DEP + T2DM
F4 ~~ 0*F1
F4 ~~ 0*F2
F4 ~~ 0*F3
F1 ~~ 0*F2
F1 ~~ 0*F3
F2 ~~ 0*F3
SMK ~~  a*SMK
a > .001'

bifactor <- usermodel(LDSCoutput_all, estimation = "DWLS", model = bifactor, 
                      CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
bifactor
```

*The path diagram of the bifactor model:*

<p align="center">
  <img width="700" src="/images/bifactor_all.png">
</p>

## Models without the APOE region

Since Alzheimer’s disease risk is known to be disproportionately
influenced by variants in the *APOE* gene region, we ran the same models
as above in all of the autosomes except for the *APOE* region on
chromosome 19. We removed SNPs in the *APOE* region +/- 100kB (chr19:
45,309,039–45,512,650) based on the position in the genome assembly that
the GWAS summary statistics had been provided in (GRCh37/hg19).

We downloaded the *APOE* snplist from the [UCSC Genome
Browser](https://genome.ucsc.edu/index.html).

``` r
## read in snplist of APOE region +/- 100kb each side from UCSC browser for genome build 37 to match GWAS sumstats - 249 SNPs

apoe <- read.table("APOE_snps.txt", header = F, sep = "", check.names = F)

## Read in munged sumstats files and remove any SNPs that match the rsIDs in the APOE snplist file

AD <- fread("AD.sumstats.gz")
AD_no_apoe <- anti_join(AD, apoe, by = c("SNP" = "V4"))
write.table(AD_no_apoe, file = "AD_no_apoe.sumstats.gz", row.names = F, quote = F)

MDD <- fread("MDD.sumstats.gz")
MDD_no_apoe <- anti_join(MDD, apoe, by = c("SNP" = "V4"))
write.table(MDD_no_apoe, file = "MDD_no_apoe.sumstats.gz", row.names = F, quote = F)

INS <- fread("insomnia.sumstats.gz")
INS_no_apoe <- anti_join(INS, apoe, by = c("SNP" = "V4"))
write.table(INS_no_apoe, file = "INS_no_apoe.sumstats.gz", row.names = F, quote = F)

LON <- fread("loneliness_UKB.sumstats.gz")
LON_no_apoe <- anti_join(LON, apoe, by = c("SNP" = "V4"))
write.table(LON_no_apoe, file = "LON_no_apoe.sumstats.gz", row.names = F, quote = F)

LSA <- fread("no_social_activity_UKB.sumstats.gz")
LSA_no_apoe <- anti_join(LSA, apoe, by = c("SNP" = "V4"))
write.table(LSA_no_apoe, file = "LSA_no_apoe.sumstats.gz", row.names = F, quote = F)

HD <- fread("hearing_difficulty_UKB.sumstats.gz")
HD_no_apoe <- anti_join(HD, apoe, by = c("SNP" = "V4"))
write.table(HD_no_apoe, file = "HD_no_apoe.sumstats.gz", row.names = F, quote = F)

LED <- fread("low_education_UKB.sumstats.gz")
LED_no_apoe <- anti_join(LED, apoe, by = c("SNP" = "V4"))
write.table(LED_no_apoe, file = "LED_no_apoe.sumstats.gz", row.names = F, quote = F)

LPA <- fread("physical_inactivity_UKB.sumstats.gz")
LPA_no_apoe <- anti_join(LPA, apoe, by = c("SNP" = "V4"))
write.table(LPA_no_apoe, file = "LPA_no_apoe.sumstats.gz", row.names = F, quote = F)

SMK <- fread("current_smoker_UKB.sumstats.gz")
SMK_no_apoe <- anti_join(SMK, apoe, by = c("SNP" = "V4"))
write.table(SMK_no_apoe, file = "SMK_no_apoe.sumstats.gz", row.names = F, quote = F)

ALC <- fread("alcohol_intake_UKB.sumstats.gz")
ALC_no_apoe <- anti_join(ALC, apoe, by = c("SNP" = "V4"))
write.table(ALC_no_apoe, file = "ALC_no_apoe.sumstats.gz", row.names = F, quote = F)

BMI <- fread("BMI_UKB.sumstats.gz")
BMI_no_apoe <- anti_join(BMI, apoe, by = c("SNP" = "V4"))
write.table(BMI_no_apoe, file = "BMI_no_apoe.sumstats.gz", row.names = F, quote = F)

DEP <- fread("deprivation_UKB.sumstats.gz")
DEP_no_apoe <- anti_join(DEP, apoe, by = c("SNP" = "V4"))
write.table(DEP_no_apoe, file = "DEP_no_apoe.sumstats.gz", row.names = F, quote = F)

SBP <- fread("systolic_BP_UKB.sumstats.gz")
SBP_no_apoe <- anti_join(SBP, apoe, by = c("SNP" = "V4"))
write.table(SBP_no_apoe, file = "SBP_no_apoe.sumstats.gz", row.names = F, quote = F)

T2DM <- fread("Diabetes.sumstats.gz")
T2DM_no_apoe <- anti_join(T2DM, apoe, by = c("SNP" = "V4"))
write.table(T2DM_no_apoe, file = "T2DM_no_apoe.sumstats.gz", row.names = F, quote = F)


# Multivariable LD score regression of data without APOE

traits.apoe <- c("AD_no_apoe.sumstats.gz", "MDD_no_apoe.sumstats.gz", "INS_no_apoe.sumstats.gz", "LON_no_apoe.sumstats.gz", 
                 "LSA_no_apoe.sumstats.gz", "HD_no_apoe.sumstats.gz", "LED_no_apoe.sumstats.gz", "LPA_no_apoe.sumstats.gz", 
                 "SMK_no_apoe.sumstats.gz", "ALC_no_apoe.sumstats.gz", "BMI_no_apoe.sumstats.gz", "DEP_no_apoe.sumstats.gz",
                 "SBP_no_apoe.sumstats.gz", "T2DM_no_apoe.sumstats.gz")

LDSCoutput_exAPOE <- ldsc(traits.apoe, sample.prev, population.prev, ld, wld, trait.names, 
                   ldsc.log="E:/GSEM_study/GSEM_AD_RISK_FACTORS_STUDY_ALL_LDSC_noAPOE")

save(LDSCoutput_exAPOE, file= "E:/GSEM_study/GSEM_AD_RISK_FACTORS_STUDY_ALL_LDSC_noAPOE.RData")

# CFA model

CFAofEFA_exAPOE <- 'F1 =~ NA*LPA + SMK + DEP + LED
F2 =~ NA*LPA + LSA + LED + ALC + LON + AD + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput_exAPOE <- usermodel(LDSCoutput_exAPOE, estimation = "DWLS", model = CFAofEFA_exAPOE, 
                            CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput_exAPOE

# Common factor

CommonFactor_exAPOE<- commonfactor(covstruc = LDSCoutput_exAPOE, estimation="DWLS")
CommonFactor_exAPOE

# Hierarchal model

hierarchal_exAPOE <- 'F1 =~ NA*LPA + SMK + DEP + LED
F2 =~ NA*LPA + LSA + LED + ALC + LON + AD + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS
F4 =~ F1 + F2 + F3'

hierarchal_exAPOE <- usermodel(LDSCoutput_exAPOE, estimation = "DWLS", model = hierarchal_exAPOE, 
                              CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
hierarchal_exAPOE

# Bifactor model

bifactor_exPOE <- 'F1 =~ NA*LPA + SMK + DEP + LED
F2 =~ NA*LPA + LSA + LED + ALC + LON + AD + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS
F4 =~ NA*AD + MDD + INS + LON + LSA + HD + LED + LPA + SMK + ALC + BMI + DEP + T2DM
F4 ~~ 0*F1
F4 ~~ 0*F2
F4 ~~ 0*F3
F1 ~~ 0*F2
F1 ~~ 0*F3
F2 ~~ 0*F3
SMK ~~  a*SMK
a > .001'

bifactor_exAPOE <- usermodel(LDSCoutput_exAPOE, estimation = "DWLS", model = bifactor_exAPOE, 
                           CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

bifactor_exAPOE
```

# Sensitivity analyses

Because the factor loadings for Alzheimer’s disease varied substantially
between the odd autosome EFA, even autosome CFA and genome-wide CFA in
comparison to the risk factor traits we ran a number of follow-up
sensitivity analyses to explore why this might be.

### LD score regression in different chromosomal groupings

We used the three outputs generated in the main study using the `ldsc`
function in the *GenomicSEM* R package to compare the SNP-based
heritability and genetic correlation estimates of the included traits
when different chromosomal groupings were measured (odd vs. even
vs. all). These estimates can be found in the .log files that are
produced during previous `ldsc` steps.

We then plotted the genetic correlations using `corrplot` so that we
could more easily visualise any differences between Alzheimer’s disease
and its risk factors.

**Note:** Because the output of the `ldsc` function is already an R data
matrix we don’t need a separate .csv file like the earlier correlation
plot since we can just use the matrix directly as the input. But we need
to convert it to a correlation matrix using the `cov2cor` function.

``` r
# Smooth and convert the odd and even LDSC matrices to correlation matrices

Ssmooth_odd<-as.matrix((nearPD(LDSCoutput_odd$S, corr = FALSE))$mat)
corr_matrix_odd <- cov2cor(Ssmooth_odd)

Ssmooth_even<-as.matrix((nearPD(LDSCoutput_even$S, corr = FALSE))$mat)
corr_matrix_even <- cov2cor(Ssmooth_even)

### Make combined LDSC matrix of odd and even chromosomes 

pdf(file = "LDSC_matrix_odd_even.pdf", height = 5, width = 5.2)
corrplot(corr_matrix_odd, method = "color", col = col(200),
         addCoef.col = "black", number.cex = .6, number.digits = 2, 
         type = "lower",  order = "original",
         tl.pos = "lt", tl.col = "black", tl.srt = 90, tl.cex = .7, 
         cl.cex = .7, cl.ratio = 0.1, 
         diag = TRUE)
corrplot(corr_matrix_even, add = TRUE, method = "color", col = col(200),
         addCoef.col = "black", number.cex = .6, number.digits = 2, 
         type = "upper",  order = "original",
         diag = FALSE, tl.pos = "n", cl.pos = "n")
dev.off()
```

Output of the genetic correlation matrix with the odd autosome results
on the lower triangle and the even autosome results on the upper
triangle:

<p align="center">
  <img width="500" src="/images/odd_even_ldsc_matrix.jpg">
</p>

## EFA comparison across chromosomal groupings

Since the correlation coefficients were substantially different for AD
between the odd and even chromosomes, we also used EFA to test the
loadings for a 3-factor model based on the even autosomal data and all
autosomes to compare the factor loadings results with our initial
results based on the odd autosomes.

``` r
# 3-factor model even autosomes

EFA_even<-factanal(covmat = Ssmooth_even, 
                   factors = 3, rotation = "promax")

print(EFA_even, sort=TRUE)
print(EFA_even,digits=3,cutoff=.20,sort=TRUE)

# 3-factor model all autosomes

Ssmooth_all<-as.matrix((nearPD(LDSCoutput_all$S, corr = FALSE))$mat)
EFA_all<-factanal(covmat = Ssmooth_all, 
                   factors = 3, rotation = "promax")

print(EFA_all, sort=TRUE)
print(EFA_all,digits=3,cutoff=.20,sort=TRUE)
```

## CFA comparison across chromosomal groupings

We also ran post-hoc sensitivity CFA analyses using the outputs from the
sensitivity EFA to compare the model fit and loadings of the best
fitting models based on the even autosome and all autosome EFA to see
how they compared to our main analysis.

``` r
# CFA models of even autosome EFA run in the odd autosomes

# Model 1: CFA of 3-factor model: including traits with >.20 loadings and cross-loadings (DWLS method),
# negative loadings excluded; factor correlations present

CFAofEFA_posthoc1 <- 'F1 =~ NA*LPA + SMK + DEP
F2 =~ NA*LSA + LED + LPA + ALC + DEP + BMI + T2DM + SBP
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1
SMK ~~  a*SMK
a > .001'

CFAoutput_posthoc1 <- usermodel(LDSCoutput_odd, estimation = "DWLS", model = CFAofEFA_posthoc1, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput_posthoc1

# Model 2: CFA of 3-factor model: including traits with >.20 loadings and no cross-loadings (DWLS method),
# negative loadings excluded; factor correlations present

CFAofEFA_posthoc2 <- 'F1 =~ NA*SMK + DEP
F2 =~ NA*LSA + LED + LPA + ALC + BMI + T2DM + SBP
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1
SMK ~~  a*SMK
a > .001'

CFAoutput_posthoc2 <- usermodel(LDSCoutput_odd, estimation = "DWLS", model = CFAofEFA_posthoc2, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput_posthoc2


# Model 3: CFA of 3-factor model: including traits with >.25 positive loadings and cross-loadings (DWLS method),
# no negative loadings; factor correlations present (best fitting parameters for odd autosomes) - had to remove cross loading for deprivation status for it to converge

CFAofEFA_posthoc3 <- 'F1 =~ NA*SMK + DEP
F2 =~ NA*LPA + LSA + LED + ALC + BMI + SBP + T2DM
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput_posthoc3 <- usermodel(LDSCoutput_odd, estimation = "DWLS", model = CFAofEFA_posthoc3, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput_posthoc3

## CFA conducted in all autosomes but based off the parameters from the even EFA (same as the one tested in odd autosomes) - sensitivity analysis

# Model 4: CFA of 3-factor model: including traits with >.20 loadings and cross-loadings (DWLS method),
# negative loadings excluded; factor correlations present

CFAofEFA_posthoc4 <- 'F1 =~ NA*LPA + SMK + DEP
F2 =~ NA*LSA + LED + LPA + ALC + DEP + BMI + T2DM + SBP
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1
SMK ~~  a*SMK
a > .001'

CFAoutput_posthoc4 <- usermodel(LDSCoutput_all, estimation = "DWLS", model = CFAofEFA_posthoc4, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput_posthoc4

# Model 5: CFA of 3-factor model: including traits with >.20 loadings and no cross-loadings (DWLS method),
# negative loadings excluded; factor correlations present

CFAofEFA_posthoc5 <- 'F1 =~ NA*SMK + DEP
F2 =~ NA*LSA + LED + LPA + ALC + BMI + T2DM + SBP
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1
SMK ~~  a*SMK
a > .001'

CFAoutput_posthoc5 <- usermodel(LDSCoutput_all, estimation = "DWLS", model = CFAofEFA_posthoc5, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput_posthoc5


# Model 6: CFA of 3-factor model: including traits with >.25 positive loadings and cross-loadings (DWLS method),
# no negative loadings; factor correlations present (best fitting parameters for odd autosomes)

CFAofEFA_posthoc6 <- 'F1 =~ NA*SMK + DEP
F2 =~ NA*LPA + LSA + LED + ALC + BMI + SBP + T2DM + DEP
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput_posthoc6 <- usermodel(LDSCoutput_all, estimation = "DWLS", model = CFAofEFA_posthoc6, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput_posthoc6

# CFA models of all autosome EFA run in all autosomes

# Model 7: CFA of 3-factor model: including traits with >.25 positive loadings and cross-loadings (DWLS method),
# no negative loadings; factor correlations present (KEPT INSOMNIA IN AS .24)

CFAofEFA_posthoc7 <- 'F1 =~ NA*LPA + SMK + DEP + LED
F2 =~ NA*LPA + LSA + LED + ALC + BMI + T2DM + SBP
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput_posthoc7 <- usermodel(LDSCoutput_all, estimation = "DWLS", model = CFAofEFA_posthoc7, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput_posthoc7

# Model 8: CFA of 3-factor model: including traits with >.20 loadings and cross-loadings (DWLS method),
# negative loadings excluded; factor correlations present

CFAofEFA_posthoc8 <- 'F1 =~ NA*LPA + SMK + DEP + LED
F2 =~ NA*LSA + LED + LPA + ALC + INS + BMI + T2DM + SBP
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput_posthoc8 <- usermodel(LDSCoutput_all, estimation = "DWLS", model = CFAofEFA_posthoc8, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput_posthoc8

# Model 9: CFA of 3-factor model: including traits with >.20 loadings and no cross-loadings (DWLS method),
# negative loadings excluded; factor correlations present

CFAofEFA_posthoc9 <- 'F1 =~ NA*SMK + DEP
F2 =~ NA*LSA + LED + LPA + ALC + BMI + T2DM + SBP
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput_posthoc9 <- usermodel(LDSCoutput_all, estimation = "DWLS", model = CFAofEFA_posthoc9, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAoutput_posthoc9
```

# Final remarks

We hope that you find this Markdown useful. If you use this code in your
own analyses please cite this Github page. If you run into any issues or
queries related to this work please don’t hesitate to contact Isabelle
Foote (<isabelle.foote@colorado.edu>).
