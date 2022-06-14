# CLONE REPO AND SET ENVIRONMENT

git clone https://github.com/bulik/ldsc.git
cd ldsc
conda env create --file LDSC_env.yml
source activate LDSC_env

cd /path/to/directory

## MUNGE DATA

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

## CALCULATE SNP-BASED HERITABILITY

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

## CALCULATE GENETIC CORRELATION

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