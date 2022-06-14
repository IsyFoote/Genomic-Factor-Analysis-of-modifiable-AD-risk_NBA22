## Code script for genomic factor analysis and sensitivity analyses from our paper (outlined in detail in the code walkthrough readme)

setwed("/path/to/directory")

# Install genomic SEM

install.packages("devtools")
library(devtools)
install_github("GenomicSEM/GenomicSEM")
library(GenomicSEM)

# Munge GWAS data

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

# Run multivariable LD score regression odd autosomes

library(Matrix)
library(stats)

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

## Run EFA
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

# Make EFA loadings plot (need to premake csv file as described in walkthrough)
library(ggplot2)
library(tidyr)

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

# Run multivariable LDSC in even autosomes and perform CFA

LDSCoutput_even <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, 
                        ldsc.log="E:/GSEM_study/GSEM_AD_RISK_FACTORS_STUDY_EVEN_LDSC", 
                        select = "EVEN") #specify even autosomes

save(LDSCoutput_even, file= "E:/GSEM_study/GSEM_AD_RISK_FACTORS_STUDY_EVEN_LDSC.RData")

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

## Run multivariable LD score regression in all autosomes and do genome-wide CFA

LDSCoutput_all <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, 
                       ldsc.log="E:/GSEM_study/GSEM_AD_RISK_FACTORS_STUDY_ALL_LDSC")

save(LDSCoutput_all, file= "E:/GSEM_study/GSEM_AD_RISK_FACTORS_STUDY_ALL_LDSC.RData")

CFAofEFAall <- 'F1 =~ NA*LPA + SMK + DEP + LED
F2 =~ NA*LPA + LSA + LED + ALC + LON + AD + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS
F1 ~~ F2
F2 ~~ F3
F3 ~~ F1'

CFAoutput_all <- usermodel(LDSCoutput_all, estimation = "DWLS", model = CFAofEFAall, 
                           CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

CFAoutput_all

# Common Factor model
CommonFactor<- commonfactor(covstruc = LDSCoutput_all, estimation="DWLS")

CommonFactor

# Hierarchal model
hierarchal <- 'F1 =~ NA*LPA + SMK + DEP + LED
F2 =~ NA*LPA + LSA + LED + ALC + LON + AD + BMI + T2DM
F3 =~ NA*MDD + LON + HD + INS
F4 =~ F1 + F2 + F3'

hierarchal <- usermodel(LDSCoutput_all, estimation = "DWLS", model = hierarchal, 
                        CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
hierarchal

# Bifactor model
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

## Models without APOE
library(dplyr)
library(data.table)

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

## SENSITIVITY ANALYSES

# Genetic correlation plot comparing odd and even autosomes

library(corrplot)

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

# Sensitivity EFA

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

# Sensitivity CFA

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