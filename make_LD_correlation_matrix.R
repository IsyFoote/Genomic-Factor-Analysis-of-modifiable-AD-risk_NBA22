## R script to create LD correlation matrix

setwd("/path/to/directory")

# Make genetic correlation matrix - assumes you have made the relevant csv files as described in the walkthrough
library(tidyverse)
library(corrplot)

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
