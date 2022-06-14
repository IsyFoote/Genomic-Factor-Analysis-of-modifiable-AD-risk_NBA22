## Example code to curate GWAS summary statistics from Neale lab UK Biobank GWAS for use with LDSC and genomicSEM software

setwd("/path/to/directory")
library(dplyr)
library(data.table)

# read in the generic variants file provided by the Neale lab

variant_data <- fread("variants.tsv")
variant_data <- as.data.table(variant_data)
head(variant_data)

variant_data_edited <- variant_data %>% select(variant, ref, alt, rsid, info)
head(variant_data_edited)

write.table(variant_data_edited, file = "variant_data_edited_UKB.txt", row.names = FALSE,
            quote = FALSE)

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