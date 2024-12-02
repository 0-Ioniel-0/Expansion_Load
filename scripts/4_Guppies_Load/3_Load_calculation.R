### R v.4.3.3
### This script calculates load for each individual and nucleotide diversity (then effective population size), please note that population names are the ones used in our study.
###
###
###

print("Starting.")
starting_time <- Sys.time()
print(starting_time)

library(tidyr)
library(dplyr)

# Load big data sets
current_dir <- getwd()
setwd(current_dir)

# This file is an output from pixy
PI <- read.delim("~/Documents/Expantion_Load/PI.txt")
# This file is an output from Data_matrix_creation script
load <- read.delim("~/Documents/Expantion_Load/load.txt")

## Number of all derived alleles
DAcount <- all_derived_alleles %>%
  group_by(Sample,Population) %>%
  summarize(TotalDerivedAlleles=sum(Derived_Alleles))

## Only deletrious!
deletrious <- all_derived_alleles %>%
  filter(Deletriousness == "deletrious")

## Sum of [deleterious alleles * CS]
deletrious_sum <- deletrious %>%
  group_by(Sample,Population) %>%
  summarize(DelSum=sum(Derived_Alleles))

### Calaculate CS-based Load [deletrious DA / all DA]
indirect_load <- cbind(DAcount, "DelSum" = deletrious_sum$DelSum)
indirect_load$Load <- indirect_load$DelSum / indirect_load$TotalDerivedAlleles
# It's good to save the df at this point
write.csv(indirect_load, file = "indirect_load.csv")

## Effective Population Size
PI$pop <- factor(PI$pop, levels = c("DR", "ROX", "ARI", "CAU", "LOP", "SCR", "LSA", "QUA_up", "QUA_lo", "ORO_up", "ORO_lo", "TUR_up", "TUR_mid_up", "TUR_mid_lo", "TUR_lo"))
nucleotide_diversity <- PI %>%
  group_by(pop) %>%
  summarize(TotalDifferences=sum(count_diffs), TotalComparisons=sum(count_comparisons))
nucleotide_diversity$Pi <- nucleotide_diversity$TotalDifferences/nucleotide_diversity$TotalComparisons
nucleotide_diversity$Ne <- nucleotide_diversity$Pi/(4*2.9e-9)

# Populations order
indirect_load$Population <- factor(indirect_load$Population, levels = c("DR", "ROX", "ARI", "CAU", "LOP", "SCR", "LSA", "QUA_up", "QUA_lo", "ORO_up", "ORO_lo", "TUR_up", "TUR_mid_up", "TUR_mid_lo", "TUR_lo"))

### Average Load across populations and add to df with Ne
mean_load <- indirect_load %>%
  group_by(Population) %>%
    summarize(Load=mean(Load))
nucleotide_diversity$Load <- mean_load$Load

## Homozygotes
homo_deletrious <- all_derived_alleles %>%
  filter(Deletriousness == "deletrious") %>%
  filter(Derived_Alleles == 2)
homo_deletrious_sum <- homo_deletrious %>%
  group_by(Sample,Population) %>%
  summarize(DelSum=sum(Conservation_Score))
homo_indirect_load <- cbind(DAcount, "DelSum" = homo_deletrious_sum$DelSum)
homo_indirect_load$Load <- homo_indirect_load$DelSum / homo_indirect_load$TotalDerivedAlleles
#
homo_indirect_load$Population <- factor(homo_indirect_load$Population, levels = c("DR", "ROX", "ARI", "CAU", "LOP", "SCR", "LSA", "QUA_up", "QUA_lo", "ORO_up", "ORO_lo", "TUR_up", "TUR_mid_up", "TUR_mid_lo", "TUR_lo"))

## Heterozygotes
hetero_deletrious <- hetero %>%
  filter(Deletriousness == "deletrious") %>%
  filter(Derived_Alleles == 1)
hetero_deletrious_sum <- hetero_deletrious %>%
  group_by(Sample,Population) %>%
  summarize(DelSum=sum(Conservation_Score))
hetero_indirect_load <- cbind(DAcount, "DelSum" = hetero_deletrious_sum$DelSum)
hetero_indirect_load$Load <- hetero_indirect_load$DelSum / hetero_indirect_load$TotalDerivedAlleles
#
hetero_indirect_load$Population <- factor(hetero_indirect_load$Population, levels = c("DR", "ROX", "ARI", "CAU", "LOP", "SCR", "LSA", "QUA_up", "QUA_lo", "ORO_up", "ORO_lo", "TUR_up", "TUR_mid_up", "TUR_mid_lo", "TUR_lo"))

### Loss of Function data
LOF <- all_derived_alleles %>%
  filter(Effect == "HIGH")
summary(LOF)
LOF$Derived_Alleles <- as.factor(LOF$Derived_Alleles)
LOF_alleles <- LOF %>%
  group_by(Sample, Population) %>%
  count(Derived_Alleles)

## LOF Homozygotes 
LOF_homo <- LOF_alleles %>%
  filter(Derived_Alleles == "2")
LOF_homo_joined <- left_join(LOF_homo, DAcount, by = c("Population", "Sample"))
LOF_homo_joined$Normalized <- LOF_homo_joined$n/LOF_homo_joined$TotalDerivedAlleles
## LOF Heterozygotes
LOF_hetero <- LOF_alleles %>%
  filter(Derived_Alleles == "1")
LOF_hetero_joined <- left_join(LOF_hetero, DAcount, by = c("Population", "Sample"))
LOF_hetero_joined$Normalized <- LOF_hetero_joined$n/LOF_hetero_joined$TotalDerivedAlleles
# Let's save
LOF_joined <- rbind(LOF_homo_joined, LOF_hetero_joined)
write.csv(LOF_joined, file = "LOF_joined.csv")
