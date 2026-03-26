########################################
## LOAD LIBRARIES
########################################
library(DHARMa)
library(glmmTMB)
library(broom)
library(broom.mixed)
library(dplyr)
library(openxlsx)

########################################
## SET WORKING DIRECTORY
########################################
setwd("/dir/")

########################################
## LOAD DATA
########################################
LOF_Pi <- read.delim("LOF_Pi.txt")
MISS_Pi <- read.delim("MISS_Pi.txt")
RGL_Pi <- read.delim("RGL_Pi.txt")

LOF_site <- read.csv("LOF_site.txt")
MISS_site <- read.csv("MISS_site.txt")
RGL_site <- read.csv("RGL_site.txt")

LOF_TUR <- read.csv("LOF_TUR.txt")
MISS_TUR <- read.csv("MISS_TUR.txt")
RGL_TUR <- read.csv("RGL_TUR.txt")

LOF_bottleneck <- read.delim("LOF_bottleneck.txt")
MISS_bottleneck <- read.csv("MISS_bottleneck.txt")
RGL_bottleneck <- read.csv("RGL_bottleneck.txt")

########################################
## SET FACTOR LEVELS (CRITICAL)
########################################
RGL_TUR$Population <- factor(RGL_TUR$Population,
                             levels = c("TUR upper","TUR middle","TUR lower"))

MISS_TUR$Population <- factor(MISS_TUR$Population,
                              levels = c("TUR upper","TUR middle","TUR lower"))

LOF_TUR$Population <- factor(LOF_TUR$Population,
                             levels = c("TUR upper","TUR middle","TUR lower"))

RGL_site$Population <- factor(RGL_site$Population, levels = c("ORO","QUA"))
RGL_site$Site <- factor(RGL_site$Site, levels = c("lower","upper"))

MISS_site$Population <- factor(MISS_site$Population, levels = c("ORO","QUA"))
MISS_site$Site <- factor(MISS_site$Site, levels = c("lower","upper"))

LOF_site$Population <- factor(LOF_site$Population, levels = c("ORO","QUA"))
LOF_site$Site <- factor(LOF_site$Site, levels = c("lower","upper"))

########################################
## CORE FUNCTION
########################################
extract_model <- function(model, dependent_name, table_name) {
  
  if (inherits(model, "lm")) {
    tidy_mod <- broom::tidy(model)
    stat_type <- "t"
    
  } else if (inherits(model, "glmmTMB")) {
    tidy_mod <- broom.mixed::tidy(model, effects = "fixed")
    stat_type <- "z"
    
  } else if (inherits(model, "glm")) {
    tidy_mod <- broom::tidy(model)
    stat_type <- "z"
    
  } else {
    stop("Unsupported model type")
  }
  
  tidy_mod %>%
    mutate(
      Table = table_name,
      Dependent = dependent_name,
      Estimate = round(estimate, 4),
      SE = round(std.error, 4),
      Statistic = round(statistic, 3),
      P = p.value,
      Stat_type = stat_type
    ) %>%
    select(Table, Dependent, term, Estimate, SE, Statistic, P, Stat_type)
}

########################################
## HELPER FUNCTIONS
########################################
format_table <- function(df) {
  df %>%
    mutate(
      P = ifelse(P < 0.001, "<0.001", sprintf("%.4f", P))
    )
}

clean_terms <- function(df) {
  df %>%
    mutate(term = case_when(
      term == "(Intercept)" ~ "intercept",
      grepl("Population", term) ~ gsub("Population", "", term),
      grepl("Site", term) ~ gsub("Site", "", term),
      TRUE ~ term
    ))
}


########################################
########################################
## TABLE 3 — SITE + POPULATION
########################################
########################################
results_T3 <- list()

### RGL
results_T3[[1]] <- extract_model(
  lm(Total ~ Population + Site + Coverage + Batch, data = RGL_site),
  "Total RGL", "Table 3")

results_T3[[2]] <- extract_model(
  lm(Heterozygotes ~ Population + Site + Coverage + Batch, data = RGL_site),
  "Masked RGL", "Table 3")

results_T3[[3]] <- extract_model(
  lm(Homozygotes ~ Population + Site + Coverage + Batch, data = RGL_site),
  "Realized RGL", "Table 3")

### RML (FULL)
results_T3[[4]] <- extract_model(
  glm(cbind(Total_MISS, SynDerivedAlleles) ~ Population + Site + Coverage + Batch,
      data = MISS_site, family = quasibinomial),
  "Total RML", "Table 3")

results_T3[[5]] <- extract_model(
  glm(cbind(Heterozygotes_MISS, SynDerivedAlleles) ~ Population + Site + Coverage + Batch,
      data = MISS_site, family = quasibinomial),
  "Masked RML", "Table 3")

results_T3[[6]] <- extract_model(
  glm(cbind(Homozygotes_MISS, SynDerivedAlleles) ~ Population + Site + Coverage + Batch,
      data = MISS_site, family = quasibinomial),
  "Realized RML", "Table 3")

##RLL should be binomial and glmmTB, I think, however stats look very similar anyway 
### RLL (FULL)
results_T3[[7]] <- extract_model(
  glm(cbind(Total_LOF, SynDerivedAlleles) ~ Population + Site + Coverage + Batch,
      data = LOF_site, family = quasibinomial),
  "Total RLL", "Table 3")

results_T3[[8]] <- extract_model(
  glm(cbind(Heterozygotes_LOF, SynDerivedAlleles) ~ Population + Site + Coverage + Batch,
      data = LOF_site, family = quasibinomial),
  "Masked RLL", "Table 3")

results_T3[[9]] <- extract_model(
  glm(cbind(Homozygotes_LOF, SynDerivedAlleles) ~ Population + Site + Coverage + Batch,
      data = LOF_site, family = quasibinomial),
  "Realized RLL", "Table 3")

table_T3 <- bind_rows(results_T3)

########################################
########################################
## TABLE 4 — BOTTLENECK
########################################
########################################
results_T4 <- list()

### RGL
results_T4[[1]] <- extract_model(
  lm(Total ~ Population + Coverage, data = RGL_bottleneck),
  "Total RGL", "Table 4")

results_T4[[2]] <- extract_model(
  lm(log2(Heterozygotes) ~ Population + Coverage, data = RGL_bottleneck),
  "Masked RGL", "Table 4")

results_T4[[3]] <- extract_model(
  lm(Homozygotes ~ Population + Coverage, data = RGL_bottleneck),
  "Realized RGL", "Table 4")

### RML (FULL)
results_T4[[4]] <- extract_model(
  glmmTMB(cbind(Total_MISS, SynDerivedAlleles) ~ Population + Coverage,
          data = MISS_bottleneck, family = binomial),
  "Total RML", "Table 4")

results_T4[[5]] <- extract_model(
  glmmTMB(cbind(Heterozygotes_MISS, SynDerivedAlleles) ~ Population + Coverage,
          data = MISS_bottleneck, family = binomial),
  "Masked RML", "Table 4")

results_T4[[6]] <- extract_model(
  glmmTMB(cbind(Homozygotes_MISS, SynDerivedAlleles) ~ Population + Coverage,
          data = MISS_bottleneck, family = binomial),
  "Realized RML", "Table 4")

### RLL (FULL)
results_T4[[7]] <- extract_model(
  glmmTMB(cbind(Total_LOF, SynDerivedAlleles) ~ Population + Coverage,
          data = LOF_bottleneck, family = binomial),
  "Total RLL", "Table 4")

results_T4[[8]] <- extract_model(
  glmmTMB(cbind(Heterozygotes_LOF, SynDerivedAlleles) ~ Population + Coverage,
          data = LOF_bottleneck, family = binomial),
  "Masked RLL", "Table 4")

results_T4[[9]] <- extract_model(
  glmmTMB(cbind(Homozygotes_LOF, SynDerivedAlleles) ~ Population + Coverage,
          data = LOF_bottleneck, family = binomial),
  "Realized RLL", "Table 4")

table_T4 <- bind_rows(results_T4)
########################################
########################################
## TABLE 5 — TURURE
########################################
########################################
results_T5 <- list()

### RGL
results_T5[[1]] <- extract_model(
  lm(Total ~ Population + Coverage + Batch, data = RGL_TUR),
  "Total RGL", "Table 5")

results_T5[[2]] <- extract_model(
  lm(Heterozygotes ~ Population + Coverage + Batch, data = RGL_TUR),
  "Masked RGL", "Table 5")

results_T5[[3]] <- extract_model(
  lm(Homozygotes ~ Population + Coverage + Batch, data = RGL_TUR),
  "Realized RGL", "Table 5")

### RML (FULL)
results_T5[[4]] <- extract_model(
  glmmTMB(cbind(Total_MISS, SynDerivedAlleles) ~ Population + Coverage + Batch,
          data = MISS_TUR, family = binomial),
  "Total RML", "Table 5")

results_T5[[5]] <- extract_model(
  glmmTMB(cbind(Heterozygotes_MISS, SynDerivedAlleles) ~ Population + Coverage + Batch,
          data = MISS_TUR, family = binomial),
  "Masked RML", "Table 5")

results_T5[[6]] <- extract_model(
  glmmTMB(cbind(Homozygotes_MISS, SynDerivedAlleles) ~ Population + Coverage + Batch,
          data = MISS_TUR, family = binomial),
  "Realized RML", "Table 5")

### RLL (FULL)
results_T5[[7]] <- extract_model(
  glmmTMB(cbind(Total_LOF, SynDerivedAlleles) ~ Population + Coverage + Batch,
          data = LOF_TUR, family = binomial),
  "Total RLL", "Table 5")

results_T5[[8]] <- extract_model(
  glmmTMB(cbind(Heterozygotes_LOF, SynDerivedAlleles) ~ Population + Coverage + Batch,
          data = LOF_TUR, family = binomial),
  "Masked RLL", "Table 5")

results_T5[[9]] <- extract_model(
  glmmTMB(cbind(Homozygotes_LOF, SynDerivedAlleles) ~ Population + Coverage + Batch,
          data = LOF_TUR, family = binomial),
  "Realized RLL", "Table 5")

table_T5 <- bind_rows(results_T5)

########################################
########################################
## EXPORT TO EXCEL
########################################
########################################
wb <- createWorkbook()

addWorksheet(wb, "Table 3")
writeData(wb, "Table 3", table_T3)

addWorksheet(wb, "Table 4")
writeData(wb, "Table 4", table_T4)

addWorksheet(wb, "Table 5")
writeData(wb, "Table 5", table_T5)

saveWorkbook(wb, "Supplementary_Tables_ALL_UPDATED2.xlsx", overwrite = TRUE)


########################################
########################################
## NEED TO DO THE SAME FOR REGIONS
########################################
########################################
##str()#### REgions
LOF_reg <- read.csv("LOF_region.txt")
MISS_reg <- read.csv("MISS_region.txt")
RGL_reg <- read.csv("RGL_region.txt")

MISS_reg$Region<-factor(MISS_reg$Region, levels=c("Oropouche","Caroni","Tobago"))
LOF_reg$Region<-factor(LOF_reg$Region, levels=c("Oropouche","Caroni","Tobago"))
RGL_reg$Region<-factor(RGL_reg$Region, levels=c("Oropouche","Caroni","Tobago"))

#
results_T1 <- list()

### RGL
results_T1[[1]] <- extract_model(
  glmmTMB(Total~Region+(1|Population)+Coverage, data=RGL_reg, family=gaussian),
  "Total RGL", "Table 1")

results_T1[[2]] <- extract_model(
  glmmTMB(log(Heterozygotes)~Region+(1|Population)+Coverage, data=RGL_reg, family=gaussian),
  "Masked RGL", "Table 1")

results_T1[[3]] <- extract_model(
  glmmTMB(Homozygotes~Region+(1|Population)+Coverage, data=RGL_reg, family=gaussian),
  "Realized RGL", "Table 1")

### RML (FULL)
results_T1[[4]] <- extract_model(
  glmmTMB(cbind(Total_MISS, SynDerivedAlleles) ~ Region+(1|Population)+Coverage,
          data = MISS_reg, family = binomial),
  "Total RML", "Table 1")

results_T1[[5]] <- extract_model(
  glmmTMB(cbind(Heterozygotes_MISS, SynDerivedAlleles) ~ Region+(1|Population)+Coverage,
          data = MISS_reg, family = binomial),
  "Masked RML", "Table 1")

results_T1[[6]] <- extract_model(
  glmmTMB(cbind(Homozygotes_MISS, SynDerivedAlleles) ~ Region+(1|Population)+Coverage,
          data = MISS_reg, family = binomial),
  "Realized RML", "Table 1")

### RLL (FULL)
results_T1[[7]] <- extract_model(
  glmmTMB(cbind(Total_LOF, SynDerivedAlleles) ~ Region+(1|Population)+Coverage,
         data = LOF_reg, family = binomial),
  "Total RLL", "Table 1")

results_T1[[8]] <- extract_model(
  glmmTMB(cbind(Heterozygotes_LOF, SynDerivedAlleles) ~ Region+(1|Population)+Coverage,
          data = LOF_reg, family = binomial),
  "Masked RLL", "Table 1")

results_T1[[9]] <- extract_model(
  glmmTMB(cbind(Homozygotes_LOF, SynDerivedAlleles) ~ Region+(1|Population)+Coverage,
          data = LOF_reg, family = binomial),
  "Realized RLL", "Table 1")

table_T1 <- bind_rows(results_T1)

wb <- createWorkbook()

addWorksheet(wb, "Table 1")
writeData(wb, "Table 1", table_T1)

saveWorkbook(wb, "Supplementary_Tables_ALL_UPDATED_TabS1_v2.xlsx", overwrite = TRUE)