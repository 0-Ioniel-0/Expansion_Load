library(ggplot2)
library(emmeans)
library(gridExtra)
### Loss Of Function Alleles Analysis ###

### Loading Files
setwd("~/Genetic Load/Data & Stats")
LOF_Ne <- read.delim("LOF_Ne.txt")
LOF_drainage <- read.delim("LOF_drainage.txt")
LOF_island <- read.delim("LOF_island.txt")
LOF_Turure <- read.delim("LOF_Turure.txt")
LOF_site <- read.delim("LOF_site.txt")
LOF_bottleneck <- read.delim("LOF_Turure_others.txt")

### PART 1 - Effective Population Size

## 1.1 LOF vs. Ne
# 1.1.1 model
model_LOF_Ne <- lm(LOF ~ Ne + Region, data = LOF_Ne)
par(mfrow= c(2,2))
plot(model_LOF_Ne)
summary(model_LOF_Ne)
# 1.1.2 plot
ggplot(data = LOF_Ne, aes(x = Ne, y = LOF, fill = Population)) + theme_classic() +
  scale_fill_manual(values = c("DR" = "#0c2f55",
                               "ROX" = "#0c2f55",
                               "ARI" = "#00a5ca",
                               "CAU" = "#00a5ca",
                               "LOP" = "#00a5ca",
                               "SCR" = "#00a5ca",
                               "LSA" = "#6fbb86",
                               "QUA" = "#6fbb86",
                               "ORO" = "#6fbb86")) +
  geom_label(aes(label = Population, size = NULL), color = "white") +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nEffective Population Size") + ylab("Normalized count of LOF alleles\n") +
  scale_x_continuous(labels = scales::scientific)

## 1.2 LOF in homozygotes vs. Ne
# 1.2.1 model
model_LOF_Ne_ho <- lm(HomozygotesLOF ~ Ne + Region, data = LOF_Ne)
par(mfrow= c(2,2))
plot(model_LOF_Ne_ho)
summary(model_LOF_Ne_ho)
# 1.2.2 plot
plot_1_2_LOF <- ggplot(data = LOF_Ne, aes(x = Ne, y = HomozygotesLOF, fill = Population)) +
  theme_classic() +
  scale_fill_manual(values = c("DR" = "#0c2f55",
                               "ROX" = "#0c2f55",
                               "ARI" = "#00a5ca",
                               "CAU" = "#00a5ca",
                               "LOP" = "#00a5ca",
                               "SCR" = "#00a5ca",
                               "LSA" = "#6fbb86",
                               "QUA" = "#6fbb86",
                               "ORO" = "#6fbb86")) +
  geom_label(aes(label = Population, size = NULL), color = "white") +
  theme(legend.position = "none", axis.title.y = element_text(size = 14, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nEffective Population Size") + ylab("Normalized count of LOF alleles\n") +
  scale_x_continuous(labels = scales::scientific)

## 1.3 LOF in heterozygotes vs. Ne
# 1.3.1 model
model_LOF_Ne_het <- lm(HeterozygotesLOF ~ Ne + Region, data = LOF_Ne)
par(mfrow= c(2,2))
plot(model_LOF_Ne_het)
summary(model_LOF_Ne_het)
# 1.3.2 plot
plot_1_3_LOF <- ggplot(data = LOF_Ne, aes(x = Ne, y = HeterozygotesLOF, fill = Population)) +
  theme_classic() +
  scale_fill_manual(values = c("DR" = "#0c2f55",
                               "ROX" = "#0c2f55",
                               "ARI" = "#00a5ca",
                               "CAU" = "#00a5ca",
                               "LOP" = "#00a5ca",
                               "SCR" = "#00a5ca",
                               "LSA" = "#6fbb86",
                               "QUA" = "#6fbb86",
                               "ORO" = "#6fbb86")) +
  geom_label(aes(label = Population, size = NULL), color = "white") +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nEffective Population Size") + ylab("Normalized count of LOF alleles\n") +
  scale_x_continuous(labels = scales::scientific)

### PART 2 - Island difference

## 2.1 LOF in islands
# 2.1.1 model
model_LOF_island <- glm(cbind(TotalLOF, TotalDerivedAlleles) ~ Island, data = LOF_island, family = "binomial")
summary(model_LOF_island)
sum(residuals(model_LOF_island, type = "pearson")^2) / model_LOF_island$df.residual
# 2.1.1 plot
plot_2_1_L <- ggplot(data = LOF_island, aes(x = Island, y = LOF, fill = Island)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("Trinidad" = "gold2", "Tobago" = "#0c2f55")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nIsland") + ylab("Normalized LOF alleles\n") #+ ("All LOF Alleles")

## 2.2 LOF in homozygotes in islands
# 2.2.1 model
model_LOF_island_ho <- glm(cbind(TotalHomo, TotalDerivedAlleles) ~ Island, data = LOF_island, family = quasibinomial())
summary(model_LOF_island_ho)
sum(residuals(model_LOF_island_ho, type = "pearson")^2) / model_LOF_island_ho$df.residual

# 2.2.2 plot
plot_2_2 <- ggplot(data = LOF_island, aes(x = Island, y = HomoLOF, fill = Island)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("Triniad" = "darkgray",
                               "Tobago" = "#0c2f55")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nIsland") + ylab("Normalized LOF alleles\n") #+ ("Homozygotes")

## 2.3 LOF in heterozygotes in islands
# 2.3.1 model
model_LOF_island_he <- glm(cbind(TotalHetero, TotalDerivedAlleles) ~ Island, data = LOF_island, family = "binomial")
summary(model_LOF_island_he)
sum(residuals(model_LOF_island_he, type = "pearson")^2) / model_LOF_island_he$df.residual

# 2.3.2 plot
plot_2_3 <- ggplot(data = LOF_island, aes(x = Island, y = HeteroLOF, fill = Island)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("Triniad" = "darkgray",
                               "Tobago" = "#0c2f55")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nIsland") + ylab("Normalized LOF alleles\n") #+ ("Heterozygotes")

### PART 3 - Drainage difference
## 3.1 LOF in drainages
# 3.1.1 model
model_LOF_drainage <- glm(cbind(TotalLOF, TotalDerivedAlleles) ~ Drainage, data = LOF_drainage, family = "binomial")
summary(model_LOF_drainage)
sum(residuals(model_LOF_drainage, type = "pearson")^2) / model_LOF_drainage$df.residual

# 3.1.2 plot
plot_3_1_L <- ggplot(data = LOF_drainage, aes(x = Drainage, y = LOF, fill = Drainage)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("Caroni" = "#00a5ca",
                               "Oropouche" = "#6fbb86")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nDrainage") + ylab("Normalized LOF alleles\n") #+ ("All LOF Alleles")

## 3.2 LOF in homozygotes in drainages
# 3.2.1 model
model_LOF_drainage_ho <- glm(cbind(TotalHomo, TotalDerivedAlleles) ~ Drainage, data = LOF_drainage, family = quasibinomial())
summary(model_LOF_drainage_ho)
sum(residuals(model_LOF_drainage_ho, type = "pearson")^2) / model_LOF_drainage_ho$df.residual

# 3.2.2 plot
plot_3_2 <- ggplot(data = LOF_drainage, aes(x = Drainage, y = HomoLOF, fill = Drainage)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("Caroni" = "#00a5ca",
                               "Oropouche" = "#6fbb86")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nDrainage") + ylab("Normalized LOF alleles\n") #+ ("Homozygotes")

## 3.3 LOF in heterozygotes in drainages
# 3.3.1 model
model_LOF_drainage_he <- glm(cbind(TotalHetero, TotalDerivedAlleles) ~ Drainage, data = LOF_drainage, family = "binomial")
summary(model_LOF_drainage_he)
sum(residuals(model_LOF_drainage_he, type = "pearson")^2) / model_LOF_drainage_he$df.residual

# 3.3.2 plot
plot_3_3 <- ggplot(data = LOF_drainage, aes(x = Drainage, y = HeteroLOF, fill = Drainage)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("Caroni" = "#00a5ca",
                               "Oropouche" = "#6fbb86")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nDrainage") + ylab("Normalized LOF alleles\n") #+ ("Heterozygotes")

### PART 4 - Site difference
## 4.1 LOF in sites
# 4.1.1 model
model_LOF_site <- glm(cbind(TotalLOF, TotalDerivedAlleles) ~ Site, data = LOF_site_O, family = quasibinomial())
summary(model_LOF_site)
sum(residuals(model_LOF_site, type = "pearson")^2) / model_LOF_site$df.residual
LOF_site_O <- LOF_site %>%
  filter(Population != "TUR")
# 4.1.2 plot
plot_4_1_L <- ggplot(data = LOF_site, aes(x = Site, y = LOF, fill = Site)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("upper" = "mediumorchid3", "lower" = "black")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nSite") + ylab("Normalized LOF alleles\n") #+ ("All LOF Alleles")

## 4.2 LOF in homozygotes in sites
# 4.2.1 model
model_LOF_site_ho <- glm(cbind(TotalHomo, TotalDerivedAlleles) ~ Site, data = LOF_site_O, family = quasibinomial())
summary(model_LOF_site_ho)
sum(residuals(model_LOF_site_ho, type = "pearson")^2) / model_LOF_site_ho$df.residual

# 4.2.2 plot
plot_4_2 <- ggplot(data = LOF_site, aes(x = Site, y = HomoLOF, fill = Site)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("upper" = "white",
                               "lower" = "black")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nSite") + ylab("Normalized LOF alleles\n") #+ ("Homozygotes")

## 4.3 LOF in heterozygotes in sites
# 4.3.1 model
model_LOF_site_he <- glm(cbind(TotalHetero, TotalDerivedAlleles) ~ Site, data = LOF_site_O, family = "binomial")
summary(model_LOF_site_he)
sum(residuals(model_LOF_site_he, type = "pearson")^2) / model_LOF_site_he$df.residual

# 4.3.2 plot
plot_4_3 <- ggplot(data = LOF_site, aes(x = Site, y = HeteroLOF, fill = Site)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("upper" = "white",
                               "lower" = "black")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nSite") + ylab("Normalized LOF alleles\n") #+ ("Heterozygotes")

### PART 5 - Turure differences
## 5.1 LOF in Turure
# 5.1.1 model
model_LOF_Turure <- glm(cbind(TotalLOF, TotalDerivedAlleles) ~ Population, data = LOF_Turure, family = quasibinomial())
summary(model_LOF_Turure)
emmeans_results <- emmeans(model_LOF_Turure, ~ Population)
pairwise_comparisons <- pairs(emmeans_results)
print(pairwise_comparisons)
sum(residuals(model_LOF_Turure, type = "pearson")^2) / model_LOF_Turure$df.residual

# 5.1.2 plot
plot_5_1_L <- ggplot(data = LOF_Turure, aes(x = Population, y = LOF)) +
  geom_boxplot(fill = "#fd1f4b", alpha = .4) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nTurure Site") + ylab("Normalized LOF alleles\n") #+ ("All LOF Alleles")

## 5.2 LOF in homozygotes in Turure
# 5.2.1 model
model_LOF_Turure_ho <- glm(cbind(TotalHomo, TotalDerivedAlleles) ~ Population, data = LOF_Turure, family = quasibinomial())
summary(model_LOF_Turure_ho)
emmeans_results <- emmeans(model_LOF_Turure_ho, ~ Population)
pairwise_comparisons <- pairs(emmeans_results)
print(pairwise_comparisons)
sum(residuals(model_LOF_Turure_ho, type = "pearson")^2) / model_LOF_Turure_ho$df.residual

# 5.2.2 plot
plot_5_2 <- ggplot(data = LOF_Turure, aes(x = Population, y = HomoLOF)) +
  geom_boxplot(fill = "#fd1f4b", alpha = .4) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nPopulation") + ylab("Normalized LOF alleles\n") #+ ("Homozygotes")

## 5.3 LOF in heterozygotes in Turure
# 5.3.1 model
model_LOF_Turure_he <- glm(cbind(TotalHetero, TotalDerivedAlleles) ~ Population, data = LOF_Turure, family = quasibinomial())
summary(model_LOF_Turure_he)
plot(model_LOF_Turure_he)
emmeans_results <- emmeans(model_LOF_Turure_he, ~ Population)
pairwise_comparisons <- pairs(emmeans_results)
print(pairwise_comparisons)
sum(residuals(model_LOF_Turure_he, type = "pearson")^2) / model_LOF_Turure_he$df.residual

# 5.3.2 plot
plot_5_3 <- ggplot(data = LOF_Turure, aes(x = Population, y = HeteroLOF)) +
  geom_boxplot(fill = "#fd1f4b", alpha = .4) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nPopulation") + ylab("Normalized LOF alleles\n") #+ ("Heterozygotes")

### PART 6 - Turure vs. lower locations
## 6.1 load in Turure vs. lower locations
# 6.1.1 model
model_6_1_L <- glm(cbind(TotalLOF, TotalDerivedAlleles) ~ State, data = LOF_bottleneck, family = quasibinomial())
summary(model_6_1_L)
plot(model_LOF_Turure)
plot(model_6_1_L)
sum(residuals(model_6_1_L, type = "pearson")^2) / model_6_1_L$df.residual
# 6.1.2 plot
plot_6_1_L <- ggplot(data = LOF_bottleneck, aes(x = State, y = LOF, color = Population)) +
  geom_boxplot(alpha = .4) +
  scale_color_manual(values = c("QUA" = "#6fbb86",
                                "ORO" = "#6fbb86",
                                "TUR" = "#fd1f4b")) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nState") + ylab("Normalized LOF alleles\n")

## 6.2 Homozygotes load in Turure vs. lower locations
# 6.2.1 model
model_6_2_L <- glm(cbind(HomoLOF, TotalDerivedAlleles) ~ State, data = LOF_bottleneck, family = quasibinomial())
summary(model_6_2_L)
sum(residuals(model_6_2_L, type = "pearson")^2) / model_6_2_L$df.residual
# 6.1.2 plot
plot_6_2_L <- ggplot(data = LOF_bottleneck, aes(x = State, y = HomoLOF, color = Population)) +
  geom_boxplot(alpha = .4) +
  scale_color_manual(values = c("QUA" = "#6fbb86",
                                "ORO" = "#6fbb86",
                                "TUR" = "#fd1f4b")) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nState") + ylab("Normalized LOF alleles\n")
## 6.3 Heterozygotes load in Turure vs. lower locations
# 6.3.1 model
model_6_3_L <- glm(cbind(HeteroLOF, TotalDerivedAlleles) ~ State, data = LOF_bottleneck, family = quasibinomial())
summary(model_6_3_L)
sum(residuals(model_6_3_L, type = "pearson")^2) / model_6_3_L$df.residual
# 6.1.2 plot
plot_6_3_L <- ggplot(data = LOF_bottleneck, aes(x = State, y = HeteroLOF, color = Population)) +
  geom_boxplot(alpha = .4) +
  scale_color_manual(values = c("QUA" = "#6fbb86",
                                "ORO" = "#6fbb86",
                                "TUR" = "#fd1f4b")) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nState") + ylab("Normalized LOF alleles\n")
### Finishing
# Print in a grid

LOF_island_long <- LOF_island %>% pivot_longer(cols=c('LOF', 'HeteroLOF', 'HomoLOF'),
                                       names_to='Genotype',
                                       values_to='Normalized_Count')
LOF_drainage_long <- LOF_drainage %>% pivot_longer(cols=c('LOF', 'HeteroLOF', 'HomoLOF'),
                                                 names_to='Genotype',
                                                 values_to='Normalized_Count')
LOF_site_long <- LOF_site_O %>% pivot_longer(cols=c('LOF', 'HeteroLOF', 'HomoLOF'),
                                             names_to='Genotype',
                                             values_to='Normalized_Count')
LOF_Turure_long <- LOF_Turure %>% pivot_longer(cols=c('LOF', 'HeteroLOF', 'HomoLOF'),
                                               names_to='Genotype',
                                               values_to='Normalized_Count')
LOF_bottleneck_long <- LOF_bottleneck %>% pivot_longer(cols=c('LOF', 'HeteroLOF', 'HomoLOF'),
                                               names_to='Genotype',
                                               values_to='Normalized_Count')
gt.labs <- c("Heterozygotes", "Homozygotes", "Total")
names(gt.labs) <- c("HeteroLOF", "HomoLOF", "LOF")

plot_I <- ggplot() + geom_boxplot(data = LOF_island_long, aes(x = Island, y = Normalized_Count, alpha = Genotype)) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nIsland") + ylab("Normalized LOF alleles\n") +
  facet_wrap(~Genotype, labeller = labeller(Genotype = gt.labs))
plot_D <- ggplot() +  geom_boxplot(data = LOF_drainage_long, aes(x = Drainage, y = Normalized_Count, alpha = Genotype)) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nDrainage") + ylab("Normalized LOF alleles\n") + facet_wrap(~Genotype, labeller = labeller(Genotype = gt.labs))
plot_S <- ggplot() +  geom_boxplot(data = LOF_site_long, aes(x = Site, y = Normalized_Count, color = Population)) +
  theme_classic() +
  scale_color_manual(values = c("QUA" = "black",
                                "ORO" = "black")) +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nLocation") + ylab("Normalized LOF alleles\n") + facet_wrap(~Genotype, labeller = labeller(Genotype = gt.labs))
plot_T <- ggplot() +  geom_boxplot(data = LOF_Turure_long, aes(x = Population, y = Normalized_Count, alpha = Genotype)) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nTurure Location") + ylab("Normalized LOF alleles\n") + facet_wrap(~Genotype, labeller = labeller(Genotype = gt.labs))
plot_B <- ggplot() +  geom_boxplot(data = LOF_bottleneck_long, aes(x = State, y = Normalized_Count, color = Population)) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nFounder effect") + ylab("Normalized LOF alleles\n") + facet_wrap(~Genotype, labeller = labeller(Genotype = gt.labs)) +
  scale_x_discrete(labels = c("yes", "no")) +
  scale_color_manual(values = c("QUA" = "black",
                                "ORO" = "black",
                                "TUR" = "black"))

grid.arrange(plot_I, plot_D, plot_S, plot_T, ncol = 2)

pdf("LOF.pdf")
