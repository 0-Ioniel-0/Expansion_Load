library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library("lme4")
library(gridExtra)
### Loss Of Function Alleles Analysis ###

### Loading Files
setwd("~/Genetic Load/Data & Stats")
load_Ne <- read.delim("load_Ne.txt")
load_drainage <- read.delim("load_drainage.txt")
load_island <- read.delim("load_island.txt")
load_Turure <- read.delim("load_Turure2.txt")
load_site <- read.delim("load_site.txt")
load_bottleneck <- read.delim("load_Turure_others.txt")

### PART 1 - Effective Population Size

## 1.1 load vs. Ne
# 1.1.1 model
model_load_Ne <- lm(MeanLoad ~ Ne + Region, data = load_Ne)
par(mfrow= c(2,2))
plot(model_load_Ne)
summary(model_load_Ne)
summary.lm(model_load_Ne)
# 1.1.2 plot
plot_1_1 <- ggplot(data = load_Ne, aes(x = Ne, y = MeanLoad, fill = Population)) + theme_classic() +
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
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nEffective Population Size") + ylab("Relative genetic load\n") +
  scale_x_continuous(labels = scales::scientific)

## 1.2 load in homozygotes vs. Ne
# 1.2.1 model
model_load_Ne_ho <- lm(HomozygotesLoad ~ Ne + Region, data = load_Ne)
par(mfrow= c(2,2))
plot(model_load_Ne_ho)
summary(model_load_Ne_ho)
# 1.2.2 plot
plot_1_2_load <- ggplot(data = load_Ne, aes(x = Ne, y = HomozygotesLoad, fill = Population)) +
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
        axis.text = element_text(size = 10)) + xlab("\nEffective Population Size") + ylab("Relative genetic load\n") +
  scale_x_continuous(labels = scales::scientific)

## 1.3 load in heterozygotes vs. Ne
# 1.3.1 model
model_load_Ne_het <- lm(HeterozygotesLoad ~ Ne + Region, data = load_Ne)
par(mfrow= c(2,2))
plot(model_load_Ne_het)
summary(model_load_Ne_het)
# 1.3.2 plot
plot_1_3_load <- ggplot(data = load_Ne, aes(x = Ne, y = HeterozygotesLoad, fill = Population)) +
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
        axis.text = element_text(size = 10)) + xlab("\nEffective Population Size") + ylab("Relative genetic load\n") +
  scale_x_continuous(labels = scales::scientific)

### PART 2 - Island difference

## 2.1 load in islands
# 2.1.1 model
model_load_island <-lmer(Load ~ Island + (1|Population), data = load_island)
model_load_island_2 <-lmer(Load ~ 1 + (1|Population), data = load_island)
anova(model_load_island, model_load_island_2, test="Chi")
summary(model_load_island)
residuals_model_island <- residuals(model_load_island)

plot(model_load_island)
qqnorm(residuals_model_island)
qqline(residuals_model_island, col = "red")
hist(residuals_model_island, breaks = 20)
influence_island <- influence(model_load_island, group = "Population")
cooks.distance(influence_island)

# 2.1.1 plot
plot_2_1 <- ggplot(data = load_island, aes(x = Island, y = Load, fill = Island)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("Trinidad" = "gold2", "Tobago" = "#0c2f55")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) + xlab("\nIsland") + ylab("Relative genetic load\n")

## 2.2 load in homozygotes in islands
# 2.2.1 model
model_load_island_ho <-lmer(Homozygotes ~ Island + (1|Population), data = load_island)
model_load_island_ho_2 <-lmer(Homozygotes ~ 1 + (1|Population), data = load_island)
anova(model_load_island_ho, model_load_island_ho_2, test="Chi")
summary(model_load_island_ho)
residuals_model_island_ho <- residuals(model_load_island_ho)

plot(model_load_island_ho)
qqnorm(residuals_model_island_ho)
qqline(residuals_model_island_ho, col = "red")
hist(residuals_model_island_ho, breaks = 20)
influence_island_ho <- influence(model_load_island_ho, group = "Population")
cooks.distance(influence_island_ho)

# 2.2.2 plot
plot_2_2 <- ggplot(data = load_island, aes(x = Island, y = Homozygotes, fill = Island)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("Trinidad" = "gold2", "Tobago" = "#0c2f55")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nIsland") + ylab("Relative genetic load\n") + ggtitle("Homozygotes")

## 2.3 load in heterozygotes in islands
# 2.3.1 model
model_load_island_he <-lmer(Heterozygotes ~ Island + (1|Population), data = load_island)
model_load_island_he_2 <-lmer(Heterozygotes ~ 1 + (1|Population), data = load_island)
anova(model_load_island_he, model_load_island_he_2, test="Chi")
summary(model_load_island_he)
residuals_model_island_he <- residuals(model_load_island_he)

plot(model_load_island_he)
qqnorm(residuals_model_island_he)
qqline(residuals_model_island_he, col = "red")
hist(residuals_model_island_he, breaks = 20)
influence_island <- influence(model_load_island_he, group = "Population")
cooks.distance(influence_island)

# 2.3.2 plot
plot_2_3 <- ggplot(data = load_island, aes(x = Island, y = Heterozygotes, fill = Island)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("Trinidad" = "gold2", "Tobago" = "#0c2f55")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nIsland") + ylab("Relative genetic load\n") + ggtitle("Heterozygotes")

### PART 3 - Drainage difference
## 3.1 load in drainages
# 3.1.1 model
model_load_drainage <-lmer(Load ~ Drainage + (1|Population), data = load_drainage)
model_load_drainage_2 <-lmer(Load ~ 1 + (1|Population), data = load_drainage)
anova(model_load_drainage, model_load_drainage_2, test="Chi")
summary(model_load_drainage)
residuals_model_drainage <- residuals(model_load_drainage)

plot(model_load_drainage)
qqnorm(residuals_model_drainage)
qqline(residuals_model_drainage, col = "red")
hist(residuals_model_drainage, breaks = 20)
influence_drainage <- influence(model_load_drainage, group = "Population")
cooks.distance(influence_drainage)

# 3.1.2 plot
plot_3_1 <- ggplot(data = load_drainage, aes(x = Drainage, y = Load, fill = Drainage)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("Caroni" = "#00a5ca",
                               "Oropouche" = "#6fbb86")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) + xlab("\nDrainage") + ylab("Relative genetic load\n")

## 3.2 load in homozygotes in drainages
# 3.2.1 model
model_load_drainage_ho <-lmer(Homozygotes ~ Drainage + (1|Population), data = load_drainage)
model_load_drainage_ho_2 <-lmer(Homozygotes ~ 1 + (1|Population), data = load_drainage)
anova(model_load_drainage_ho, model_load_drainage_ho_2, test="Chi")
summary(model_load_drainage_ho)
residuals_model_drainage_ho <- residuals(model_load_drainage_ho)

plot(model_load_drainage_ho)
qqnorm(residuals_model_drainage_ho)
qqline(residuals_model_drainage, col = "red")
hist(residuals_model_drainage, breaks = 20)
influence_island <- influence(model_load_island, group = "Population")
cooks.distance(influence_island)

# 3.2.2 plot
plot_3_2 <- ggplot(data = load_drainage, aes(x = Drainage, y = HomozygotesLoad, fill = Drainage)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("Caroni" = "#00a5ca",
                               "Oropouche" = "#6fbb86")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nDrainage") + ylab("Relative genetic load\n") + ggtitle("Homozygotes")

## 3.3 load in heterozygotes in drainages
# 3.3.1 model
model_load_drainage_he <-lmer(Heterozygotes ~ Drainage + (1|Population), data = load_drainage)
model_load_drainage_he_2 <-lmer(Heterozygotes ~ 1 + (1|Population), data = load_drainage)
anova(model_load_drainage_he, model_load_drainage_he_2, test="Chi")
summary(model_load_drainage_he)
residuals_model_drainage_he <- residuals(model_load_drainage_he)

plot(model_load_drainage_he)
qqnorm(residuals_model_drainage_he)
qqline(residuals_model_drainage, col = "red")
hist(residuals_model_drainage, breaks = 20)
influence_island_he <- influence(model_load_island_he, group = "Population")
cooks.distance(influence_island_he)

# 3.3.2 plot
plot_3_3 <- ggplot(data = load_drainage, aes(x = Drainage, y = HeterozygotesLoad, fill = Drainage)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("Caroni" = "#00a5ca",
                               "Oropouche" = "#6fbb86")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nDrainage") + ylab("Relative genetic load\n") + ggtitle("Heterozygotes")

### PART 4 - Site difference
## 4.1 load in sites
# 4.1.1 model
load_site_O <- load_site %>%
  filter(Population != "TUR")
model_load_site <- aov(Load ~ Site * Population, data = load_site_O)
summary(model_load_site)
coef(model_load_site)
summary.lm(model_load_site)
par(mfrow=c(2,2))
plot(model_load_site)
# 4.1.2 plot
plot_4_1 <- ggplot(data = load_site, aes(x = Site, y = Load, fill = Site)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("upper" = "mediumorchid3", "lower" = "black")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) + xlab("\nSite") + ylab("Relative genetic load\n")

## 4.2 load in homozygotes in sites
# 4.2.1 model
model_load_site_ho <- aov(Homozygotes ~ Site * Population, data = load_site_O)
summary(model_load_site_ho)
summary.lm(model_load_site_ho)
coef(model_load_site_ho)
par(mfrow=c(2,2))
plot(model_load_site_ho)
# 4.2.2 plot
plot_4_2 <- ggplot(data = load_site, aes(x = Site, y = Homozygotes, fill = Site)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("upper" = "mediumorchid3", "lower" = "black")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nSite") + ylab("Relative genetic load\n") + ggtitle("Homozygotes")

## 4.3 load in heterozygotes in sites
# 4.3.1 model
model_load_site_he <- aov(Heterozygotes ~ Site * Population, data = load_site_O)
summary(model_load_site_he)
summary.lm(model_load_site_he)
coef(model_load_site_he)
par(mfrow=c(2,2))
plot(model_load_site_he)
# 4.3.2 plot
plot_4_3 <- ggplot(data = load_site, aes(x = Site, y = Heterozygotes, fill = Site)) +
  geom_boxplot(alpha = .4) +
  scale_fill_manual(values = c("upper" = "mediumorchid3", "lower" = "black")) + theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nSite") + ylab("Relative genetic load\n") + ggtitle("Heterozygotes")

### PART 5 - Turure differences
## 5.1 load in Turure
# 5.1.1 model
model_load_Turure <- aov(Load ~ Population, data = load_Turure)
summary(model_load_Turure)
summary.lm(model_load_Turure)
coef(model_load_Turure)
pairwise.t.test(load_Turure$Load, load_Turure$Population, p.adjust.method = "bonferroni")
par(mfrow=c(2,2))
plot(model_load_Turure)
emmeans_load_Turure <- emmeans(model_load_Turure, ~ Population)
custom_contrasts <- list(
  "1 vs 2" = c(-1, 1, 0),   # Contrast for 4 vs 6 cylinders
  "1 vs 3" = c(-1, 0, 1),   # Contrast for 4 vs 8 cylinders
  "2 vs 3" = c(0, -1, 1)    # Contrast for 6 vs 8 cylinders
)
contrast(emmeans_load_Turure, custom_contrasts)

# 5.1.2 plot
plot_5_1 <- ggplot(data = load_Turure, aes(x = Population, y = Load)) +
  geom_boxplot(fill = "#fd1f4b", alpha = .4) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) + xlab("\nTurure Site") + ylab("Relative genetic load\n")

## 5.2 load in homozygotes in Turure
# 5.2.1 model
model_load_Turure_ho <- aov(sqrt(Homozygotes) ~ Population, data = load_Turure)
summary(model_load_Turure_ho)
summary.lm(model_load_Turure_ho)
coef(model_load_Turure_ho)
pairwise.t.test(load_Turure$Homozygotes, load_Turure$Population, p.adjust.method = "bonferroni")
par(mfrow=c(2,2))
plot(model_load_Turure_ho)
emmeans_load_Turure_ho <- emmeans(model_load_Turure_ho, ~ Population)
custom_contrasts <- list(
  "TUR_lo vs TUR_mid_up" = c(-1, 1, 0),   # Contrast for 4 vs 6 cylinders
  "TUR_lo vs TUR_up" = c(-1, 0, 1),   # Contrast for 4 vs 8 cylinders
  "TUR_mid_up vs TUR_up" = c(0, -1, 1)    # Contrast for 6 vs 8 cylinders
)
contrast(emmeans_load_Turure_ho, custom_contrasts)

# 5.2.2 plot
plot_5_2 <- ggplot(data = load_Turure, aes(x = Population, y = Homozygotes)) +
  geom_boxplot(fill = "#fd1f4b", alpha = .4) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nPopulation") + ylab("Relative genetic load\n") + ggtitle("Homozygotes")
## 5.3 load in heterozygotes in Turure
# 5.3.1 model
model_load_Turure_he <- aov(1/Heterozygotes ~ Population, data = load_Turure)
summary(model_load_Turure_he)
summary.lm(model_load_Turure_he)
coef(model_load_Turure_he)
par(mfrow=c(2,2))
plot(model_load_Turure_he)
emmeans_load_Turure_he <- emmeans(model_load_Turure_he, ~ Population)
custom_contrasts <- list(
  "TUR_lo vs TUR_mid_up" = c(-1, 1, 0),   # Contrast for 4 vs 6 cylinders
  "TUR_lo vs TUR_up" = c(-1, 0, 1),   # Contrast for 4 vs 8 cylinders
  "TUR_mid_up vs TUR_up" = c(0, -1, 1)    # Contrast for 6 vs 8 cylinders
)
contrast(emmeans_load_Turure_he, custom_contrasts)
# 5.3.2 plot
plot_5_3 <- ggplot(data = load_Turure, aes(x = Population, y = Heterozygotes)) +
  geom_boxplot(fill = "#fd1f4b", alpha = .4) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) + xlab("\nPopulation") + ylab("Relative genetic load\n") + ggtitle("Heterozygotes")

### PART 6 - Turure vs. lower locations
## 6.1 load in Turure vs. lower locations
# 6.1.1 model
model_6_1 <- aov((Load) ~ State * Population, data = load_bottleneck)
summary(model_6_1)
summary.lm(model_6_1)
coef(model_6_1)
par(mfrow=c(2,2))
plot(model_6_1)
# 6.1.2 plot
plot_6_1 <- ggplot(data = load_bottleneck, aes(x = State, y = Load, color = Population)) +
  geom_boxplot(alpha = .4) +
  scale_color_manual(values = c("QUA" = "#6fbb86",
                                "ORO" = "#6fbb86",
                                "TUR" = "#fd1f4b")) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nState") + ylab("Relative genetic load\n")

## 6.2 Homozygotes load in Turure vs. lower locations
# 6.2.1 model
model_6_2 <- aov(sqrt(Homozygotes) ~ State * Population, data = load_bottleneck)
summary(model_6_2)
summary.lm(model_6_2)
coef(model_6_2)
par(mfrow=c(2,2))
plot(model_6_2)
# 6.2.2 plot
plot_6_2 <- ggplot(data = load_bottleneck, aes(x = State, y = Homozygotes, color = Population)) +
  geom_boxplot(alpha = .4) +
  scale_color_manual(values = c("QUA" = "#6fbb86",
                                "ORO" = "#6fbb86",
                                "TUR" = "#fd1f4b")) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nState") + ylab("Relative genetic load\n")

## 6.3 Heterozygotes load in Turure vs. lower locations
# 6.3.1 model
model_6_3 <- aov((Heterozygotes)^2 ~ State * Population, data = load_bottleneck)
summary(model_6_3)
summary.lm(model_6_3)
coef(model_6_3)
par(mfrow=c(2,2))
plot(model_6_3)
# 6.3.2 plot
plot_6_3 <- ggplot(data = load_bottleneck, aes(x = State, y = Heterozygotes, color = Population)) +
  geom_boxplot(alpha = .4) +
  scale_color_manual(values = c("QUA" = "#6fbb86",
                                "ORO" = "#6fbb86",
                                "TUR" = "#fd1f4b")) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nState") + ylab("Relative genetic load\n")

### Finishing
### Finishing
# Print in a grid

load_island_long <- load_island %>% pivot_longer(cols=c('Load', 'Heterozygotes', 'Homozygotes'),
                                               names_to='Genotype',
                                               values_to='Relative_Load')
load_drainage_long <- load_drainage %>% pivot_longer(cols=c('Load', 'Heterozygotes', 'Homozygotes'),
                                               names_to='Genotype',
                                               values_to='Relative_Load')
load_site_long <- load_site_O %>% pivot_longer(cols=c('Load', 'Heterozygotes', 'Homozygotes'),
                                             names_to='Genotype',
                                             values_to='Relative_Load')  
load_Turure_long <- load_Turure %>% pivot_longer(cols=c('Load', 'Heterozygotes', 'Homozygotes'),
                                                 names_to='Genotype',
                                                 values_to='Relative_Load')  
load_bottleneck_long <- load_bottleneck %>% pivot_longer(cols=c('Load', 'Heterozygotes', 'Homozygotes'),
                                                 names_to='Genotype',
                                                 values_to='Relative_Load')  

lgt.labs <- c("Heterozygotes", "Homozygotes", "Total")
names(lgt.labs) <- c("Heterozygotes", "Homozygotes", "Load")
blgt.labs <- c("no", "yes")
names(blgt.labs) <- c("pre-bottleneck", "post-bottleneck")
plot_I_lo <- ggplot() + geom_boxplot(data = load_island_long, aes(x = Island, y = Relative_Load, alpha = Genotype)) +
  facet_wrap(~Genotype, labeller = labeller(Genotype = lgt.labs)) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nIsland") + ylab("Relative genetic load\n")

plot_D_lo <- ggplot() +  geom_boxplot(data = load_drainage_long, aes(x = Drainage, y = Relative_Load, alpha = Genotype)) +
  facet_wrap(~Genotype, labeller = labeller(Genotype = lgt.labs)) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nDrainage") + ylab("Relative genetic load\n")

plot_S_lo <- ggplot() +  geom_boxplot(data = load_site_long, aes(x = Site, y = Relative_Load, color = Population)) +
  facet_wrap(~Genotype, labeller = labeller(Genotype = lgt.labs)) +
  theme_classic() +
  scale_color_manual(values = c("QUA" = "black",
                               "ORO" = "black")) +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nLocation") + ylab("Relative genetic load\n")
plot_T_lo <- ggplot() +  geom_boxplot(data = load_Turure_long, aes(x = Population, y = Relative_Load, alpha = Genotype)) +
  facet_wrap(~Genotype, labeller = labeller(Genotype = lgt.labs)) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nTurure Location") + ylab("Relative genetic load\n")
plot_B_lo <- ggplot() +  geom_boxplot(data = load_bottleneck_long, aes(x = State, y = Relative_Load, color = Population)) +
  facet_wrap(~Genotype, labeller = labeller(Genotype = lgt.labs)) +
  theme_classic() +
  scale_color_manual(values = c("QUA" = "black",
                                "ORO" = "black",
                                "TUR" = "black")) +
  theme(legend.position = "none", axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + xlab("\nFounder effect") + ylab("Relative genetic load\n") +
  scale_x_discrete(labels = c("yes", "no"))

grid.arrange(plot_I_lo, plot_I,
             plot_S_lo, plot_S,
             plot_T_lo, plot_T,
             plot_B_lo, plot_B,
             ncol = 2, nrow = 4)

pdf("CS-load.pdf")
