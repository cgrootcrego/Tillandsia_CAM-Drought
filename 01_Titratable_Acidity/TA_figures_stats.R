# 31.08.2023
# This script generates figures and statistics for the titratable acidity section of the CAM drought paper
pacman::p_load("ggplot2", "dplyr", "forcats","grid", "gridExtra", "broom.mixed", "tidyr", "coin", "grDevices", "car", "lme4", "lmtest", "ggsignif", "Cairo")

dat <- read.delim('Titratable_acidity_Tio_Tlei.txt', sep = '\t', header = T)
# Add order for plotting timepoints correctly
dat$order <- rep(c(rep(4,3),rep(1,3),rep(2,3),rep(3,3)),4)

dat$group <- factor(
  interaction(dat$species, dat$watering_regime),
  levels = unique(interaction(dat$species, dat$watering_regime))
)

# Calculate significance: are any of the timepoints significantly different from each other within a species and watering regime?
sign_within <- data.frame()
# Iterate over each group and perform pairwise Wilcoxon signed-rank tests
for (group_df in split(dat, f = interaction(dat$species, dat$watering_regime))) {
  group <- unique(group_df$group)
  test_result <- pairwise.wilcox.test(
    group_df$titratable_acid, group_df$timepoint,
    paired = TRUE,
    p.adjust.method = "none"  # Adjust p-values for multiple comparisons
  )
  pval <- broom::tidy(test_result)
  pval$group <- group
  pval <- as.data.frame(pval)
  sign_within <- rbind(sign_within, pval)
}

medians <- dat %>%
  group_by(timepoint, watering_regime, species)  %>%
  summarize(median_value = median(titratable_acid))

# Calculata significance using a mixed-effects linear model
dat$timepoint <- as.factor(dat$timepoint)
dat$watering_regime <- as.factor(dat$watering_regime)
dat$species <- as.factor(dat$species)
dat$replicate <- as.factor(dat$replicate)

model <- lmer(titratable_acid ~ species * watering_regime * timepoint + (1 | replicate), data = dat)
summary(model)

# Check if a model without interactions fits better. However, the anova results give a significantly lower AIC and BIC and a significant very low p-value, so the model with interations fits much better.
model_simplified <- lmer(titratable_acid ~ species + watering_regime + timepoint + (1 | replicate), data = dat)
anova(model, model_simplified)
# Data: dat
# Models:
#   model_simplified: titratable_acid ~ species + watering_regime + timepoint + (1 | replicate)
# model: titratable_acid ~ species * watering_regime * timepoint + (1 | replicate)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# model_simplified    8 323.82 338.79 -153.91   307.82
# model              18 297.72 331.40 -130.86   261.72 46.101 10  1.375e-06 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
par(mfrow=c(1,2))

# Check linearity of the model
plot(fitted(model), resid(model),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red")
# The plot shows clusters, which shouldn't be.

# Check normality of residuals
qqnorm(resid(model))
qqline(resid(model), col = "red")

shapiro.test(resid(model))
# at the extremes the datapoints diverge from the red line, suggesting there is an issue with normality

# Try general linear models with gamma distribution
model_gamma <- glmer(titratable_acid ~ species + watering_regime + timepoint + (1 | replicate), data = dat, family = Gamma(link = "log"))

resids <- residuals(model_gamma, type = "response")  # response residuals
# Plot residuals
plot(predict(model_gamma, type = "response"), resids, xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
# works well but shows no interactions, which are necessary
qqnorm(resids)
qqline(resids, col = "red")

summary(model_gamma)

# Trying with interactions causes the model to not converge, likely because we have too few data points per group (only 3)
model_gamma_interact <- glmer(titratable_acid ~ species * watering_regime * timepoint + (1 | replicate),family = Gamma(link = "log"), data = dat)

# Try with a simpler model by doing it separate for each species:
VAN <- subset(dat, dat$species=="Tillandsia vanhyningii")
LEI <- subset(dat, dat$species=="Tillandsia leiboldiana")

shapiro.test(VAN$titratable_acid) # W = 0.95047, p-value = 0.2774
shapiro.test(LEI$titratable_acid) # W = 0.96584, p-value = 0.5662

# Since the TA values within species are normally distributed, we use LME
library(lmerTest)
model_VAN <- lmer(titratable_acid ~ watering_regime * timepoint + (1 | replicate), data = VAN)
model_LEI <- lmer(titratable_acid ~ watering_regime * timepoint + (1 | replicate), data = LEI)

resids_VAN <- residuals(model_VAN, type = "response")  # response residuals
resids_LEI <- residuals(model_LEI, type = "response")  # response residuals


# Plot residuals
pdf("LM_Diagnostic_Resid-to-Fitted_VAN.pdf", width = 5, height = 5)
plot(predict(model_VAN, type = "response"), resids_VAN, xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
dev.off()

pdf("LM_Diagnostic_Resid-to-Fitted_LEI.pdf", width = 5, height = 5)
plot(predict(model_LEI, type = "response"), resids_LEI, xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
dev.off()

pdf("LM_Diagnostic_QQplot_VAN.pdf", width = 5, height = 5)
qqnorm(resids_VAN)
qqline(resids_VAN, col = "red")
dev.off()

pdf("LM_Diagnostic_QQplot_LEI.pdf", width = 5, height = 5)
qqnorm(resids_LEI)
qqline(resids_LEI, col = "red")
dev.off()

shapiro.test(resids_VAN) # W = 0.97571, p-value = 0.8056
shapiro.test(resids_LEI) # W = 0.94208, p-value = 0.1815

# post-hoc test
library(emmeans)
emm_VAN <- emmeans(model_VAN, ~ watering_regime * timepoint)
pairs_emm_VAN <- pairs(emm_VAN)
summary(pairs_emm_VAN, adjust = "tukey")  # applying Tukey adjustment for multiple comparisons
emm_LEI <- emmeans(model_LEI, ~ watering_regime * timepoint)
pairs_emm_LEI <- pairs(emm_LEI)
summary(pairs_emm_LEI, adjust = "tukey")  # applying Tukey adjustment for multiple comparisons

# Export fixed effects
model_tidy_VAN <- broom.mixed::tidy(model_VAN)
write.csv(model_tidy_VAN, "LM_VAN_fixed_effects.csv")
# Export model summary statistics
model_glance_VAN <- broom.mixed::glance(model_VAN)
write.csv(model_glance_VAN, "LM_VAN_overall_summary.csv")
# Export random effects
model_random_VAN <- broom.mixed::tidy(model_VAN, effects = "ran_pars")
write.csv(model_random_VAN, "LM_VAN_random_effects.csv")
# Export post-hoc Tukey tests
post_hoc_VAN <- summary(pairs_emm_VAN, adjust = "tukey")
write.csv(post_hoc_VAN, "LM_VAN_post-hoc_Tukey.csv")

# Export fixed effects
model_tidy_LEI <- broom.mixed::tidy(model_LEI)
write.csv(model_tidy_LEI, "LM_LEI_fixed_effects.csv")
# Export model summary statistics
model_glance_LEI <- broom.mixed::glance(model_LEI)
write.csv(model_glance_LEI, "LM_LEI_overall_summary.csv")
# Export random effects
model_random_LEI <- broom.mixed::tidy(model_LEI, effects = "ran_pars")
write.csv(model_random_LEI, "LM_LEI_random_effects.csv")
# Export post-hoc Tukey tests
post_hoc_LEI <- summary(pairs_emm_LEI, adjust = "tukey")
write.csv(post_hoc_LEI, "LM_LEI_post-hoc_Tukey.csv")

#### FIGURES ####

# Make a separate figure for each species and place them next to eachother
VAN <- subset(dat, dat$species=="Tillandsia vanhyningii")
LEI <- subset(dat, dat$species=="Tillandsia leiboldiana")

p_van <- ggplot(VAN, aes(
  x=fct_reorder(timepoint, order),
  y=titratable_acid,
  fill=watering_regime)) +
  geom_boxplot(color = "black")+
  theme_bw() +
  labs(title="Tillandsia vanhyningii",
       x ="Sampling time", y = "Titratable acid (μequivalents / g FW)")+
  theme(plot.title = element_text(size = 14, face = "italic")) +
  theme(axis.title = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour="black",size=10,face="plain"),
        axis.text.y = element_text(colour="black",size=10,face="plain")) +
  scale_fill_manual(values = alpha(c("control" = "#ffa400", "drought" = "#ffa400"), c(1,.2))) +
  scale_y_continuous(breaks=seq(0,55,10), limits=c(0,55))+
  labs(fill = "Watering regime")

p_lei <- ggplot(LEI, aes(
  x=fct_reorder(timepoint, order),
  y=titratable_acid,
  fill=watering_regime)) +
  geom_boxplot(color = "black")+
  theme_bw() +
  labs(title="Tillandsia leiboldiana",
       x ="Sampling time", y = "Titratable acid (μequivalents / g FW)")+
  theme(plot.title = element_text(size = 14, face = "italic")) +
  theme(axis.title = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour="black",size=10,face="plain"),
        axis.text.y = element_text(colour="black",size=10,face="plain")) +
  scale_fill_manual(values = alpha(c("control" = "#2a9d8f", "drought" = "#2a9d8f"), c(1,.2))) +
  scale_y_continuous(breaks=seq(0,55,10), limits=c(0,55)) +
  labs(fill = "Watering regime")

pdf("Figure2_TitratableAcidity_TwoPlots.pdf", height = 8, width = 14)
grid.arrange(p_van, p_lei, nrow = 1)
dev.off()

# Make one figure where both species are combined
p_all <- ggplot(dat, aes(
  x=fct_reorder(timepoint, order),
  y=titratable_acid,
  fill=group)) +
  geom_boxplot(color = "black")+
  theme_bw() +
  labs(title="",
       x ="Sampling time", y = "Titratable acid (μequivalents / g FW)")+
  theme(axis.title = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour="black",size=10,face="plain"),
        axis.text.y = element_text(colour="black",size=10,face="plain")) +
  scale_fill_manual(labels = c("Tillandsia vanhyningii.control" = "T. vanhyningii, control", "Tillandsia vanhyningii.drought"="T. vanhyningii, drought", "Tillandsia leiboldiana.control" = "T. leiboldiana, control", "Tillandsia leiboldiana.drought"="T. leiboldiana, drought"),
                    values = c("Tillandsia vanhyningii.control" = "#ffa400", "Tillandsia vanhyningii.drought"="#FFE4B4", "Tillandsia leiboldiana.control" = "#2a9d8f", "Tillandsia leiboldiana.drought"="#BFEEE8")) +
  scale_y_continuous(breaks=seq(0,55,10), limits=c(0,55))

cairo_pdf("Figure2_TitratableAcidity_OnePlots.pdf", height = 8, width = 8)
p_all
dev.off()

# Take the delta titratable acidity between D10 and N10
deltaTA <- dat %>%
  filter(timepoint %in% c("D+10", "N+10")) %>%
  group_by(replicate, watering_regime, species) %>%
  summarise(Titratable_acidity_diff = -diff(titratable_acid))

deltaTA_avg <- deltaTA %>%
  group_by(species, watering_regime) %>%
  summarise(mean_diff = mean(Titratable_acidity_diff),
            sd_diff = sd(Titratable_acidity_diff))

deltaTA_avg$group <- factor(
  interaction(deltaTA_avg$species, deltaTA_avg$watering_regime),
  levels = unique(interaction(deltaTA_avg$species, deltaTA_avg$watering_regime))
)
write.table(deltaTA, file = "deltaTA_N10_D10.txt", row.names = FALSE, sep = "\t", quote = FALSE)

p_delta <- ggplot(deltaTA_avg, aes(x=group, y=mean_diff, fill=group)) +
  theme_bw()+
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_diff-sd_diff, ymax=mean_diff+sd_diff), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(labels = c("Tillandsia vanhyningii.control" = "T. vanhyningii, control", "Tillandsia vanhyningii.drought"="T. vanhyningii, drought", "Tillandsia leiboldiana.control" = "T. leiboldiana, control", "Tillandsia leiboldiana.drought"="T. leiboldiana, drought"),
                    values = c("Tillandsia vanhyningii.control" = "#ffa400", "Tillandsia vanhyningii.drought"="#FFE4B4", "Tillandsia leiboldiana.control" = "#2a9d8f", "Tillandsia leiboldiana.drought"="#BFEEE8")) +
  labs(title="",
       x ="",
       y = "ΔTA (μequivalents / g FW)",
       fill="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=12),
        legend.position = "bottom",
        legend.text = element_text(size=12))+
  ylim(c(-3,45))

pdf("Figure2_TitratableAcidity_Delta.05052024.pdf", height = 7, width = 9)
p_delta
dev.off()

# Perform tests on deltaTA
# I will first test if we can perform ANOVA on this data. For that, I need to test if the residuals are normally distributed.
# Fit the model
model <- aov(Titratable_acidity_diff ~ watering_regime * species, data = deltaTA)
model$
# Shapiro-Wilk test, test of normality
shapiro.test(residuals(model)) # p-value = 0.3574
# Residuals are normally distributed

# Levene's Test, homogeneity of variance
leveneTest(Titratable_acidity_diff ~ watering_regime * species, data = deltaTA) # p-value = 0.315
# Variances are homogeneous.
# Conclusion, ANOVA can be performed.
summary(model)

# Perform post-hoc test
posthoc_result <- TukeyHSD(model)
posthoc_result

tukey_results <- as.data.frame(posthoc_result$`watering_regime:species`)  # Replace <factor_name> with the name of the relevant factor

# Optionally, clean up the data frame to make it more understandable
names(tukey_results) <- c("Estimate", "Conf.Low", "Conf.High", "p.value")
tukey_results$Comparison <- rownames(tukey_results)
rownames(tukey_results) <- NULL

write.csv(tukey_results, file = "DeltaTA_PostHoc_TukeyHSD_interactionterm.csv", quote = F)
