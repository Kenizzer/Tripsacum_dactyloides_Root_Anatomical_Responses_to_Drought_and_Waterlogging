### Tripsacum greenhouse experiment
### Phenotypic analysis
### Code by Joel Swift 04-2024

library(tidyverse)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)

# Theme set and Color Palettes
theme_set(theme_bw())
treatment_palette <- c("Waterlogged 24" = "#26008E",
                       "Waterlogged 48" = "#134B99",
                       "Waterlogged 72" = "#0096A4",
                       "Control" = "#5CB270",
                       "Drought" = "#DAA520")

# Function from the R cookbook
# From: http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

#### 1) load data ####
Trip_df_full <- readxl::read_xlsx("Tripsacum_greenhouse_metadata.xlsx")
summary(Trip_df_full)
# Remove plants that were not assigned to a treatment group
Trip_df_full <- Trip_df_full[Trip_df_full$treatment != "NA", ]
summary(Trip_df_full)

#### 2) Height plot ####
Trip_df_full
# Gather data in long dataframe and reorder
long_df_ht <- gather(Trip_df_full, key = "day", value = "Height", 4:10, -treatment, -germination_cohort)
long_df_ht$day <- factor(long_df_ht$day, c("height_DAP24_mm", "height_DAP28_mm", "height_DAP30_mm",
                                           "height_DAP32_mm", "height_DAP34_mm", "height_DAP36_mm",
                                           "final_height_mm"))
# Plot
height_plot <- ggplot(long_df_ht, aes(x = day, y = Height, fill = treatment)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05)) +
  geom_boxplot(alpha = 0.8, outlier.colour = NA, color = "black") +
  scale_fill_manual(name = "Treatment", values = treatment_palette) +
  ylab("Shoot Height (mm)") +
  scale_x_discrete(labels=c("DAP24","DAP28","DAP30", "DAP32", "DAP34", "DAP36", "Final")) +
  theme(axis.title.x=element_blank(),
        legend.position = 'right')
#save
#ggsave("../figures/Shoot_height_by_treatment_timeline.svg", height_plot, height = 6, width = 8)
#ggsave("../figures/Shoot_height_by_treatment_timeline.png", height_plot, height = 6, width = 8)

#### 3) Dry weight Plot ####
weight_plot <- ggplot(Trip_df_full, aes(x = treatment, y = shoot_dry_weight_mg, fill = treatment)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
  geom_boxplot(alpha = 0.8, outlier.colour = NA, color = "black") +
  scale_fill_manual(name = "Treatment", values = treatment_palette) +
  ylab("Shoot dry Weight (mg)") + theme(axis.title.x=element_blank(), legend.position = 'right')
# Save
ggsave("../figures/Shoot_fresh_weight_by_treatment.svg", weight_plot, height = 85, width = 85, units = "mm")
ggsave("../figures/Shoot_fresh_weight_by_treatment.png", weight_plot, height = 85, width = 85, units = "mm")

#### 4) Combo plot final height and weight ####
weight_plot <- ggplot(Trip_df_full, aes(x = treatment, y = shoot_dry_weight_mg, fill = treatment)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
  geom_boxplot(alpha = 0.8, outlier.colour = NA, color = "black") +
  scale_fill_manual(name = "Treatment", values = treatment_palette) +
  ylab("Shoot dry Weight (mg)") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

final_height_plot <- ggplot(Trip_df_full, aes(x = treatment, y = final_height_mm, fill = treatment)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
  geom_boxplot(alpha = 0.8, outlier.colour = NA, color = "black") +
  scale_fill_manual(name = "Treatment", values = treatment_palette) +
  ylab("Shoot Height (mm)") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
# arrange
both_plots <- ggarrange(weight_plot, final_height_plot, common.legend = TRUE, legend = 'right')
# Save
ggsave("../figures/Figure1_Shoot_height_and_weight_treatment.svg", both_plots, height = 85, width = 170, units = "mm")
ggsave("../figures/Figure1_Shoot_height_and_weight_treatment.png", both_plots, height = 85, width = 170, units = "mm")

# Stats
summarySE(Trip_df_full, measurevar = "shoot_dry_weight_mg")
summarySE(Trip_df_full, measurevar = "final_height_mm")

#### 5) Linear models STATS
ht_mod <- lmer(final_height_mm ~ treatment + (1|germination_cohort), data = Trip_df_full)
anova(ht_mod)
rand(ht_mod)
mod.emm <- emmeans(ht_mod,~treatment)
contrast(mod.emm, "trt.vs.ctrl", ref = 1)

wt_mod <- lmer(shoot_dry_weight_mg ~ treatment + (1|germination_cohort), data = Trip_df_full)
anova(wt_mod)
rand(wt_mod)
mod.emm <- emmeans(wt_mod,~treatment)
contrast(mod.emm, "trt.vs.ctrl", ref = 1)
