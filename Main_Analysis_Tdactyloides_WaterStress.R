### Tripsacum greenhouse experiment
### Root section analysis
### Code by Joel Swift 04-2024

library(tidyverse)
library(readxl)
library(emmeans)
library(ggpubr)
library(lmerTest)
library(viridis)

# Set theme
theme_set(theme_bw())
treatment_palette <- c("Waterlogged 24" = "#26008E",
                       "Waterlogged 48" = "#134B99",
                       "Waterlogged 72" = "#0096A4",
                       "Control" = "#5CB270",
                       "Drought" = "#DAA520")


#### 1) Load data and calculate derived areas and percentages ####
section_data <- read_xlsx("Section_data_plus_meta.xlsx", sheet = 1)
# Convert couple of columns with NAs to numeric and segment to factor
numeric <- section_data %>% select(`Root Perimeter`:`Root Hair Number`) %>% mutate_if(is.character, as.numeric)
meta <- section_data %>% select(Plant_id:`Image Name`)
section_data <- bind_cols(meta, numeric)
rm(numeric, meta)
# convert all character columns to factor
section_data[sapply(section_data, is.character)] <- lapply(section_data[sapply(section_data, is.character)], as.factor)
section_data$Segment <- as.factor(section_data$Segment)
# Sanity check
summary(section_data)

# Calculate areas (um2)
section_data$`Cortex Area` <- (section_data$`Root Area` - section_data$`Stele Area`)
section_data$`Living Cortical Area` <- section_data$`Cortex Area` - section_data$`Aerenchyma Area`

# Calculate percentages (%)
section_data$`Living Cortical Percent` <- ((section_data$`Cortex Area` - section_data$`Aerenchyma Area`) / section_data$`Root Area`) * 100
section_data$`Aerenchyma Percent` <- (section_data$`Aerenchyma Area` / section_data$`Cortex Area`) * 100
section_data$`Aerenchyma Percent Whole Root` <- (section_data$`Aerenchyma Area` / section_data$`Root Area`) * 100

section_data$`Stele Percent` <- (section_data$`Stele Area` / section_data$`Root Area`) * 100
section_data$`Metaxylem Percent` <- (section_data$`Metaxylem Vessel Area` / section_data$`Stele Area`) * 100

# Ratios (X:Y)
section_data$`Stele:Cortex Ratio` <- section_data$`Stele Area` / section_data$`Cortex Area`
section_data$`Stele:Root Ratio` <- section_data$`Stele Area` / section_data$`Root Area`

# Scale root hair number by perimeter measurement 
section_data$`Root Hair Number` <- section_data$`Root Hair Number` / (section_data$`Root Perimeter`/1000) # Root hairs per mm of root perimeter

# Sanity check
str(section_data)


#### 2) Check the variance between sections of the same segment of root ####
# Select all plants that have three sections per segment (i.e, have a replicate A, B, and C)
section_data_index <- section_data %>% filter(Replicate %in% c("B", "C"))
section_data_replicated <- section_data %>% filter(Plant_id %in% unique(section_data_index$Plant_id))
rm(section_data_index)

# Calculate coefficient of variation by Replicate, Segment, and Plant ID.
cv <- function(x, na.rm = FALSE) (sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm))*100

#  CV for Replicate
CVs <- section_data_replicated %>%
  select(Plant_id, Segment, Replicate, Treatment, `Root Area`:`Stele:Root Ratio`) %>%
  group_by(Plant_id, Segment, Treatment) %>%
  mutate(across(`Root Area`:`Stele:Root Ratio`, cv)) %>%
  ungroup() 

Rep_CV_table <- CVs %>%
  group_by(Treatment) %>%
  select(Treatment, `Root Area`:`Stele:Root Ratio`) %>%
  filter_all(all_vars(!is.infinite(.))) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  pivot_longer(names_to = "Traits", values_to = "CV", -Treatment) %>%
  mutate(across(where(is.numeric), round, 0))

# Segment CV
rep_means <- section_data_replicated %>%
  select(Plant_id, Segment, Replicate, Treatment, `Root Area`:`Stele:Root Ratio`) %>%
  group_by(Plant_id, Segment, Treatment) %>%
  summarise(across(`Root Area`:`Stele:Root Ratio`, mean)) %>%
  ungroup() 

CVs <- rep_means %>%
  select(Plant_id, Segment, Treatment, `Root Area`:`Stele:Root Ratio`) %>%
  group_by(Plant_id, Treatment) %>%
  mutate(across(`Root Area`:`Stele:Root Ratio`, cv)) %>%
  ungroup() 

Seg_CV_table <- CVs %>%
  group_by(Treatment) %>%
  select(Treatment, `Root Area`:`Stele:Root Ratio`) %>%
  filter_all(all_vars(!is.infinite(.))) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  pivot_longer(names_to = "Traits", values_to = "CV", -Treatment) %>%
  mutate(across(where(is.numeric), round, 0))

# Plant ID CV
rep_means <- section_data_replicated %>%
  select(Plant_id, Segment, Replicate, Treatment, `Root Area`:`Stele:Root Ratio`) %>%
  group_by(Plant_id, Segment, Treatment) %>%
  summarise(across(`Root Area`:`Stele:Root Ratio`, mean)) %>%
  ungroup() 

CVs <- rep_means %>%
  select(Plant_id, Segment, Treatment, `Root Area`:`Stele:Root Ratio`) %>%
  group_by(Treatment, Segment) %>%
  mutate(across(`Root Area`:`Stele:Root Ratio`, cv)) %>%
  ungroup()

PID_CV_table <- CVs %>%
  group_by(Treatment) %>%
  select(Treatment, `Root Area`:`Stele:Root Ratio`) %>%
  filter_all(all_vars(!is.infinite(.))) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  pivot_longer(names_to = "Traits", values_to = "CV", -Treatment) %>%
  mutate(across(where(is.numeric), round, 0))

# Reorder factor levels to be more sensible
# Root surface > Cortex > Stele > Metaxylem Vessels
# Table 1
Rep_CV_table$Traits <- factor(Rep_CV_table$Traits, levels = rev(c("Root Area", "Epidermis Depth", "Root Hair Number",
                                                              "Cortex Area", "Cortical Cell File Number", "Living Cortical Percent", "Aerenchyma Percent",
                                                              "Stele Area", "Stele:Cortex Ratio", "Stele:Root Ratio",
                                                              "Metaxylem Vessel Area", "Metaxylem Vessel Number", "Metaxylem Vessel Mean Area")))
Rep_CV_table <- Rep_CV_table[!is.na(Rep_CV_table$Traits), ]
# Table 2
Seg_CV_table$Traits <- factor(Seg_CV_table$Traits, levels = rev(c("Root Area", "Epidermis Depth", "Root Hair Number",
                                                                  "Cortex Area", "Cortical Cell File Number", "Living Cortical Percent", "Aerenchyma Percent",
                                                                  "Stele Area", "Stele:Cortex Ratio", "Stele:Root Ratio",
                                                                  "Metaxylem Vessel Area", "Metaxylem Vessel Number", "Metaxylem Vessel Mean Area")))
Seg_CV_table <- Seg_CV_table[!is.na(Seg_CV_table$Traits), ]
# Table 3
PID_CV_table$Traits <- factor(PID_CV_table$Traits, levels = rev(c("Root Area", "Epidermis Depth", "Root Hair Number",
                                                                  "Cortex Area", "Cortical Cell File Number", "Living Cortical Percent", "Aerenchyma Percent",
                                                                  "Stele Area", "Stele:Cortex Ratio", "Stele:Root Ratio",
                                                                  "Metaxylem Vessel Area", "Metaxylem Vessel Number", "Metaxylem Vessel Mean Area")))
PID_CV_table <- PID_CV_table[!is.na(PID_CV_table$Traits), ]

# Tile plots with text labels
a <- ggplot(Rep_CV_table, aes(x = Treatment, y = Traits, fill = CV)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = CV), color = "black", size = 4) +
  scale_fill_gradient2(low='firebrick1', high='dodgerblue1', mid='white', name="Coefficent of Variation (%)", limits = c(0,177)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = "Among Replicate Sections")

b <- ggplot(Seg_CV_table, aes(x = Treatment, y = Traits, fill = CV)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = CV), color = "black", size = 4) +
  scale_fill_gradient2(low='firebrick1', high='dodgerblue1', mid='white', name="Coefficent of Variation (%)", limits = c(0,177)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(title = "Among Root Segments")


c <- ggplot(PID_CV_table, aes(x = Treatment, y = Traits, fill = CV)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = CV), color = "black", size = 4) +
  scale_fill_gradient2(low='firebrick1', high='dodgerblue1', mid='white', name="Coefficent of Variation (%)", limits = c(0,177)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank())  +
  labs(title = "Among Plants")

CV_plots <- ggarrange(a,b,c, common.legend = TRUE, nrow = 1, legend = "right", labels = "AUTO", align = 'hv')

ggsave("figures/Figure3_CV_plots.png", CV_plots, width = 16, height = 8, units = "in")
ggsave("figures/Figure3_CV_plots.svg", CV_plots, width = 16, height = 8, units = "in")

#### 3) Analysis with A sections, expanded sampling ####
# We use only single section from each of the samples that had 3 sections (labeled A,B,C) imaged to avoid
# including technical replicates into the main analysis.
## Select only A sections
A_section_data <- section_data %>% filter(Replicate == "A")

# flip order from 8 == root tip, 1 == root base
# to: 1 == root tip, 8 == Root base
# This is more sensible for presentation.
A_section_data <- A_section_data %>%
  group_by(Plant_id) %>%
  mutate(Segment = rev(Segment))

#### 4) Root length plot ####
# Root lengths by treatment
Root_length_count <- A_section_data %>%
  group_by(Treatment) %>%
  count(Segment)

Root_length_count_plot <- ggplot(Root_length_count, aes(x = Treatment, y = n, fill = Treatment)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks = seq(1, 12, by = 1)) +
  ylab("Number of Samples") +
  scale_fill_manual(values = treatment_palette, name = "Treatment") +
  theme(legend.text.align = 0,
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = 'right') + 
  xlab("Root Segment (cm) | Root tip to root base") +
  facet_grid(~Segment)

ggsave("figures/FigureS1_Root_length_count.svg", Root_length_count_plot, width = 10, height = 4, units = "in")
ggsave("figures/FigureS1_Root_length_count.png", Root_length_count_plot, width = 10, height = 4, units = "in")

# Sanity check
str(A_section_data)

#### 4.1) Root aerenchmya plot, in relation to whole root to compare to Gong et. al 2019 ####
# Gong et al. 2019 used the aerenchyma as a percentage of the whole root, here I make a similar plot and summarise the stats.
A_section_data %>%
  group_by(Treatment, Segment) %>%
  select(`Aerenchyma Percent Whole Root`) %>%
  summarise(mean = mean(`Aerenchyma Percent Whole Root`, na.rm = TRUE), sd = sd(`Aerenchyma Percent Whole Root`, na.rm = TRUE)) %>%
  ggplot(., aes(y = Segment, x = mean, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, color = "black") +
  scale_fill_manual(values = treatment_palette) +
  geom_errorbar( aes(y=Segment, xmin=mean, xmax=mean+sd),stat = "identity", position = "dodge", width=0.8, colour="black", alpha=0.8, size=0.5)

stats <- A_section_data %>%
  group_by(Treatment, Segment) %>%
  select(`Aerenchyma Percent Whole Root`) %>%
  summarise(mean = mean(`Aerenchyma Percent Whole Root`, na.rm = TRUE), sd = sd(`Aerenchyma Percent Whole Root`, na.rm = TRUE))

#### 5) Modeling and treatment to control plot ####
# Models
RA_mod <- lmer(`Root Area` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)
EP_mod <- lmer(`Epidermis Depth` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)
RH_mod <- lmer(`Root Hair Number` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)
CA_mod <- lmer(`Cortex Area` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)
CF_mod <- lmer(`Cortical Cell File Number` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)
LC_mod <- lmer(`Living Cortical Percent` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)
AR_mod <- lmer(`Aerenchyma Percent` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)
SA_mod <- lmer(`Stele Area` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)
SC_mod <- lmer(`Stele:Cortex Ratio` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)
SR_mod <- lmer(`Stele:Root Ratio` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)
MA_mod <- lmer(`Metaxylem Vessel Area` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)
MN_mod <- lmer(`Metaxylem Vessel Number` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)
MV_mod <- lmer(`Metaxylem Vessel Mean Area` ~ Treatment*Segment + (1|Plant_id), data = A_section_data)

# Anovas
anova(RA_mod)
anova(EP_mod)
anova(RH_mod)
anova(CA_mod)
anova(CF_mod)
anova(LC_mod)
anova(AR_mod)
anova(SA_mod)
anova(SC_mod)
anova(SR_mod)
anova(MA_mod)
anova(MN_mod)
anova(MV_mod)

# Contrasts T*S
contrast(emmeans(RA_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(EP_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(RH_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(CA_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(CF_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(LC_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(AR_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(SA_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(SC_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(SR_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(MA_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(MN_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(MV_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)

test <- A_section_data %>%
  filter(Treatment == "Waterlogged 48")


# Histograms of residuals
hist(resid(RA_mod))
hist(resid(EP_mod))
hist(resid(RH_mod))
hist(resid(CA_mod))
hist(resid(CF_mod))
hist(resid(LC_mod))
hist(resid(AR_mod))
hist(resid(SA_mod))
hist(resid(SC_mod))
hist(resid(SR_mod))
hist(resid(MA_mod))
hist(resid(MN_mod))
hist(resid(MV_mod))
# Residuals look mostly fine in the histograms, often have long tails.

# Posthoc dunnet tests
Dunn_emmeans <- function(mod){
  anova(mod)
  mod.emm <- emmeans(mod,~Treatment|Segment)
  con <- contrast(mod.emm, "trt.vs.ctrl", ref = 1)
  return(con)
}

plot(Dunn_emmeans(RA_mod))
plot(Dunn_emmeans(EP_mod))
plot(Dunn_emmeans(RH_mod))
plot(Dunn_emmeans(CA_mod))
plot(Dunn_emmeans(CF_mod))
plot(Dunn_emmeans(LC_mod))
plot(Dunn_emmeans(AR_mod))
plot(Dunn_emmeans(SA_mod))
plot(Dunn_emmeans(SC_mod))
plot(Dunn_emmeans(SR_mod))
plot(Dunn_emmeans(MA_mod))
plot(Dunn_emmeans(MN_mod))
plot(Dunn_emmeans(MV_mod))

# Dunnett Figure df function
Dunn_emmeans_DF <- function(mod){
  anova(mod)
  mod.emm <- emmeans(mod,~Treatment|Segment)
  con <- contrast(mod.emm, "trt.vs.ctrl", ref = 1)
  # These indexes break the emmgrid into pieces to make a dataframe (easier for later plotting)
  df <- rbind(data.frame(con[1:4]),
              data.frame(con[5:8]),
              data.frame(con[9:12]),
              data.frame(con[13:16]),
              data.frame(con[17:20]),
              data.frame(con[21:24]),
              data.frame(con[25:28]),
              data.frame(con[29:32]))
  return(df)
}

RA_df <- data.frame(Dunn_emmeans_DF(RA_mod), "Trait" = replicate(32, "Root Area"))
EP_df <- data.frame(Dunn_emmeans_DF(EP_mod), "Trait" = replicate(32, "Epidermis Depth"))
RH_df <- data.frame(Dunn_emmeans_DF(RH_mod), "Trait" = replicate(32, "Root Hair Number"))
CA_df <- data.frame(Dunn_emmeans_DF(CA_mod), "Trait" = replicate(32, "Cortex Area"))
CF_df <- data.frame(Dunn_emmeans_DF(CF_mod), "Trait" = replicate(32, "Cortical Cell File Number"))
LC_df <- data.frame(Dunn_emmeans_DF(LC_mod), "Trait" = replicate(32, "Living Cortical Percent"))
AR_df <- data.frame(Dunn_emmeans_DF(AR_mod), "Trait" = replicate(32, "Aerenchyma Percent"))
SA_df <- data.frame(Dunn_emmeans_DF(SA_mod), "Trait" = replicate(32, "Stele Area"))
SC_df <- data.frame(Dunn_emmeans_DF(SC_mod), "Trait" = replicate(32, "Stele:Cortex Ratio"))
SR_df <- data.frame(Dunn_emmeans_DF(SR_mod), "Trait" = replicate(32, "Stele:Root Ratio"))
MA_df <- data.frame(Dunn_emmeans_DF(MA_mod), "Trait" = replicate(32, "Metaxylem Vessel Area"))
MN_df <- data.frame(Dunn_emmeans_DF(MN_mod), "Trait" = replicate(32, "Metaxylem Vessel Number"))
MV_df <- data.frame(Dunn_emmeans_DF(MV_mod), "Trait" = replicate(32, "Metaxylem Vessel Mean Area"))

# Create long dataframe
full_df <- rbind(RA_df, EP_df, RH_df, CA_df,
                 CF_df, LC_df, AR_df, SA_df, 
                 SC_df, SR_df, MA_df, MN_df,
                 MV_df)

# Rearrange trait factor order
full_df$Trait <- factor(full_df$Trait, levels = rev(c("Root Area", "Epidermis Depth", "Root Hair Number",
                                                      "Cortex Area", "Cortical Cell File Number", "Living Cortical Percent", "Aerenchyma Percent",
                                                      "Stele Area", "Stele:Cortex Ratio", "Stele:Root Ratio",
                                                      "Metaxylem Vessel Area", "Metaxylem Vessel Number", "Metaxylem Vessel Mean Area")))

# Create dataframes for plotting by treatment
Drought_DF <- full_df %>%
  filter(contrast == "Drought - Control") %>%
  group_by(Trait) %>%
  mutate(`Relative Change`= round(((estimate/ sum(abs(estimate)))*100), digits = 0),
         "Significance" = case_when(p.value < 0.05 ~ "*",
                                    p.value < 0.1 ~ "-",
                                    p.value > 0.05 ~ ""))

Flood24_DF <- full_df %>%
  filter(contrast == "Waterlogged 24 - Control") %>%
  group_by(Trait) %>%
  mutate(`Relative Change`= round(((estimate/ sum(abs(estimate)))*100), digits = 0),
         "Significance" = case_when(p.value < 0.05 ~ "*",
                                    p.value < 0.1 ~ "-",
                                    p.value > 0.05 ~ ""))

Flood48_DF <- full_df %>%
  filter(contrast == "Waterlogged 48 - Control") %>%
  group_by(Trait) %>%
  mutate(`Relative Change`= round(((estimate/ sum(abs(estimate)))*100), digits = 0),
         "Significance" = case_when(p.value < 0.05 ~ "*",
                                    p.value < 0.1 ~ "-",
                                    p.value > 0.05 ~ ""))

Flood72_DF <- full_df %>%
  filter(contrast == "Waterlogged 72 - Control") %>%
  group_by(Trait) %>%
  mutate(`Relative Change`= round(((estimate/ sum(abs(estimate)))*100), digits = 0),
         "Significance" = case_when(p.value < 0.05 ~ "*",
                                    p.value < 0.1 ~ "-",
                                    p.value > 0.05 ~ ""))

# get limits for plotting at same gradient scale
max(Flood72_DF$`Relative Change`, Flood24_DF$`Relative Change`,
    Drought_DF$`Relative Change`, Flood48_DF$`Relative Change`)

min(Flood72_DF$`Relative Change`, Flood24_DF$`Relative Change`,
    Drought_DF$`Relative Change`, Flood48_DF$`Relative Change`)

# plots
a <- ggplot(Drought_DF, aes(x = Segment , y = Trait, fill = `Relative Change`)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = `Significance`), color = "black", size = 8) +
  scale_fill_gradient2(low='firebrick1', high='dodgerblue1', mid='white', name="Relative Change from Control (%)", limits = c(-45,40)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = "Drought")

b <- ggplot(Flood24_DF, aes(x = Segment , y = Trait, fill = `Relative Change`)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = `Significance`), color = "black", size = 8) +
  scale_fill_gradient2(low='firebrick1', high='dodgerblue1', mid='white', name="Relative Change from Control (%)", limits = c(-45,40)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(title = "Flood 24")

c <- ggplot(Flood48_DF, aes(x = Segment , y = Trait, fill = `Relative Change`)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = `Significance`), color = "black", size = 8) +
  scale_fill_gradient2(low='firebrick1', high='dodgerblue1', mid='white', name="Relative Change from Control (%)", limits = c(-45,40)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(title = "Flood 48")

d <- ggplot(Flood72_DF, aes(x = Segment , y = Trait, fill = `Relative Change`)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = `Significance`), color = "black", size = 8) +
  scale_fill_gradient2(low='firebrick1', high='dodgerblue1', mid='white', name="Relative Change from Control (%)", limits = c(-40,40)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(title = "Flood 72")

#arrange and save
Segment_by_treatment <- ggarrange(a,b,c,d, labels = "AUTO", nrow = 1, legend = "right", common.legend = TRUE, align = "hv")
ggsave("figures/Figure4_heatmap_bytreatment_plots.png", Segment_by_treatment, width = 16, height = 8, units = "in")
ggsave("figures/Figure4_heatmap_bytreatment_plots.svg", Segment_by_treatment, width = 16, height = 8, units = "in")

#### 7) Emmeans Cortex Plot ####
contrast(emmeans(CF_mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)

Cortical_file_plot <- ggplot(plot(emmeans(CF_mod,~Treatment|Segment), plotit = FALSE), aes(x = Segment, y = the.emmean, color = Treatment)) +
  ylab("Cortical Cell File Number") +
  scale_y_continuous(limits = c(10, 20), breaks = seq(10, 20, by = 2)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, fill = Treatment), size = 0.50,
                  linewidth = 1.75, fill = "white", shape = 22, position = position_dodge(width = 0.75)) +
  scale_color_manual(values = treatment_palette, name = "Treatment") +
  xlab("Root Segment (cm) | Root tip to root base") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank())

LCA_plot <- ggplot(plot(emmeans(LC_mod,~Treatment|Segment), plotit = FALSE), aes(x = Segment, y = the.emmean, color = Treatment)) +
  ylab("Living Cortical Percent") +
  scale_y_continuous(limits = c(50, 90), breaks = seq(50, 90, by = 10), labels = function(x) paste0(x, "%")) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, fill = Treatment), size = 0.50,
                  linewidth = 1.75, fill = "white", shape = 22, position = position_dodge(width = 0.75)) +
  scale_color_manual(values = treatment_palette, name = "Treatment") +
  xlab("Root Segment (cm) | Root tip to root base")

CO_LCA_plot <- ggarrange(Cortical_file_plot, LCA_plot, nrow = 2, common.legend = TRUE,
                         legend = "right", labels = "AUTO", align = 'hv')
ggsave("figures/Figure5_Cell_file_LCA_plot.svg", CO_LCA_plot, width = 170, height = 135, units = "mm")
ggsave("figures/Figure5_Cell_file_LCA_plot.png", CO_LCA_plot, width = 170, height = 135, units = "mm")

#### 8) Emmeans Metaxylem Vessel Plot ####
contrast(emmeans(SA_mod,~Treatment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(MN_mod,~Treatment), "trt.vs.ctrl", ref = 1)
contrast(emmeans(MV_mod,~Treatment), "trt.vs.ctrl", ref = 1)

a <- ggplot(plot(emmeans(SA_mod,~Treatment), plotit = FALSE), aes(x = Treatment, y = the.emmean, color = Treatment)) +
  ylab("Stele Area (um^2)") +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, fill = Treatment),
                  size = 0.75, linewidth = 1.75, fill = "white", shape = 22) +
  scale_color_manual(values = treatment_palette, name = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank())

b <- ggplot(plot(emmeans(MN_mod,~Treatment), plotit = FALSE), aes(x = Treatment, y = the.emmean, color = Treatment)) +
  ylab("Metaxylem Vessel Count") +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, fill = Treatment),
                  size = 0.75, linewidth = 1.75, fill = "white", shape = 22,
                  position = position_dodge(width = 0.5)) +
  scale_color_manual(values = treatment_palette, name = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank())

c <- ggplot(plot(emmeans(MV_mod,~Treatment), plotit = FALSE), aes(x = Treatment, y = the.emmean, color = Treatment)) +
  ylab("Metaxylem Vessel Mean Area (um^2)") +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, fill = Treatment),
                  size = 0.75, linewidth = 1.75, fill = "white", shape = 22) +
  scale_color_manual(values = treatment_palette, name = "Treatment") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank())


A_section_data$`Metaxylem Vessel Area`

A_section_data$`Metaxylem Vessel Mean Area`

Metaxylem_plot <- ggarrange(a,b,c, nrow = 1, common.legend = TRUE, legend = 'right', labels = "AUTO")
ggsave("figures/Figure6_metaxylem_plot.svg", Metaxylem_plot, width = 170, height = 95, units = "mm")
ggsave("figures/Figure6_metaxylem_plot.png", Metaxylem_plot, width = 170, height = 95, units = "mm")


#### 9) Drought vs Flood 72 Comparison | S:R Ratio plot####
# filter to select stress classes
Only_F72_D <- A_section_data %>%
  filter(Treatment %in% c("Drought", "Waterlogged 72", "Control"))

# Sanity check
str(Only_F72_D)

# STATS
anova(lmer(`Stele:Root Ratio` ~ Treatment*Segment + (1|Plant_id), data = Only_F72_D))
mod <- lmer(`Stele:Root Ratio` ~ Treatment*Segment + (1|Plant_id), data = Only_F72_D)
contrast(emmeans(mod,~Treatment|Segment), "trt.vs.ctrl", ref = 1)
pairs(emmeans(mod,~Treatment|Segment))

SR_ratio_plot <- ggplot(plot(emmeans(mod,~Treatment*Segment), plotit = FALSE), aes(x = Segment, y = the.emmean, color = Treatment)) +
  ylab("Stele:Root ratio") +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, color = Treatment),
                  size = 0.75, linewidth = 1.75, fill = "white", shape = 22,
                  position = position_dodge(width = 0.5)) +
  scale_color_manual(values = treatment_palette, name = "Treatment") +
  xlab("Root Segment (cm) | Root tip to root base")

ggsave("figures/Figure7_Stele-Root_Ratio_plot.svg", SR_ratio_plot, width = 170, height = 60, unit = "mm")
ggsave("figures/Figure7_Stele-Root_Ratio_plot.png", SR_ratio_plot, width = 170, height = 60, unit = "mm")

