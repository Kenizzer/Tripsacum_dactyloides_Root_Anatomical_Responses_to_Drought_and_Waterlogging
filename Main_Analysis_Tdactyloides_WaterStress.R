### Tripsacum greenhouse experiment
### Root section analysis
### Code by Joel Swift 04-2024

library(tidyverse)
library(readxl)
library(emmeans)
library(ggpubr)
library(lmerTest)
library(viridis)
library(scales)

# Set theme
theme_set(theme_bw())
treatment_palette <- c("Waterlogged 24" = "#26008E",
                       "Waterlogged 48" = "#134B99",
                       "Waterlogged 72" = "#0096A4",
                       "Control" = "#5CB270",
                       "Drought" = "#DAA520")


#### 1) Load data and calculate derived areas and percentages ####
section_data <- read_xlsx("Section_data_plus_metadata.xlsx", sheet = 1)
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
section_data$`Living Cortical Percent` <- (section_data$`Living Cortical Area` / section_data$`Root Area`) * 100
section_data$`Aerenchyma Percent` <- (section_data$`Aerenchyma Area` / section_data$`Root Area`) * 100

section_data$`Stele Percent` <- (section_data$`Stele Area` / section_data$`Root Area`) * 100
section_data$`Metaxylem Percent` <- (section_data$`Metaxylem Vessel Area` / section_data$`Stele Area`) * 100

# Ratios (X:Y)
section_data$`Stele:Cortex Ratio` <- section_data$`Stele Area` / section_data$`Cortex Area`
section_data$`Stele:Root Ratio` <- section_data$`Stele Area` / section_data$`Root Area`

# Scale root hair number by perimeter measurement 
section_data$`Root Hair Number` <- section_data$`Root Hair Number` / (section_data$`Root Perimeter`/1000) # Root hairs per mm of root perimeter

# Sanity check
str(section_data)

### Diameter range of roots ###
section_data$root_diameter <- (sqrt(section_data$`Root Area` / pi) * 2)

diameter_df <- section_data %>%
  select(root_diameter, Plant_id, Treatment) %>%
  group_by(Plant_id, Treatment) %>%
  summarize(mean_diameter_mm = (mean(root_diameter, na.rm = TRUE)/1000),
            count = n())

ggplot(diameter_df, aes(x = Treatment, y = mean_diameter_mm, fill = Treatment)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(alpha = 0.8, outlier.colour = NA, color = "black") +
  scale_fill_manual(name = "Treatment", values = treatment_palette) +
  ylab("Root Diameter (mm)") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

x <- section_data %>%
  select(root_diameter, Plant_id, Treatment) %>%
  group_by(Plant_id) %>%
  summarize(mean_diameter_mm = (mean(root_diameter, na.rm = TRUE)/1000),
            count = n())

# mean and standard error
summary(x$mean_diameter_mm)
sd(x$mean_diameter_mm) / sqrt(length(x$mean_diameter_mm))

#### 2) Check the variance between sections of the same segment of root ####
# Select all plants that have three sections per segment (i.e, have a replicate A, B, and C)
section_data_index <- section_data %>% filter(Replicate %in% c("B", "C"))
section_data_replicated <- section_data %>% filter(Plant_id %in% unique(section_data_index$Plant_id))
rm(section_data_index)

##### Percent variance explained  #####
pve_by_t <- function(trait, treatment_var = "Treatment") {
  # Trait with back ticks
  trait_sym <- as.name(trait)
  # Treatment  levels
  treatments <- unique(section_data_replicated[[treatment_var]])
  # Initialize output list
  output <- list()
  
  # Loop through treatment levels
  for (t in treatments) {
    # Subset data
    data_subset <- section_data_replicated %>% filter(.data[[treatment_var]] == t)
    # Conditional Metaxylem Vessel Number
    formula <- if (trait == "Metaxylem Vessel Number") {
      as.formula(bquote(.(trait_sym) ~ (1 | Plant_id)))
    } else {
      as.formula(bquote(.(trait_sym) ~ (1 | Plant_id) + (1 | Plant_id:Segment)))
    }
    # Fit model
    model <- lmer(formula, data = data_subset)
    
    # Get variance components
    var_comps <- as.data.frame(VarCorr(model))
    total_var <- sum(var_comps$vcov)
    var_comps$PVE <- (var_comps$vcov / total_var)
    
    # Tidy results
    result <- var_comps[, c("grp", "vcov", "PVE")]
    result$Trait <- trait
    result$Treatment <- t
    names(result) <- c("Group", "Variance", "Percent_Variance_Explained", "Trait", "Treatment")
    
    # Rename groups
    result$Group <- dplyr::recode(result$Group,
                                  "Residual" = "Among Replicates",
                                  "Plant_id:Segment" = "Among Segments within Plants",
                                  "Plant_id" = "Among Plants")
    #loop result
    output[[t]] <- result
  }
  
  # Combine and return
  temp <- bind_rows(output)
  return( temp %>%
    select(-Variance) %>%
    pivot_longer(cols = Percent_Variance_Explained, values_to = "PVE"))
}

#Modeling Strategy 
# Variance partitioned by Sources
# (1|Plant_id) = Among plants (regardless of segment) [All Plant_id]
# (1|Plant_id:Segment) = Among Segments [1,2,3,4,5,6,7,8]
# (1|Residual) = Among Replicates [A,B,C]
PVE_df <- rbind(pve_by_t("Root Area"),
                pve_by_t("Epidermis Depth"),
                pve_by_t("Root Hair Number"),
                pve_by_t("Cortex Area"),
                pve_by_t("Cortical Cell File Number"),
                pve_by_t("Living Cortical Percent"),
                pve_by_t("Aerenchyma Percent"),
                pve_by_t("Stele Area"),
                pve_by_t("Stele:Cortex Ratio"),
                pve_by_t("Stele:Root Ratio"),
                pve_by_t("Metaxylem Vessel Area"),
                pve_by_t("Metaxylem Vessel Number"),
                pve_by_t("Metaxylem Vessel Mean Area"))

# There is no variance within replicates for metaxylem number
# so the "residual term" is re-coded for that trait here
PVE_df <- PVE_df %>%
  mutate(Group = if_else(
    Trait == "Metaxylem Vessel Number" & Group == "Among Replicates",
    "Among Segments within Plants",
    Group
  ))

# set factor levels for plotting
PVE_df$Group <- factor(PVE_df$Group, levels = c("Among Replicates",
                                                "Among Segments within Plants",
                                                "Among Plants"))

# Get mean PVE across traits by treatment
mean_pve_df <- PVE_df %>%
  group_by(Treatment, Group) %>%
  summarise(PVE = mean(PVE), .groups = "drop") %>%
  mutate(Trait = "Mean PVE", name = "Percent_Variance_Explained")
PVE_df <- bind_rows(PVE_df, mean_pve_df)

# Factor level for plot right to left. Rename root hair number.
PVE_df$Trait <- factor(PVE_df$Trait, levels = rev(c("Mean PVE","Root Area", "Epidermis Depth", "Root Hair Number",
                                                    "Cortex Area", "Cortical Cell File Number", "Living Cortical Percent", "Aerenchyma Percent",
                                                    "Stele Area", "Stele:Cortex Ratio", "Stele:Root Ratio",
                                                    "Metaxylem Vessel Mean Area", "Metaxylem Vessel Area", "Metaxylem Vessel Number")))
levels(PVE_df$Trait)[levels(PVE_df$Trait) == "Root Hair Number"] <- "Root Hair Density"
# PLOTTING 
Barplot_PVE <- ggplot(PVE_df, aes(x = Trait, y = PVE, fill = Group)) +
  annotate("rect",xmin = 13.5, xmax = 14.5, ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.8) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.1) +
  facet_wrap(~ Treatment, ncol = 1) +
  labs(y = "Percent Variance Explained", fill = "Group") +
  scale_fill_brewer(palette = "Greys", name = "Source") +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, by = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), axis.title.x = element_blank(),
        strip.text = element_text(size = 6, margin = margin(t = 1, b = 1)), legend.position = "top",
        legend.text = element_text(size = 6), legend.title = element_text(size = 10),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 2),
        legend.margin = margin(b = 0.5, t = 0.5)) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))
Barplot_PVE

ggsave("figures/Figure3_PVE_Plot.png", Barplot_PVE, width = 3.25, height = 6, units = "in")
ggsave("figures/Figure3_PVE_Plot.svg", Barplot_PVE, width = 3.25, height = 6, units = "in")

# Numbers
section_data_replicated %>%
  distinct(Plant_id, Treatment) %>%
  count(Treatment, name = "N_treatment") # Number plants by treatment

section_data_replicated %>%
  distinct(Plant_id, Treatment, Segment) %>%
  count(Treatment, name = "N_treatment") # Number Segments by treatment

section_data_replicated %>%
  distinct(Plant_id, Segment, Treatment, Replicate) %>%
  count(Treatment, name = "N_treatment") # Number of replicates by treatment

# Stats, dunnett comparing rep to plants/segments
lm_by_traits <- lm(PVE ~ Group*Trait, data = PVE_df %>% filter(Trait != "Mean PVE"))
contrast(emmeans(lm_by_traits,~Group|Trait), "trt.vs.ctrl", ref = 1)

PVE_df %>%
  filter(Trait == "Mean PVE") %>%
  select(Group, Treatment, PVE) %>%
  filter(Group == "Among Replicates")

PVE_df %>%
  filter(Trait == "Mean PVE") %>%
  select(Group, Treatment, PVE) %>%
  filter(Group == "Among Segments within Plants")

PVE_df %>%
  filter(Trait == "Mean PVE") %>%
  select(Group, Treatment, PVE) %>%
  filter(Group == "Among Plants")

section_data_replicated %>%
  summarize(CV = sd(`Epidermis Depth`, na.rm = TRUE) / mean(`Epidermis Depth`, na.rm = TRUE))
section_data_replicated %>%
  summarize(CV = sd(`Root Area`, na.rm = TRUE) / mean(`Root Area`, na.rm = TRUE))
section_data_replicated %>%
  summarize(CV = sd(`Metaxylem Vessel Number`, na.rm = TRUE) / mean(`Metaxylem Vessel Number`, na.rm = TRUE))


section_data_replicated %>%
  pivot_longer(
    cols = where(is.numeric), # or specify trait columns directly
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  group_by(Trait) %>%
  summarize(
    CV = sd(Value, na.rm = TRUE) / mean(Value, na.rm = TRUE)
  ) %>%
  arrange(desc(CV))


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
  select(`Aerenchyma Percent`) %>%
  summarise(mean = mean(`Aerenchyma Percent`, na.rm = TRUE), sd = sd(`Aerenchyma Percent`, na.rm = TRUE)) %>%
  ggplot(., aes(y = Segment, x = mean, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, color = "black") +
  scale_fill_manual(values = treatment_palette) +
  geom_errorbar( aes(y=Segment, xmin=mean, xmax=mean+sd),stat = "identity", position = "dodge", width=0.8, colour="black", alpha=0.8, size=0.5)

stats <- A_section_data %>%
  group_by(Treatment, Segment) %>%
  select(`Aerenchyma Percent`) %>%
  summarise(mean = mean(`Aerenchyma Percent`, na.rm = TRUE), sd = sd(`Aerenchyma Percent`, na.rm = TRUE))

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
  xlab("Root Segment (cm)") +
  scale_fill_gradient2(low='dodgerblue1', high='firebrick1', mid='white', name="Relative Change from Control (%)", limits = c(-45,40)) +
  theme(axis.title.y = element_blank()) +
  labs(title = "Drought")

b <- ggplot(Flood24_DF, aes(x = Segment , y = Trait, fill = `Relative Change`)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = `Significance`), color = "black", size = 8) +
  xlab("Root Segment (cm)") +
  scale_fill_gradient2(low='dodgerblue1', high='firebrick1', mid='white', name="Relative Change from Control (%)", limits = c(-45,40)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Flood 24")

c <- ggplot(Flood48_DF, aes(x = Segment , y = Trait, fill = `Relative Change`)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = `Significance`), color = "black", size = 8) +
  xlab("Root Segment (cm)") +

  scale_fill_gradient2(low='dodgerblue1', high='firebrick1', mid='white', name="Relative Change from Control (%)", limits = c(-45,40)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Flood 48")

d <- ggplot(Flood72_DF, aes(x = Segment , y = Trait, fill = `Relative Change`)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = `Significance`), color = "black", size = 8) +
  xlab("Root Segment (cm)") +
  scale_fill_gradient2(low='dodgerblue1', high='firebrick1', mid='white', name="Relative Change from Control (%)", limits = c(-40,40)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Flood 72")

#arrange and save
Segment_by_treatment <- ggarrange(a,b,c,d, labels = "AUTO", nrow = 1, legend = "right", common.legend = TRUE, align = "hv")
ggsave("figures/Figure4_heatmap_bytreatment_plots.png", Segment_by_treatment, width = 16, height = 8, units = "in")
ggsave("figures/Figure4_heatmap_bytreatment_plots.svg", Segment_by_treatment, width = 16, height = 8, units = "in")

#### 7) Emmeans Cortex Plot ####
RH_plot_df <- plot(emmeans(RH_mod,~Treatment|Segment), plotit = FALSE) %>% select(Segment, Treatment, upper.CL)
Root_hair_density_plot <- ggplot(plot(emmeans(RH_mod,~Treatment|Segment), plotit = FALSE), aes(x = Segment, y = the.emmean, color = Treatment)) +
  ylab("Root Hair Density") +
  scale_y_continuous(limits = c(-7, 25), breaks = seq(-5, 25, by = 5)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, fill = Treatment), size = 0.50,
                  linewidth = 1.75, fill = "white", shape = 22, position = position_dodge(width = 0.75)) +
  scale_color_manual(values = treatment_palette, name = "Treatment") +
  theme(axis.title.x = element_blank()) +
  annotate("text", x = 0.85, y = RH_plot_df %>% filter(Treatment == "Drought" & Segment == 1) %>% pull(upper.CL) * 1.1 , label = "*", color = treatment_palette["Drought"], size = 6) +
  annotate("text", x = 1.85, y = RH_plot_df %>% filter(Treatment == "Drought" & Segment == 2) %>% pull(upper.CL) * 1.1 , label = "*", color = treatment_palette["Drought"], size = 6)

CF_plot_df <- plot(emmeans(CF_mod,~Treatment|Segment), plotit = FALSE) %>% select(Segment, Treatment, upper.CL)
Cortical_file_plot <- ggplot(plot(emmeans(CF_mod,~Treatment|Segment), plotit = FALSE), aes(x = Segment, y = the.emmean, color = Treatment)) +
  ylab("Cortical Cell File Number") +
  scale_y_continuous(limits = c(10, 20), breaks = seq(10, 20, by = 2)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, fill = Treatment), size = 0.50,
                  linewidth = 1.75, fill = "white", shape = 22, position = position_dodge(width = 0.75)) +
  scale_color_manual(values = treatment_palette, name = "Treatment") +
  theme(axis.title.x = element_blank()) +
  annotate("text", x = 2.85, y = CF_plot_df %>% filter(Treatment == "Drought" & Segment == 3) %>% pull(upper.CL) * 1.1, label = "-", color = treatment_palette["Drought"], size = 6) +
  annotate("text", x = 4.85, y = CF_plot_df %>% filter(Treatment == "Drought" & Segment == 5) %>% pull(upper.CL) * 1.1, label = "-", color = treatment_palette["Drought"], size = 6) +
  annotate("text", x = 5.85, y = CF_plot_df %>% filter(Treatment == "Drought" & Segment == 6) %>% pull(upper.CL) * 1.1, label = "*", color = treatment_palette["Drought"], size = 6) +
  annotate("text", x = 6.85, y = CF_plot_df %>% filter(Treatment == "Drought" & Segment == 7) %>% pull(upper.CL) * 1.1, label = "-", color = treatment_palette["Drought"], size = 6) +
  annotate("text", x = 6, y = CF_plot_df %>% filter(Treatment == "Waterlogged 24" & Segment == 6) %>% pull(upper.CL) * 1.1, label = "*", color = treatment_palette["Waterlogged 24"], size = 6) +
  annotate("text", x = 3.3, y = CF_plot_df %>% filter(Treatment == "Waterlogged 72" & Segment == 3) %>% pull(upper.CL) * 1.1, label = "-", color = treatment_palette["Waterlogged 72"], size = 6) +
  annotate("text", x = 5.3, y = CF_plot_df %>% filter(Treatment == "Waterlogged 72" & Segment == 5) %>% pull(upper.CL) * 1.1, label = "*", color = treatment_palette["Waterlogged 72"], size = 6) +
  annotate("text", x = 6.3, y = CF_plot_df %>% filter(Treatment == "Waterlogged 72" & Segment == 6) %>% pull(upper.CL) * 1.1, label = "*", color = treatment_palette["Waterlogged 72"], size = 6) +
  annotate("text", x = 7.3, y = CF_plot_df %>% filter(Treatment == "Waterlogged 72" & Segment == 7) %>% pull(upper.CL) * 1.1, label = "*", color = treatment_palette["Waterlogged 72"], size = 6) +
  annotate("text", x = 8.3, y = CF_plot_df %>% filter(Treatment == "Waterlogged 72" & Segment == 8) %>% pull(upper.CL) * 1.1, label = "-", color = treatment_palette["Waterlogged 72"], size = 6)

LC_plot_df <- plot(emmeans(LC_mod,~Treatment|Segment), plotit = FALSE) %>% select(Segment, Treatment, upper.CL)
LCA_plot <- ggplot(plot(emmeans(LC_mod,~Treatment|Segment), plotit = FALSE), aes(x = Segment, y = the.emmean, color = Treatment)) +
  ylab("Living Cortical Percent") +
  scale_y_continuous(limits = c(50, 90), breaks = seq(50, 90, by = 10), labels = function(x) paste0(x, "%")) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, fill = Treatment), size = 0.50,
                  linewidth = 1.75, fill = "white", shape = 22, position = position_dodge(width = 0.75)) +
  scale_color_manual(values = treatment_palette, name = "Treatment") +
  xlab("Root Segment (cm) | Root tip to root base") +
  annotate("text", x = 5.3, y = LC_plot_df %>% filter(Treatment == "Waterlogged 72" & Segment == 5) %>% pull(upper.CL) * 1.1, label = "*", color = treatment_palette["Waterlogged 72"], size = 6) +
  annotate("text", x = 2.85, y = LC_plot_df %>% filter(Treatment == "Drought" & Segment == 3) %>% pull(upper.CL) * 1.1, label = "*", color = treatment_palette["Drought"], size = 6) 

RH_CN_LCA_plot <- ggarrange(Root_hair_density_plot, Cortical_file_plot, LCA_plot, nrow = 3, common.legend = TRUE,
                         legend = "right", labels = "AUTO", align = 'hv')


ggsave("figures/Figure5_RH_Cell_file_LCA_plot.svg", RH_CN_LCA_plot, width = 170, height = 200, units = "mm")
ggsave("figures/Figure5_RH_Cell_file_LCA_plot.png", RH_CN_LCA_plot, width = 170, height = 200, units = "mm")

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

