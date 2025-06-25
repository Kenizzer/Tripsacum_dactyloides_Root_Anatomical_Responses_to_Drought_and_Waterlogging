### Section Randomization testing
### Code by Joel Swift 06/2025

library(tidyverse)
library(readxl)
library(emmeans)
library(ggpubr)
library(lmerTest)

# Set theme
theme_set(theme_bw())

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

# Ratios (X:Y)
section_data$`Stele:Root Ratio` <- section_data$`Stele Area` / section_data$`Root Area`
section_data$`Metaxylem:Stele Ratio` <- section_data$`Metaxylem Vessel Area` /section_data$`Stele Area`

# Scale root hair number by perimeter measurement 
section_data$`Root Hair Density` <- section_data$`Root Hair Number` / (section_data$`Root Perimeter`/1000) # Root hairs per mm of root perimeter

# Root Diameter
section_data$Root_diameter <- (sqrt(section_data$`Root Area` / pi) * 2)
section_data$Stele_diameter <- (sqrt(section_data$`Stele Area` / pi) * 2)

# Cortical Cell file size normalized by cortex depth
section_data$`CCFN-Normalized` <- section_data$`Cortical Cell File Number` / ((section_data$Root_diameter - section_data$Stele_diameter) / 2)

# Sanity check
str(section_data)

# Select all plants that have three sections per segment (i.e, have a replicate A, B, and C)
section_data_index <- section_data %>% filter(Replicate %in% c("B", "C"))
section_data_replicated <- section_data %>% filter(Plant_id %in% unique(section_data_index$Plant_id))
rm(section_data_index)


#### Randomization testing R2 ####
section_data_replicated_clean <- section_data_replicated %>% clean_names()
averaged_sections_clean <- averaged_sections %>% clean_names()

# Create a function for one iteration
resample_model_R2 <- function(trait, df) {
  
  # Stitch together formula
  formula <- as.formula(paste0("`", trait, "` ~ treatment*segment + (1|plant_id)"))
  
  # Sample one replicate section per Plant_id + Segment
  sampled_data <- df %>%
    group_by(plant_id, segment) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  # Fit model (adjust the formula to fit your data)
  model <- lmer(formula, data = sampled_data)
  
  # Extract marginal R² (fixed effects only)
  r2 <- performance::r2(model)$R2_marginal
  
  return(r2)
}

# Define number of iterations
n_iter <- 1000
set.seed(08123850)  # for reproducibility
RA_R2 <- map_dbl(1:n_iter, ~ resample_model_R2("root_area", section_data_replicated_clean))
ED_R2 <- map_dbl(1:n_iter, ~ resample_model_R2("epidermis_depth", section_data_replicated_clean))
RP_R2 <- map_dbl(1:n_iter, ~ resample_model_R2("root_perimeter", section_data_replicated_clean))
RD_R2 <- map_dbl(1:n_iter, ~ resample_model_R2("root_hair_density", section_data_replicated_clean))
CC_R2 <- map_dbl(1:n_iter, ~ resample_model_R2("cortical_cell_file_number", section_data_replicated_clean))
CN_R2 <- map_dbl(1:n_iter, ~ resample_model_R2("ccfn_normalized", section_data_replicated_clean))
LC_R2 <- map_dbl(1:n_iter, ~ resample_model_R2("living_cortical_percent", section_data_replicated_clean))
AR_R2 <- map_dbl(1:n_iter, ~ resample_model_R2("aerenchyma_percent", section_data_replicated_clean))
SR_R2 <- map_dbl(1:n_iter, ~ resample_model_R2("stele_root_ratio", section_data_replicated_clean))
MS_R2 <- map_dbl(1:n_iter, ~ resample_model_R2("metaxylem_stele_ratio", section_data_replicated_clean))
MV_R2 <- map_dbl(1:n_iter, ~ resample_model_R2("metaxylem_vessel_number", section_data_replicated_clean))
MA_R2 <- map_dbl(1:n_iter, ~ resample_model_R2("metaxylem_vessel_mean_area", section_data_replicated_clean))

# Observed Values
RA_OBS <- lmer(root_area ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)
RP_OBS <- lmer(root_perimeter ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)
SA_OBS <- lmer(stele_area ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)
CC_OBS <- lmer(cortical_cell_file_number ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)
ED_OBS <- lmer(epidermis_depth ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)
AR_OBS <- lmer(aerenchyma_percent ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)
LC_OBS <- lmer(living_cortical_percent ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)
SR_OBS <- lmer(stele_root_ratio ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)
MS_OBS <- lmer(metaxylem_stele_ratio ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)
RD_OBS <- lmer(root_hair_density ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)
CN_OBS <- lmer(ccfn_normalized ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)
MV_OBS <- lmer(metaxylem_vessel_number ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)
MA_OBS <- lmer(metaxylem_vessel_mean_area ~ treatment*segment + (1|plant_id), data = averaged_sections_clean)

colnames(section_data_replicated)

mean(RA_R2 >= performance::r2(RA_OBS)$R2_marginal)
a <- ggplot(data.frame(R2 = RA_R2), aes(x = R2)) +
  geom_histogram(bins = 30, fill = "grey80", color = "black") +
  geom_vline(xintercept = performance::r2(RA_OBS)$R2_marginal, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Root Area", x = "R2", y = "Count") +
  annotate("text", x = Inf, y = Inf, label = paste0("p = ", signif(mean(RA_R2 >= performance::r2(RA_OBS)$R2_marginal), 2)),
           vjust = 2, hjust = 1.1, size = 5, color = "black")


mean(ED_R2 >= performance::r2(ED_OBS)$R2_marginal)
b <- ggplot(data.frame(R2 = ED_R2), aes(x = R2)) +
  geom_histogram(bins = 30, fill = "grey80", color = "black") +
  geom_vline(xintercept = performance::r2(ED_OBS)$R2_marginal, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Epidermis Depth", x = "R2", y = "Count") +
  annotate("text", x = Inf, y = Inf, label = paste0("p = ", signif(mean(ED_R2 >= performance::r2(ED_OBS)$R2_marginal), 2)),
           vjust = 2, hjust = 1.1, size = 5, color = "black")


mean(RD_R2 >= performance::r2(RD_OBS)$R2_marginal)
c <- ggplot(data.frame(R2 = RD_R2), aes(x = R2)) +
  geom_histogram(bins = 30, fill = "grey80", color = "black") +
  geom_vline(xintercept = performance::r2(RD_OBS)$R2_marginal, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Root Hair Density", x = "R2", y = "Count") +
  annotate("text", x = Inf, y = Inf, label = paste0("p = ", signif(mean(RD_R2 >= performance::r2(RD_OBS)$R2_marginal), 2)),
           vjust = 2, hjust = 1.1, size = 5, color = "black")


mean(CC_R2 >= performance::r2(CC_OBS)$R2_marginal)
d <- ggplot(data.frame(R2 = CC_R2), aes(x = R2)) +
  geom_histogram(bins = 30, fill = "grey80", color = "black") +
  geom_vline(xintercept = performance::r2(CC_OBS)$R2_marginal, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Cortical Cell File Number", x = "R2", y = "Count") +
  annotate("text", x = Inf, y = Inf, label = paste0("p = ", signif(mean(CC_R2 >= performance::r2(CC_OBS)$R2_marginal), 2)),
           vjust = 2, hjust = 1.1, size = 5, color = "black")


mean(CN_R2 >= performance::r2(CN_OBS)$R2_marginal)
e <- ggplot(data.frame(R2 = CN_R2), aes(x = R2)) +
  geom_histogram(bins = 30, fill = "grey80", color = "black") +
  geom_vline(xintercept = performance::r2(CN_OBS)$R2_marginal, color = "red", linetype = "dashed", size = 1) +
  labs(title = "CCFN-Normalized", x = "R2", y = "Count") +
  annotate("text", x = Inf, y = Inf, label = paste0("p = ", signif(mean(CN_R2 >= performance::r2(CN_OBS)$R2_marginal), 2)),
           vjust = 2, hjust = 1.1, size = 5, color = "black")


mean(LC_R2 >= performance::r2(LC_OBS)$R2_marginal)
f <- ggplot(data.frame(R2 = LC_R2), aes(x = R2)) +
  geom_histogram(bins = 30, fill = "grey80", color = "black") +
  geom_vline(xintercept = performance::r2(LC_OBS)$R2_marginal, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Living Cortical Percent", x = "R2", y = "Count") +
  annotate("text", x = Inf, y = Inf, label = paste0("p = ", signif(mean(LC_R2 >= performance::r2(LC_OBS)$R2_marginal), 2)),
           vjust = 2, hjust = 1.1, size = 5, color = "black")


mean(AR_R2 >= performance::r2(AR_OBS)$R2_marginal)
g <- ggplot(data.frame(R2 = AR_R2), aes(x = R2)) +
  geom_histogram(bins = 30, fill = "grey80", color = "black") +
  geom_vline(xintercept = performance::r2(AR_OBS)$R2_marginal, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Aerenchyma Percent", x = "R2", y = "Count") +
  annotate("text", x = Inf, y = Inf, label = paste0("p = ", signif(mean(AR_R2 >= performance::r2(AR_OBS)$R2_marginal), 2)),
           vjust = 2, hjust = 1.1, size = 5, color = "black")

mean(SR_R2 >= performance::r2(SR_OBS)$R2_marginal)
h <- ggplot(data.frame(R2 = SR_R2), aes(x = R2)) +
  geom_histogram(bins = 30, fill = "grey80", color = "black") +
  geom_vline(xintercept = performance::r2(SR_OBS)$R2_marginal, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Stele:Root Ratio", x = "R2", y = "Count") +
  annotate("text", x = Inf, y = Inf, label = paste0("p = ", signif(mean(SR_R2 >= performance::r2(SR_OBS)$R2_marginal), 2)),
           vjust = 2, hjust = 1.1, size = 5, color = "black")

mean(MS_R2 >= performance::r2(MS_OBS)$R2_marginal)
i <- ggplot(data.frame(R2 = MS_R2), aes(x = R2)) +
  geom_histogram(bins = 30, fill = "grey80", color = "black") +
  geom_vline(xintercept = performance::r2(MS_OBS)$R2_marginal, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Metaxylem:Stele Ratio", x = "R2", y = "Count") +
  annotate("text", x = Inf, y = Inf, label = paste0("p = ", signif(mean(MS_R2 >= performance::r2(MS_OBS)$R2_marginal), 2)),
           vjust = 2, hjust = 1.1, size = 5, color = "black") 

mean(MA_R2 >= performance::r2(MA_OBS)$R2_marginal)
j <- ggplot(data.frame(R2 = MA_R2), aes(x = R2)) +
  geom_histogram(bins = 30, fill = "grey80", color = "black") +
  geom_vline(xintercept = performance::r2(MA_OBS)$R2_marginal, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Metaxylem Vessel Mean Area", x = "R2", y = "Count") +
  annotate("text", x = Inf, y = Inf, label = paste0("p = ", signif(mean(MA_R2 >= performance::r2(MA_OBS)$R2_marginal), 2)),
           vjust = 2, hjust = 1.1, size = 5, color = "black") 


R2_plots <- ggarrange(a,b,c,d,e,f,g,h,i,j, align = "hv", nrow = 2, ncol = 5, labels = "AUTO")

ggsave("figures/Figure_S3_R2_plots.svg", R2_plots, width = 16, height = 8, units = "in")
ggsave("figures/Figure_S3_R2_plots.png", R2_plots, width = 16, height = 8, units = "in")
