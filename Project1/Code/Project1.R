################################################################################
###                                 Project 0                                ###
################################################################################
###                         Initial Data Exploration                         ###
################################################################################

# Load libraries
library(readr)
library(knitr)
library(ggplot2)
library(MASS)
library(emmeans)
library(multcomp)
library(kableExtra)
library(stats)
library(dplyr)
library(tibble)
library(labelled)
library(tableone)
library(lme4)
library(lmerTest)
library(nlme)
library(gtsummary)
library(gt)
library(ggeffects)
library(UpSetR)
library(visdat)
library(lubridate)
library(see)
library(data.table)
library(gridExtra)
library(kableExtra)

##### Read in Data
proj1_dat <- read.csv("../Data/hiv_6624_final.csv")

################################ Missingness ###################################

# Cumulative sum of missingness for each variable
naniar::miss_var_summary(proj1_dat)

# How many variables 0 - 5 missing values
naniar::miss_case_table(proj1_dat)

# Visualize the missingness
vis_dat(proj1_dat)

################################### Table 1 ####################################

# Adjust labels for collection sample
# data <- wake_data %>%
#   mutate(
#     `Collection Sample` = factor(
#       `Collection Sample`,
#       levels = c(1, 2, 3, 4),
#       labels = c("Wake Time", "+30 Min", "Pre Lunch", "+600 Min")
#     )
#   )

# Subset proj1 data to years 0 - 2
data <- proj1_dat %>% 
  filter(years %in% c(0, 1, 2))

# Create subject level summaries - NOT observation level
data <- data %>%
  group_by(newid, years) %>%
  summarise(
    mean_cesd = mean(CESD, na.rm = TRUE),
    mean_vload = mean(VLOAD, na.rm = TRUE),
    mean_bmi = mean(BMI, na.rm = TRUE),
    mean_agg_ment = mean(AGG_MENT, na.rm = TRUE),
    mean_hbp = mean(HBP, na.rm = TRUE),
    mean_tchol = mean(TCHOL, na.rmm = TRUE),
    mean_agg_phys = mean(AGG_PHYS, na.rm = TRUE),
    mean_diab = mean(DIAB, na.rm = TRUE),
    mean_trig= mean(TRIG, na.rm = TRUE),
    mean_ldl = mean(LDL, na.rm = TRUE),
    mean_leu3n = mean(LEU3N, na.rm = TRUE),
    mean_age = mean(age, na.rm = TRUE),
    .groups = "drop")

# Define variables for the table one
variables <- colnames(data)

# Separate by collection time
strata <- "years"

# modify labels
# var_label(data$mean_booklet) <- "Booklet Minutes"
# var_label(data$mean_mem) <- "MEM Minutes"
# var_label(data$mean_wake) <- "Minutes Since Wake"
# var_label(data$mean_dhea) <- "DHEA (nmol/L)"
# var_label(data$mean_cort) <- "Cortisol (nmol/L)"

# Create the table one
studypop_tab1 <- CreateTableOne(
  vars = variables,
  strata = strata,
  test = FALSE,
  data = data
)

studypop_tab_forkbl <- print(
  studypop_tab1,
  varLabels = TRUE,
  printToggle = FALSE
)

# Display with kbl
kbl(
  studypop_tab_forkbl,
  caption = "Summary Per Subject For Each Year",
  booktabs = TRUE,
  align = "c"
) %>%
  kable_styling(latex_options = "hold_position", font_size = 10)

# HARD DRUGS AT BASELINE
# Subset data to just baseline aka year 1 and one observat
baseline <- proj1_dat[proj1_dat$years == 0, ]

# How many were on hard drugs at baseline?
table(baseline$hard_drugs)
# 0   1 
# 649  66

baseline %>% 
  summarise(mean(hard_drugs, na.rm = TRUE))
# Approximately 10% were on hard drugs on baseline

# Pretty unbalanced that could be problematic with analysis

# Subset data to the first two years
analysis_data <- proj1_dat %>% 
  filter(years %in% c(0, 1, 2))

################################## VLOAD #######################################

# What is the distribution of VLOAD
summary(analysis_data$VLOAD)

# Plot VLOAD trajectories
ggplot(analysis_data, aes(y = VLOAD, 
                          x = years, 
                          colour = factor(newid))) + 
  geom_line() + 
  geom_point() +
  theme_lucid() + 
  theme(legend.position = "none")

# Notice some outliers here - what are unrealistic values of VLOAD?

# Filter to those that have VLOAD higher than the 3rd quartile
analysis_data %>% 
  filter(VLOAD >= 30929) %>% 
  reframe(newid, years, VLOAD) %>% 
  arrange(desc(VLOAD))

# Plot without these outliers for a bit
filtered <- analysis_data %>% 
  filter(VLOAD < 2520009)

# Plot
ggplot(filtered,
       aes(x = years, y = VLOAD, colour = factor(newid))) + 
  geom_line() + 
  geom_point() +
  theme_lucid() + 
  theme(legend.position = "none")

################################################################################

################################## LEU3N #######################################

# What is the distribution of VLOAD
summary(analysis_data$LEU3N)

# Plot VLOAD trajectories
ggplot(analysis_data, aes(y = LEU3N, 
                          x = years, 
                          colour = factor(newid))) + 
  geom_path(aes(group = newid)) + #spaghetti plot
  geom_point() +
  theme_lucid() + 
  theme(legend.position = "none") + 
  facet_grid( ~ hard_drugs)
  
## Look into plotting mean trajectories
ggplot(analysis_data, aes(y = LEU3N, 
                          x = years, 
                          colour = factor(newid))) + 
  geom_path(aes(group = newid)) + #spaghetti plot
  geom_point() +
  theme_lucid() + 
  theme(legend.position = "none") + 
  facet_grid( ~ hard_drugs) +
  stat_summary(fun.y = mean,
               geom = "line",
               lwd = 1, aes(group = 1))


############## Randomly sample like 20 subjects
set.seed(645)
plot_ids <- analysis_data %>% 
  distinct(newid) %>% 
  sample_n(30) %>% 
  # Obtain just the ids
  pull(newid)

# Subset analysis data frame
plot_data <- analysis_data %>% 
  filter(newid %in% plot_ids) %>% 
  filter(!is.na(LEU3N))

# Plot
# Plot VLOAD trajectories
ggplot(plot_data, aes(y = LEU3N, 
                          x = years, 
                          colour = factor(newid))) + 
  geom_line() + 
  geom_point() +
  theme_lucid() + 
  theme(legend.position = "none") + 
  facet_grid( ~ hard_drugs)
# On initial glance, looking like CD4 counts increasing








