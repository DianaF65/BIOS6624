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
library(cmdstanr)
library(bayesplot)
library(posterior)
library(bayestestR)
library(mcmcse)
library(loo)
library(finalfit)
color_scheme_set("brightblue")


##### Read in Data
proj1_dat <- read.csv("../Data/hiv_6624_final.csv")

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
  filter(years %in% c(0, 2))

# Create subject level summaries - NOT observation level
data <- data %>%
  group_by(newid, years) %>%
  summarise(
    mean_vload = mean(VLOAD, na.rm = TRUE),
    mean_bmi = mean(BMI, na.rm = TRUE),
    mean_agg_ment = mean(AGG_MENT, na.rm = TRUE),
    mean_agg_phys = mean(AGG_PHYS, na.rm = TRUE),
    mean_leu3n = mean(LEU3N, na.rm = TRUE),
    mean_agg_ment = mean(AGG_MENT, na.rm = TRUE),
    mean_agg_phys = mean(AGG_PHYS, na.rm = TRUE),
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


################################## Analysis DF #################################
# Subset data to the first two years
analysis_data <- proj1_dat %>% 
  filter(years %in% c(0, 2))

# Subset to cols of interest
cols <- c("newid", 
          "AGG_MENT", "AGG_PHYS", 
          "VLOAD", "LEU3N",
          "SMOKE", "RACE", "EDUCBAS", "age",
          "ART", "everART", "years",  "hard_drugs")

# Subset analysis data to these cols
analysis_data <- analysis_data[ ,colnames(analysis_data) %in% cols ]

# Ensure that data types are correct
str(analysis_data)

### Convert to correct data types
## Factors
# new id
analysis_data$newid <- as.factor(analysis_data$newid)
analysis_data$SMOKE <- as.factor(analysis_data$SMOKE)
analysis_data$RACE <- as.factor(analysis_data$RACE)
analysis_data$EDUCBAS <- as.factor(analysis_data$EDUCBAS)
analysis_data$ART <- as.factor(analysis_data$ART)
analysis_data$everART <- as.factor(analysis_data$everART)

##### Table1

################################ Missingness ###################################

#### Really nice tutorial: https://cran.r-project.org/web/packages/finalfit/vignettes/missing.html

# Cumulative sum of missingness for each variable
naniar::miss_var_summary(analysis_data)

# How many variables 0 - 5 missing values
naniar::miss_case_table(analysis_data)

# Visualize the missingness
vis_dat(analysis_data)

# Another way to visualize missingness
missing_plot(analysis_data)

## Looking for patterns of missingness
#### VLOAD
explanatory <- c("AGG_MENT", "AGG_PHYS", "SMOKE", "RACE",   
                 "EDUCBAS", "age", "ART", "everART", "years",
                 "hard_drugs")
# VLOAD
vload <- c("VLOAD")

analysis_data %>% 
missing_pattern(vload, explanatory)

# This plot lets us see which variables are missing together. 
# There does not seem to be a relationship with
# missingness between any pair of predictors or predictor and outcome.
# Each row represents a unique combination of observed vs. missing values. 
# In other words, a row is a pattern
# shared by many people.
# The numbers on the right are the number of missing variables.
# The numbers on the bottom are the number of missing observations missing for each variable.
# The numbers on the left shows the number of cases with a specific pattern.
# For example, looking at the bottom row, there is 1 case where 3 variables are missing,
# years, SMOKE, and VLOAD
# With the first row, there are 149 cases with the “0” missing pattern or no variables missing.
# Overall, we don’t see any alarmingly large numbers on the left for particular missing patterns

### LEU3N
# VLOAD
leu3n <- c("LEU3N")

analysis_data %>% 
  missing_pattern(leu3n, explanatory)

# No alarming patterns here

### AGG MENT
explanatory <- c("AGG_PHYS", "SMOKE", "RACE",   
                 "EDUCBAS", "age", "ART", "everART", "years",
                 "hard_drugs")

agg_ment <- c("AGG_MENT")

analysis_data %>% 
  missing_pattern(agg_ment, explanatory)

# Nothing concerning here

### AGG PHYS
explanatory <- c("AGG_MENT", "SMOKE", "RACE",   
                 "EDUCBAS", "age", "ART", "everART", "years",
                 "hard_drugs")
agg_phys <- c("AGG_PHYS")

analysis_data %>% 
  missing_pattern(agg_phys, explanatory)

# No concerning patterns here either

## Based on this brief analyssis, we will not continue with missingness investigation


################################# Baseline #####################################
# Subset data to just baseline aka year 1 and one observat
baseline <- analysis_data[analysis_data$years == 0, ]
# There are 280 subjects at baseline

# How many were on hard drugs at baseline?
table(baseline$hard_drugs)
# 0   1 
# 252  28 

baseline %>% 
  summarise(mean(hard_drugs, na.rm = TRUE))
# Approximately 10% were on hard drugs on baseline

# Pretty unbalanced that could be problematic with analysis

############################### Data for Year 0 & 2 ############################

# Subset data to the first two years
analysis_data <- proj1_dat %>% 
  filter(years %in% c(0, 2))

# Subset to variables of interest
cols <- c("newid", "AGG_MENT", "AGG_PHYS", "income",
          "SMOKE", "LEU3N", "VLOAD",
          "RACE", "EDUCBAS", "age",
          "ART", "everART", "years", "hard_drugs")
analysis_data <- analysis_data[ , colnames(analysis_data) %in% cols]

### Missingness here
# Cumulative sum of missingness for each variable
naniar::miss_var_summary(analysis_data)

# How many variables 0 - 5 missing values
naniar::miss_case_table(analysis_data)

# Visualize the missingness
vis_dat(analysis_data)

#################################### Year 2 ####################################
year2 <- analysis_data[analysis_data$years == 2, ]
# There are 185 individuals at year 2
## Suppose we will have to subset analysis data to these 185 individuals..?

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

# Histogram of the distribution of VLOAD
ggplot(analysis_data, aes(x = VLOAD)) + 
  geom_histogram(bins = 20) +
  theme_lucid() + 
  theme(legend.position = "none")

# Log transform VLOAD
ggplot(analysis_data, aes(x = log(VLOAD))) + 
  geom_histogram(bins = 20) +
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

# Histogram of the distribution of VLOAD
ggplot(analysis_data, aes(x = LEU3N)) + 
  geom_histogram(bins = 20) +
  theme_lucid() + 
  theme(legend.position = "none")

# Log transform VLOAD
ggplot(analysis_data, aes(x = log(LEU3N))) + 
  geom_histogram(bins = 20) +
  theme_lucid() + 
  theme(legend.position = "none")

# Sqrt transform VLOAD
ggplot(analysis_data, aes(x = sqrt(LEU3N))) + 
  geom_histogram(bins = 20) +
  theme_lucid() + 
  theme(legend.position = "none")

################################################################################

################################## AGG PHYS ####################################

# What is the distribution of Aggregate phsycial quality of life score
summary(analysis_data$AGG_PHYS)

# Histogram of distribution
ggplot(analysis_data, aes(x = na.omit(AGG_PHYS))) + 
  geom_histogram() + 
  theme_lucid() + 
  theme(legend.position = "none") + 
  facet_grid( ~ hard_drugs)



################################## AGG MENT ####################################

summary(analysis_data$AGG_MENT)

################################################################################

################################# Covariates ###################################

################################## Hard Drugs  #################################

table(baseline$hard_drugs)

ggplot(baseline, aes(x = factor(na.omit(hard_drugs)),
                          fill = factor(na.omit(hard_drugs)))) + 
  geom_bar() + 
  theme_lucid() + 
  theme(legend.position = "none")

################################### Ever ART ###################################

table(analysis_data$everART)

####################################### Age  ###################################

# Summary of age at baseline
summary(baseline$age)


###################################### Race  ###################################

# 1= White, non-Hispanic
# 2= White, Hispanic
# 3= Black, non-Hispanic
# 4= Black, Hispanic
# 5= American Indian or Alaskan
# Native
# 6= Asian or Pacific Islander
# 7= Other
# 8= Other Hispanic (created
#                    for 2001-03 new recruits)
# Blank= Missing

# Table of distribution of race
table(factor(baseline$RACE))

# Barplot
ggplot(baseline, aes(x = factor(na.omit(RACE)),
                     fill = factor(RACE))) + 
  geom_bar() + 
  theme_lucid() + 
  theme(legend.position = "none")
# Pretty unbalanced


<<<<<<< HEAD
################################## Education  #################################




################################# Smoking Status ##############################
=======
##################################### SMOKE  ###################################
table(factor(baseline$SMOKE))
# 1   2   3 
# 71  81 128 
table(factor(analysis_data$SMOKE))
# 1   2   3 
# 120 142 202 

##################################### EDUCBAS  #################################

# 1= 8th grade or less
# 2= 9, 10, or 11th grade
# 3= 12th grade
# 4= At least one year college
# but no degree
# 5= Four years college / got
# degree
# 6= Some graduate work
# 7= Post-graduate degree
# Blank= Missing


table(factor(baseline$EDUCBAS))
# 1   2   3   4   5   6   7 
# 2  20  49 104  54  18  33 

# Barplot
ggplot(baseline, aes(x = factor(na.omit(EDUCBAS)),
                     fill = factor(na.omit(EDUCBAS)))) + 
  geom_bar() + 
  theme_lucid() + 
  theme(legend.position = "none")
# Pretty unbalanced
>>>>>>> 08e6a7d73543f005202622615e0ab35d4926c367



################################ STAN INSTALLATION #############################

### Checking if its connectivity from inside R
# Basic connectivity check
capabilities("libcurl")
getOption("download.file.method")

# Try to hit GitHub API directly
u <- "https://api.github.com/repos/stan-dev/cmdstan/releases/latest"
try(readLines(u, n = 5), silent = FALSE)

# Specify a certain version
cmdstanr::install_cmdstan(version = "2.37.0", cores = 2)
### This worked

############################# CMDSTAN Vignette #################################

file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file)

# names correspond to the data block in the Stan program
data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)

