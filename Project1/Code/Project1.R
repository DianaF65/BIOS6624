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
library(naniar)
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

############################### Data for Year 0 & 2 ############################

# Subset data to the first two years
analysis_data <- proj1_dat %>% 
  filter(years %in% c(0, 2))

# Subset to cols of interest
cols <- c("newid", 
          "AGG_MENT", "AGG_PHYS", 
          "VLOAD", "LEU3N",
          "SMOKE", "RACE", "EDUCBAS", "age",
          "ART", "everART", "years",  "hard_drugs",
          "ADH")

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
analysis_data$ADH <- as.factor(analysis_data$ADH)

#### Filter to complete observations - OBS. for BOTH Year 0 and year 2
# Filter to these observations
complete_analysis_data <- analysis_data %>%
  group_by(newid) %>%
  filter(n_distinct(years) == 2) %>%
  ungroup()
# From 715 patients to 506 patients

##### Table1

################################ Missingness ###################################

#### Really nice tutorial: https://cran.r-project.org/web/packages/finalfit/vignettes/missing.html

# Cumulative sum of missingness for each variable
naniar::miss_var_summary(complete_analysis_data)

# How many variables 0 - 5 missing values
naniar::miss_case_table(complete_analysis_data)

# Visualize the missingness
vis_dat(complete_analysis_data)

# Another way to visualize missingness
missing_plot(complete_analysis_data)

## Looking for patterns of missingness

###################################### VLOAD ###################################

explanatory <- c("AGG_MENT", "AGG_PHYS", "SMOKE", "RACE",   
                 "EDUCBAS", "age", "ART", "everART", "years",
                 "hard_drugs", "ADH")
# VLOAD
vload <- c("VLOAD")

complete_analysis_data %>% 
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


### Compare Not Missing and Missing
complete_analysis_data %>% 
  missing_compare(vload, explanatory) %>% 
  knitr::kable(row.names=FALSE, align = c("l", "l", "r", "r", "r")) 

# Based on this, missingness does not appear to differ by any covariates
## Data is not MAR or MCAR

####################################### LEU3N ##################################
# VLOAD
leu3n <- c("LEU3N")

complete_analysis_data %>% 
  missing_pattern(leu3n, explanatory)

# No alarming patterns here

#### Compare Not Missing and Missing
complete_analysis_data %>% 
  missing_compare(leu3n, explanatory) %>% 
  knitr::kable(row.names=FALSE, align = c("l", "l", "r", "r", "r")) 

# No significant differences between missing and not missing values
# Missingness for LEU3N does not differ by covariates
# Not MAR or MCAR

#################################### AGG MENT ##################################
explanatory <- c("AGG_PHYS", "SMOKE", "RACE",   
                 "EDUCBAS", "age", "ART", "everART", "years",
                 "hard_drugs", "ADH")

agg_ment <- c("AGG_MENT")

complete_analysis_data %>% 
  missing_pattern(agg_ment, explanatory)

complete_analysis_data %>% 
  group_by(years) %>% 
  summarise(means = mean(AGG_MENT, na.rm = TRUE))
# 1     0  45.4
# 2     2  47.5
## Not a big difference in means between both years

# Values do not significantly differ between missing and not missing 
# Nothing concerning here

################################### AGG PHYS ###################################
explanatory <- c("AGG_MENT", "SMOKE", "RACE",   
                 "EDUCBAS", "age", "ART", "everART", "years",
                 "hard_drugs", "ADH")
agg_phys <- c("AGG_PHYS")

complete_analysis_data %>% 
  missing_pattern(agg_phys, explanatory)

# No concerning patterns here either
complete_analysis_data %>% 
  group_by(years) %>% 
  summarise(means = mean(AGG_PHYS, na.rm = TRUE))
# 1     0  50.1
# 2     2  49.4
## Not much difference in means between the two years

# Values do not significantly differ between missing and not missing 
# Nothing concerning here



################################# Baseline #####################################
# Subset data to just baseline aka year 1 and one observat
baseline <- complete_analysis_data[complete_analysis_data$years == 0, ]
# There are 280 subjects at baseline

# How many were on hard drugs at baseline?
baseline %>% 
  summarise(prop_hd = mean(hard_drugs, na.rm = TRUE))
# 0.0771 not on hard drugs

table(baseline$hard_drugs)
# 0   1 
# 467  39
# So around 8% were on hard drugs at baseline

################################## VLOAD #######################################

# What is the distribution of VLOAD
summary(complete_analysis_data$VLOAD)

# Plot VLOAD trajectories
ggplot(complete_analysis_data, aes(y = VLOAD, 
                          x = years, 
                          colour = factor(newid))) + 
  geom_line() + 
  geom_point() +
  theme_lucid() + 
  theme(legend.position = "none")

# Notice some outliers here - what are unrealistic values of VLOAD?

# Histogram of the distribution of VLOAD
ggplot(complete_analysis_data, aes(x = VLOAD)) + 
  geom_histogram(bins = 20) +
  theme_lucid() + 
  theme(legend.position = "none")

# Log transform VLOAD
ggplot(complete_analysis_data, aes(x = log(VLOAD))) + 
  geom_histogram(bins = 20,
                 fill = "blue",
                 col = "black") +
  geom_density() + 
  theme_lucid() + 
  theme(legend.position = "none")

## Based on this, might consider mixture modeling
## Frequentist
# https://jef.works/blog/2017/08/05/a-practical-introduction-to-finite-mixture-models/
## Bayesian
# https://rpubs.com/jensroes/mixture-models-tutorial


# Square root
ggplot(complete_analysis_data, aes(x = sqrt(VLOAD))) + 
  geom_histogram() + 
  theme_lucid() + 
  theme(legend.position = "none")

################################################################################

################################## LEU3N #######################################

# What is the distribution of LEU3n
summary(complete_analysis_data$LEU3N)

# Plot VLOAD trajectories
ggplot(complete_analysis_data, aes(y = LEU3N, 
                          x = years, 
                          colour = factor(newid))) + 
  geom_path(aes(group = newid)) + #spaghetti plot
  geom_point() +
  theme_lucid() + 
  theme(legend.position = "none") + 
  facet_grid( ~ hard_drugs)
  
## Look into plotting mean trajectories
ggplot(complete_analysis_data, aes(y = LEU3N, 
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
plot_ids <- complete_analysis_data %>% 
  distinct(newid) %>% 
  sample_n(30) %>% 
  # Obtain just the ids
  pull(newid)

# Subset analysis data frame
plot_data <- complete_analysis_data %>% 
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
ggplot(complete_analysis_data, aes(x = LEU3N)) + 
  geom_histogram(bins = 20,
                 fill = "seagreen3",
                 col = "black") +
  theme_lucid() + 
  theme(legend.position = "none")

# Log transform VLOAD
ggplot(complete_analysis_data, aes(x = log(LEU3N))) + 
  geom_histogram(bins = 20) +
  theme_lucid() + 
  theme(legend.position = "none")
## Log transformation too strong - conduct a test of normality to determine
shapiro.test(complete_analysis_data$LEU3N) 
# We reject the H0 and data is not normally distributed


# Sqrt transform VLOAD
ggplot(complete_analysis_data, aes(x = sqrt(LEU3N))) + 
  geom_histogram(bins = 20) +
  theme_lucid() + 
  theme(legend.position = "none")

################################################################################

################################## AGG MENT ####################################

summary(complete_analysis_data$AGG_MENT)

# Histogram of distribution
ggplot(complete_analysis_data, aes(x = (AGG_MENT))) + 
  geom_histogram(fill = "purple", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none") 

# Log transform
ggplot(complete_analysis_data, aes(x = sqrt(AGG_MENT))) + 
  geom_histogram(fill = "purple", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none") 

### Based on this consider beta regression
# Would have to revert percentages to decimal values
complete_analysis_data$prop_AGG_MENT <- complete_analysis_data$AGG_MENT / 100

# Histogram of this
ggplot(complete_analysis_data, aes(x = (prop_AGG_MENT))) + 
  geom_histogram(fill = "purple", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none") 

################################## AGG PHYS ####################################

# What is the distribution of Aggregate phsycial quality of life score
summary(complete_analysis_data$AGG_PHYS)

# Histogram of distribution
ggplot(complete_analysis_data, aes(x = (AGG_PHYS))) + 
  geom_histogram(fill = "orange", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none")

# We also observe a left skew here
### We will also consider beta regression with this
# Convert percentages to decimal values
complete_analysis_data$prop_AGG_PHYS <- complete_analysis_data$AGG_PHYS / 100

# Histogram of this
ggplot(complete_analysis_data, aes(x = (prop_AGG_PHYS))) + 
  geom_histogram(fill = "orange", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none")

summary(complete_analysis_data$prop_AGG_PHYS)


################################################################################

################################# Covariates ###################################

################################## Hard Drugs  #################################

table(complete_analysis_data$hard_drugs)
# 0   1 
# 649  66
# > table(baseline$hard_drugs)

0   1 
467  39 

## Also include for people that have follow up to 2 years
ggplot(baseline, aes(x = factor((hard_drugs)),
                          fill = factor((hard_drugs)))) + 
  geom_bar() + 
  theme_lucid() + 
  theme(legend.position = "none")

################################### Ever ART ###################################

table(complete_analysis_data$everART)
table(baseline$everART)

####################################### Age  ###################################

# Summary of age at baseline
summary(baseline$age)

# Histogram of distribution
ggplot(baseline, aes(x = age)) + 
  geom_histogram(col = "black", fill = "yellow") + 
  theme_lucid()


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

# Collapse into these categories
# White
# Black
# American Indian / Alaska Native
# Asian / Pacific Islander
# Other


##################################### SMOKE  ###################################
table(factor(baseline$SMOKE))
# 1   2   3 
# 71  81 128 
table(factor(complete_analysis_data$SMOKE))
# 1   2   3 
# 120 142 202 

ggplot(baseline, aes(x = SMOKE)) + 
  geom_bar() + 
  theme_lucid()

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


################################## Adherence ###################################

table(complete_analysis_data$ADH)

# Barplot
ggplot(complete_analysis_data, aes(x = ADH)) + 
  geom_bar() + 
  theme_lucid()

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

