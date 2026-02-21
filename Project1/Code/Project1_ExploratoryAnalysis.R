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

################################################################################
###                        Data Exploration and Cleaning                     ###
################################################################################

############################### Data for Year 0 & 2 ############################

# Subset data to the first two years
analysis_data <- proj1_dat %>% 
  filter(years %in% c(0, 2))

# Subset to cols of interest
cols <- c("newid", 
          "AGG_MENT", "AGG_PHYS", 
          "VLOAD", "LEU3N",
          "SMOKE", "RACE", "EDUCBAS", "age",
          "BMI",
          "ART", "everART", "years",  "hard_drugs",
          "ADH")

# Subset analysis data to these cols
analysis_data <- analysis_data[ ,colnames(analysis_data) %in% cols ]

# Ensure that data types are correct
str(analysis_data)

##### Table1

#### First glimpse of the data
glimpse(analysis_data)

# Summary of the variables in the data frame
summary(analysis_data)

# We note some concerning values for BMI: -1 and 999
summary(analysis_data$BMI)

# What is an unreasonable BMI?
# Based on this google and this article: 
# https://pmc.ncbi.nlm.nih.gov/articles/PMC2930234/#:~:text=The%20BMI%20is%20a%20surrogate,upper%20limit%20of%20'starvation'.
# We will filter out individuals that have a BMI within the ranges of (10, 60)
# Anything less than 10 indicates severe starvation and > 60 is extreme obesity

# Filter out unrealistic bmi values from data set
# We only care about baseline values
# Create data frame to year 1 then filter out these BMI values
year0 <- analysis_data %>% 
  # Filter to baseline
  filter(years == 0) %>% # 715 obs and 715 subjects
  # Filter BMI values
  filter(BMI > 10 & BMI < 60) # 671 obs and 671 subjects
# Double check this worked
summary(year0$BMI) 

# Create data frame with year 2 data
year2 <- analysis_data %>% 
  filter(years == 2)
# 506 obs and 506 subjects

# Merge these two data frames by id
merged_data <- rbind(year0, year2)

# Organize df so subjects are together
merged_data <- merged_data %>% 
  # group_by(newid) %>% # This does something weird with the df...
  arrange(newid)

### Clean environment
rm(year0)
rm(year2)

# Merged data will have BMIs out of this range
summary(merged_data$BMI)

### From 1221 obs to 1177 obs
### From 715 subjects to 689 subjects

################################ Missingness ###################################

#### Really nice tutorial: https://cran.r-project.org/web/packages/finalfit/vignettes/missing.html

# Cumulative sum of missingness for each variable
naniar::miss_var_summary(merged_data)
# Looks like ADH has a lot of missingness but this is expected since adh is
# measured at second visit

# How many variables 0 - 5 missing values
naniar::miss_case_table(merged_data)

# Visualize the missingness
vis_dat(merged_data)

# Another way to visualize missingness
missing_plot(merged_data)

###### Investigating any patterns of missingness below

########################## Assessing patterns of missingness ###################

# Explore patterns of missingness to determine if data is MNAR, MAR, MCAR

###################################### VLOAD ###################################

# Explore patterns of missingness for VLOAD

explanatory <- c("AGG_PHYS", "AGG_MENT", "SMOKE", "RACE",   
                 "EDUCBAS", "age", "years",
                 "BMI",
                 "hard_drugs", "ADH")
# VLOAD
vload <- c("VLOAD")

merged_data %>% 
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


### Compare Not Missing and Missing for 
merged_data %>% 
  missing_compare(vload, explanatory) %>% 
  knitr::kable(row.names=FALSE, align = c("l", "l", "r", "r", "r")) 

## Looking at some categorical variables since there were some issues
# Usingn simulated pvalues from chisquare test
# SMOKE 
merged_data %>% 
  dplyr::summarise(
    pval_smoke = chisq.test(SMOKE, VLOAD, simulate.p.value = TRUE)$p.value,
    pval_race = chisq.test(RACE, VLOAD, simulate.p.value = TRUE)$p.value,
    pval_hdrugs = chisq.test(hard_drugs, VLOAD, simulate.p.value = TRUE)$p.value,
    pval_adh = chisq.test(ADH, VLOAD, simulate.p.value = TRUE)$p.value
  )
# There appears to be a significant difference between expected and obs values
# For SMOKE and ADH
# However, we will not delve too much more into this. 
# The sample sizes for some categorical variables are extremley unbalanced.
# Based on this, missingness does not appear to differ by any covariates
## Data is not MAR or MCAR

####################################### LEU3N ##################################
# VLOAD
leu3n <- c("LEU3N")

merged_data %>% 
  missing_pattern(leu3n, explanatory)

# No alarming patterns here

#### Compare Not Missing and Missing
merged_data %>% 
  missing_compare(leu3n, explanatory) %>% 
  knitr::kable(row.names=FALSE, align = c("l", "l", "r", "r", "r")) 

# No significant differences between missing and not missing values
# Missingness for LEU3N does not differ by covariates
# Not MAR or MCAR

## Looking at some categorical variables since there were some issues
# Usingn simulated pvalues from chisquare test
# SMOKE 
merged_data %>% 
  dplyr::summarise(
    pval_smoke = chisq.test(SMOKE, LEU3N, simulate.p.value = TRUE)$p.value,
    pval_race = chisq.test(RACE, LEU3N, simulate.p.value = TRUE)$p.value,
    pval_hdrugs = chisq.test(hard_drugs, LEU3N, simulate.p.value = TRUE)$p.value,
    pval_adh = chisq.test(ADH, LEU3N, simulate.p.value = TRUE)$p.value
  )
# There appears to be a significant difference between expected and obs values
# For hard drugs and ADH
# However, we will not delve too much more into this. 
# The sample sizes for some categorical variables are extremley unbalanced.
# Based on this, missingness does not appear to differ by any covariates
## Data is not MAR or MCAR

#################################### AGG MENT ##################################
explanatory <- c("SMOKE", "RACE",   
                 "EDUCBAS", "age", "ART", "everART", "years",
                 "BMI",
                 "hard_drugs", "ADH")

agg_ment <- c("AGG_MENT")

merged_data %>% 
  missing_pattern(agg_ment, explanatory)

merged_data %>% 
  group_by(years) %>% 
  summarise(means = mean(AGG_MENT, na.rm = TRUE))
# 1     0  45.4
# 2     2  47.5
## Not a big difference in means between both years

# Values do not significantly differ between missing and not missing 
# Nothing concerning here

################################### AGG PHYS ###################################
explanatory <- c("SMOKE", "RACE",   
                 "EDUCBAS", "age", "years", "BMI",
                 "hard_drugs", "ADH")
agg_phys <- c("AGG_PHYS")

merged_data %>% 
  missing_pattern(agg_phys, explanatory)

# No concerning patterns here either
merged_data %>% 
  group_by(years) %>% 
  summarise(means = mean(AGG_PHYS, na.rm = TRUE))
# 1     0  50.1
# 2     2  49.4
## Not much difference in means between the two years

#### Compare Not Missing and Missing
merged_data %>% 
  missing_compare(agg_phys, explanatory) %>% 
  knitr::kable(row.names=FALSE, align = c("l", "l", "r", "r", "r")) 

# Values do not significantly differ between missing and not missing 
# Nothing concerning here


######################### Complete Case Data frame (-ADH) #####################
# Create data frame with complete observations - with exception of ADH
exclude_col <- "ADH"

# Get logical vector of rows that are complete except for ADH
complete_rows <- complete.cases(merged_data[ ,
        setdiff(names(merged_data), exclude_col)])

# Create data frame
complete_analysis_data <- merged_data[complete_rows, ]
## We will proceed with this data frame
# From 1177 obs to 1085 obs
# From 689 obs. to 660 subjects

# Assess missingness once again
# Cumulative sum of missingness for each variable
naniar::miss_var_summary(complete_analysis_data)

# How many variables 0 - 5 missing values
naniar::miss_case_table(complete_analysis_data)
# All just missing ADH values

# Visualize the missingness
vis_dat(complete_analysis_data)
# Only see missingness for ADH variable

# Another way to visualize missingness
missing_plot(complete_analysis_data)

########################## Subjects with Year1 & Year 2 ########################

####### Filter to complete observations - OBS. for BOTH Year 0 and year 2 ######
# Filter to these subjects
both_yrs_data <- complete_analysis_data %>%
  group_by(newid) %>%
  filter(n_distinct(years) == 2) %>%
  ungroup()
# From 1085 obs to 850 obs
# From 689 subjects to 425 obs

# Assess missingness once again
# Cumulative sum of missingness for each variable
naniar::miss_var_summary(both_yrs_data)

# How many variables 0 - 5 missing values
naniar::miss_case_table(both_yrs_data)
# All just missing ADH values

# Visualize the missingness
vis_dat(both_yrs_data)
# Only see missingness for ADH variable

# Another way to visualize missingness
missing_plot(both_yrs_data)


################################# Baseline DF ##################################

# Data frame of only baseline information
# Subset data to just baseline aka year 1 and one observat
baseline <- both_yrs_data[both_yrs_data$years == 0, ]
# There are 425 subjects at baseline

# How many were on hard drugs at baseline?
baseline %>% 
  summarise(prop_hd = mean(hard_drugs, na.rm = TRUE))
# 0.0824 not on hard drugs

table(baseline$hard_drugs)
# 0   1 
# 390  35
# So around 9% were on hard drugs at baseline

################################################################################
###                        Models for Primary RQs                           ###
################################################################################

### For our outcomes, we will not just be using each variable as they are
# VLOAD gave us a weird bimodal distribution, 
# LEU3N, gave us a heavily right skewed distribution,
# And the QOL scores had heavily left skewed distributions.

### Instead, our models will look like for the primary questions of interest: 
# VLOAD_yr2 ~ VLOAD_yr1 + hard_drug + covariates
# This answers the question of the VLOAD diff at year 2 for nonhard and hard drug
# LEU3N_yr2 ~ LEU3N_yr1 + hard_drug + covariates
# AGG_MENT_yr2 ~ AGG_MENT_yr1 + hard_drug + covariates
# AGG_PHYS_yr2 ~ AGG_PHYS_yr1 + hard_drug + covariates

########## Create these variables in the data frames

## Creating a wide data frame - we need one row for each subject
# because we are NOT doing a longitudinal analysis D:

# Create vector of covariates for data manipulation



wide_data <- both_yrs_data %>%
  group_by(newid) %>%
  summarise(
    AGG_MENT_yr0 = AGG_MENT[years == 0],
    AGG_MENT_yr2 = AGG_MENT[years == 2],
    VLOAD_yr0    = VLOAD[years == 0],
    VLOAD_yr2    = VLOAD[years == 2],
    LEU3N_yr0    = LEU3N[years == 0],
    LEU3N_yr2    = LEU3N[years == 2],
    hard_drugs      = hard_drugs[years == 0],
    across(all_of(covariates), ~ .x[years == 0])
  )

# Assessing changes
both_yrs_data %>% 
  group_by(newid) %>% 
  reframe(newid, years, AGG_MENT, AGG_MENT_yr2, AGG_MENT_yr0)


################################################################################
###                             Distn. of Outcomes                           ###
################################################################################

################################## VLOAD #######################################

# What is the distribution of VLOAD
summary(both_yrs_data$VLOAD_yr2)

# Histogram of the distribution of VLOAD
ggplot(both_yrs_data, aes(x = VLOAD_yr2)) + 
  geom_histogram(bin2 = 40,
                 fill = "blue",
                 col = "black") +
  theme_lucid() + 
  theme(legend.position = "none")

# Log transform VLOAD
ggplot(both_yrs_data, aes(x = log10(VLOAD_yr2))) + 
  geom_histogram(bins = 30,
                 fill = "blue",
                 col = "black") +
  geom_density() + 
  theme_lucid() + 
  theme(legend.position = "none")

## Formal assessment of normality
shapiro.test(both_yrs_data$VLOAD_yr2) 
# We reject the H0 and data is not normally distributed...
# However, Camille said log10 distribution is sufficient


################################################################################

################################## LEU3N #######################################

# What is the distribution of LEU3n
summary(both_yrs_data$LEU3N_yr2)

# Histogram of the distribution of LEU3N
ggplot(both_yrs_data, aes(x = LEU3N_yr2)) + 
  geom_histogram(bins = 30,
                 fill = "seagreen3",
                 col = "black") +
  theme_lucid() + 
  theme(legend.position = "none")

# Log transform LEU3N - do not log transform
ggplot(both_yrs_data, aes(x = log10(LEU3N_yr2))) + 
  geom_histogram(bins = 30,
                 fill = "seagreen3",
                 col = "black") +
  theme_lucid() + 
  theme(legend.position = "none")

# Sqrt transform VLOAD 
ggplot(both_yrs_data, aes(x = sqrt(LEU3N_yr2))) + 
  geom_histogram(bins = 30,
                 fill = "seagreen3",
                 col = "black") +
  theme_lucid() + 
  theme(legend.position = "none")
# However, due to interpretations with square root we will not proceed with this

## Formal assessment of normality
shapiro.test(both_yrs_data$LEU3N_yr2) 
# We reject the H0 and data is not normally distributed but we proceed anyways


################################################################################

################################## AGG MENT ####################################

summary(both_yrs_data$AGG_MENT_yr2)

# Histogram of distribution
ggplot(both_yrs_data, aes(x = (AGG_MENT_yr2))) + 
  geom_histogram(fill = "purple", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none") 

# Log transform
ggplot(both_yrs_data, aes(x = sqrt(AGG_MENT_yr2))) + 
  geom_histogram(fill = "purple", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none") 

# Histogram of this
ggplot(both_yrs_data, aes(x = (prop_AGG_MENT))) + 
  geom_histogram(fill = "purple", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none") 

################################## AGG PHYS ####################################

# What is the distribution of Aggregate phsycial quality of life score
summary(both_yrs_data$AGG_PHYS)

# Histogram of distribution
ggplot(both_yrs_data, aes(x = (AGG_PHYS))) + 
  geom_histogram(fill = "orange", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none")

# We also observe a left skew here
### We will also consider beta regression with this
# Convert percentages to decimal values
both_yrs_data$prop_AGG_PHYS <- both_yrs_data$AGG_PHYS / 100

# Histogram of this
ggplot(both_yrs_data, aes(x = (prop_AGG_PHYS))) + 
  geom_histogram(fill = "orange", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none")

summary(both_yrs_data$prop_AGG_PHYS)


################################################################################

################################# Covariates ###################################

################################## Hard Drugs  #################################

table(both_yrs_data$hard_drugs)
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

#################################### BMI ######################################

summary(both_yrs_data$BMI)
# There are some values that have 999 and even -1

# Plot the distribution
ggplot(both_yrs_data, aes(x = BMI)) + 
  geom_histogram(bins = 10) + 
  theme_lucid()


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
table(factor(both_yrs_data$SMOKE))
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

table(both_yrs_data$ADH)

# Barplot
ggplot(both_yrs_data, aes(x = ADH)) + 
  geom_bar() + 
  theme_lucid()



################################### Table 1 ####################################

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

