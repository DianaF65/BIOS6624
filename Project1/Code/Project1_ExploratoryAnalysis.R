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
          "ART", "years",  "hard_drugs",
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

################################ Missingness ###################################

#### Really nice tutorial: https://cran.r-project.org/web/packages/finalfit/vignettes/missing.html

# Cumulative sum of missingness for each variable
naniar::miss_var_summary(analysis_data)
# Looks like ADH has a lot of missingness but this is expected since adh is
# measured at second visit

# How many variables 0 - 5 missing values
naniar::miss_case_table(analysis_data)

# Visualize the missingness
vis_dat(analysis_data)

# Another way to visualize missingness
missing_plot(analysis_data)

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


### Compare Not Missing and Missing for 
analysis_data %>% 
  missing_compare(vload, explanatory) %>% 
  knitr::kable(row.names=FALSE, align = c("l", "l", "r", "r", "r")) 

## Looking at some categorical variables since there were some issues
# Usingn simulated pvalues from chisquare test
# SMOKE 
analysis_data %>% 
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

analysis_data %>% 
  missing_pattern(leu3n, explanatory)

# No alarming patterns here

#### Compare Not Missing and Missing
analysis_data %>% 
  missing_compare(leu3n, explanatory) %>% 
  knitr::kable(row.names=FALSE, align = c("l", "l", "r", "r", "r")) 

# No significant differences between missing and not missing values
# Missingness for LEU3N does not differ by covariates
# Not MAR or MCAR

## Looking at some categorical variables since there were some issues
# Usingn simulated pvalues from chisquare test
# SMOKE 
analysis_data %>% 
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

analysis_data %>% 
  missing_pattern(agg_ment, explanatory)

analysis_data %>% 
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

analysis_data %>% 
  missing_pattern(agg_phys, explanatory)

# No concerning patterns here either
analysis_data %>% 
  group_by(years) %>% 
  summarise(means = mean(AGG_PHYS, na.rm = TRUE))
# 1     0  50.1
# 2     2  49.4
## Not much difference in means between the two years

#### Compare Not Missing and Missing
analysis_data %>% 
  missing_compare(agg_phys, explanatory) %>% 
  knitr::kable(row.names=FALSE, align = c("l", "l", "r", "r", "r")) 

# Values do not significantly differ between missing and not missing 
# Nothing concerning here


######################### Complete Case Data frame (-ADH) #####################
# Create data frame with complete observations - with exception of ADH
exclude_col <- "ADH"

# Get logical vector of rows that are complete except for ADH
complete_rows <- complete.cases(analysis_data[ ,
        setdiff(names(analysis_data), exclude_col)])

# Create data frame
complete_analysis_data <- analysis_data[complete_rows, ]
## We will proceed with this data frame
# From 1221 obs to 1101 obs.
# From 715 subjecs to 663 subjects

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
# From 1101 obs to 876 obs
# From 663 subjects to 438 subjects

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
# There are 438 subjects at baseline

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
covariates <- c("BMI", "SMOKE", "RACE",
               "EDUCBAS",  "age",
               "hard_drugs")


wide_bth_yrs <- both_yrs_data %>%
  group_by(newid) %>%
  summarise(
    # AGG MENT
    AGG_MENT_yr0 = AGG_MENT[years == 0],
    AGG_MENT_yr2 = AGG_MENT[years == 2],
    # AGG PHYS
    AGG_PHYS_yr0 = AGG_PHYS[years == 0],
    AGG_PHYS_yr2 = AGG_PHYS[years == 2],
    # VLOAD
    VLOAD_yr0    = VLOAD[years == 0],
    VLOAD_yr2    = VLOAD[years == 2],
    # LEU3N
    LEU3N_yr0    = LEU3N[years == 0],
    LEU3N_yr2    = LEU3N[years == 2],
    # Only need 2nd year of ADH
    ADH_yr2      = ADH[years == 2],
    # Baseline for all other covariates
    across(all_of(covariates), ~ .x[years == 0])
  )

# Examining data frame
summary(wide_bth_yrs)

#### Filtering out abnormal BMI values

# Looking into BMI
summary(wide_bth_yrs$BMI)
# Need to filter out some unrealisitc BMI values.

# What is an unreasonable BMI?
# Based on this google and this article: 
# https://pmc.ncbi.nlm.nih.gov/articles/PMC2930234/#:~:text=The%20BMI%20is%20a%20surrogate,upper%20limit%20of%20'starvation'.
# We will filter out individuals that have a BMI within the ranges of (10, 60)
# Anything less than 10 indicates severe starvation and > 60 is extreme obesity

# Filter out unrealistic bmi values from data set
# We only care about baseline values
# Create data frame to year 1 then filter out these BMI values
wide_bth_yrs <- wide_bth_yrs %>% 
  # Filter BMI values
  filter(BMI > 10 & BMI < 60) 
# From 438 obs. to 425 obs.
# Same subjects

# Double check this worked
summary(wide_bth_yrs$BMI) 


################################################################################
###                             Distn. of Outcomes                           ###
################################################################################

################################## VLOAD #######################################

# What is the distribution of VLOAD
summary(wide_bth_yrs$VLOAD_yr2)

# Histogram of the distribution of VLOAD
ggplot(wide_bth_yrs, aes(x = VLOAD_yr2)) + 
  geom_histogram(bin2 = 40,
                 fill = "blue",
                 col = "black") +
  theme_lucid() + 
  theme(legend.position = "none")

# Log transform VLOAD
ggplot(wide_bth_yrs, aes(x = log10(VLOAD_yr2))) + 
  geom_histogram(bins = 30,
                 fill = "blue",
                 col = "black") +
  geom_density() + 
  theme_lucid() + 
  theme(legend.position = "none")

## Formal assessment of normality
shapiro.test(wide_bth_yrs$VLOAD_yr2) 
# We reject the H0 and data is not normally distributed...
# However, Camille said log10 distribution is sufficient


################################################################################

################################## LEU3N #######################################

# What is the distribution of LEU3n
summary(wide_bth_yrs$LEU3N_yr2)

# Histogram of the distribution of LEU3N
ggplot(wide_bth_yrs, aes(x = LEU3N_yr2)) + 
  geom_histogram(bins = 30,
                 fill = "seagreen3",
                 col = "black") +
  theme_lucid() + 
  theme(legend.position = "none")

# Log transform LEU3N - do not log transform
ggplot(wide_bth_yrs, aes(x = log10(LEU3N_yr2))) + 
  geom_histogram(bins = 30,
                 fill = "seagreen3",
                 col = "black") +
  theme_lucid() + 
  theme(legend.position = "none")

# Sqrt transform VLOAD 
ggplot(wide_bth_yrs, aes(x = sqrt(LEU3N_yr2))) + 
  geom_histogram(bins = 30,
                 fill = "seagreen3",
                 col = "black") +
  theme_lucid() + 
  theme(legend.position = "none")
# However, due to interpretations with square root we will not proceed with this

## Formal assessment of normality
shapiro.test(wide_bth_yrs$LEU3N_yr2) 
# We reject the H0 and data is not normally distributed but we proceed anyways


################################################################################

################################## AGG MENT ####################################

summary(wide_bth_yrs$AGG_MENT_yr2)

# Histogram of distribution
ggplot(wide_bth_yrs, aes(x = (AGG_MENT_yr2))) + 
  geom_histogram(fill = "purple", col = "black") + 
  theme_lucid() + 
  scale_x_continuous(breaks = seq(0, 100, by = 10),
                     limits = c(0, 100)) + 
  theme(legend.position = "none") 

# Log transform
ggplot(wide_bth_yrs, aes(x = sqrt(AGG_MENT_yr2))) + 
  geom_histogram(fill = "purple", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none") 

# Histogram of this
ggplot(wide_bth_yrs, aes(x = (prop_AGG_MENT))) + 
  geom_histogram(fill = "purple", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none") 

################################## AGG PHYS ####################################

# What is the distribution of Aggregate phsycial quality of life score
summary(wide_bth_yrs$AGG_PHYS_yr2)

# Histogram of distribution
ggplot(wide_bth_yrs, aes(x = (AGG_PHYS_yr2))) + 
  geom_histogram(bins = 40,
                 fill = "orange", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none")

# We also observe a left skew here
### We will also consider beta regression with this
# Convert percentages to decimal values
wide_bth_yrs$prop_AGG_PHYS <- wide_bth_yrs$AGG_PHYS / 100

# Histogram of this
ggplot(wide_bth_yrs, aes(x = (prop_AGG_PHYS))) + 
  geom_histogram(fill = "orange", col = "black") + 
  theme_lucid() + 
  theme(legend.position = "none")

summary(wide_bth_yrs$prop_AGG_PHYS)


################################################################################

################################# Covariates ###################################

################################## Hard Drugs  #################################

table(wide_bth_yrs$hard_drugs)
# 0   1 
# 390  35 

## Also include for people that have follow up to 2 years
ggplot(wide_bth_yrs, aes(x = factor((hard_drugs)),
                          fill = factor((hard_drugs)))) + 
  geom_bar() + 
  theme_lucid() + 
  theme(legend.position = "none")

#################################### BMI ######################################

summary(wide_bth_yrs$BMI)
# Looks good

# Plot the distribution
ggplot(wide_bth_yrs, aes(x = BMI)) + 
  geom_histogram(bins = 30,
                 fill = "brown",
                 col = "black") + 
  theme_lucid()


####################################### Age  ###################################

# Summary of age at baseline
summary(both_yrs_data$age)

# Histogram of distribution
ggplot(both_yrs_data, aes(x = age)) + 
  geom_histogram(col = "black", fill = "yellow") + 
  theme_lucid()


###################################### Race  ###################################

### The initial race categories are: 
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
table(factor(wide_bth_yrs$RACE))


# Change to White vs. non-white
wide_bth_yrs <- wide_bth_yrs %>% 
  mutate(collapse_RACE = factor(ifelse(
    RACE == 1,
    "White_NH", "Other"
  ))
  )

## Make sure it is a factor

# Barplot
ggplot(wide_bth_yrs, aes(x = collapse_RACE,
                     fill = collapse_RACE)) + 
  geom_bar() + 
  theme_lucid() + 
  theme(legend.position = "none")


##################################### SMOKE  ###################################
table(factor(wide_bth_yrs$SMOKE))
# 1   2   3 
# 237 287 326

# Visualize the distribution
ggplot(wide_bth_yrs, aes(x = factor(SMOKE),
                         fill = factor(SMOKE))) + 
  geom_bar() + 
  theme_lucid(legend.position = "none") 

##################################### EDUCBAS  #################################

# Initial categories for education are: 
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

# We will collapse this to Complete college or higher and did not complete college
wide_bth_yrs <- wide_bth_yrs %>% 
  mutate(collapsed_EDUCBAS = ifelse(
    EDUCBAS %in% c(5, 6, 7),
    "Comp_Coll_Higher",
    "Less_Than_College"
  ))

table(factor(wide_bth_yrs$collapsed_EDUCBAS))
# Comp_Coll_Higher Less_Than_College 
# 182               243  

# Barplot
ggplot(wide_bth_yrs, aes(x = factor(collapsed_EDUCBAS),
                     fill = factor(collapsed_EDUCBAS))) + 
  geom_bar() + 
  theme_lucid() + 
  theme(legend.position = "none")
# Pretty unbalanced


################################## Adherence ###################################

table(wide_bth_yrs$ADH_yr2)

# Barplot
ggplot(wide_bth_yrs, aes(x = factor(ADH_yr2),
                         fill = factor(ADH_yr2))) + 
  geom_bar() + 
  theme_lucid(legend.position = "none")

################################################################################
###                               Final Data frame                           ###
################################################################################

# Remove the original education and race variables
wide_bth_yrs <- subset(wide_bth_yrs, 
                       select = -c(RACE, EDUCBAS))

# Save the final data frame for analysis
write.csv(wide_bth_yrs,
          "../Data/final_hiv_data.csv")


################################################################################
###                                   Table 1                                ###
################################################################################

# Create a Table 1

### Convert to correct data types
## Factors
# new id
wide_bth_yrs$newid <- as.factor(wide_bth_yrs$newid)
wide_bth_yrs$SMOKE <- as.factor(wide_bth_yrs$SMOKE)
wide_bth_yrs$collapse_RACE <- as.factor(wide_bth_yrs$collapse_RACE)
wide_bth_yrs$collapsed_EDUCBAS <- as.factor(wide_bth_yrs$collapsed_EDUCBAS)
wide_bth_yrs$ADH_yr2 <- as.factor(wide_bth_yrs$ADH_yr2)
wide_bth_yrs$hard_drugs <- as.factor(wide_bth_yrs$hard_drugs)


### Attempt 2
library(table1)

# Create table 1
table1(~ VLOAD_yr2 + LEU3N_yr2 + 
         AGG_MENT_yr2 + AGG_PHYS_yr2 +
         hard_drugs + ADH_yr2 + SMOKE + collapse_RACE + collapsed_EDUCBAS
       + BMI + age, 
       data = wide_bth_yrs)

# Separate by hard drug exposure at baseline
table1(~ VLOAD_yr2 + LEU3N_yr2 + 
         AGG_MENT_yr2 + AGG_PHYS_yr2 +
         ADH_yr2 + SMOKE + collapse_RACE + collapsed_EDUCBAS
       + BMI + age | hard_drugs, 
       data = wide_bth_yrs)

################################################################################
###                                   STAN                                  ###
################################################################################


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

