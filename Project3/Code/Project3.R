################################################################################
###                                  Project 3                               ###
################################################################################

# Load libraries
library(readr)
library(here)
library(dplyr)
library(ggplot2)
library(naniar)
library(see)
library(pillar)
library(visdat)
library(finalfit)
library(table1)

# Read in dataset
frame_data <- read_csv(here("Project3", "Data", "frmgham2.csv"))

################################################################################
###                        Data Exploration and Cleaning                     ###
################################################################################

#### First glimpse of the data
glimpse(frame_data)

# Summary of vars in data
summary(frame_data)

# Which types of variables
str(frame_data)

# Make sure vars are what they should be
## Factors
# RANDID
frame_data$RANDID <- factor(frame_data$RANDID)
# Sex
frame_data$SEX <- factor(frame_data$SEX)
# Cursmoke
frame_data$CURSMOKE <- factor(frame_data$CURSMOKE)
# Diabetes
frame_data$DIABETES <- factor(frame_data$DIABETES)
# BP meds
frame_data$BPMEDS <- factor(frame_data$BPMEDS)
# eEducation
frame_data$educ <- factor(frame_data$educ)
# PREVCHD
frame_data$PREVCHD <- factor(frame_data$PREVCHD)
# PREVAP 
frame_data$PREVAP <- factor(frame_data$PREVAP)
# PREVMI
frame_data$PREVMI <- factor(frame_data$PREVMI)
# PREVSTRK 
frame_data$PREVMI <- factor(frame_data$PREVMI)
# PREVHYP
frame_data$PREVHYP <- factor(frame_data$PREVHYP)
# Death
frame_data$DEATH <- factor(frame_data$DEATH)
# Angine
frame_data$ANGINA <- factor(frame_data$ANGINA)
# HOSPI
frame_data$HOSPMI <- factor(frame_data$HOSPMI)
# mi fchd
frame_data$MI_FCHD <- factor(frame_data$MI_FCHD)
# PREVCHD
frame_data$PREVCHD <- factor(frame_data$PREVCHD)
# stroke
frame_data$STROKE <- factor(frame_data$STROKE)
# CVD
frame_data$CVD <- factor(frame_data$CVD)
# Hyptertension
frame_data$HYPERTEN <- factor(frame_data$HYPERTEN)

############################### Survival Analysis DF ###########################

## Data exploring
# Subjects that had a stroke at time 0
stroke_before <- frame_data %>% 
  group_by(RANDID) %>% 
  filter(STROKE == 1 & TIMESTRK == 0) %>% 
  # filter(STROKE == 1 & TIME == 0) %>% 
  reframe(RANDID, STROKE, TIMESTRK, TIME,
          # Create new variables
          STROKE10 = 1,
         TIMESTRK10 = 0 ) 
# Using STROKE == 1 & TIMESTROKE == 0 because TIME = 0 removes too much
# There were 32 subjects that had strokes before 

# Create new STROKE and STRKE10 years relative to 10 years
# For subjects that had strokes before, we will adjust the time and stroke == 0
# Indicator for whether data is within 10 years
pre_surv_df <- frame_data %>% 
  # group_by(RANDID) %>% 
  mutate(# Stroke at TIME = 0 indicator
         strokeBL = ifelse((TIMESTRK == 0 & STROKE == 1),
                           1, 0),
         # Death within 10 years
         deathwithin10 = ifelse(TIMEDTH <= 3650, 1, 0),
         # Stroke within 10 years indicator
         strokewithin10 = ifelse((strokeBL == 0 |
                                   deathwithin10 == 0),
                            ifelse(STROKE == 1 & between(TIMESTRK, 1, 3650),
                                   1, 0),
                            0))%>% 
  relocate(strokewithin10, .after = TIMESTRK) %>% 
  relocate(strokeBL, .after = STROKE) 

# Double checking stroke within 10 years
check1 <- pre_surv_df %>% 
  filter(PERIOD == 1 & strokewithin10 == 1)

# Looking at vars of interest
check1 <- subset(check1, select = c("RANDID", "TIME", "PERIOD",
                                    "TIMESTRK", "STROKE",
                                    "strokeBL", "strokewithin10",
                                    "DEATH", "TIMEDTH",
                                    "deathwithin10"))
# Seems to be good

# Investigate subjects who died
## Did this after and went back and added deathwithin 10 var
mort <- subset(pre_surv_df, 
               select = c("RANDID", "TIME", "PERIOD",
                          "TIMESTRK", "STROKE",
                          "strokeBL", "strokewithin10", 
                          "DEATH", "TIMEDTH"))



# Create baby data frame to mess around and figure out stroke variables
baby <- pre_surv_df %>% 
  reframe(RANDID, STROKE, PREVSTRK, PERIOD, TIME,
          TIMESTRK, strokewithin10)

# Create a new stroke and time to stroke varible in regard to 10 years
# We also want to adjust for those that already had a stroke before start of study
pre_surv_df <- pre_surv_df %>% 
  # For subjects who did not have stroke at time 0
  mutate(STROKE10 = factor(ifelse(strokewithin10 == 1,
                           1, 0)),
         # The time to stroke for strokes within 10 years
        STROKE10TIME = ifelse(strokewithin10 == 1,
                              TIMESTRK, 0))
# There are currently 11627 obs for 4434 subjects
# Looking at vars of interest
b <- subset(pre_surv_df, select = c("RANDID", "TIME", "PERIOD",
                                    "TIMESTRK", "STROKE",
                                    "strokeBL", "STROKE10",
                                    "STROKE10TIME"))


################################ Missingness ###################################
#### Really nice tutorial: https://cran.r-project.org/web/packages/finalfit/vignettes/missing.html

# Missinginess with all information
# Primary covariates of interest are: Age, Diabetes, Blood Pressure
# Additional covariates: CHD, BP Meds, Smoke status, CHOL, BMI, SYS BP

# Cumulative sum of missingness for each variable
m1 <- pre_surv_df[,-1] %>% 
  naniar::miss_var_summary(order = TRUE)
# variable n_miss pct_miss
# 1 LDLC       8601  74.0   
# 2 HDLC       8600  74.0   
# 3 GLUCOSE    1440  12.4   
# 4 BPMEDS      593   5.10  #####
# 5 TOTCHOL     409   3.52  #####
# 6 educ        295   2.54  
# 7 CIGPDAY      79   0.679 
# 8 BMI          52   0.447 #####
# 9 HEARTRTE      6   0.0516
# 10 SEX           0   0 
## LDLC, HDLC, and GLUCOSE have highest missingness

# How many variables 0 - 5 missing values
m2 <- frame_data[, -1] %>% 
  naniar::miss_case_table()
# n_miss_in_case n_cases pct_cases
# 1              0    2236   19.2   
# 2              1     614    5.28  
# 3              2    7240   62.3   
# 4              3    1196   10.3   
# 5              4     314    2.70  
# 6              5      25    0.215 
# 7              6       2    0.0172
## Missingness highest for 2 variables

# Visualize the missingness - 40 cols
# First 20
vis_dat(frame_data[1:20])
# Next 20
vis_dat(frame_data[21:39]) # HDLC very missing

# Another way to visualize missingness
# Break up into chunks - 40 cols
# First 20
missing_plot(frame_data[1:20])
# Next 20
missing_plot(frame_data[21:39])
# Whole data set
missing_plot(frame_data[2:39])

## Overall, no missingness for survival variables: STROKE and TIMESTRK
## No Missingness for primrary covariates
## Additional covariates: BPMEDS, TOTCHOL, BMI


########################## Assessing patterns of missingness ###################
### Really good package for survival analysis: 
# https://sentinelinitiative.org/sites/default/files/documents/smdi_r_pharma2023_vF.pdf

# Explore patterns of missingness to determine if data is MNAR, MAR, MCAR

##################################### STROKE 10 ################################

# Explore patterns of missingness for stroke - yes/no
explanatory <- c(# Primary vars
  "SEX", "AGE", "SYSBP", "DIABETES", "STROKE10TIME", 
  # Additional covariates
  "PREVCHD", "TIMECHD", "BPMEDS", "CURSMOKE", "TOTCHOL", "BMI")
# STROKE
stroke <- c("STROKE10")

# Plot for patterns 
pre_surv_df %>% 
  missing_pattern(stroke, explanatory)

# There does not seem to be a relationship with
# missingness between any pair of predictors or predictor and STROKE

### Compare Not Missing and Missing values 
pre_surv_df[, -1] %>% 
  ungroup() %>% 
  missing_compare(stroke, explanatory) %>% 
  knitr::kable(row.names=FALSE) 
# Cannot run this for some reason

# Errors with above
pre_surv_df[, -1] %>% 
  dplyr::summarise(
    pval_diabetews = chisq.test(DIABETES, STROKE10, simulate.p.value = TRUE)$p.value,
  )

######## Missinginess with SMDI package
library(smdi)
library(sm)

pre_surv_df %>% 
  gg_miss_upset()
# We see that there are 6 observations that are missing all covariates
# Does not appear to be any subjects that are missing both TOTCHOL and BPMEDs
# Good amount of TOTCHOL missing with HDLC and LDLC though = 270
# BPMEDS missing with HDLC and LDLC = 136

# Create a missingness df without the RANDIDs
miss_df <- pre_surv_df[, -1]

# smdi function
# This function outputs means for all observed and missing obs for all vars
ah <- smdi_asmd(miss_df, includeNA = TRUE)

# Look at TOTCHOL
ah$TOTCHOL$asmd_table1

############################### TIME TO STROKE #################################
# Explore patterns of missingness for time to stroke - yes/no
# Explanatory vars
explanatory2 <- c(# Primary vars
  "SEX", "AGE", "SYSBP", "DIABETES", 
  # Additional covariates
  "PREVCHD", "TIMECHD", "BPMEDS", "CURSMOKE", "TOTCHOL", "BMI")

# TIME TO STROKE
timestroke <- c("STROKE10TIME")

pre_surv_df %>% 
  missing_pattern(timestroke, explanatory2)

# Mean with and without missing vals
pre_surv_df %>% 
  ungroup() %>% 
  missing_compare(timestroke, explanatory2) %>% 
  knitr::kable(row.names=FALSE) 

# No relationship with missingness and time to stroke

############################### Time Varying Covariates ########################

# Time varying Covariates of interest: Age, Diabetes, SYS BP
time_covs <- c("RANDID", "AGE", "SEX",
               "DIABETES", "SYSBP", "TIME", "PERIOD",
               "STROKE10", "STROKE10TIME")

# Subset df
time_df <- subset(pre_surv_df, 
                  select = time_covs)

# Create a Table 1 with these variables
table1(~ factor(STROKE10) + STROKE10TIME + 
         AGE + SYSBP + factor(DIABETES) + factor(SEX)
       | factor(PERIOD), 
       data = time_df)

# Pehaps make some plots eventually to visualize the changes in these vars 
# over time


##################################### PERIOD 1 ################################

# Create one observation for each subject
# Filter to PERIOD 1
per1_surv_df <- pre_surv_df %>% 
  group_by(RANDID) %>% 
  filter(PERIOD == 1) %>% 
  ungroup(RANDID)
# There are 4434 obs for 4434 subjects

# Filter to columns of interest for analysis
cols <- c("RANDID", 
          # Primary vars
          "SEX", "AGE", "SYSBP", "DIABETES",
          "STROKE10", "STROKE10TIME", "deathwithin10",
          # Additional covariates
          "PREVCHD", "TIMECHD", "BPMEDS", "CURSMOKE", "TOTCHOL", "BMI")

# Filter
per1_surv_df <- subset(per1_surv_df, 
                       select = cols)

# Check how many subjects had stroke within 10 years
check3 <- per1_surv_df %>% 
  filter(STROKE10 == 1)
# Verify female and male count
table(check3$SEX)
# 1  2 
# 49 62 

#### Data frame for Analysis
# Look at missingness one more time

# Missinginess with all information
# Primary covariates of interest are: Age, Diabetes, Blood Pressure
# Additional covariates: CHD, BP Meds, Smoke status, CHOL, BMI, SYS BP

# Cumulative sum of missingness for each variable
m3 <- per1_surv_df[,-1] %>% 
  miss_var_summary(order = TRUE)
# variable     n_miss pct_miss
# 1 BPMEDS           61    1.38 
# 2 TOTCHOL          52    1.17 
# 3 BMI              19    0.429
# 4 SEX               0    0    
# 5 AGE               0    0  

# How many variables 0 - 5 missing values
m4 <- per1_surv_df[, -1] %>% 
  naniar::miss_case_table()
# n_miss_in_case n_cases pct_cases
# 1              0    4304   97.1   
# 2              1     128    2.89  
# 3              2       2    0.0451

# Another way to visualize missingness
# Whole data set
missing_plot(per1_surv_df[, -1])



######### Save DF
write.csv(per1_surv_df, here("Project3",
                             "Data",
                             "frmg_survival_df.csv"),
          row.names = FALSE)

################################ Sex Specific DFs ##############################

# Males
males_df <- per1_surv_df %>% 
  filter(SEX == 1)
# There are 1944 obs for 1944 subjects

# Females
females_df <- per1_surv_df %>% 
  filter(SEX == 2)
# There are 2490 obs for 2490 subjects


############################## Data Distributions ##############################

############ Overall DF
# Sex
table(per1_surv_df$SEX)
# 1    2 
# 1944 2490

###### Outcomes 

# Summary stats of stroke and time for 10 years
table(per1_surv_df$STROKE10)

# Plot of stroke within 10 years
ggplot(per1_surv_df, aes(x = factor(STROKE10),
                         fill = factor(STROKE10))) + 
  geom_bar() + 
  theme_lucid()

#### Males
table(males_df$STROKE10)
# 0    1 
# 1895   49 

# Plot of stroke within 10 years
ggplot(males_df, aes(x = STROKE10,
                         fill = STROKE10)) + 
  geom_bar() + 
  theme_lucid()

#### Females
table(females_df$STROKE10)
# 0    1 
# 2428   62

# Plot of stroke within 10 years
ggplot(females_df, aes(x = STROKE10,
                         fill = STROKE10)) + 
  geom_bar() + 
  theme_lucid()


# Distribution times of strokes within 10 years

### Males
ggplot(males_df, aes(x = STROKE10TIME)) + 
  geom_histogram(color = "magenta4", fill = "magenta4") + 
  theme_lucid() 

# Summary
summary(males_df$STROKE10TIME)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     0.0     0.0   261.8     0.0  3595.0 

table(males_df$STROKE10TIME)
# 0   26   45   87  133  266  267  287  294  305  346  350  378  424  430  442 
# 1895    1    1    1    1    1    1    1    1    1    1    1    1    1    2    1 
# Pretty heavy zero inflated

### Females
ggplot(females_df, aes(x = STROKE10TIME)) + 
  geom_histogram(color = "magenta4", fill = "magenta4") + 
  theme_lucid() 

# Summary
summary(females_df$STROKE10TIME)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     0.0     0.0   177.9     0.0  3619.0 

table(females_df$STROKE10TIME)
# 0   22   47   58   73  101  110  126  145  146  150  168  178  182  184  234 
# 2448    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
# Pretty heavy zero inflated

##################################### Covariates ###############################

######################################## Age ###################################
## Males
summary(males_df$AGE)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 33.00   42.00   49.00   49.79   57.00   69.00

# Plot
ggplot(males_df, aes(x = AGE)) + 
  geom_histogram(fill = "red", color = "black") +
  theme_lucid()

## Females
summary(females_df$AGE)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 32.00   43.00   49.00   50.03   57.00   70.00 

# Plot
ggplot(females_df, aes(x = AGE)) + 
  geom_histogram(fill = "magenta2", color = "black") +
  theme_lucid()
 
##################################### Diabetes #################################
## Males
table(males_df$DIABETES)
# 0    1 
# 1885   59 

# Plot
ggplot(males_df, aes(x = DIABETES, fill = factor(DIABETES))) + 
  geom_bar() +
  theme_lucid()

## Females
table(females_df$DIABETES)
# 0    1 
# 2428   62

# Plot
ggplot(females_df, aes(x = DIABETES, fill = factor(DIABETES))) + 
  geom_bar() +
  theme_lucid()

## Overal diabetes is pretty unbalanced

###################################### SYSBP ###################################
## Males
summary(males_df$SYSBP)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 83.5   118.0   129.0   131.7   141.5   235.0 

# Plot
ggplot(males_df, aes(x = SYSBP)) + 
  geom_histogram(color = "black", fill = "orange") +
  theme_lucid()

# Log transformed
ggplot(males_df, aes(x = log(SYSBP))) + 
  geom_histogram(color = "black", fill = "orange") +
  theme_lucid()
##### Will log transform

## Females
summary(females_df$SYSBP)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 83.5   116.0   128.5   133.8   146.5   295.0 

# Plot
ggplot(females_df, aes(x = SYSBP)) + 
  geom_histogram(color = "black", fill = "orange") +
  theme_lucid()

# Log transformed
ggplot(females_df, aes(x = log(SYSBP))) + 
  geom_histogram(color = "black", fill = "orange") +
  theme_lucid()
##### Will log transform

#################################### PREV CHD ##################################
## Males
summary(males_df$PREVCHD)
# 0    1 
# 1234  710 

# Plot
ggplot(males_df, aes(x = PREVCHD)) + 
  geom_bar() +
  theme_lucid()

## Feales
summary(females_df$PREVCHD)
# 0    1 
# 1960  530

# Plot
ggplot(females_df, aes(x = PREVCHD)) + 
  geom_bar() +
  theme_lucid()

##################################### BP MEDS ##################################

## Males
summary(males_df$BPMEDS)
# 0    1 NA's 
# 1880   42   22 

# Plot
ggplot(males_df, aes(x = BPMEDS)) + 
  geom_bar() +
  theme_lucid()

## Females
summary(females_df$BPMEDS)
# 0    1 NA's 
# 2349  102   39 

# Plot
ggplot(females_df, aes(x = BPMEDS)) + 
  geom_bar() +
  theme_lucid()

################################## SMOKE STATUS ################################

## Males
summary(males_df$CURSMOKE)
# 0    1 
# 769 1175

# Plot
ggplot(males_df, aes(x = CURSMOKE)) + 
  geom_bar() +
  theme_lucid()

## Females
summary(females_df$CURSMOKE)
# 0    1 
# 1484 1006 

# Plot
ggplot(females_df, aes(x = CURSMOKE)) + 
  geom_bar() +
  theme_lucid()


###################################### TOTCHOL #################################

## Males
summary(males_df$TOTCHOL)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   113.0   206.0   231.0   233.6   259.0   696.0       7

# Plot
ggplot(males_df, aes(x = TOTCHOL)) + 
  geom_histogram(color = "black", fill = "orangered") +
  theme_lucid()

# Log transformed
ggplot(males_df, aes(x = log(TOTCHOL))) + 
  geom_histogram(color = "black", fill = "orangered") +
  theme_lucid()
#### Log transform 

## Females
summary(females_df$TOTCHOL)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   107.0   206.0   237.0   239.7   269.0   600.0      45 

# Plot
ggplot(females_df, aes(x = TOTCHOL)) + 
  geom_histogram(color = "black", fill = "orangered") +
  theme_lucid()

# Log transformed
ggplot(females_df, aes(x = log(TOTCHOL))) + 
  geom_histogram(color = "black", fill = "orangered") +
  theme_lucid()
### Log transform

######################################## BMI ###################################
## Males
summary(males_df$BMI)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 15.54   23.97   26.08   26.17   28.32   40.38       5 

# Plot
ggplot(males_df, aes(x = BMI)) + 
  geom_histogram(color = "black", fill = "blue3") +
  theme_lucid()

# Log transformed
ggplot(males_df, aes(x = log(BMI))) + 
  geom_histogram(color = "black", fill = "blue3") +
  theme_lucid()
# Potentially log transform

## Females
summary(females_df$BMI)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 15.96   22.54   24.83   25.59   27.82   56.80      14 

# Plot
ggplot(females_df, aes(x = BMI)) + 
  geom_histogram(color = "black", fill = "blue3") +
  theme_lucid()

# Log transformed
ggplot(females_df, aes(x = log(BMI))) + 
  geom_histogram(color = "black", fill = "blue3") +
  theme_lucid()
# Actually definiitely log transform


##################################### Table1 ###################################

# Create table 1
table1(~ factor(STROKE10) + STROKE10TIME + 
         AGE + SYSBP + factor(DIABETES) | factor(SEX), 
       data = per1_surv_df)


################################################################################
###                              Survival Analysis                           ###
################################################################################

############################ Variable Selection ################################

# Perform some variable selection prevchd, bpmeds, smokecur, total chol, and BMI









