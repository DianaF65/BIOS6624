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

# Read in dataset
frame_data <- read_csv(here("Project3", "Data", "frmgham2.csv"))

################################################################################
###                        Data Exploration and Cleaning                     ###
################################################################################

#### First glimpse of the data
glimpse(frame_data)

# Summary of vars in data
summary(frame_data)

################################ Missingness ###################################
#### Really nice tutorial: https://cran.r-project.org/web/packages/finalfit/vignettes/missing.html

# Missinginess with all information
# Primary covariates of interest are: Age, Diabetes, Blood Pressure
# Additional covariates: CHD, BP Meds, Smoke status, CHOL, BMI, SYS BP

# Cumulative sum of missingness for each variable
m1 <- frame_data %>% 
  miss_var_summary(order = TRUE)
# variable n_miss pct_miss
# 1 LDLC       8601  74.0   
# 2 HDLC       8600  74.0   
# 3 GLUCOSE    1440  12.4   
# 4 BPMEDS      593   5.10  
# 5 TOTCHOL     409   3.52  
# 6 educ        295   2.54  
# 7 CIGPDAY      79   0.679 
# 8 BMI          52   0.447 
# 9 HEARTRTE      6   0.0516
# 10 RANDID        0   0  
## LDLC, HDLC, and GLUCOSE have highest missingness

# How many variables 0 - 5 missing values
m2 <- frame_data %>% 
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
missing_plot(frame_data)

## Overall, no missingness for survival variables: STROKE and TIMESTRK
## No Missingness for primrary covariates
## Additional covariates: BPMEDS, TOTCHOL, BMI


########################## Assessing patterns of missingness ###################
### Really good package for survival analysis: 
# https://sentinelinitiative.org/sites/default/files/documents/smdi_r_pharma2023_vF.pdf

# Explore patterns of missingness to determine if data is MNAR, MAR, MCAR

###################################### STROKE ###################################

# Explore patterns of missingness for stroke - yes/no
explanatory <- c("AGE", "SYSBP", "DIABETES")
# STROKE
stroke <- c("STROKE")

frame_data %>% 
  missing_pattern(stroke, explanatory)

# There does not seem to be a relationship with
# missingness between any pair of predictors or predictor and STROKE

############################### TIME TO STROKE #################################

# Explore patterns of missingness for time to stroke - yes/no
explanatory <- c("AGE", "SYSBP", "DIABETES")
# TIME TO STROKE
timestroke <- c("TIMESTRK")

frame_data %>% 
  missing_pattern(timestroke, explanatory)

# No relationship with missingness and time to stroke

############################ Time Changing Variables ###########################
### Description of time changing variables (Diabetes and BP specfically)
# Should these be considered for future time varying analyses?



############################### Survival Analysis DF ###########################


### Data exploration
# Subjects that had a stroke at time 0
stroke_before <- frame_data %>% 
  group_by(RANDID) %>% 
  filter(STROKE == 1 & TIMESTRK == 0) %>% 
  # filter(STROKE == 1 & TIME == 0) %>% 
  reframe(RANDID, STROKE, TIMESTRK, TIME,
          # Create new variables
          STROKE10 = 1,
         TIMESTRK10 = 0 ) 
# There were 32 subjects that had strokes before 

# Create new STROKE and STRKE10 years relative to 10 years
# For subjects that had strokes before, we will adjust the time and stroke == 0
# Indicator for whether data is within 10 years
frame_data <- frame_data %>% 
  group_by(RANDID) %>% 
  mutate(strokewithin10 = ifelse(TIMESTRK > 0 & TIMESTRK <= 3600, 1, 0)) %>% 
  relocate(c(strokewithin10), .after = TIMESTRK) 

# Create baby data frame to mess around and figure out
baby <- frame_data %>% 
  reframe(RANDID, STROKE, PREVSTRK, PERIOD, TIME,
          TIMESTRK, strokewithin10)

# Create a new stroke and time to stroke varible in regard to 10 years
# We also want to adjust for those that already had a stroke before start of study
baby2 <- baby %>% 
  mutate(STROKE10 = ifelse(PREVSTRK == 0 & TIME == 0,
                           ifelse(strokewithin10 == 1,
                                  1,0),
                           0),
        STROKE10TIME = ifelse(PREVSTRK == 0 & TIME == 0,
                                ifelse(STROKE10 == 1,
                                       TIMESTRK, 0),
                              0))
# There are currently 11627 obs for 4434

# Looking into those subjects that did have stroke within 10 years


# Create one observation for each subject
pre_surv_df <- frame_data %>% 
  group_by(RANDID) %>% 
  filter(PERIOD == 1)
  
  







