# Project 0 Code Analysis

# To make commits, in terminal type: 
#   Git add .
# Git commit -m “message”
# Git push
# Check that changes are online


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

# Read in the dataset
proj0_dat <- read.csv("RawData/Project0_Clean_v2.csv",
                      header = TRUE,
                      check.names = FALSE)

### Glimpse at the data
glimpse(proj0_dat)

### Missingness
# How many NAs
# Cumulative sum of missingness for each variable
naniar::miss_var_cumsum(proj0_dat)

# A summary of missingness for each variable
naniar::miss_var_summary(proj0_dat)
# Looks MEMs: Sample interval Decimal Time is missing completely for some subj.

# Which variables have missing values
naniar::miss_var_which(proj0_dat)

# How many variables 0 - 5 missing values
naniar::miss_case_table(proj0_dat)

# Visualize the missingness
naniar::vis_miss(proj0_dat)

### Summarize the subject's and MEMs recording of sampling times
proj0_dat %>% 
  group_by(SubjectID, "Collection Sample") %>% 
  summarise(`MEMs: Clock Time`, `Booket: Clock Time`)
# Look into how to visualize time data


