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

# Read in dataset
frame_data <- read_csv(here("Project3", "Data", "frmgham2.csv"))

################################################################################
###                        Data Exploration and Cleaning                     ###
################################################################################

#### First glimpse of the data
glimpse(frame_data)

# Summary of vars in data
summary(frame_data)

# Interest is time to stroke within first 10 years of study


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
baby3 <- baby2 %>% 
  group_by(RANDI) %>% 
  # Pull first obs for each subject 
  slice_head(n = 1)
  
  
  








