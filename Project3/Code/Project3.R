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
  mutate(within10 = ifelse(TIME > 0 & TIME <= 3600, 1, 0),
         strokewithin10 = ifelse(TIMESTRK > 0 & TIMESTRK <= 3600, 1, 0)) %>% 
  relocate(c(within10, strokewithin10), .after = TIME) 

# Create baby data frame to mess around and figure out
baby <- frame_data %>% 
  reframe(RANDID, within10, STROKE, PREVSTRK, TIMESTRK, strokewithin10)

# Create a new stroke and time to stroke varible in regard to 10 years
baby <- baby %>% 
  mutate(STROKE10 = )
  
  








