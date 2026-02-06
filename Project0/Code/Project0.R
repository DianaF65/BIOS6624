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


### Read in the data
# The data is located in folder called RawData
proj0_dat <- read.csv("../RawData/Project0_Clean_v2.csv",
  header = TRUE,
  check.names = FALSE)

### Modify data for Person 3029
proj0_dat <- proj0_dat %>%
  mutate(DAYNUMB = case_when(
    SubjectID == 3029 & `Collection Date` == "10/7/2018" ~ 1,
    SubjectID == 3029 & `Collection Date` == "10/10/2018" ~ 2,
    SubjectID == 3029 & `Collection Date` == "10/12/2018" ~ 3,
    TRUE ~ DAYNUMB
  )) %>%
  ungroup()

############################ Data Description ##################################
####### Table1

# Fix Days to be a factor and show nicely for table1
data <- wake_data %>%
  mutate(
    DAYNUMB = factor(
      DAYNUMB,
      levels = c(1, 2, 3),
      labels = c("Day 1", "Day 2", "Day 3")),
    `Collection Sample` = factor(
      `Collection Sample`,
      levels = c(1, 2, 3, 4),
      labels = c("Wake Time", "+30 Min", "Pre Lunch", "+600 Min")
    )
  )



# for creating Table1
variables <- c("minutes_since_wake",
  "Booklet_Minutes",
  "MEM_Minutes",
  "DHEA (nmol/L)",
  "Cortisol (nmol/L)")

strata <- "Collection Sample"

# modify labels
var_label(data$Booklet_Minutes) <- "Booklet Time (Min Since Midnight)"
var_label(data$minutes_since_wake) <- "Minutes Since Waking (Midnight)"
var_label(data$MEM_Minutes) <- "MEM Time (Min Since Midnight)"


# create Table1
studypop_tab1 <- CreateTableOne(
  vars = variables,
  strata = strata,
  test = FALSE,
  data = data,
  addOverall = T
)

studypop_tab_forkbl <- print(studypop_tab1,
  varLabels = T,
  printToggle = F)

# Create table of results
kbl(studypop_tab_forkbl,
  caption =
    "Characteristics of Study Cohort",
  booktabs = T,
  align = "c") %>%
  kable_styling(latex_options = "HOLD_position", font_size = 9) %>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle")

###### Convert some of these back to hh:mm
# Create function to convert minutes to hh:mm
min_to_clock <- function(x) {
  h <- floor(x / 60)
  m <- round(x %% 60)
  sprintf("%02d:%02d", h, m)
}


### Get the means
# Mean minutes since wake
mean_minwake <- mean(data$`Sleep Diary Wake Minutes`,
  na.rm = TRUE)
# min_to_clock(mean_minwake)
# 06:52

# Mean booklet minutes
mean_bookmin <- mean(data$Booklet_Minutes,
  na.rm = TRUE) # 654.3472
# min_to_clock(mean_bookmin)
# 10:54

# Mean MEM minutes
mean_memmin <- mean(data$MEM_Minutes,
  na.rm = TRUE) # 684.2379
# min_to_clock(mean_memmin)
# 11:24

################################################################################
###                               Question                                  ###
################################################################################

############################## Question 1 Data #################################

## Make a subset of the data with the relevant cols for Q1 - Q3

# Vector of the cols of interest
cols <- c("SubjectID",
  "Collection Date", "DAYNUMB",
  "Collection Sample",
  "Booket: Clock Time", "MEMs: Clock Time",
  "Sleep Diary reported wake time",
  "Cortisol (nmol/L)", "DHEA (nmol/L)")

# Subset the original data frame
wake_data <- proj0_dat[, colnames(proj0_dat) %in% cols]

# Convert the Sleep diary Wake time to minutes
wake_data$`Sleep Diary Wake Minutes` <-
  # Convert to date objects and convert to minutes
  lubridate::hour(hm(wake_data$`Sleep Diary reported wake time`)) * 60 +
  lubridate::minute(hm(wake_data$`Sleep Diary reported wake time`))

# Move this to after Collection Sample
wake_data <- wake_data %>%
  relocate(`Sleep Diary Wake Minutes`,
    .after = `Sleep Diary reported wake time`)

# Convert the Booklet Time to minutes
wake_data$Booklet_Minutes <-
  # Convert to date objects and convert to minutes
  lubridate::hour(hm(wake_data$`Booket: Clock Time`)) * 60 +
  lubridate::minute(hm(wake_data$`Booket: Clock Time`))

# Move this to after Booklet: Clock Time
wake_data <- wake_data %>%
  relocate(Booklet_Minutes, .after = `Booket: Clock Time`)

# Convert the MEM Time to minutes
wake_data$MEM_Minutes <-
  # Convert to date objects and convert to minutes
  lubridate::hour(hm(wake_data$`MEMs: Clock Time`)) * 60 +
  lubridate::minute(hm(wake_data$`MEMs: Clock Time`))

# Move this to after MEMs: Clock Time
wake_data <- wake_data %>%
  relocate(MEM_Minutes, .after = `MEMs: Clock Time`)

# Center MEM time
wake_data$MEM_C <-
  # Convert to date objects and convert to minutes
  lubridate::hour(hm(wake_data$`MEMs: Clock Time`)) * 60 +
  lubridate::minute(hm(wake_data$`MEMs: Clock Time`))

# Move this to after MEMs: Clock Time
wake_data <- wake_data %>%
  relocate(MEM_C, .after = `MEMs: Clock Time`)

# Create a Minutes SINCE WAKE variable
wake_data <- wake_data %>%
  mutate(
    clock_min = lubridate::hour(hm(`Booket: Clock Time`)) * 60 +
      lubridate::minute(hm(`Booket: Clock Time`))) %>%
  group_by(SubjectID, DAYNUMB) %>%
  mutate(wake_min = clock_min[`Collection Sample` == 1][1],
    minutes_since_wake = clock_min - wake_min) %>%
  ungroup() %>%
  relocate(minutes_since_wake, .after = `Sleep Diary Wake Minutes`)


### Modify data for Person 3029
wake_data <- wake_data %>%
  mutate(DAYNUMB = case_when(
    SubjectID == 3029 & `Collection Date` == "10/7/2018" ~ 1,
    SubjectID == 3029 & `Collection Date` == "10/10/2018" ~ 2,
    SubjectID == 3029 & `Collection Date` == "10/12/2018" ~ 3,
    TRUE ~ DAYNUMB
  )) %>%
  ungroup()

# Create a subset of the data frame that only has Wake times
# AKA Collection Sample == 1
wakes <- wake_data[wake_data$`Collection Sample` == 1, ]


######### Model
# LMM with random intercept for subject - Day added
model_q1 <- lme(Booklet_Minutes ~ MEM_C + factor(DAYNUMB),
  random = ~ 1 | SubjectID,
  na.action = na.omit,
  method = "REML",
  data = wakes)


# Display model summary with kable
# Extract t and p-values from model for data frame
summary_df1 <- as.data.frame(summary(model_q1)$tTable)

# kable
kable(
  summary_df1,
  format = "latex",
  booktabs = TRUE,
  digits = 3,
  caption = "Model Summary"
) |>
  kable_styling(latex_options = c("striped", "hold_position"))

# intervals(model_q1)

#  Fixed effects:
#                       lower       est.      upper
# (Intercept)      -5.0891734  3.6513329 12.3918392
# MEM_C             0.9667782  0.9859424  1.0051067
# factor(DAYNUMB)2 -6.6349646 -3.4589031 -0.2828415
# factor(DAYNUMB)3 -5.0521491 -1.8810020  1.2901450

######## Assess bias
# Create a variable that is the difference between MEM and Booklet times
wake_data <- wake_data %>%
  mutate(diff_minutes = Booklet_Minutes - MEM_Minutes)

# Model to assess bias were difference is the outcome
bias_model <- lme(
  diff_minutes ~ 1 + factor(DAYNUMB),
  random = ~ 1 | SubjectID,
  data = wake_data,
  na.action = na.omit,
  method = "REML"
)

# summary(bias_model)
# intervals(bias_model)

######### Visualization

# Relationship Between booklet and MEM Wake times

ggplot(wakes, aes(x = MEM_Minutes, y = Booklet_Minutes,
  color = factor(DAYNUMB))) +
  geom_point(size = 3, alpha = 0.7) +
  theme_lucid() +
  geom_abline(intercept = 0, slope = 1, linetype = "solid") +
  labs(color = "Day") +
  labs(title = "Relationship Between Booklet and MEM Wake Up Times",
    x = "MEM Wake Up Time (Minutes)",
    y = "Booklet Wake Up Time (Minutes)")


################################################################################
###                               Question 2                                 ###
################################################################################

############################## Question 2 Data #################################



# Df with Wake time, 30 min post wake, and 600 min post wake
wake_check <- wake_data %>%
  # filter(`Collection Sample` %in% c(1, 2)) %>%
  group_by(SubjectID, DAYNUMB) %>%
  reframe(
    # Wake time
    wake_time = `Sleep Diary Wake Minutes`[`Collection Sample` == 1],
    # The 30 minute time with Booklet
    Book_30time = Booklet_Minutes[`Collection Sample` == 2],
    # The 600 miniute time with booklet
    Book_600time = Booklet_Minutes[`Collection Sample` == 4],
    # The 30 minute time with MEM
    MEM_30time = MEM_Minutes[`Collection Sample` == 2],
    # The 600 minute time with MEM
    MEM_600time = MEM_Minutes[`Collection Sample` == 4]
  )


# Define adherence for BOOK
book_wake_check <- wake_check %>%
  mutate(
    # 30 min after waking target time
    target_30time = wake_time + 30,
    # 600 min after waking target time
    target_600time = wake_time + 600,
    # Difference between 30 min later and the recorded time
    diff_min30  = Book_30time - target_30time,
    # Difference between 600 min later and the recorded time
    diff_min600 = Book_600time - target_600time,
    # For 30 min: Define good adherence - within 7.5 minutes
    good_adherent30 = abs(diff_min30) <= 7.5,
    # For 600 min: Define good adherence
    good_adherent600 = abs(diff_min600) <= 7.5,
    # For 30 min: Define adequate adherence - within 15 minutes
    adequate_ad30 = abs(diff_min30) <= 15,
    # For 600 min: Define adequate adherence
    adequate_ad600 = abs(diff_min600) <= 15
  ) %>%
  relocate(c(target_30time, diff_min30,
    good_adherent30, adequate_ad30),
  .after = Book_30time) %>%
  relocate(c(target_600time, diff_min600,
    good_adherent600, adequate_ad600),
  .after = Book_600time)

# Define adherence for MEM
MEM_wake_check <- wake_check %>%
  mutate(
    # 30 min after waking target time
    target_30time = wake_time + 30,
    # 600 min after waking target time
    target_600time = wake_time + 600,
    # Difference between 30 min later and the recorded time
    diff_min30  = MEM_30time - target_30time,
    # Difference between 600 min later and the recorded time
    diff_min600 = MEM_600time - target_600time,
    # For 30 min: Define good adherence - within 7.5 minutes
    good_adherent30 = abs(diff_min30) <= 7.5,
    # For 600 min: Define good adherence
    good_adherent600 = abs(diff_min600) <= 7.5,
    # For 30 min: Define adequate adherence - within 15 minutes
    adequate_ad30 = abs(diff_min30) <= 15,
    # For 600 min: Define adequate adherence
    adequate_ad600 = abs(diff_min600) <= 15
  ) %>%
  relocate(c(target_30time, diff_min30,
    good_adherent30, adequate_ad30),
  .after = MEM_30time) %>%
  relocate(c(target_600time, diff_min600,
    good_adherent600, adequate_ad600),
  .after = MEM_600time)


################################## Good Adherence ##############################

############################## Adherence by Subject ############################

################################### +30 Minutes ################################

######## Booklet
subject_good_30adherence <- book_wake_check %>%
  group_by(SubjectID) %>%
  summarise(
    # How many days out of 3 were they adherent?
    days_good_adherent = sum(good_adherent30),
    # What is the proportion
    prop_good_adherent = mean(good_adherent30),
    .groups = "drop"
  )
# subject_good_30adherence

# How many subjects adhered all 3 days
# sum(subject_good_30adherence$prop_good_adherent == 1, na.rm = T)
# Proportion and Percentage of subjects
# 14/31 = 0.45 or 45%

######## MEM
MEM_subject_good_30adherence <- MEM_wake_check %>%
  group_by(SubjectID) %>%
  summarise(
    # How many days out of 3 were they adherent?
    days_good_adherent = sum(good_adherent30),
    # What is the proportion
    prop_good_adherent = mean(good_adherent30),
    .groups = "drop"
  )
# MEM_subject_good_30adherence

# How many subjects adhered all 3 days
# sum(MEM_subject_good_30adherence$prop_good_adherent == 1, na.rm = T)
# Proportion and Percentage of subjects
# 6/31




################################# Adequate Adherence ###########################

############################## Adherence by Subject ############################

################################### +30 Minutes ################################

###### Book

subject_ad_30adherence <- book_wake_check %>%
  group_by(SubjectID) %>%
  summarise(
    # How many days out of 3 were they adherent?
    days_ad_adherent = sum(adequate_ad30),
    # What is the proportion
    prop_ad_adherent = mean(adequate_ad30),
    .groups = "drop"
  )

# subject_ad_30adherence

# How many subjects adhered all 3 days
# sum(subject_ad_30adherence$prop_ad_adherent == 1, na.rm = T)
# Proportion and Percentage of subjects
# 20/31 = 0.65 or 65%

###### MEM
MEM_subject_ad_30adherence <- MEM_wake_check %>%
  group_by(SubjectID) %>%
  summarise(
    # How many days out of 3 were they adherent?
    days_ad_adherent = sum(adequate_ad30),
    # What is the proportion
    prop_ad_adherent = mean(adequate_ad30),
    .groups = "drop"
  )

# MEM_subject_ad_30adherence

# How many subjects adhered all 3 days
# sum(MEM_subject_ad_30adherence$prop_ad_adherent == 1, na.rm = T)
# Proportion and Percentage of subjects
# 10/31 = 0.32 or 32%

################################## Good Adherence ##############################

################################ Adherence Overall #############################

################################### +30 Minutes ################################

######### Book
# Adherence overall
good_book <- book_wake_check %>%
  # group_by(SubjectID) %>%
  summarise(
    # What is the proportion
    prop_good_adherent = round(mean(good_adherent30, na.rm = TRUE), 4),
    # Percentage
    per_good_adherence = (prop_good_adherent * 100)
  )
# good_book
# A tibble: 1 × 2
#   prop_good_adherent per_good_adherence
#                <dbl>              <dbl>
# 1              0.782               78.2


######### MEM
# Adherence overall
good_MEM <- MEM_wake_check %>%
  # group_by(SubjectID) %>%
  summarise(
    # What is the proportion
    prop_good_adherent = round(mean(good_adherent30, na.rm = TRUE), 4),
    # Percentage
    per_good_adherence = (prop_good_adherent * 100)
  )
# good_MEM
# A tibble: 1 × 2
#   prop_good_adherent per_good_adherence
#                <dbl>              <dbl>
# 1              0.529               52.9

################################# Adequate Adherence ###########################

################################ Adherence Overall #############################

################################### +30 Minutes ################################

######## Book

# Adherence overall
adequate_book <- book_wake_check %>%
  summarise(
    # What is the proportion
    prop_adeq_adherent = round(mean(adequate_ad30, na.rm = TRUE), 4),
    # Percentage
    per_adeq_adherence = (prop_adeq_adherent * 100)
  )
# adequate_book
#   prop_adeq_adherent per_adeq_adherence
#                <dbl>              <dbl>
# 1              0.897               89.7

######## MEM
# Adherence overall
adequate_MEM <- MEM_wake_check %>%
  summarise(
    # What is the proportion
    prop_adeq_adherent = round(mean(adequate_ad30, na.rm = TRUE), 4),
    # Percentage
    per_adeq_adherence = (prop_adeq_adherent * 100)
  )
# adequate_MEM
#   prop_adeq_adherent per_adeq_adherence
#                <dbl>              <dbl>
# 1              0.714               71.4

###### Create Data frame to display with this
# Combine all the proportions into a table to display nicely
props <- data.frame(Good = c(good_book$per_good_adherence,
  good_MEM$per_good_adherence),
Adequate = c(adequate_book$per_adeq_adherence,
  adequate_MEM$per_adeq_adherence))

# Add rownames
rownames(props) <- c("Booklet", "MEM")

################################## Good Adherence ##############################

################################ Adherence Overall #############################

#################################### +600 Minutes ##############################

##### Book
# Adherence overall
good_book600 <- book_wake_check %>%
  # group_by(SubjectID) %>%
  summarise(
    # What is the proportion
    prop_good_adherent = round(mean(good_adherent600, na.rm = TRUE), 4),
    # Percentage
    per_good_adherence = (prop_good_adherent * 100)
  )
#   prop_good_adherent per_good_adherence
#                <dbl>              <dbl>
# 1              0.462               46.2

##### MEM
# Adherence overall
good_MEM600 <- MEM_wake_check %>%
  # group_by(SubjectID) %>%
  summarise(
    # What is the proportion
    prop_good_adherent = round(mean(good_adherent600, na.rm = TRUE), 4),
    # Percentage
    per_good_adherence = (prop_good_adherent * 100)
  )
#   prop_good_adherent per_good_adherence
#                <dbl>              <dbl>
# 1              0.325               32.5

############################### Adequate Adherence #############################

################################ Adherence Overall #############################

#################################### +600 Minutes ##############################

# Adherence overall
adequate_book600 <- book_wake_check %>%
  summarise(
    # What is the proportion
    prop_adeq_adherent = round(mean(adequate_ad600, na.rm = TRUE), 4),
    # Percentage
    per_adeq_adherence = (prop_adeq_adherent * 100)
  )
#   prop_adeq_adherent per_adeq_adherence
#                <dbl>              <dbl>
# 1              0.562               56.2

# Adherence overall
adequate_MEM600 <- MEM_wake_check %>%
  summarise(
    # What is the proportion
    prop_adeq_adherent = round(mean(adequate_ad600, na.rm = TRUE), 4),
    # Percentage
    per_adeq_adherence = (prop_adeq_adherent * 100)
  )

#   prop_adeq_adherent per_adeq_adherence
#                <dbl>              <dbl>
# 1              0.398               39.8

###### Create Data frame to display with this
# Combine all the proportions into a table to display nicely
props600 <- data.frame(Good = c(good_book600$per_good_adherence,
  good_MEM600$per_good_adherence),
Adequate = c(adequate_book600$per_adeq_adherence,
  adequate_MEM600$per_adeq_adherence))

# Add rownames
rownames(props600) <- c("Booklet", "MEM")


####### Combine information in a more cohesive way

props30_long <- props %>%
  as.data.frame() %>%
  rownames_to_column("Method") %>%
  mutate(Window = "+30 minutes")

props60_long <- props600 %>%
  as.data.frame() %>%
  rownames_to_column("Method") %>%
  mutate(Window = "+60 minutes")

final_props <- bind_rows(props30_long, props60_long) %>%
  relocate(Window, Method)

kable(
  final_props,
  format = "latex",
  booktabs = TRUE,
  digits = 2,
  caption = "Percentages for Good and Adequate Adherence by Collection Window"
) |>
  kable_styling(latex_options = c("striped", "hold_position")) |>
  collapse_rows(columns = 1, latex_hline = "major")



#################################### Visualize #################################

# Visualize the distribution of minutes from the sampling times

# Book +30 min
p1 <- ggplot(book_wake_check, aes(x = diff_min30)) +
  geom_histogram(binwidth = 3, col = "black", fill = "turquoise") +
  geom_vline(xintercept = c(-7.5, 7.5), col = "red",
    linetype = "dashed") +
  geom_vline(xintercept = c(-15, 15), color = "hotpink",
    linetype = "dashed") +
  theme_lucid() +
  scale_y_continuous(breaks = seq(0, 40, by = 10)) +
  labs(title = "Minutes from +30 Minute Sampling Time",
    subtitle = "Booklet",
    x = "Minutes from target (+30)",
    y = "Count")

# MEM + 30 min
p2 <- ggplot(MEM_wake_check, aes(x = diff_min30)) +
  geom_histogram(binwidth = 10, col = "black", fill = "turquoise") +
  geom_vline(xintercept = c(-7.5, 7.5), col = "red",
    linetype = "dashed") +
  geom_vline(xintercept = c(-15, 15), color = "hotpink",
    linetype = "dashed") +
  theme_lucid() +
  scale_y_continuous(breaks = seq(0, 40, by = 10),
    limits = c(0, 40)) +
  labs(title = "Minutes from +30 Minute Sampling Time",
    subtitle = "MEM",
    x = "Minutes from target (+30)",
    y = "Count")

# MEM +30 Min
p3 <- ggplot(book_wake_check, aes(x = diff_min600)) +
  geom_histogram(binwidth = 20, col = "black", fill = "turquoise") +
  geom_vline(xintercept = c(-7.5, 7.5), col = "red",
    linetype = "dashed") +
  geom_vline(xintercept = c(-15, 15), color = "hotpink",
    linetype = "dashed") +
  theme_lucid() +
  scale_x_continuous(breaks = seq(-50, 360, by = 50)) +
  scale_y_continuous(breaks = seq(0, 40, by = 10),
    limits = c(0, 40)) +
  labs(title = "Minutes from +600 Sampling Time",
    subtitle = "Booklet",
    x = "Minutes from target (+600)",
    y = "Count")

# MEM +600 min
p4 <- ggplot(MEM_wake_check, aes(x = diff_min600)) +
  geom_histogram(binwidth = 20, col = "black", fill = "turquoise") +
  geom_vline(xintercept = c(-7.5, 7.5), col = "red",
    linetype = "dashed") +
  geom_vline(xintercept = c(-15, 15), color = "hotpink",
    linetype = "dashed") +
  theme_lucid() +
  scale_x_continuous(breaks = seq(-50, 360, by = 50)) +
  scale_y_continuous(breaks = seq(0, 40, by = 10),
    limits = c(0, 40)) +
  labs(title = "Minutes from +600 Minute Sampling Time",
    subtitle = "MEM",
    x = "Minutes from target (+600)",
    y = "Count")


# Arrange
grid.arrange(p1, p3, p2, p4, nrow = 2)


################################################################################
###                               Question 3                                 ###
################################################################################

##### Create data frame for question 3
## Make a subset of the data with the relevant cols for Q3

# Vector of the cols of interest
cols <- c("SubjectID",
  "Collection Date", "DAYNUMB",
  "Collection Sample",
  "Booket: Clock Time", "MEMs: Clock Time",
  "Sleep Diary reported wake time",
  "Cortisol (nmol/L)", "DHEA (nmol/L)")

# Subset the original data frame
wake_data_q3 <- proj0_dat[, colnames(proj0_dat) %in% cols]

### Modify data for Person 3029
wake_data_q3 <- wake_data_q3 %>%
  mutate(DAYNUMB = case_when(
    SubjectID == 3029 & `Collection Date` == "10/7/2018" ~ 1,
    SubjectID == 3029 & `Collection Date` == "10/10/2018" ~ 2,
    SubjectID == 3029 & `Collection Date` == "10/12/2018" ~ 3,
    TRUE ~ DAYNUMB
  )) %>%
  ungroup()

# Convert the Sleep diary Wake time to minutes
wake_data_q3$`Sleep Diary Wake Minutes` <-
  # Convert to date objects and convert to minutes
  lubridate::hour(hm(wake_data_q3$`Sleep Diary reported wake time`)) * 60 +
  lubridate::minute(hm(wake_data_q3$`Sleep Diary reported wake time`))

# Move this to after Collection Sample
wake_data_q3 <- wake_data_q3 %>%
  relocate(`Sleep Diary Wake Minutes`,
    .after = `Sleep Diary reported wake time`)

# Convert the Booklet Time to minutes
wake_data_q3$Booklet_Minutes <-
  # Convert to date objects and convert to minutes
  lubridate::hour(hm(wake_data_q3$`Booket: Clock Time`)) * 60 +
  lubridate::minute(hm(wake_data_q3$`Booket: Clock Time`))

# Move this to after Booklet: Clock Time
wake_data_q3 <- wake_data_q3 %>%
  relocate(Booklet_Minutes, .after = `Booket: Clock Time`)

# Convert the MEM Time to minutes
wake_data_q3$MEM_Minutes <-
  # Convert to date objects and convert to minutes
  lubridate::hour(hm(wake_data_q3$`MEMs: Clock Time`)) * 60 +
  lubridate::minute(hm(wake_data_q3$`MEMs: Clock Time`))

# Move this to after MEMs: Clock Time
wake_data_q3 <- wake_data_q3 %>%
  relocate(MEM_Minutes, .after = `MEMs: Clock Time`)

# Create a Minutes SINCE WAKE variable
wake_data_q3 <- wake_data_q3 %>%
  group_by(SubjectID, DAYNUMB) %>%
  mutate(wake_min = Booklet_Minutes[`Collection Sample` == 1][1],
    minutes_since_wake = Booklet_Minutes - wake_min) %>%
  ungroup() %>%
  relocate(minutes_since_wake, .after = `Sleep Diary Wake Minutes`)

# Rename cortisol col
wake_data_q3 <- wake_data_q3 %>%
  rename(Cortisol = `Cortisol (nmol/L)`)

#################################### DHEA ######################################

### Remove DHEA values that equaled 5.205
wake_data_q3 <- wake_data_q3 %>%
  filter(
    is.na(`DHEA (nmol/L)`) |
      !(dplyr::near(`DHEA (nmol/L)`, 5.205) | SubjectID == "3037"))

# 14 observations removed - 358 observations left

## Visualize DHEA over time
ggplot(wake_data_q3, aes(x = `Collection Sample`,
  y = `DHEA (nmol/L)`,
  color = factor(SubjectID))) +
  geom_point() +
  geom_line() +
  facet_wrap(~DAYNUMB)

##### Data checks

# A plot of the distribution of the outcome
ggplot(wake_data_q3, aes(x = `DHEA (nmol/L)`)) +
  geom_histogram(bins = 20, fill = "palevioletred2",
    col = "black") +
  theme_lucid()

# A plot of the distribution of the outcome
ggplot(wake_data_q3, aes(x = log(`DHEA (nmol/L)`))) +
  geom_histogram(bins = 20, fill = "palevioletred2",
    col = "black") +
  theme_lucid()

# Log transform DHEA


######## Model

#### Fit the LMM

# Create a logged DHEA variable
wake_data_q3$log_DHEA <- log(wake_data_q3$`DHEA (nmol/L)`)

# Knot at 30 minutes - Collection Sample 2
wake_data_q3 <- wake_data_q3 %>%
  mutate(
    # Slope from waking to 30 miutes
    time_pre30  = pmin(minutes_since_wake, 30),
    # Slope from after 30 minutes
    time_post30 = pmax(minutes_since_wake - 30,
      0)
  )

# Center Sleep Diary wake minutes
wake_data_q3 <- wake_data_q3 %>%
  mutate(Sleep_Centered = `Sleep Diary Wake Minutes` -
    mean(`Sleep Diary Wake Minutes`, na.rm = TRUE),
  Book_MinC = Booklet_Minutes - mean(Booklet_Minutes,
    na.rm = TRUE)) %>%
  relocate(Sleep_Centered, .after = `Sleep Diary Wake Minutes`)
# Mean wake up time is - 409.4545 minutes


##### Model
model_q3 <- lme(log_DHEA ~ time_pre30 + time_post30,
  random = ~ 1 | SubjectID,
  method = "REML",
  na.action = na.omit,
  data = wake_data_q3)

# summary(model_q3)
# confint(model_q3)

# Display model summary with kable
# Extract t and p-values from model for data frame
summary_df2 <- as.data.frame(summary(model_q3)$tTable)

# kable
kable(
  summary_df2,
  format = "latex",
  booktabs = TRUE,
  digits = 5,
  caption = "Model Summary"
) |>
  kable_styling(latex_options = c("striped", "hold_position"))

#### Back transform betas
# Pre 30
a <- -0.01724 * 30
a_orig <- (exp(a) - 1) * 100 # -38.52871

# Confidence interval
a1 <- (exp(-0.02313 * 30) - 1) * 100 # -71.6346
a2 <- (exp(-0.0113 * 30) - 1) * 100 # -25.92

# Post 30
b <- -0.00153 * 600
b_org <- (exp(b) - 1) * 100 # -37.07332
b1 <- (exp(-0.00183 * 600) - 1) * 100
b2 <- (exp(-0.00122 * 600) - 1) * 100


#### Visualize DHEA trajectories with knot at 30 minutes

# For plotting purposes here, we use a model with lmer()!!
model_plot_dhea <- lmer(log_DHEA ~ time_pre30 + time_post30
  + (1 | SubjectID),
REML = TRUE,
data = wake_data_q3)

t_max <- quantile(wake_data_q3$minutes_since_wake, 0.95, na.rm = TRUE)

pred_grid <- tibble(
  minutes_since_wake = seq(0, as.numeric(t_max), by = 5)) %>%
  mutate(
    time_pre30  = pmin(minutes_since_wake, 30),
    time_post30 = pmax(minutes_since_wake - 30, 0)
  )

# Get fixed effects predictions
pred_grid$fit <- predict(model_plot_dhea, newdata = pred_grid, re.form = NA)

# Add a 95% CI
X <- model.matrix(~ time_pre30 + time_post30, data = pred_grid)
V <- vcov(model_plot_dhea)
se_fit <- sqrt(diag(X %*% V %*% t(X)))

pred_grid <- pred_grid %>%
  mutate(
    lwr = fit - 1.96 * se_fit,
    upr = fit + 1.96 * se_fit
  )

# Plot
ggplot() +
  geom_point(data = wake_data_q3,
    aes(x = minutes_since_wake, y = log_DHEA),
    alpha = 0.25, size = 1) +
  geom_ribbon(data = pred_grid,
    aes(x = minutes_since_wake,
      ymin = lwr,
      ymax = upr),
    alpha = 0.2,
    fill = "blue") +
  geom_line(data = pred_grid,
    aes(x = minutes_since_wake, y = fit),
    linewidth = 1,
    col = "blue") +
  geom_vline(xintercept = 30, linetype = 2,
    col = "red") +
  labs(x = "Minutes since waking", y = "log(DHEA)",
    title = "DHEA Pattern Over Time",
    subtitle = "Knot at 30 minutes") +
  theme_lucid()


###### Model diagnostics
# Obtain data used to fit model
model_dataq3 <- model_q3@frame

# Obtain residuals of model
model_dataq3$resids <- resid(model_q3)

# Obtain fitted values
model_dataq3$fitted <- fitted(model_q3)

# Residual histogram
ggplot(model_dataq3, aes(x = resids)) +
  geom_histogram(bins = 25,
    fill = "cadetblue2",
    color = "black") +
  theme_lucid()

# residual vs. fitted
ggplot(model_dataq3, aes(x = fitted, y = resids)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_classic()

################################## Cortisol ##################################

# Determine how many measurements are above 80
# wake_data_q3 %>%
#   filter(`Cortisol (nmol/L)` >= 26) %>%
#   reframe(SubjectID, `Collection Sample`,
#           DAYNUMB,`Cortisol (nmol/L)`)

# Exclude the observation with 89 cortisol
wake_data_q3 <- wake_data_q3 %>%
  filter(!dplyr::near(Cortisol, 89.55714))

# Visualize cortisol over time
ggplot(wake_data_q3, aes(x = `Collection Sample`,
  y = Cortisol,
  color = factor(SubjectID))) +
  geom_point() +
  geom_line() +
  facet_wrap(~DAYNUMB) +
  theme_lucid()

######### Model


# Model
model_cortisol <- lme(Cortisol ~ time_pre30 + time_post30,
  random = ~ 1 | SubjectID,
  method = "REML",
  na.action = na.omit,
  data = wake_data_q3)

# summary(model_cortisol)
# confint(model_cortisol)

# Extract t and p-values from model for data frame
summary_df3 <- as.data.frame(summary(model_cortisol)$tTable)

# kable
kable(
  summary_df3,
  format = "latex",
  booktabs = TRUE,
  digits = 5,
  caption = "Model Summary"
) |>
  kable_styling(latex_options = c("striped", ""))


###### Model diagnostics
# Obtain data used to fit model
model_data_cort <- model_cortisol@frame

# Obtain residuals of model
model_data_cort$resids <- resid(model_cortisol)

# Obtain fitted values
model_data_cort$fitted <- fitted(model_cortisol)

# Residual histogram
ggplot(model_data_cort, aes(x = resids)) +
  geom_histogram(bins = 25,
    fill = "cadetblue2",
    color = "black")

# residual vs. fitted
ggplot(model_data_cort, aes(x = fitted, y = resids)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_classic()

#### Visualize
# For plotting purposes here, we use a model with lmer()!!
model_plot_cort <- lmer(Cortisol ~ time_pre30 + time_post30
  + (1 | SubjectID),
REML = TRUE,
data = wake_data_q3)

t_max <- quantile(wake_data_q3$minutes_since_wake, 0.95, na.rm = TRUE)

pred_grid <- tibble(
  minutes_since_wake = seq(0, as.numeric(t_max), by = 5)) %>%
  mutate(
    time_pre30  = pmin(minutes_since_wake, 30),
    time_post30 = pmax(minutes_since_wake - 30, 0))

# Get fixed effects predictions
pred_grid$fit_cort <- predict(model_plot_cort,
  newdata = pred_grid,
  re.form = NA)

# Add a 95% CI
X_cortisol <- model.matrix(~ time_pre30 + time_post30,
  data = pred_grid)
V_cortisol <- vcov(model_plot_cort)

se_fit_cort <- sqrt(diag(X_cortisol %*% V_cortisol %*% t(X_cortisol)))

pred_grid <- pred_grid %>%
  mutate(
    lwr_cort = fit_cort - 1.96 * se_fit_cort,
    upr_cort = fit_cort + 1.96 * se_fit_cort
  )

# Plot
ggplot() +
  geom_point(data = wake_data_q3,
    aes(x = minutes_since_wake,
      y = Cortisol),
    alpha = 0.25, size = 1) +
  geom_ribbon(data = pred_grid,
    aes(x = minutes_since_wake,
      ymin = lwr_cort,
      ymax = upr_cort),
    alpha = 0.2,
    fill = "blue") +
  geom_line(data = pred_grid,
    aes(x = minutes_since_wake,
      y = fit_cort),
    linewidth = 1,
    col = "blue") +
  geom_vline(xintercept = 30, linetype = 2,
    col = "red") +
  labs(x = "Minutes since waking", y = "Cortisol",
    title = "Cortisol Pattern Over Time",
    subtitle = "Knot at 30 minutes") +
  theme_lucid()
