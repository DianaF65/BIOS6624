# Code for Project0

# Libraries
library(powertools)
library(here)
library(readr)

# Read in the data
dat <- read_csv(here("Project2", "Data", "PrelimData.csv"))

### Calculating correlations for power analyses
# Correlations are used in power calculations when measuring the strength
# and direction of a relationship between two continuous variables instead of
# comparing group means

colnames(dat)
# CVLT_CNG3 - Measure of episodic memory
#   - Recall/Recognition cognitive outcome variable
# CORT_CNG3 - Measure of corticol thickness
# IL_6 - A novel Cytokine/chemokine of interest
# MCP_1 - A tested cytokine/chemokine of interest

#### Aim 1
# Correlation between IL6 and episodic memory
cor(dat$IL_6, dat$CVLT_CNG3)
# -0.2585617

# Correlation between MCP1 and episodic memory
cor(dat$MCP_1, dat$CVLT_CNG3)
# -0.3183585

# Correlation between IL6 and corticol thickness
cor(dat$IL_6, dat$CORT_CNG3)
# -0.5993366

# Correlation between MCP1 and corticol thickness
cor(dat$MCP_1, dat$CORT_CNG3)
# -0.685273

#### Correlation between cytokines/chemokines
# Correlation between IL6 and MCP1
cor(dat$IL_6, dat$MCP_1)
# 0.9335558
# - Pretty highly correlated so perhaps we can reduce the number of chemokines/
#cytokines to just 1 instead of examining all 4 which will also address the
# multiple testing issue

### OR JUST DO THIS
cor(dat)

### Pre power analysis calculations
# For Aim 1, correlations range from -0.26 to -0.69
# We will use the corr.1samp() function
# NOTE: The correlation coefficient becomes an effect size when squared
# AKA the r-squared value
test1 <- corr.1samp(N = 175,
                    rho0 = 0,
                    rhoA = (0.26)^2,
                    alpha = 0.05,
                    power = NULL,
                    sides = 1)

test1

test2 <- corr.1samp(N = 175,
                    rho0 = 0,
                    rhoA = 0.7,
                    alpha = 0.05,
                    power = NULL,
                    sides = 1)

test2



