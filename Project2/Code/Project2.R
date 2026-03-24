# Code for Project0

# Libraries
library(powertools)
library(here)
library(readr)

# Read in the data
dat <- read_csv(here("Project2", "Data", "PrelimData.csv"))

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

### Pre power analysis calculations


