################################## Worksheet 1 #################################

# Load libraries
library(ggplot2)

# The purpose of this work is to write code to simulate from a multiple regression
# model and to evaluate some common statistical properties of the simulation.

# 1.) 
# Write a program to create 20 simulated datasets from the following linear 
# regression model.

# a.	N=100, X1= 0 or 1 with ½ in each group, 
# X2 = age, mean = 39 years and SD = 5, 
# X3 = the interaction between X1 and X2.  
# Set beta0 = 1, beta1 = 1, beta2 = 0.5, beta3 = 0.25. 
# Variance of model error = 1


# Set a seed
set.seed(7)

# Parameters
n <- 100
beta0 <- 1
beta1 <- 1
beta2 <- 0.5
beta3 <- 0.25
sigma <- 1   # because variance = 1, so SD = 1

# Create empty list to store the 20 datasets
dfs_list <- vector("list", 20)

# Generate 20 simulated datasets
for (i in 1:20) {
  
  X1 <- rbinom(n, size = 1, prob = 0.5)
  X2 <- rnorm(n, mean = 39, sd = 5)
  X3 <- X1 * X2
  error <- rnorm(n, mean = 0, sd = sigma)
  
  Y <- beta0 + beta1 * X1 + beta2 * X2 + beta3 * X3 + error
  
  dfs_list[[i]] <- data.frame(
    sim = i,
    X1 = X1,
    X2 = X2,
    X3 = X3,
    Y = Y
  )
}

# Look at first data set
head(dfs_list[[1]])

# 2.) 
# Plot of the simulated datasets with X2 on the X axis and the simulated
# outcome on the Y axis.  
# Should we make any changes to the parameters?  
# Is there a reasonable amount of error and variation across the simulated 
# datasets?

# Get all datasets into one data frame
all_sims <- do.call(rbind, dfs_list)

# Plot each simulated dataset
ggplot(all_sims, aes(x = X2, y = Y)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ sim) +
  labs(
    title = "Simulated Datasets",
    x = "X2 (Age)",
    y = "Simulated Outcome Y"
  ) +
  theme_minimal()

# 3.)







