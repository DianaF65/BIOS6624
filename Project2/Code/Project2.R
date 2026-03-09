# Code for Project0

# Libraries
library(powertools)

##### Worksheet 3 - Study Design Module

# Question 1 part b
W3_Q1_partb <- ttest.2samp(n1 = 65,
            n.ratio = 35/65,
            delta = 0.5,
            sd1 = 1,
            sd.ratio = 1,
            df.method = c("welch"),
            alpha = 0.05,
            power = NULL,
            sides = 2,
            v = FALSE)
W3_Q1_partb

# Question 2 part b
W3_Q2_partb <- ttest.2samp(n1 = 65,
                           # exposed / unexposed
                           n.ratio = 35/65,
                           delta = 0.5,
                           # exposed / unexposed
                           sd1 = 1.5,
                           sd.ratio = 1.5,
                           df.method = c("welch"),
                           alpha = 0.05,
                           power = NULL,
                           sides = 2,
                           v = FALSE)
W3_Q2_partb

W3_Q3_partb <- ttest.2samp(n1 = 65,
                           # exposed / unexposed
                           n.ratio = 35/65,
                           delta = 0.5,
                           # exposed / unexposed
                           sd1 = 1.5,
                           sd.ratio = 1/1.5,
                           alpha = 0.05,
                           power = NULL,
                           sides = 2,
                           v = FALSE)
W3_Q2_partb

#### Worksheet 4

# Question 5
W4_Q5 <- ttest.2samp(n1 = NULL,
                     n.ratio = 1,
                     delta = 8.7-5.3,
                     sd1 = 8.6,
                     sd.ratio = 1,
                      alpha = 0.05,
                      power = 0.8,
                      sides = 2
                      )
W4_Q5

# Question 6
W4_Q6 <- ttest.2samp(n1 = 25,
                     n.ratio = 75/25,
                     delta = 8.7-5.3,
                     sd1 = 8.6,
                     sd.ratio = 1,
                     alpha = 0.05,
                     power = NULL,
                     sides = 2
)
W4_Q6




