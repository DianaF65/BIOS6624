####### Simulation Worksheet 2

# Load libraries
library(hdrm)


# In worksheet 1 we simulated data from a regression model with two covariates.
# For the project we will work with a much larger set of covariates and investigate
# how well various model selection approaches work at identifying the variables 
# that are indeed associated with the outcome.  The purpose of this work is to 
# investigate the hdrm (high dimensional regression modeling package) for 
# simulating multiple regression models with a larger number of X variables.

#### Question 4

# Investigating the gen_data function from the package
testdata <- gen_data(100, 10, 5)

# a.	What is the sample size of this dataset?
## N = 100

# b.	How many X variables are there in the dataset?
## There are 10 variables

# c.	How many of the X variables are associated with the Y (outcome)? 
# Verify by checking the beta coefficients.
## We see that the first 5 are nonzero or are associated with Y

# d.	What are the means and SD’s of each X variable?
## 
# Set seed
set.seed(123)
# Means: 
colMeans(testdata$X)
# V01         V02         V03         V04         V05         V06 
# -0.10754680  0.12046511 -0.03622291  0.10585093 -0.04229996 -0.14964414 
# V07         V08         V09         V10 
# 0.10587355  0.09358971 -0.01919274  0.11994058 
# SDs: 
apply(testdata$X, MARGIN = 2, FUN = sd, na.rm=T)
# V01       V02       V03       V04       V05       V06       V07       V08 
# 0.9669866 0.9498790 1.0387812 0.9893458 0.9387282 1.0282366 1.0100100 1.0518066 
# V09       V10 
# 1.0203317 1.0429982


# e.	Are the X’s correlated with each other or independent?
## 
apply(testdata$X, MARGIN = 2, FUN = cor, testdata$X)
# V01          V02         V03          V04         V05           V06
# [1,]  1.00000000  0.030579031  0.04383271 -0.130622187  0.11448792  0.0781128798
# [2,]  0.03057903  1.000000000 -0.04486571 -0.024848379  0.01821008  0.0086851768
# [3,]  0.04383271 -0.044865707  1.00000000 -0.019259855 -0.08991276 -0.0637873268
# [4,] -0.13062219 -0.024848379 -0.01925986  1.000000000  0.20661771 -0.0065025161
# [5,]  0.11448792  0.018210080 -0.08991276  0.206617706  1.00000000 -0.0646926020
# [6,]  0.07811288  0.008685177 -0.06378733 -0.006502516 -0.06469260  1.0000000000
# [7,] -0.03312350 -0.115029463  0.16789328 -0.140653506  0.09415014  0.0002653148
# [8,] -0.04532832 -0.053281912 -0.16506751 -0.039583467  0.07436348 -0.1289100521
# [9,] -0.09149035 -0.014442042  0.24521164 -0.016043315 -0.02505532 -0.0204695736
# [10,] -0.02372887 -0.160557106  0.19574314  0.034231649 -0.04127698  0.0849030583
# V07           V08         V09         V10
# [1,] -0.0331234961 -0.0453283184 -0.09149035 -0.02372887
# [2,] -0.1150294630 -0.0532819116 -0.01444204 -0.16055711
# [3,]  0.1678932800 -0.1650675124  0.24521164  0.19574314
# [4,] -0.1406535056 -0.0395834667 -0.01604332  0.03423165
# [5,]  0.0941501419  0.0743634830 -0.02505532 -0.04127698
# [6,]  0.0002653148 -0.1289100521 -0.02046957  0.08490306
# [7,]  1.0000000000 -0.0008120057  0.02425240  0.14588393
# [8,] -0.0008120057  1.0000000000 -0.02272065 -0.08666105
# [9,]  0.0242524030 -0.0227206515  1.00000000  0.17228854
# [10,]  0.1458839328 -0.0866610453  0.17228854  1.00000000
## The Xs are not correlated with each other the correlation values are very
## small. 


# f.	What is the model error (sigma^2_e)? This will take some investigating. 
# Hint: Can you figure out yhat and get an estimate of sigma^2_e?
## 
model <- lm(testdata$y ~ .,
            data = as.data.frame(testdata$X))
yhat <- fitted(model)
residuals <- testdata$y - yhat
mean(residuals^2)
# ] 0.9102451

# g.	What is the SNR for test data?  SNR=R^2/(1-R^2), I think.
## 
# Get adjusted r-squared from model
r <- 0.438

# Formula for SNR
snr <- r/(1-(r))
#   0.7793594

##### Question 5
# 5.	Now set the beta coefficients to range between 0.5 and 2 for the X variables
# with an association and 0 otherwise.

# Create beta2 for the new beta coefficients
beta_nonzero <- runif(5, min = 0.5, max = 2)
beta <- c(beta_nonzero, rep(0, 5))
testdata_new <- gen_data(n = 100,
                         p = 10,
                         p1 = 5,
                         beta = beta)

# a.	What are the means and SD’s of each X variable?
# Set seed
set.seed(123)
# Means: 
colMeans(testdata_new$X)
# V01         V02         V03         V04         V05         V06 
# -0.04042018 -0.04307763  0.10498690 -0.15090051 -0.05498782 -0.00249087 
# V07         V08         V09         V10 
# -0.09704134 -0.09324871  0.10382732 -0.06696130 

# SDs: 
apply(testdata_new$X, MARGIN = 2, FUN = sd, na.rm=T)
# V01       V02       V03       V04       V05       V06       V07       V08 
# 0.9225508 0.8959588 1.0697922 0.8743863 1.0199163 1.0334833 0.9902554 0.9144675 
# V09       V10 
# 0.8733836 1.0121417 

# b.	Are the X’s correlated with each other or independent?
# apply(testdata_new$X, MARGIN = 2, FUN = cor, testdata$X)
# V01         V02         V03         V04          V05         V06
# [1,] -0.14160603 -0.10157362  0.10639085 -0.05447619 -0.026845176  0.07090806
# [2,] -0.01945526  0.08458584 -0.03721682  0.03022851 -0.002509401 -0.02165745
# [3,] -0.10899321 -0.01067434 -0.16598678 -0.16364230 -0.209390200 -0.04048147
# [4,] -0.05268314 -0.10647093  0.11843658 -0.07236968  0.072606586  0.04797665
# [5,] -0.15776267 -0.09128442  0.08142274 -0.05897889  0.120074839 -0.01257004
# [6,]  0.15578368  0.04694443  0.06023847 -0.10847234 -0.049566839 -0.10743750
# [7,] -0.09554410  0.09091271  0.12263666 -0.05773639  0.008745663  0.04457701
# [8,] -0.07304593  0.11968374  0.07997393  0.01449302 -0.082997780  0.29428390
# [9,]  0.05402104 -0.01915568 -0.21120226 -0.04791543 -0.201021922 -0.04888372
# [10,]  0.07428896 -0.07374662 -0.14324305 -0.20833156  0.228329440 -0.17167302
# V07         V08          V09         V10
# [1,]  0.01926931  0.04962184  0.094316708  0.13762323
# [2,]  0.06904835 -0.02309141 -0.085954791 -0.03813380
# [3,]  0.14921151  0.14579098  0.113951663  0.05445518
# [4,] -0.11471551 -0.01824112 -0.071452187 -0.01754340
# [5,]  0.12974362  0.05939958 -0.005565728  0.20026421
# [6,]  0.06757797  0.02928102  0.121477406 -0.04920250
# [7,]  0.16433411  0.11768799  0.044346161 -0.01202248
# [8,] -0.10225983 -0.03637480 -0.013968102  0.14246851
# [9,]  0.30278170 -0.02992681 -0.025385471  0.10050944
# [10,]  0.12825243  0.18138904  0.001061304 -0.15013264


# c.	What is the model error (sigma^2_e)? This will take some investigating. 
# Hint: Can you figure out yhat and get an estimate of sigma^2_e? 
## 
model_new <- lm(testdata_new$y ~ .,
            data = as.data.frame(testdata_new$X))
yhat_new <- fitted(model_new)
residuals_new <- testdata_new$y - yhat_new
mean(residuals_new^2)
# ] 0.9130782


# d.	What is the SNR for test data?  SNR=R^2/(1-R^2), I think.
# Get adjusted r-squared from model
r_new <- 0.9207

# Formula for SNR
snr_new <- r_new /(1-(r_new))
# 11.61034

##### Question 6
# Now set the beta coefficients to range between 0.5 and 2 for the X variables 
# with an association and 0 otherwise and introduce an exchangeable correlation of
# 0.5 between the X variables. 

# a.	What are the means and SD’s of each X variable?


# b.	What is the var-cov matrix of the X’s?


# c.	What is the model error (sigma^2_e)? 
# This will take some investigating. 
# Hint: Can you figure out yhat and get an estimate of sigma^2_e? 
  
  
# d.	What is the SNR for test data?  SNR=R^2/(1-R^2), I think.








