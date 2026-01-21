#Implementing the automatic democratic method

library(dplyr)
library(lpSolve)
library(nnls)


#Normalization function
normalize <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))  # constant column safeguard
  (x - rng[1]) / (rng[2] - rng[1])
}


# Normalize data assuming that the data has columns arranged as: Species, Probability_Ratio, Rel_ab, Prevalence
data_norm <- data_norm_olm %>%
  filter(!is.na(Probability_Ratio),
         !is.na(Rel_ab),
         !is.na(Prevalence))

data_norm$Probability_Ratio <- normalize(data_norm$Probability_Ratio)
data_norm$Rel_ab            <- normalize(data_norm$Rel_ab)
data_norm$Prevalence        <- normalize(data_norm$Prevalence)


#Computing the optimal score of each species
compute_optimal <- function(traits) {
  if(any(is.na(traits)) || any(is.infinite(traits))) return(NA)
  
  n <- length(traits)
  f.obj <- traits
  f.con <- matrix(rep(1, n), nrow = 1)
  f.dir <- "="
  f.rhs <- 1
  
  lp_solution <- lp("max", f.obj, f.con, f.dir, f.rhs)
  
  if(lp_solution$status != 0) return(NA)
  
  sum(lp_solution$solution * traits)
}

data_norm$OptimalScore <- apply(
  data_norm[, c("Probability_Ratio","Log_ab","Prevalence")],
  1, compute_optimal
)

#Regression step using Non-Negative Least Regression (NNLS)
X <- as.matrix(data_norm[, c("Probability_Ratio","Log_ab","Prevalence")])
y <- data_norm$OptimalScore

nn <- nnls(X, y)                 
coefs <- coef(nn)                
weights_log <- coefs / sum(coefs)    # this involves normalizing the sum to 1
weights_log



fit <- lm.fit(X, y)
coefs <- fit$coefficients
weights_log_lm<- coefs / sum(coefs)
weights_log_lm

#Generated weights
weights_df <- data.frame(
  Feature = c("Probability_Ratio","Rel_ab","Prevalence"),
  Weight  = round(weights, 4)
)

print(weights_df)
