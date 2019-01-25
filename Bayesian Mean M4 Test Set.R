# ------------------------------------------------------
# This file generates the Bayesian ensemble forecasts
# for the 100,000 time series of the M4 competition.
# -----------------------------------------------------

# --------------------------
# Install and load packages.
# --------------------------

# install.packages("doSNOW")
# install.packages("foreach")
# install.packages("doParallel")
# install.packages("Rmpfr")
# install.packages("permute")
# install.packages("gtools")
# install.packages("Matrix")
# install.packages("Metrics")

# Packages for parallel processing
library(foreach)
library(doSNOW)
library(doParallel)

# Other packages
library(Rmpfr) # For creating high precision numbers
library(permute) # For calculating permutations
library(gtools) # For finding all the combinations
library(Matrix) 
library(Metrics)

# --------------------------
# Load the data
# --------------------------

load("M4.RData") # Load in the M4 time series data
library(M4metaresults) # Load in the base model predictions in the test set

# -------------------------------
# Set up for parallel processing
# -------------------------------

cores = detectCores()
cl <- makeCluster(cores)
registerDoSNOW(cl)

# -------------------------------
# Rearrange the M4 time series
# -------------------------------

yearly_submission_M4 <- Filter(function(l) l$period == "Yearly", submission_M4)
y_len <- length(yearly_submission_M4)  # 23,000
quarterly_submission_M4 <- Filter(function(l) l$period == "Quarterly", submission_M4)
q_len <-length(quarterly_submission_M4)  # 24,000
monthly_submission_M4 <- Filter(function(l) l$period == "Monthly", submission_M4)
m_len <-length(monthly_submission_M4)  # 48,000
weekly_submission_M4 <- Filter(function(l) l$period == "Weekly", submission_M4)
w_len <-length(weekly_submission_M4)  # 359
daily_submission_M4 <- Filter(function(l) l$period == "Daily", submission_M4)
d_len <-length(daily_submission_M4)  # 4,227
hourly_submission_M4 <- Filter(function(l) l$period == "Hourly", submission_M4)
h_len <-length(hourly_submission_M4)  # 414

rearranged_submission_M4 <- append(daily_submission_M4, hourly_submission_M4)
rearranged_submission_M4 <- append(rearranged_submission_M4, monthly_submission_M4)
rearranged_submission_M4 <- append(rearranged_submission_M4, quarterly_submission_M4)
rearranged_submission_M4 <- append(rearranged_submission_M4, weekly_submission_M4)
rearranged_submission_M4 <- append(rearranged_submission_M4, yearly_submission_M4)

submission_M4 <- rearranged_submission_M4

# -------------------
# Performance metrics
# -------------------

smape_cal <- function(outsample, forecasts){
  # Used to estimate sMAPE
  outsample <- as.numeric(outsample) ; forecasts<-as.numeric(forecasts)
  smape <- (abs(outsample-forecasts)*200)/(abs(outsample)+abs(forecasts))
  return(smape)
}

mase_cal <- function(insample, outsample, forecasts){
  # Used to estimate MASE
  frq <- frequency(insample)
  forecastsNaiveSD <- rep(NA,frq)
  for (j in (frq+1):length(insample)){
    forecastsNaiveSD <- c(forecastsNaiveSD, insample[j-frq])
  }
  masep<-mean(abs(insample-forecastsNaiveSD),na.rm = TRUE)
  outsample <- as.numeric(outsample) ; forecasts <- as.numeric(forecasts)
  mase <- (abs(outsample-forecasts))/masep
  return(mase)
}

# -----------------------------------------------
# Set hyper-parameters for the Bayesian ensemble
# -----------------------------------------------

# Set the number of clusters and the number of experts in each cluster
perm_comb <- c(2,1,2)

# Define n0, n and delta
n0 <- 1
n1 <- 0.5
n2 <- 5
n3 <- 85
# Assume that n2 = optimal n, which is tau^(-0.5) - n0
tau <- (1/(n0+n2))^2 
n <- c(rep(n1,perm_comb[1]), rep(n2, perm_comb[2]), rep(n3, perm_comb[3]))
delta <- tau*n

# -------------------------------------
# Find all the combinations to evaluate
# -------------------------------------

# Find all the combinations and store them into a matrix
n_model <- rep(1:length(perm_comb), perm_comb)
k <- length(n_model)
perms_matrix <- permutations(length(n_model),length(n_model),n_model,set=FALSE)
perms_matrix_unique <- unique(perms_matrix) # In this matrix, each row represents a possible combination. 1, 2 and 3 in the matrix denote which cluster the model belongs to. For example, (1,1,2,3,3) means that the first two models are in cluster 1, the third model is in cluster 2, and the fourth and fifth models are in cluster 3.

# For the convenience of calculating the Bayesian weights later in this file, rewrite the combination matrix "perms_matrix_unique" using number 1 to 5. For example, if perm_comb = c(2,1,2), now in matrix "perms_matrix_unique" 1 and 2 denote cluster 1, 3 denote cluster 2, and 4 and 5 denote cluster 3. 
for (l in 1:dim(perms_matrix_unique)[1]) { # For each possible combination
  for (j in length(perm_comb):2) { # For each cluster
    perms_matrix_unique[l,perms_matrix_unique[l,]==j] <- c(1:5)[(sum(perm_comb[1:(j-1)])+1) : (sum(perm_comb[1:(j-1)]) + perm_comb[j])]
  }
  perms_matrix_unique[l,perms_matrix_unique[l,]==1] <- c(1:perm_comb[1])
}
num_combs <- nrow(perms_matrix_unique) # Get the number of combinations

# -------------------------------------
# Define matrix A, B, G and Delta
# -------------------------------------

A <- diag(1/n)
b_vec <- t(t(n/(n0+n)))
B <- diag(as.vector(b_vec))
Delta <- diag(delta)
g_vec <- delta + n/(n0 + n)^2
G <- diag(g_vec)

# -----------------------------------------------------
# Define a function to calculate the Bayesian mean
# -----------------------------------------------------

# This function returns the Bayesian ensemble prediction for series i
parallelFunction <- function(i) {
  # Select forecasts only from the arima, ets, tbats, random walk with drift, and the theta model
  M4_baseForecasts_1 <- submission_M4[[i]]$ff[c(1,2,4,6,7),] 
  # Replace negative forecasts with zeros
  M4_baseForecasts_1[M4_baseForecasts_1 < 0] <- 0 
  h <- submission_M4[[i]]$h
  
  bayesian_mean <- rep(NA,h)
  bayesian_weights <- matrix(0, k, h)
  for (t in 1:h) {
    baseForecasts_t <- M4_baseForecasts_1[,t]
    # Estimate theta_0 as the trimmed mean of the forecasts
    # theta_0 <- mean(baseForecasts_t, trim=0.3) 
    # Estimate theta_0 as the naive forecast
    theta_0 <- as.numeric(tail(submission_M4[[i]]$x,1))
    # Estimate sigma^2 (Eq. 54)
    s2 <- (1/k)*sum((baseForecasts_t-mean(baseForecasts_t))^2) 
    e <- t(t(rep(1,k)))
    sigma2 <- as.numeric(s2*solve(1/k*sum(g_vec+b_vec^2/n0) - (t(e)/k)%*%(G + b_vec%*%t(b_vec)/n0)%*%(e/k)))
    
    ordered_mu <- sort(baseForecasts_t) # Find the ordered forecasts
    
    # Calculate the Bayesian mean of the t-step ahead forecast
    numerators <- rep(NA, num_combs)
    weight_1 <- matrix(NA, k, num_combs)
    for (j in 1:num_combs) { # For each possible combination
      perm_j <- as.numeric(perms_matrix_unique[j,])
      
      # Define the permutation matrix of the jth combination. 
      perm_matrix <- matrix(0,k,k)
      for (r in 1: k) {
        perm_matrix[r,perm_j[r]] = 1
      }
      
      # Calculate the numerator of w_pi (the weight of the jth combination, Eq. 56)
      G_j <- perm_matrix%*%G%*%solve(perm_matrix)
      b_j <- b_vec[perm_j]
      numerators[j] <- exp(-1/2*t(ordered_mu - e*theta_0)%*%solve(sigma2*(G_j+b_j%*%t(b_j)/n0))%*%(ordered_mu - e*theta_0))
      
      # Calculate the numerator of the final weights (Eq. 47)
      weight_1[,j] <- numerators[j]*perm_matrix%*%solve(G)%*%b_vec
    }
    
    # Calculate the final weights of the base forecasts
    bayesian_weights[,t] <- rowSums(weight_1)/as.vector((n0+t(b_vec)%*%solve(G)%*%b_vec))/sum(numerators)
    
    # Calculate the Bayesian mean (Eq. 46)
    bayesian_mean[t] <- (1 - sum(bayesian_weights[,t]))*theta_0 + ordered_mu%*%t(t(bayesian_weights[,t])) # Because of Eq. 48, the weight put on theta_0 can be calculated as (1 - sum(bayesian_weights[,t]))
    
    # Sometimes "numerators" is a zero vector because function "exp(x)" returns 0 if x < -745. 
    # We use the Rmpfr package to solve this problem
    if (is.na(bayesian_mean[t])) {
      weight_1 <- 0
      numerators <- 0
      for (j in 1:num_combs) {
        perm_j <- as.numeric(perms_matrix_unique[j,])
        perm_matrix <- matrix(0,k,k)
        for (r in 1: k) {
          perm_matrix[r,perm_j[r]] = 1
        }
        G_j <- perm_matrix%*%G%*%solve(perm_matrix)
        b_j <- b_vec[perm_j]
        numerators <-  numerators + exp(as(-1/2*t(ordered_mu - e*theta_0)%*%solve(sigma2*(G_j+b_j%*%t(b_j)/n0))%*%(ordered_mu - e*theta_0), "mpfr")) # Calculate the sum of the numerators of w_pi
        weight_1 <- weight_1+(exp(as(-1/2*t(ordered_mu - e*theta_0)%*%solve(sigma2*(G_j+b_j%*%t(b_j)/n0))%*%(ordered_mu - e*theta_0), "mpfr"))[1]*perm_matrix%*%solve(G)%*%b_vec) # Calculate the numerator of the final weights (Eq. 47)
      }
      bayesian_weights[,t] <- as.numeric(weight_1/t(numerators%*%t(rep(1,k)))/as.vector((n0+t(b_vec)%*%solve(G)%*%b_vec)))
      bayesian_mean[t] <- (1 - sum(bayesian_weights[,t]))*theta_0 + ordered_mu%*%t(t(bayesian_weights[,t]))
    }
  }
  
  # Return the final weights and the bayesian ensemble forecasts
  outcome <- list(weights = bayesian_weights, prediction = bayesian_mean) 
  return(outcome)
}

# ----------------------------------------------
# Calculate the Bayesian mean for the M4 series 
# ----------------------------------------------

tempMatrix <- c()
resultsMatrix <- c()
nlength <- length(submission_M4)
set.seed(201)
resultsMatrix <- finalMatrix <- foreach(i=1:nlength, .combine=c) %dopar% {
  library(Rmpfr)
  tempMatrix = parallelFunction(i)
  c(resultsMatrix, tempMatrix)
}

# ----------------------------------------------
# Calculate SMAPE and MASE
# ----------------------------------------------

smape_all <- rep(NA,nlength) 
mase_all <- rep(NA,nlength) 
for (i in 1:nlength) {
  actual_i <- M4[[i]]$xx
  prediction_i <- resultsMatrix[i*2]$prediction
  smape_all[i] <- mean(smape_cal(actual_i,prediction_i))
  mase_all[i] <- mean(mase_cal(M4[[i]]$x, actual_i,prediction_i))
}

# --------------------------------------------------------
# Print out SMAPE, MASE and OWA
# --------------------------------------------------------

# Find the indices of the yearly, quarterly, monthly and others series
nyearly <- which(sapply(submission_M4, function(x) x$period == 'Yearly')==TRUE)
nquarterly <- which(sapply(submission_M4, function(x) x$period == 'Quarterly')==TRUE)
nmonthly <- which(sapply(submission_M4, function(x) x$period == 'Monthly')==TRUE)
nothers <- which(! c(1:length(submission_M4)) %in% c(nyearly,nquarterly,nmonthly))

# Calculate the average smape, mase and owa for yearly, quarterly, monthly and other series
smape_vec <- c(mean(smape_all[nyearly]),mean(smape_all[nquarterly]),mean(smape_all[nmonthly]),mean(smape_all[nothers]), mean(smape_all))
mase_vec <- c(mean(mase_all[nyearly]),mean(mase_all[nquarterly]),mean(mase_all[nmonthly]),mean(mase_all[nothers]), mean(mase_all))
owa_vec <- c((mean(smape_all[nyearly])/16.34217 + mean(mase_all[nyearly])/3.97436)/2,
             (mean(smape_all[nquarterly])/11.01157 + mean(mase_all[nquarterly])/1.37138)/2,
             (mean(smape_all[nmonthly])/14.42683 + mean(mase_all[nmonthly])/1.06336)/2,
             (mean(smape_all[nothers])/4.75434 + mean(mase_all[nothers])/3.16930)/2,
             (mean(smape_all)/13.56407 + mean(mase_all)/1.91211)/2)
accuracy_matrix <- rbind(smape_vec,mase_vec,owa_vec)

rownames(accuracy_matrix) <- c('smape', 'mase', 'owa')
colnames(accuracy_matrix) <- c('Yearly', 'Quarterly', 'Monthly', 'Others', 'Average')

accuracy_matrix
