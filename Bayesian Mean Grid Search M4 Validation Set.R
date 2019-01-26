# -----------------------------------------------------------
# This file finds the optimal hyperparameters of the Bayesian 
# ensemble method using grid search in the validation set.
# -----------------------------------------------------------

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

# Load in the base model predictions in the validation set
load("M4_valid_forec_monthly_1.RData")
load("M4_valid_forec_monthly_2.RData")
load("M4_valid_forec_monthly_3.RData")
load("M4_valid_forec_monthly_4.RData")
load("M4_valid_forec_yearly.RData")
load("M4_valid_forec_quarterly.RData")
load("M4_valid_forec_daily_1.RData")
load("M4_valid_forec_daily_2.RData")
load("M4_valid_forec_hourly.RData")
load("M4_valid_forec_weekly.RData")

baseForecasts_valid <- c(baseForecasts_valid_monthly_1, 
                         baseForecasts_valid_monthly_2, 
                         baseForecasts_valid_monthly_3, 
                         baseForecasts_valid_monthly_4, 
                         baseForecasts_valid_yearly,
                         baseForecasts_valid_quarterly,
                         baseForecasts_valid_daily_1,
                         baseForecasts_valid_daily_2,
                         baseForecasts_valid_hourly,
                         baseForecasts_valid_weekly)

rm(baseForecasts_valid_monthly_1)
rm(baseForecasts_valid_monthly_2)
rm(baseForecasts_valid_monthly_3)
rm(baseForecasts_valid_monthly_4)
rm(baseForecasts_valid_yearly)
rm(baseForecasts_valid_quarterly)
rm(baseForecasts_valid_daily_1)
rm(baseForecasts_valid_daily_2)
rm(baseForecasts_valid_hourly)
rm(baseForecasts_valid_weekly)

# -------------------------------
# Set up for parallel processing
# -------------------------------

cores = detectCores()
cl <- makeCluster(cores)
registerDoSNOW(cl)

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

# --------------------------------------------------------------------
# Set the number of clusters and the number of experts in each cluster
# ---------------------------------------------------------------------

perm_comb <- c(2,1,2)

# -------------------------------------
# Find all the combinations to evaluate
# -------------------------------------

# Find all the combinations and store them into a matrix
n_model <- rep(1:length(perm_comb), perm_comb)
k <- length(n_model)
perms_matrix <- permutations(length(n_model),length(n_model),n_model,set=FALSE)
perms_matrix_unique <- unique(perms_matrix) # In this matrix, each row represents a possible combination. 1, 2 and 3 in the matrix denote which cluster the model belongs to. For example, (1,1,2,3,3) means that the first two models are in cluster 1, the third model is in cluster 2, and the fourth and fifth models are in cluster 3.

# For the convenience of calculating the Bayesian weights later in this file, rewrite the combination matrix "perms_matrix_unique" using number 1 to 5. For example, if perm_comb = c(2,1,2), now in matrix "perms_matrix_unique" 1 and 2 denote cluster 1, 3 denote cluster 2, and 4 and 5 denote cluster 3. 
for (l in 1:dim(perms_matrix_unique)[1]) { # for each possible combination
  for (j in length(perm_comb):2) { # for each cluster
    perms_matrix_unique[l,perms_matrix_unique[l,]==j] <- c(1:5)[(sum(perm_comb[1:(j-1)])+1) : (sum(perm_comb[1:(j-1)]) + perm_comb[j])]
  }
  perms_matrix_unique[l,perms_matrix_unique[l,]==1] <- c(1:perm_comb[1])
}
num_combs <- nrow(perms_matrix_unique) # get the number of combinations

# -------------------------------------------
# Set possible values of the hyperparameters
# -------------------------------------------

n0 <- 1 # Fix n0 to 1
n1_grid <- c(0.5, 1, 2, 5)
n2_grid <- c(1, 2, 5, 10, 25, 50)
n3_grid <- c(2, 5, 10, 25, 50)

# Set up matrices to store SMAPE, MASE and OWA values
smape_grid <- array(NA, c(length(n2_grid), length(n3_grid), length(n1_grid)))
mase_grid <- array(NA, c(length(n2_grid), length(n3_grid), length(n1_grid)))
owa_grid <- array(NA, c(length(n2_grid), length(n3_grid), length(n1_grid)))

dimnames(smape_grid)[[3]] <- paste('n1=', n1_grid, sep = '') 
dimnames(smape_grid)[[1]] <- paste('n2=', n2_grid, sep = '') 
dimnames(smape_grid)[[2]] <- paste('n3=', n3_grid, sep = '') 
dimnames(mase_grid)[[3]] <- paste('n1=', n1_grid, sep = '') 
dimnames(mase_grid)[[1]] <- paste('n2=', n2_grid, sep = '') 
dimnames(mase_grid)[[2]] <- paste('n3=', n3_grid, sep = '') 
dimnames(owa_grid)[[3]] <- paste('n1=', n1_grid, sep = '') 
dimnames(owa_grid)[[1]] <- paste('n2=', n2_grid, sep = '') 
dimnames(owa_grid)[[2]] <- paste('n3=', n3_grid, sep = '') 

# -----------------------------------------------------------------
# Evaluate Bayesian ensembles with different hyperparameters
# -----------------------------------------------------------------
for (b in 1:length(n1_grid)) { # For each possible value of n2
  n1 <- n1_grid[b]
  for (m in 1:length(n2_grid)) { # For each possible value of n2
    n2 <- n2_grid[m]
    for (d in 1:length(n3_grid)) { # For each possible value of n3
      n3 <- n3_grid[d]
      if (n1 >= n2) next # In the paper, we assume n1 < n2 < n3
      if (n2 >= n3) next
      
      # Define n0, n and delta
      n <- n <- c(rep(n1,perm_comb[1]), rep(n2, perm_comb[2]), rep(n3, perm_comb[3]))
      tau <- (1/(n0+n2))^2 # assume that n2 = optimal n, which is tau^(-0.5) - n0
      delta <- tau*n
      
      # Define matrix A, B, G and Delta
      A <- diag(1/n)
      b_vec <- t(t(n/(n0+n)))
      B <- diag(as.vector(b_vec))
      g_vec <- delta + n/(n0 + n)^2
      G <- diag(g_vec)
      Delta <- diag(delta)
      
      # This function returns the Bayesian ensemble prediction for series i
      parallelFunction <- function(i) {
        # select forecasts only from the arima, ets, tbats, random walk with drift, and the theta model
        M4_baseForecasts_1 <- baseForecasts_valid[[i]]$ff[c(1,2,4,6,7),]
        # replace negative forecasts with zeros
        M4_baseForecasts_1[M4_baseForecasts_1 < 0] <- 0
        h <- baseForecasts_valid[[i]]$h
        
        bayesian_mean <- rep(NA,h)
        bayesian_weights <- matrix(0, k, h)
        for (t in 1:h) {
          baseForecasts_t <- M4_baseForecasts_1[,t]
          # estimate theta_0 as the trimmed mean of the forecasts
          # theta_0 <- mean(baseForecasts_t, trim=0.3) 
          # estimate theta_0 as the naive forecast
          theta_0 <- as.numeric(baseForecasts_valid[[i]]$naive)
          # estimate sigma^2 (Eq. 54)
          s2 <- (1/k)*sum((baseForecasts_t-mean(baseForecasts_t))^2)
          # set s2 to a small value in the case that forecasts from all five methods are equal
          if (s2 == 0) s2 <- 0.1
          e <- t(t(rep(1,k)))
          sigma2 <- as.numeric(s2*solve(1/k*sum(g_vec+b_vec^2/n0) - (t(e)/k)%*%(G + b_vec%*%t(b_vec)/n0)%*%(e/k)))
          
          ordered_mu <- sort(baseForecasts_t) # find the ordered forecasts
          
          # Calculate the Bayesian mean of the t-step ahead forecasts
          numerators <- rep(NA, num_combs)
          weight_1 <- matrix(NA, k, num_combs)
          for (j in 1:num_combs) { # for each possible combination
            perm_j <- as.numeric(perms_matrix_unique[j,])
            
            # Calculate the permutation matrix of the jth combination. 
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
          
          # Calculate the final weights for the base forecasts
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
              numerators <- numerators + exp(as(-1/2*t(ordered_mu - e*theta_0)%*%solve(sigma2*(G_j+b_j%*%t(b_j)/n0))%*%(ordered_mu - e*theta_0), "mpfr")) # Calculate the sum of the numerators of w_pi
              weight_1 <- weight_1+(exp(as(-1/2*t(ordered_mu - e*theta_0)%*%solve(sigma2*(G_j+b_j%*%t(b_j)/n0))%*%(ordered_mu - e*theta_0), "mpfr"))[1]*perm_matrix%*%solve(G)%*%b_vec) # Calculate the numerator of the final weights (Eq. 47)
            }
            bayesian_weights[,t] <- as.numeric(weight_1/t(numerators%*%t(rep(1,k)))/as.vector((n0+t(b_vec)%*%solve(G)%*%b_vec)))
            bayesian_mean[t] <- (1 - sum(bayesian_weights[,t]))*theta_0 + ordered_mu%*%t(t(bayesian_weights[,t]))
          }
        }
        
        # return the final weights and the bayesian ensemble forecasts
        outcome <- list(weights = bayesian_weights, prediction = bayesian_mean) 
        return(outcome)
      }
      
      # Calculate the Bayesian mean in the validation set
      tempMatrix <- c()
      resultsMatrix <- c()
      set.seed(201)
      resultsMatrix <- finalMatrix <- foreach(i=1:length(baseForecasts_valid), .combine=c) %dopar% {
        library(Rmpfr)
        tempMatrix = parallelFunction(i)
        c(resultsMatrix, tempMatrix)
      }
      
      # Calculate SMAPE, MASE, and OWA
      smape_all <- matrix(NA,length(baseForecasts_valid),48) 
      mase_all <- matrix(NA,length(baseForecasts_valid),48) 
      for (i in 1:length(baseForecasts_valid)) {
        actual_i <- baseForecasts_valid[[i]]$xx
        prediction_i <- resultsMatrix[i*2]$prediction
        smape_all[i,1:length(actual_i)] <- smape_cal(actual_i,prediction_i)
        mase_all[i,1:length(actual_i)] <- mase_cal(baseForecasts_valid[[i]]$x, actual_i,prediction_i)
      }
      smape_grid[m,d,b] <- mean(rowMeans(smape_all, na.rm=TRUE),na.rm=TRUE)
      mase_grid[m,d,b] <- mean(rowMeans(mase_all,na.rm=TRUE)[rowMeans(mase_all,na.rm=TRUE) < Inf], na.rm=TRUE)
      owa_grid[m,d,b] <- (smape_grid[m,d,b]/13.64627914 + mase_grid[m,d,b]/2.070712887)/2
    }
  }
}

# print out accuracy results
smape_grid
mase_grid
owa_grid

# -------------------------------------------
# Find and print the optimal hyperparameters 
# -------------------------------------------

n_hyper <- which(owa_grid == min(owa_grid, na.rm=TRUE), arr.ind = TRUE)
cat(paste(' optimal n1:', n1_grid[n_hyper[3]],"\n",
          'optimal n2:', n2_grid[n_hyper[1]],"\n",
          'optimal n3:', n3_grid[n_hyper[2]],"\n",
          'lowest owa:', min(owa_grid, na.rm=TRUE)))
