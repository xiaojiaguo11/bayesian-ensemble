# -----------------------------------------------------------
# This file generates trimmed means and metrically trimmed
# means for the 100,000 time series of the M4 competition.
# -----------------------------------------------------------

# --------------------------
# Load the data
# --------------------------
load("M4_valid_forec_100k.RData") # Load in the base model predictions in the validation set
load("M4.RData") # Load in the M4 time series data
library(M4metaresults) # Load in the base model predictions in the test set

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

# --------------
# Trimmed mean
# --------------

# --------------------------------------------------------
# Find optimal trimming level using the validation set
# --------------------------------------------------------

# Define possible trimming levels

trim_value <- c(0, 0.2, 0.5) # 0: simple mean; 0.2: trim away 1 model on each side; 0.5: median

# Set up matrices for the predictions and accuricies

trimmed_mean <- matrix(NA,length(baseForecasts_valid),48)
smape_grid <- rep(NA, length(trim_value))
mase_grid <- rep(NA, length(trim_value))
owa_grid <- rep(NA, length(trim_value))

# Generate forecasts using each trimming level defined in trim_value and calculate the accuracy

for (j in 1:length(trim_value)) {
  trim_j <- trim_value[j]
  for (i in 1:length(baseForecasts_valid)) {
    # Select forecasts only from the arima, ets, tbats, random walk with drift, and the theta model
    M4_baseForecasts_1 <- baseForecasts_valid[[i]]$ff[c(1,2,4,6,7),] 
    # Replace negative forecasts with zeros
    M4_baseForecasts_1[M4_baseForecasts_1 < 0] <-  0
    h <- baseForecasts_valid[[i]]$h
    
    # Calculate trimmed mean
    for (k in 1:h) {
      x <- M4_baseForecasts_1[,k]
      trimmed_mean[i,k] <- mean(x,trim = trim_j) 
    }
  }
  
  # Calculate SMAPE and MASE for each series
  smape_all <- rep(NA,length(baseForecasts_valid)) 
  mase_all <- rep(NA,length(baseForecasts_valid))
  for (i in 1:length(baseForecasts_valid)) {
    ts_i <- baseForecasts_valid[i]
    actual_i <- ts_i[[1]]$xx
    smape_all[i] <- mean(smape_cal(actual_i,trimmed_mean[i,1:length(actual_i)]))
    mase_all[i] <- mean(mase_cal(baseForecasts_valid[[i]]$x, actual_i,trimmed_mean[i,1:length(actual_i)]))
  }
  
  # Calculate the average SMAPE, MASE and OWA across all series
  smape_grid[j] <- mean(smape_all,na.rm = TRUE)
  mase_grid[j] <- mean(mase_all[mase_all<Inf],na.rm = TRUE) # Note that mase equals Inf for some series
  owa_grid[j] <- (smape_grid[j]/13.64627914 + mase_grid[j]/2.070712887)/2
}

# Print out accuracy results

smape_grid
mase_grid
owa_grid

# Find the optimal trimming level

trim_opt <- trim_value[which.min(owa_grid)] # the optimal trimming level is 0.2
trim_opt

# --------------------------------------------------------
# Generate trimmed mean forecasts in the test set 
# --------------------------------------------------------

trimmed_mean <- matrix(NA,length(submission_M4),48)
for (i in 1:length(submission_M4)) {
  # Select forecasts only from the arima, ets, tbats, random walk with drift, and the theta model
  M4_baseForecasts_1 <- submission_M4[[i]]$ff[c(1,2,4,6,7),]
  # Replace negative forecasts with zeros
  M4_baseForecasts_1[M4_baseForecasts_1 < 0] <- 0
  h <- submission_M4[[i]]$h
  
  # Calculate trimmed mean values
  for (k in 1:h) {
    x <- M4_baseForecasts_1[,k]
    trimmed_mean[i,k] <- mean(x,trim = trim_opt) 
  }
}

# Calculate SMAPE and MASE
smape_all <- rep(NA,length(submission_M4)) 
mase_all <- rep(NA,length(submission_M4))
for (i in 1:length(submission_M4)) {
  actual_i <- as.numeric(M4[[i]]$xx)
  smape_all[i] <- mean(smape_cal(actual_i,trimmed_mean[i,1:length(actual_i)]))
  mase_all[i] <- mean(mase_cal(M4[[i]]$x, actual_i,trimmed_mean[i,1:length(actual_i)]))
}

# --------------------------------------------------------
# Print out SMAPE, MASE and OWA
# --------------------------------------------------------

# Find the indices of yearly, quarterly, monthly and others series
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

# --------------------------
# Metrically trimmed mean
# --------------------------

# --------------------------------------------------------
# Find optimal z_star using the validation set
# --------------------------------------------------------

# Define possible z_star values

z_star_value <- c(0,0.3,0.5,1,3,5)

# Set up matrices for the predictions and accuricies

metrically_trimmed_mean <- matrix(NA,length(baseForecasts_valid),48)
smape_grid <- rep(NA, length(z_star_value))
mase_grid <- rep(NA, length(z_star_value))
owa_grid <- rep(NA, length(z_star_value))

# Generate forecasts using each z_star and calculate the accuracy

for (m in 1:length(z_star_value)) {
  z_star <- z_star_value[m]
  for (i in 1:length(baseForecasts_valid)) {
    # Select forecasts only from the arima, ets, tbats, random walk with drift, and the theta model
    M4_baseForecasts_1 <- baseForecasts_valid[[i]]$ff[c(1,2,4,6,7),] 
    # Replace negative forecasts with zeros
    M4_baseForecasts_1[M4_baseForecasts_1 < 0] <-  0
    h <- baseForecasts_valid[[i]]$h
    
    # Calculate metrically trimmed means
    for (j in 1:h) {
      x <- M4_baseForecasts_1[,j]
      simple_mean <- mean(x)
      # Find forecasts that are between the upper and lower thresholds 
      metrically_trimmed_x <- x[x > simple_mean - z_star * sd(x) 
                                & x < simple_mean + z_star * sd(x)]
      # Calculte the mean of the selected forecasts. 
      # If none of the forecasts are selected, use median as the ensemble forecast.
      metrically_trimmed_mean[i,j] <- ifelse(length(metrically_trimmed_x) > 0, 
                            mean(metrically_trimmed_x), 
                            median(x))
    }
  }
    
  # Calculate SMAPE and MASE for each series
  smape_all <- rep(NA,length(baseForecasts_valid)) 
  mase_all <- rep(NA,length(baseForecasts_valid))
  for (i in 1:length(baseForecasts_valid)) {
    ts_i <- baseForecasts_valid[i]
    actual_i <- ts_i[[1]]$xx
    smape_all[i] <- mean(smape_cal(actual_i,metrically_trimmed_mean[i,1:length(actual_i)]))
    mase_all[i] <- mean(mase_cal(baseForecasts_valid[[i]]$x, actual_i,metrically_trimmed_mean[i,1:length(actual_i)]))
  }

  # Calculate the average SMAPE, MASE and OWA across all series
  smape_grid[m] <- mean(smape_all,na.rm = TRUE)
  mase_grid[m] <- mean(mase_all[mase_all<Inf],na.rm = TRUE)
  owa_grid[m] <- (smape_grid[m]/13.64627914 + mase_grid[m]/2.070712887)/2
}

# Print out accuracy results

smape_grid
mase_grid
owa_grid

# Find the optimal z_star

# The optimal z_star is 0, which means taking the median as the ensemble forecast gives us the lowest owa.
z_star_opt <- z_star_value[which.min(owa_grid)] 
z_star_opt

# -----------------------------------------------------------
# Generate metrically trimmed mean forecasts in the test set 
# -----------------------------------------------------------

metrically_trimmed_mean <- matrix(NA,length(submission_M4),48)

for (i in 1:length(submission_M4)) {
  # Select forecasts only from the arima, ets, tbats, random walk with drift, and the theta model
  M4_baseForecasts_1 <- submission_M4[[i]]$ff[c(1,2,4,6,7),]
  # Replace negative forecasts with zeros
  M4_baseForecasts_1[M4_baseForecasts_1 < 0] <- 0
  h <- submission_M4[[i]]$h
  
  # Calculate metrically trimmed mean
  for (j in 1:h) {
    x <- M4_baseForecasts_1[,j]
    simple_mean <- mean(x)
    # Find forecasts that are between the upper and lower thresholds 
    metrically_trimmed_x <- x[x > simple_mean - z_star_opt * sd(x) 
                              & x < simple_mean + z_star_opt * sd(x)]
    # Calculte the mean of the selected forecasts 
    # If none of the forecasts are selected, use the median as the ensemble forecast
    metrically_trimmed_mean[i,j] <- ifelse(length(metrically_trimmed_x) > 0, 
                                           mean(metrically_trimmed_x), 
                                           median(x))
  }
}

# Calculate SMAPE and MASE
smape_all <- rep(NA,length(submission_M4)) 
mase_all <- rep(NA,length(submission_M4))
owa_all <- rep(NA,length(submission_M4))
for (i in 1:length(submission_M4)) {
  actual_i <- as.numeric(M4[[i]]$xx)
  smape_all[i] <- mean(smape_cal(actual_i,metrically_trimmed_mean[i,1:length(actual_i)]))
  mase_all[i] <- mean(mase_cal(M4[[i]]$x, actual_i,metrically_trimmed_mean[i,1:length(actual_i)]))
}

# --------------------------------------------------------
# Print out SMAPE, MASE and OWA
# --------------------------------------------------------

# Find the indices of yearly, quarterly, monthly and others series
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
