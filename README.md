# Bayesian Ensembles of Unlabeled Forecasts

This code is used in our paper "Bayesian Ensembles of Unlabeled Forecasts". 

The data file "M4_valid_forec_100k.RData" contains forecasts generated by nine base models in the temporal holdout dataset (validation set). These nine base models (ARIMA, TBATS, Random Walk with a Drift, ETS, THETA, STL Decomposition with AR residuals, Naive, Seasonal Naive, and Neural Network Autoregressive) are used by Montero-Manso et al. in the M4 competition. For the test set, "M4.Rdata" provides the actuals, and the base forecasts from the nine models can be obtained by loading the "M4metalearning" package.

The file "Bayesian Mean Grid Search M4 Validation Set.R" finds the optimal hyperparameters of the Bayesian ensemble method using grid search in the validation set. Running the file "Bayesian Mean M4 Test Set.R" with the optimal hyperparameters will generate Bayesian ensemble forecasts in the test set, and calculate the average SMAPE, MASE and OWA of these forecasts. 

The file "Heuristic Means.R" evaluates the accuracy of the trimmed means and metrically trimmed means in the test set. The optimal trimming level for these two methods are obtained by grid search in the validation set.
