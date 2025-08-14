
library(quantmod)
library(rugarch)
library(Metrics)
library(GAS)
library(forecast)
library(gasmodel)

getSymbols(c("BZ=F", "NG=F"), src = "yahoo", from = "2017-01-01", to = "2024-09-10")
getSymbols(c("CL=F"), src = "yahoo", from = "2017-01-01", to = "2024-09-14")

data_WTI <- `CL=F`
data_BRENT <- `BZ=F`
data_NG <- `NG=F`

price_WTI <- Cl(data_WTI)
price_BRENT <- Cl(data_BRENT)
price_NG <- Cl(data_NG)

returns_WTI <- na.omit(100 * diff(log(price_WTI)))
returns_BRENT <- na.omit(100 * diff(log(price_BRENT)))
returns_NG <- na.omit(100 * diff(log(price_NG)))

length(returns_BRENT)
length(returns_WTI)
length(returns_NG)

spec_GARCH <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0)),
                         distribution.model = "std")

spec_EGARCH <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                          mean.model = list(armaOrder = c(0, 0)),
                          distribution.model = "std")
# Define the function forecast_and_evaluate
forecast_and_evaluate <- function(returns, actual_returns) {
  pred_vol_GARCH <- c()
  pred_vol_EGARCH <- c()
  pred_vol_GAS <- c()
  pred_RETURN_GARCH <- c()
  pred_RETURN_EGARCH <- c()
  pred_RETURN_GAS <- c()
  actual_vol1 <-  c()
  actual_vol <-  c()
  
  for (i in 1:500) {
    window_size <- 1433 
    start_train <- i
    end_train <- i + window_size - 1
    train_data <- returns[start_train:end_train]
    
    # Fit GARCH
    fit_GARCH <- ugarchfit(spec = spec_GARCH, data = train_data, solver = "hybrid")
    forecast_GARCH <- ugarchforecast(fit_GARCH, n.ahead =1 )
    pred_val_GARCH <- as.numeric(sigma(forecast_GARCH))
    pred_return_GARCH <- as.numeric(fitted(forecast_GARCH))
    
    # Fit EGARCH
    fit_EGARCH <- ugarchfit(spec = spec_EGARCH, data = train_data, solver = "hybrid")
    forecast_EGARCH <- ugarchforecast(fit_EGARCH, n.ahead = 1)
    pred_val_EGARCH <- as.numeric(sigma(forecast_EGARCH))
    pred_return_EGARCH <- as.numeric(fitted(forecast_EGARCH))
    
    # Fit GAS model
    est_GAS <- gas(y = train_data, distr = "t", regress = "sep",
                   optim_arguments = list(opts = list(algorithm = "NLOPT_LN_NELDERMEAD", maxeval = 1000)))
    forecast_GAS <- gas_forecast(est_GAS, t_ahead = 1,method = "simulated_paths")
    pred_val_GAS <- sqrt(forecast_GAS$forecast$par_tv_ahead_mean[,"var"])
    pred_return_GAS <- forecast_GAS$forecast$y_ahead_mean
    
    # Combine forecasts
    pred_vol_GARCH <- c(pred_vol_GARCH, pred_val_GARCH)
    pred_vol_EGARCH <- c(pred_vol_EGARCH, pred_val_EGARCH)
    pred_vol_GAS <- c(pred_vol_GAS, pred_val_GAS)
    
    pred_RETURN_GARCH <- c(pred_RETURN_GARCH, pred_return_GARCH)
    pred_RETURN_EGARCH <- c(pred_RETURN_EGARCH, pred_return_EGARCH)
    pred_RETURN_GAS <- c(pred_RETURN_GAS, pred_return_GAS)
    
    train_mean <- mean(train_data)
    NORMreturns_train <- returns[(end_train+1) ] - train_mean
    actual_vol1 <- NORMreturns_train^2
    actual_vol <- c(actual_vol, actual_vol1)
    
  }
  
  # Trim to 500 values
  pred_vol_GARCH <- head(pred_vol_GARCH, 500)
  pred_vol_EGARCH <- head(pred_vol_EGARCH, 500)
  pred_vol_GAS <- head(pred_vol_GAS, 500)
  
  pred_RETURN_GARCH <- head(pred_RETURN_GARCH, 500)
  pred_RETURN_EGARCH <- head(pred_RETURN_EGARCH, 500)
  pred_RETURN_GAS <- head(pred_RETURN_GAS, 500)
  
  # Calculate actual volatility
  actual_return <- as.numeric(returns[1434:1933])
  
  train_mean <- mean(train_data)
  NORMreturns_train <- returns[(end_train+1) ] - train_mean
  actual_vol1 <- NORMreturns_train^2
  actual_vol <- c(actual_vol, actual_vol1)
  
  abs_vol_garch <- abs(actual_vol - pred_vol_GARCH)
  abs_vol_egarch <- abs(actual_vol - pred_vol_EGARCH)
  abs_vol_gas <- abs(actual_vol - pred_vol_GAS)
  
  abs_return_garch <- abs(actual_return - pred_RETURN_GARCH)
  abs_return_egarch <- abs(actual_return - pred_RETURN_EGARCH)
  abs_return_gas <- abs(actual_return - pred_RETURN_GAS)
  
  
  score_vol_garch <- 0
  score_vol_egarch <- 0
  score_vol_gas <- 0
  
  score_return_garch <- 0
  score_return_egarch <- 0
  score_return_gas <- 0
  
  for (i in 1:500) {
    min_vol <- min(abs_vol_garch[i], abs_vol_egarch[i], abs_vol_gas[i])
    min_return <- min(abs_return_garch[i], abs_return_egarch[i], abs_return_gas[i])
    
    if (abs_vol_garch[i] == min_vol) score_vol_garch <- score_vol_garch + 1
    if (abs_vol_egarch[i] == min_vol) score_vol_egarch <- score_vol_egarch + 1
    if (abs_vol_gas[i] == min_vol) score_vol_gas <- score_vol_gas + 1
    
    if (abs_return_garch[i] == min_return) score_return_garch <- score_return_garch + 1
    if (abs_return_egarch[i] == min_return) score_return_egarch <- score_return_egarch + 1
    if (abs_return_gas[i] == min_return) score_return_gas <- score_return_gas + 1
    win_rate_vol_garch <- score_vol_garch / 500
    win_rate_vol_egarch <- score_vol_egarch / 500
    win_rate_vol_gas <- score_vol_gas / 500
    
    win_rate_return_garch <- score_return_garch / 500
    win_rate_return_egarch <- score_return_egarch / 500
    win_rate_return_gas <- score_return_gas / 500
    
    cssfed_return_garch_egarch<-cumsum(abs_return_garch^2-abs_return_egarch^2)
    cssfed_return_gas_egarch<-cumsum(abs_return_gas^2-abs_return_egarch^2)
    cssfed_return_gas_garch<-cumsum(abs_return_gas^2-abs_return_garch^2)
    
    cssfed_vol_garch_egarch<-cumsum(abs_vol_garch^2-abs_vol_egarch^2)
    cssfed_vol_gas_egarch<-cumsum(abs_vol_gas^2-abs_vol_egarch^2)
    cssfed_vol_gas_garch<-cumsum(abs_vol_gas^2-abs_vol_garch^2)
    # محاسبه RMSE تجمعی برای 500 گام
    cum_rmse_return_garch <- numeric(500)  
    cum_rmse_return_egarch <- numeric(500) 
    cum_rmse_return_gas <- numeric(500)   
    
    cum_rmse_vol_garch <- numeric(500)  
    cum_rmse_vol_egarch <- numeric(500)
    cum_rmse_vol_gas <- numeric(500)  
    
    
    cum_rmse_return_garch <- sqrt(mean(abs_return_garch^2) )
    cum_rmse_return_egarch <- sqrt(mean(abs_return_egarch^2) )
    cum_rmse_return_gas <- sqrt(mean(abs_return_gas^2) )
    
    cum_rmse_vol_garch <- sqrt(mean(abs_vol_garch^2) )
    cum_rmse_vol_egarch <- sqrt(mean(abs_vol_egarch^2))
    cum_rmse_vol_gas <- sqrt(mean(abs_vol_gas^2) )
    
    
    error_return_garch <- actual_return - pred_RETURN_GARCH
    error_return_egarch <- actual_return - pred_RETURN_EGARCH
    error_return_gas <- actual_return - pred_RETURN_GAS
    
    error_vol_garch <- actual_vol - pred_vol_GARCH
    error_vol_egarch <- actual_vol - pred_vol_EGARCH
    error_vol_gas <- actual_vol - pred_vol_GAS
    
    dm_return_garch_egarch <- dm.test(error_return_garch, error_return_egarch, h = 1, power = 2)
    dm_return_garch_gas <- dm.test(error_return_garch, error_return_gas, h = 1, power = 2)
    dm_return_egarch_gas <- dm.test(error_return_egarch, error_return_gas, h = 1, power = 2)
    
    dm_vol_garch_egarch <- dm.test(error_vol_garch, error_vol_egarch, h = 1, power = 2)
    dm_vol_garch_gas <- dm.test(error_vol_garch, error_vol_gas, h = 1, power = 2)
    dm_vol_egarch_gas <- dm.test(error_vol_egarch, error_vol_gas, h = 1, power = 2)
    
  }
  list(
    dm_return_garch_egarch=dm_return_garch_egarch,dm_return_garch_gas=dm_return_garch_gas,dm_return_egarch_gas=dm_return_egarch_gas,
    dm_vol_garch_egarch=dm_vol_garch_egarch,dm_vol_garch_gas=dm_vol_garch_gas,dm_vol_egarch_gas=dm_vol_egarch_gas,
    cum_rmse_vol_garch=cum_rmse_vol_garch,cum_rmse_vol_egarch=cum_rmse_vol_egarch,cum_rmse_vol_gas=cum_rmse_vol_gas,
    cum_rmse_return_garch=cum_rmse_return_garch,cum_rmse_return_egarch=cum_rmse_return_egarch,cum_rmse_return_gas=cum_rmse_return_gas,
    cssfed_vol_garch_egarch=cssfed_vol_garch_egarch,cssfed_return_gas_egarch=cssfed_return_gas_egarch,cssfed_return_gas_garch=cssfed_return_gas_garch,
    cssfed_return_garch_egarch=cssfed_return_garch_egarch,cssfed_vol_gas_egarch=cssfed_vol_gas_egarch,cssfed_vol_gas_garch=cssfed_vol_gas_garch,
    win_rate_vol_gas=win_rate_vol_gas,win_rate_vol_egarch=win_rate_vol_egarch,win_rate_vol_garch=win_rate_vol_garch,
    win_rate_return_gas=win_rate_return_gas,win_rate_return_egarch=win_rate_return_egarch,win_rate_return_garch=win_rate_return_garch,
    pred_vol_GARCH = pred_vol_GARCH, pred_vol_EGARCH = pred_vol_EGARCH, pred_vol_GAS = pred_vol_GAS,
    pred_RETURN_GARCH = pred_RETURN_GARCH, pred_RETURN_EGARCH = pred_RETURN_EGARCH,
    pred_RETURN_GAS = pred_RETURN_GAS, actual_vol = actual_vol,
    abs_return_garch=abs_return_garch,abs_return_egarch=abs_return_egarch,abs_return_gas=abs_return_gas,
    abs_vol_gas=abs_vol_gas,abs_vol_egarch=abs_vol_egarch,abs_vol_garch=abs_vol_garch
  )
}

# Call the function for BRENT returns
results_BRENT <- forecast_and_evaluate(returns_BRENT, actual_returns = returns_BRENT)
results_WTI <- forecast_and_evaluate(returns_WTI,actual_returns = returns_WTI)
results_NG <- forecast_and_evaluate(returns_NG, actual_returns=returns_NG)
# Verify length of predicted returns
print(results_BRENT$win_rate_vol_garch)
print(results_BRENT$win_rate_vol_egarch)
print(results_BRENT$win_rate_vol_gas)

print(results_NG$win_rate_return_garch)
print(results_NG$win_rate_return_egarch)
print(results_NG$win_rate_return_gas)


print(results_BRENT$win_rate_return_garch)
print(results_BRENT$win_rate_return_egarch)
print(results_BRENT$win_rate_return_gas)

print(results_WTI$win_rate_return_garch)
print(results_WTI$win_rate_return_egarch)
print(results_WTI$win_rate_return_gas)

ylim_range <- range(results_BRENT$cum_rmse_return_garch,results_BRENT$cum_rmse_return_egarch,results_BRENT$cum_rmse_return_gas)

plot(results_BRENT$cum_rmse_return_garch, type = "l", col = "blue", lwd = 2, xlab = "Time", ylab = "CSSFED", main = "CSSFED for Different Models WTI 5-steps ", ylim = ylim_range)
lines(results_BRENT$cum_rmse_return_egarch, col = "red", lwd = 2, lty = 2)
lines(results_BRENT$pred_vol_GAS, col = "green", lwd = 2, lty = 3)
lines(results_BRENT$cum_rmse_return_gas, col = "purple", lwd = 2, lty = 3)
abline(h = 0, col = "black", lty = 4) # اضافه کردن خط صفر


results_WTI$win_rate_vol_garch
results_WTI$win_rate_vol_egarch
results_WTI$win_rate_vol_gas

results_BRENT$win_rate_vol_garch
results_BRENT$win_rate_vol_egarch
results_BRENT$win_rate_vol_gas

results_NG$win_rate_vol_garch
results_NG$win_rate_vol_egarch
results_NG$win_rate_vol_gas

results_BRENT$cum_rmse_vol_garch
results_BRENT$cum_rmse_vol_egarch
results_BRENT$cum_rmse_vol_gas

results_BRENT$cum_rmse_return_garch
results_BRENT$cum_rmse_return_egarch
results_BRENT$cum_rmse_return_gas

results_NG$cum_rmse_return_garch
results_NG$cum_rmse_return_egarch
results_NG$cum_rmse_return_gas

results_WTI$cum_rmse_return_garch
results_WTI$cum_rmse_return_egarch
results_WTI$cum_rmse_return_gas




results_NG$dm_return_garch_egarch
results_NG$dm_return_garch_gas
results_NG$dm_return_egarch_gas

results_BRENT$dm_return_garch_egarch
results_BRENT$dm_return_garch_gas
results_BRENT$dm_return_egarch_gas

results_WTI$dm_return_garch_egarch
results_WTI$dm_return_garch_gas
results_WTI$dm_return_egarch_gas


results_BRENT$dm_vol_garch_egarch
results_BRENT$dm_vol_garch_gas
results_BRENT$dm_vol_egarch_gas

results_WTI$dm_vol_garch_egarch
results_WTI$dm_vol_garch_gas
results_WTI$dm_vol_egarch_gas

results_NG$dm_vol_garch_egarch
results_NG$dm_vol_garch_gas
results_NG$dm_vol_egarch_gas


cssfed_GAS_GARCH <-results_WTI$cssfed_vol_gas_garch
cssfed_GAS_eGARCH <- results_WTI$cssfed_vol_gas_egarch
cssfed_GARCH_eGARCH <-results_WTI$cssfed_return_garch_egarch
ylim_range <-range(cssfed_GAS_GARCH,cssfed_GAS_eGARCH,cssfed_GARCH_eGARCH) 

plot(cssfed_GAS_GARCH, type = "l", col = "blue", lwd = 2, xlab = "Observation", ylab = "CSSFED", main = "CSSFED BRENT 1-Step return Forecasts" , ylim = ylim_range)
lines(cssfed_GAS_eGARCH, col = "red", lwd = 2, lty = 2)
lines(cssfed_GARCH_eGARCH, col = "green", lwd = 2, lty = 3)
abline(h = 0, col = "black", lty = 4) # اضافه کردن خط صفر
legend("topleft", legend = c("GAS-GARCH", "GAS-EGARCH", "GARCH-EGARCH"), col = c("blue", "red", "green"), lty = c(1, 2, 3), lwd = 2, cex = 0.43) # کوچک کردن بیشتر لگند

results_WTI$actual_vol

results_NG$pred_vol_GAS
results_WTI$pred_vol_GAS