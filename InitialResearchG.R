require(rugarch)
require(xts)
require(PerformanceAnalytics)
require(tidyquant)
require(GAS)
require(ggplot2)
setwd("~/")
#Preliminary Data preparation
SPX <- tq_get("^GSPC",from = "1990-01-01",end = Sys.Date())
SPX$return <- c(0, diff(log(SPX$adjusted)))
SPX <- SPX[SPX$date>"2015-01-01",]
GOLD <- read.csv("Documents/Bollege/FINE 460/Project/Data/GOLD.csv")
GOLD$return <- c(0,diff(log(GOLD$close)))
GOLD$date <- as.Date(GOLD$date)
ETH <- read.csv(DATA)
BTC <- read.csv(DATA)
ADA <- read.csv(DATA)
BTC_Metrics <- read.csv(DATA)
ETH_Metrics <- read.csv(DATA)
ADA_Metrics <- read.csv(DATA)

coins <- list(BTC, ETH, ADA)
coinNames <- list("BTC", "ETH", "ADA")
names(coins) <- coinNames
metrics <- list(BTC_Metrics, ETH_Metrics, ADA_Metrics)
names(metrics) <- coinNames

i <- 0
# Calculate log returns, sort dataframes, correct date format, calculate returns squared
for(coin in coins){
  i<- i+1
  coins[[i]] <- coins[[i]][coins[[i]]$close!=0,]
  coins[[i]]$return <- c(0, diff(log(coins[[i]]$close)))
  #coins[[i]]$return <- coins[[i]]$return-mean(coins[[i]]$return)
  coins[[i]] <- coins[[i]][order(coins[[i]]$date),]
  coins[[i]]$date <- as.Date(coins[[i]]$date)
  metrics[[i]]$date <- as.Date(metrics[[i]]$date)
  coins[[i]] <- merge(coins[[i]], metrics[[i]], by="date")
  coins[[i]]$returnsq <- (coins[[i]]$return)^2
}
#End of mandatory set up, rest is optional

# Make sure that each coin has same dates
dates_intersect <- intersect(coins[[1]]$date, coins[[2]]$date)
for(coin in coins){
  dates_intersect <- intersect(coin$date, dates_intersect)
}
i <- 0
for(coin in coins){
  i <- i+1
  coins[[i]] <- coin[coin$date %in% dates_intersect,]
}
i <- 0
for(metric in metrics){
  i <- i+1
  metrics[[i]] <- metric[metric$date %in% intersect(metric$date, coins[[i]]$date),]
}

#Functions to load into session memory
#Opmimizing functions fo rfindign orders with lowest lags
optimizesGARCH <- function(i, maxOrd=2){
  maxAIC = 1000
  maxModel = ugarchfit(spec = ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(0,0), include.mean = T), distribution.model = "norm"), data=xts(coins[[i]]$return, order.by = coins[[i]]$date))
  for(gp in 1:maxOrd){
    for(gq in 1:maxOrd){
      for(ap in 0:maxOrd){
        for(aq in 0:maxOrd){
          curSpec = ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(gp,gq)),mean.model = list(armaOrder=c(ap,aq), include.mean = T), distribution.model = "norm")
          model = ugarchfit(spec = curSpec, data=xts(coins[[i]]$return, order.by = coins[[i]]$date))
          AIC = 2*(ap+aq+gp+gq)-2*model@fit$LLH
          if(AIC<maxAIC){
            maxAIC=AIC
            maxSpec = curSpec
            maxModel=model
          }
        }      
      }
    }
  }
  results = list(maxModel, maxSpec, xts(coins[[i]]$return, order.by = coins[[i]]$date))
  names(results) =list("model", "spec", "data")
  return(results)
}
optimizegjrGARCH <- function(i, maxOrd=2){
  maxAIC = 1000
  maxModel = ugarchfit(spec = ugarchspec(variance.model = list(model="gjrGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(0,0), include.mean = T), distribution.model = "norm"), data=xts(coins[[i]]$return, order.by = coins[[i]]$date))
  for(gp in 1:maxOrd){
    for(gq in 1:maxOrd){
      for(ap in 0:maxOrd){
        for(aq in 0:maxOrd){
          curSpec = ugarchspec(variance.model = list(model="gjrGARCH",garchOrder=c(gp,gq)),mean.model = list(armaOrder=c(ap,aq), include.mean = T), distribution.model = "norm")
          model = ugarchfit(spec = curSpec, data=xts(coins[[i]]$return, order.by = coins[[i]]$date))
          AIC = 2*(ap+aq+gp+gq)-2*model@fit$LLH
          if(AIC<maxAIC){
            maxAIC=AIC
            maxSpec = curSpec
            maxModel=model
          }
        }      
      }
    }
  }
  results = list(maxModel, maxSpec, xts(coins[[i]]$return, order.by = coins[[i]]$date))
  names(results) =list("model", "spec", "data")
  return(results)
}

optimizeeGARCHx <- function(i, maxOrd=2){
  maxAIC = 1000
  maxModel = ugarchfit(spec = ugarchspec(variance.model = list(model="eGARCH",garchOrder=c(1,1), external.regressors=matrix(c(0,diff(log(coins[[i]]$X1dayactive/coins[[i]]$X7dayactive))))),mean.model = list(armaOrder=c(0,0), include.mean = T), distribution.model = "norm"), data=xts(coins[[i]]$return, order.by = coins[[i]]$date))
  for(gp in 1:maxOrd){
    for(gq in 1:maxOrd){
      for(ap in 0:maxOrd){
        for(aq in 0:maxOrd){
          curspec=ugarchspec(variance.model = list(model="eGARCH",garchOrder=c(gp,gq), external.regressors=matrix(c(0,diff(log(coins[[i]]$X1dayactive/coins[[i]]$X7dayactive))))),mean.model = list(armaOrder=c(ap,aq), include.mean = T), distribution.model = "norm")
          model = ugarchfit(spec = curspec, data=xts(coins[[i]]$return, order.by = coins[[i]]$date))
          AIC = 2*(ap+aq+2*gp+gq+1)-2*model@fit$LLH
          if(AIC<maxAIC){
            maxAIC=AIC
            maxModel=model
            maxSpec=curspec
          }
        }      
      }
    }
  }
  results = list(maxModel, maxSpec, xts(coins[[i]]$return, order.by = coins[[i]]$date))
  names(results) =list("model", "spec", "data")
  return(results)
}
#Diagnostics on forecasting
forecastStats <- function(aSpec, aData, fore=F){
  print(ugarchfit(spec=aSpec, data=aData, out.sample = 500))
  forecast = ugarchroll(aSpec, aData)
  if(fore){
    return(forecast)
  }
  print(report(forecast, VaR.alpha=.05))
  fpm = fpm(forecast)
  fpm[1,1] = sqrt(fpm[1,1])
  fpm[4,1] = mean(fpm(forecast, summary = F)[,5])
  fpm[5,1] = RMSE_R2(forecast)[[1]]
  print(fpm)
  return(fpm)
}
RMSE_R2 <- function(model){ 
  RMSE = sqrt(mean(
    (tail(model@model$data^2, 500) - as.numeric(model@forecast$density[,2])^2)^2
  ))
R2 = 1 - sum((tail(model@model$data^2, 500) - as.numeric(model@forecast$density[,2])^2)^2)/
  sum((tail(model@model$data^2, 500) - mean(tail(model@model$data^2, 500)))^2)
return(c(RMSE, R2)) 
}
RMSE_rolling <- function(data){
  RMSE=sqrt(mean(data^2-rollapply(data, 30, var)^2)^2)
  return(RMSE)
}
logRatioTest <- function(umodel, rmodel){
  LR = 2*(umodel@fit$LLH-rmodel@fit$LLH)
  return(pchisq(q = LR, df=1))
}
#VaR backtesting functions
indvVaR <- function(x){
  return(VaR(x, method="historical"))
}
#VaR Benchmark results
historicalVaR <- function(i){
  data = coins[[i]][coins[[i]]$date>coins[[i]]$date[[length(coins[[i]]$date)-530]],]
  date = index(xts(data$return, order.by = data$date)[29,])
  return(BacktestVaR(xts(data$return, order.by = data$date)[data$date>date],rollapply(data$return, 30, indvVaR), .05))
}
#Function to backtest hedged portolio performance
backtestPortfolio <- function(model, l_forecast){
  totalReturnman=1
  managedSeries = rep()
  unManagedSeries = rep()
  totalReturn=1
  l_training=length(model$model@model$modeldata$data)-l_forecast
  fit = ugarchfit(model$spec, model$data[1:l_training])
  sigma = ugarchforecast(fitORspec = fit, n.ahead=1)@forecast$sigmaFor[[1]]
  for(i in l_training:length(model$model@model$modeldata$data)){
    if(i %%25==0){
      done = (i-500)/(length(model$model@model$modeldata$data)-500)*100
      print(paste(done, "%"))
    }
    fit=ugarchfit(model$spec, model$data[1:i],solver = "hybrid")
    fore = ugarchforecast(fit)
    sigma=fore@forecast$sigmaFor[[1]]
    #If predicted volatility is not above historical average, do not hedge
    if(sigma<sd(model$data[1:l_training])){
      sigma=0
    }
    #Ensures no negative hedging
    if((1-5*sigma)<0){
      sigma=.2
    }
    totalReturnman=totalReturnman*(1+model$model@model$modeldata$data[i])*(1-5*sigma)+5*sigma*totalReturnman
    managedSeries <- append(managedSeries,totalReturnman)
    totalReturn=totalReturn*(1+model$model@model$modeldata$data[i])
    unManagedSeries <- append(unManagedSeries,totalReturn)
  }
  print(totalReturnman/totalReturn)
  return(list(managedSeries, unManagedSeries))
}


#Create Benchmarks
BTC_bench <- historicalVaR(1)
ETH_bench <- historicalVaR(2)
ADA_bench <- historicalVaR(3)
#Create Models
BTC_garch <- optimizesGARCH(1)
ETH_garch <- optimizesGARCH(2)
ADA_garch <- optimizesGARCH(3)
BTC_egarchx <- optimizeeGARCHx(1)
ETH_egarchx <- optimizeeGARCHx(2)
ADA_egarchx <- optimizeeGARCHx(3)
BTC_igarch <- createiGARCH(1)
ETH_igarch <- createiGARCH(2)
ADA_igarch <- createiGARCH(3)
BTC_gjrgarch <- optimizegjrGARCH(1)
ETH_gjrgarch <- optimizegjrGARCH(2)
ADA_gjrgarch <- optimizegjrGARCH(3)
#Test Models
BTC_garchTest <- forecastStats(BTC_garch$spec, BTC_garch$data)
ETH_garchTest <- forecastStats(ETH_garch$spec, ETH_garch$data)
ADA_garchTest <- forecastStats(ADA_garch$spec, ADA_garch$data)
BTC_egarchxTest <- forecastStats(BTC_egarchx$spec, BTC_egarchx$data)
ETH_egarchxTest <- forecastStats(ETH_egarchx$spec, ETH_egarchx$data)
ADA_egarchxTest <- forecastStats(ADA_egarchx$spec, ADA_egarchx$data)
BTC_iGarchTest <- forecastStats(BTC_igarch$spec, BTC_igarch$data)
ETH_iGarchTest <- forecastStats(ETH_igarch$spec, ETH_igarch$data)
ADA_iGarchTest <- forecastStats(ADA_igarch$spec, ADA_igarch$data)
BTC_gjrGarchTest <- forecastStats(BTC_gjrgarch$spec, BTC_gjrgarch$data)
ETH_gjrGarchTest <- forecastStats(ETH_gjrgarch$spec, ETH_gjrgarch$data)
ADA_gjrGarchTest <- forecastStats(ADA_gjrgarch$spec, ADA_gjrgarch$data)
BTC_egarchTest <- forecastStats(ugarchspec(variance.model = list(model="eGARCH",garchOrder=c(2,2)),mean.model = list(armaOrder=c(1,1), include.mean = T), distribution.model = "norm"), xts(coins[[1]]$return, order.by = coins[[1]]$date))
ETH_egarchTest <- forecastStats(ugarchspec(variance.model = list(model="eGARCH",garchOrder=c(2,2)),mean.model = list(armaOrder=c(1,1), include.mean = T), distribution.model = "norm"), xts(coins[[2]]$return, order.by = coins[[2]]$date))
ADA_egarchTest <- forecastStats(ugarchspec(variance.model = list(model="eGARCH",garchOrder=c(2,2)),mean.model = list(armaOrder=c(1,1), include.mean = T), distribution.model = "norm"), xts(coins[[3]]$return, order.by = coins[[3]]$date))
#Backtest Portfolio
BTC_garch_backtest <- backtestPortfolio(BTC_garch, 500)
BTC_igarch_backtest <- backtestPortfolio(BTC_igarch, 500)
BTC_gjrgarch_backtest <- backtestPortfolio(BTC_gjrgarch, 500)
BTC_egarchx_backtest <- backtestPortfolio(BTC_egarchx, 500)
ETH_garch_backtest <- backtestPortfolio(ETH_garch, 500)
ETH_igarch_backtest <- backtestPortfolio(ETH_igarch, 500)
ETH_gjrgarch_backtest <- backtestPortfolio(ETH_gjrgarch, 500)
ETH_egarchx_backtest <- backtestPortfolio(ETH_egarchx, 500)
ADA_garch_backtest <- backtestPortfolio(ADA_garch, 500)
ADA_igarch_backtest <- backtestPortfolio(ADA_igarch, 500)
ADA_gjrgarch_backtest <- backtestPortfolio(ADA_gjrgarch, 500)
ADA_egarchx_backtest <- backtestPortfolio(ADA_egarchx, 500)

#Plotting Functions
backtestGraph <- function(seriesList){
  l_training=length(coins$BTC$date)-length(seriesList)
  ggplot()+ geom_line(aes(x=tail(coins$BTC$date, length(seriesList[[1]])), y=as.numeric(seriesList[[1]]), color="Managed"))+
    geom_line(aes(x=tail(coins$BTC$date, length(seriesList[[1]])), y=as.numeric(seriesList[[2]]), color="Unmanaged"))+
    labs(color="Legend")+ylab("Total Return")+xlab("Date")+ggtitle("Total Return of $1")
}
createGarchPlot <- function(model, model2){
  ggplot()+
    geom_line(aes(model2@model$modeldata$index, y=model2@fit$sigma, color="IGARCH"))+
    geom_line(aes(model@model$modeldata$index, y=model@fit$sigma, color="GARCH"))+
    theme_bw()+ylab("ETH Volatility")+xlab("Date")+ggtitle("Daily Conditional Volatility")+labs(color="Legend")+
    theme(legend.position = c(.99, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))
}

createoverlayGarchPlot <- function(model, model2){
  ggplot()+geom_line(aes(x=model@model$modeldata$index, y=model@model$modeldata$data, color="Returns"))+
    geom_line(aes(model2@model$modeldata$index, y=2*model2@fit$sigma, color="eGARCHX"))+
    geom_line(aes(model2@model$modeldata$index, y=-2*model2@fit$sigma, color="eGARCHX"))+
    geom_line(aes(model@model$modeldata$index, y=2*model@fit$sigma, color="sGARCH"))+
    geom_line(aes(model@model$modeldata$index, y=-2*model@fit$sigma, color="sGARCH"))+
    theme_bw()+ylab("ADA log returns")+xlab("Date")+ggtitle("Daily Returns with 2 SD")+scale_color_manual(values=c('blue','#000000', 'red'))+labs(color="Legend")+
    theme(legend.position = c(.99, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))
}





#Appendix
#This includes code for initial research that was done before choosing the topic of GARCH models, still informative about behavior of crypto currencies. Additionally, it has code that create figures found in initial project report.
require(tidyquant)
require(tidyverse)
require(forecast)
require(tseries)
require(fGarch)
require(rmgarch)

#Create daily returns graphs
dev.off()
for(i in 1:length(coins)){
  jpeg(file=paste0("Documents/Bollege/FINE 460/Project/Figures/",coinNames[[i]], "Daily Returns.jpeg"))
  plot(as.Date(coins[[i]]$date), coins[[i]]$return, type = "l",  xlab="Date", ylab="Returns")
  title(coinNames[[i]])
  dev.off()
}

createDailyReturns <- function(i){
  plot(as.Date(coins[[i]]$date), coins[[i]]$return, type = "l",  xlab="Date", ylab="Returns")
  title(i)
}

#PCA analysis
returns <- data.frame(sapply(coins, function(x){c(x$return)}))
coinPCA <- princomp(returns)

#PLOT
#plot loadings
par(mfrow=c(1,2))
barplot(coinPCA$loadings[,1], names.arg = coinNames)
title("1st Factor Loadings")

barplot(coinPCA$loadings[,2], names.arg = coinNames)
title("2nd Factor Loadings")
#plot factor importance
dev.off()
plot(100*coinPCA$sdev[1:5]^2/sum(coinPCA$sdev^2), ylab="% variance explained", xlab = "Factor", pch =21,
     bg = "skyblue", cex = 1.2, panel.first= grid(col = "grey60", lwd = 1.5), type = "b", main = "Factor Importance")
par(new = T)
plot(cumsum(100*coinPCA$sdev[1:5]^2/sum(coinPCA$sdev^2)), xaxt = "n", yaxt="n", type = "b",
     pch = 24, bg = "orange", ylim = c(0,100), xlab = "", ylab = "") 
axis(side = 4)
legend("right", c("individual (left)", "cumulated (right)"),
       pch = c(21,24), pt.bg = c("skyblue", "orange"), bty = "n")   

#ACF of different cryptos
createACFPACF <- function(i){
  par(mfrow=c(1,2))
  ACF <- Acf(coins[[i]]$return)
  PACF <- Pacf(coins[[i]]$return)
  plot(ACF$lag[-1], ACF$acf[-1], type="p", pch=21, bg="blue", xlab = "Lag", ylab="ACF")
  segments(x0=ACF$lag[-1], y0=0, x1=ACF$lag[-1], y1=ACF$acf[-1], col=c("blue"))
  abline(0,0)
  abline(qnorm(1 - 0.05 / 2) / sqrt(length(coins[[1]]$return)), 0, lty="dashed")
  abline(-qnorm(1 - 0.05 / 2) / sqrt(length(coins[[1]]$return)), 0, lty="dashed")
  title(coinNames[[i]])
  plot(PACF$lag, PACF$acf, type="p", pch=21, bg="blue", xlab = "Lag", ylab="PACF")
  segments(x0=PACF$lag, y0=0, x1=PACF$lag, y1=PACF$acf, col=c("blue"))
  abline(0,0)
  abline(qnorm(1 - 0.05 / 2) / sqrt(length(coins[[1]]$return)), 0, lty="dashed")
  abline(-qnorm(1 - 0.05 / 2) / sqrt(length(coins[[1]]$return)), 0, lty="dashed")
  title(coinNames[[i]])
} 
createACFPACF(1)

# CAPM Model
reg_results <- list()
reg_coefficients <- list()
reg_betas <- vector()
for(i in 1:(length(coins)-1))
{
  ind_reg <- lm(coins[[i]]$return ~ coins[[i]]$Volume.USD)
  reg_results[[i]] <- ind_reg
  reg_coefficients[[i]] <- ind_reg$coefficients
  reg_betas[[i]] <- ind_reg$coefficients[[2]]
}

barplot(reg_betas, names.arg = c("XRP", "BCH", "ETH", "LTC", "BTC"))

garchspec <-  ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(0,0), include.mean = T), distribution.model = "norm")
egarchspec <-  ugarchspec(variance.model = list(model="eGARCH",garchOrder=c(1,1), external.regressors = matrix(coins$BTC$X1dayactive/coins$BTC$X7dayactive)),mean.model = list(armaOrder=c(0,0), include.mean = T), distribution.model = "norm")
BTC_garchspec <-  ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(0,1), include.mean = T), distribution.model = "norm")
ETH_garchspec <-  ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(2,0), include.mean = T), distribution.model = "norm")
ADA_garchspec <-  ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,2), include.mean = T), distribution.model = "norm")
BTC_garch1 <- ugarchfit(spec = garchspec, data = (coins[["BTC"]]$return))
ugarchBTC_garch2 <- ugarchfit(spec = BTC_garchspec, data = (coins[["BTC"]]$return))
ETH_garch1 <- ugarchfit(spec = garchspec, data = coins[["ETH"]]$return)
ETH_garch2 <- ugarchfit(spec = garchspec, data = coins[["ETH"]]$return)
ADA_garch1 <- ugarchfit(spec = garchspec, data = (coins[["ADA"]]$return))
ADA_garch2 <- ugarchfit(spec = garchspec, data = (coins[["ADA"]]$return))
BTC_egarchx <- ugarchfit(spec=egarchspec, data = coins[["BTC"]]$return)

SPX_garch <- ugarchfit(spec = garchspec, data = SPX$returns)