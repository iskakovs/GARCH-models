# Install packages
install.packages("rugarch")
install.packages("forecast")
install.packages("lubridate")
install.packages("xts")
install.packages("urca")
install.packages("tseries")
install.packages("ggplot2")
install.packages("fUnitRoots")
install.packages("lmtest")
install.packages("rio")
install.packages("car")
install.packages("timeSeries")
install.packages("tidyverse")
install.packages("reshape")
install.packages("nutshell")
install.packages("gridExtra")
install.packages("Hmisc")

#Load libraries
library(rugarch)
library(forecast)
library(lubridate)
library(xts)
library(urca)
library(tseries)
library(ggplot2)
library(fUnitRoots)
library(lmtest)
library(rio)
library(car)
library(timeSeries)
library(tidyverse)
library(fpp3)
library(patchwork)
library(reshape)
library(gridExtra)
library(Hmisc)

#Load the data 
magnit = read.csv("C:\\Users\\Admin\\Desktop\\MGNT_180101_220223.csv", header = TRUE)

#Select variables we manipulate
data <- data.frame(magnit$X.DATE., magnit$X.CLOSE.)

#Select variables we manipulate
data <- data.frame(magnit$X.DATE., magnit$X.CLOSE.)

# Rename columns
colnames(data) <- c('date', 'p') 

# Insert trading days column
data = mutate(data, trading_day = row_number())
data = as_tsibble(data, index = trading_day)

# Plot the graph and since we have time series ACF and PACF
gg_tsdisplay(data, p, plot_type = 'partial')

#Check for stationarity (ADF-test)
summary(ur.df(data$p, type = "drift", selectlags = "AIC"))
#Value of t-statistics = -2.302 > -3.43 = Critical value 
# => H0 hypothesis of non-stationarity is failed to be rejected

#Take log of the prices
r = diff(log(data$p))

#Plot the graph
plot(r, type = "l", xlab = "")
# Looks pretty stationary

#Check for stationarity again for transformed series(ADF-test)
summary(ur.df(r, type = "drift", selectlags = "AIC"))
#Value of t-statistics = -18.2858 < -3.43 = Critical value 
# => H0 hypothesis of non-stationarity od the prices is rejected

#ACF and PACF estimated values
acf(r, plot = FALSE)
pacf(r, plot = FALSE)

#Visualize ACF and PACF
acf(r) #MA(2)
pacf(r) #AR(1) or may be AR(2)

#Check if the residuals are in the form of random noise
series <- arima.sim(list(order = c(0,0,0)), n = 1000)
grid.arrange(ggAcf(series), ggPacf(series), nrow = 2)

# Model 1: ARIMA (1, 0, 2)
model1 <- arima(r, order = c(1, 0, 2))
summary(model1)

# Model 2: ARIMA (2, 0, 2)
model2 <- arima(r, order = c(2, 0, 2))
summary(model2)

# Model 3: ARIMA (2, 0, 3)
model3 <- arima(r, order = c(2, 0, 3))
summary(model3)

# Model 4: ARIMA (2, 0, 3)
model4 <- arima(r, order = c(0, 0, 0))
summary(model4)

# CONCLUSION: The smallest AIC is for the  ARIMA (2, 0, 2), so we conclude that we have ARMA(2,2) model

#Check for residuals
Box.test(model2$residuals, type = "Ljung-Box")
# So, from Ljung-Box test we accepted the null hypothesis that the residuals are white noise

#Save the residuals of the model
e = model2$residuals

#Plot the residuals
ggtsdisplay(e)

#Check for GARCH-effect
#Generate lagged variables
e_1 <- Lag(e, -1)
e_2 <- Lag(e, -2)

#Create the model from square of the residuals and check for significance 
g_model <- lm(e^2 ~ e_1^2 + e_2^2)
summary(g_model)

# All our coefficients are significant, and p value is 2.726e-16. 
# So we can detect the GARCH-effect

# Check for GARCH(1,1) model
require(rugarch)
GARCHspec11 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)))
M_garch11 = ugarchfit(GARCHspec11, data = r)
M_garch11
# All our coefficients are significant at 5% level, but not at 1%.
# Let's check other models - GARCH(1,2) and GARCH(2,1)

# Check for GARCH(1,2) model
GARCHspec12 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,2)))
M_garch12 = ugarchfit(GARCHspec12, data = r)
M_garch12

# Check for GARCH(2,1) model
GARCHspec21 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2,1)))
M_garch21 = ugarchfit(GARCHspec21, data = r)
M_garch21

# So the optimal model is GARCH(1,1)

# Since both additions of sigma_2 and gamma_2 turned out to be insignificant, 
# the GARCH(1,1) model will also be preferable to the GARCH(2,2) model

#Let's check for volatilities of the GARCH(1,1) model
# Estimate the volatility
vol = M_garch11@fit$sigma

#Annualize the volatilities
vol_ann = vol*(250)^(1/2)*100

#Plot the graph of the volatilities
plot(vol_ann, type = "l", xlab = "")

#Plot the graph of the volatilities
plot(vol_ann, type = "l", xlab = "")

# Forecast the volatilities for T+1 and T+2 periods
# Forecast for 2 periods
volforec = ugarchforecast(M_garch11, r, n.ahead = 2) 
volforec

# Convert the volatility to annual volatility
vol_ann_forec = volforec@forecast$sigmaFor*(250)^(1/2)*100
vol_ann_forec
#So we have a high volatility 222.47% and 219.85%

#Create the GARCH-t(1,1) model
tGARCHspec11 <- ugarchspec(variance.model = list(model = "sGARCH", 
                                            garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0,0), 
                                            include.mean = TRUE), distribution.model = "std")
t_garch11 = ugarchfit(tGARCHspec11, data = r)
t_garch11

# According to parameters estimation this model is better than the simple GARCH(1,1) model
# So, for now tGARCH(1,1) model shows better coefficients estimates (see Pr(>|t|) - all coefficients are significant)
# All our coefficients are significant at 1% level and the Akaike information criterion is better (lower)

#Create the EGARCH(1,1) model and check for lambda asymmetry coefficient 
eGARCHspec11 <- ugarchspec(variance.model = list(model = "eGARCH", 
                          garchOrder = c(1,1)),
                          distribution.model = "norm")
e_garch11 = ugarchfit(eGARCHspec11, data = r)
e_garch11

# asymmetry coefficient is represented by aplha1 which is < 0, but insignificant
# the negativeness says that positive shock has less effect, than the negative shocks.

#Finally, let's compare the models
M_garch11
t_garch11
e_garch11

#We will use AIC
# For GARCH(1,1)  : -5.2395
# For GRACH-t     : -5.3745
# For EGRACH(1,1) : -5.2481

# The smallest value for AIC is in GRACH-t, so we conclude that this model is the best one

# Now let's calculate the volatility for this model:
vol_t = t_garch11@fit$sigma

#Annualize the volatilities
vol_t_ann = vol_t*(250)^(1/2)*100

#Plot the graph of the volatilities
plot(vol_t_ann, type = "l", xlab = "")

# And forecast the volatilities for GARCH-t model
# Forecast the volatilities for T+1 and T+2 periods
# Forecast for 2 periods
volforec_t = ugarchforecast(e_garch11, r, n.ahead = 2) 
volforec_t

# Convert the volatility to annual volatility
volt_ann_forec = volforec_t@forecast$sigmaFor*(250)^(1/2)*100
volt_ann_forec
#Still we have a high volatility for period T+1: 221.09% and period T+2: 200.81%
