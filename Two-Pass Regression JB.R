#WORKING IN PERCENTAGES and EXCESS RETURNS
#Part 1: Build required data Matrices---------------------------
#Read in Data for the USA
library(lmtest)
setwd("C:/Users/jackb/MF825/MF825-Final-Project")
FF <- read.csv(file='FF-Monthly-USA.csv', header=TRUE)
indus <- read.csv(file='49-Industries-ValueWeight-USA-Monthly.csv', header=TRUE)
return20yr <- read.csv('Bond-20y-return.csv', header=TRUE)
macro_uncert_data <- read.csv(file='MacroUncertaintyToCirculate.csv', header=TRUE)
financial_uncert_data <- read.csv(file='FinancialUncertaintyToCirculate.csv', header=TRUE)

#Drop NA rows from 20yr bond return data
return20yr <- return20yr[complete.cases(return20yr), ]

#Standardize date formats to FF monthly format
financial_uncert_data$Date <- sub("-", "", financial_uncert_data$Date)
financial_uncert_data$Date <- sub("-", "", financial_uncert_data$Date)
year <- substr(financial_uncert_data$Date, 1, 4)  # Extract characters representing year
month <- substr(financial_uncert_data$Date, 5, 6)  # Extract characters representing month
financial_uncert_data$Date <- paste0(year, month)
macro_uncert_data$Date <- sub("-", "", macro_uncert_data$Date)
macro_uncert_data$Date <- sub("-", "", macro_uncert_data$Date)
year <- substr(macro_uncert_data$Date, 1, 4)  # Extract characters representing year
month <- substr(macro_uncert_data$Date, 5, 6)  # Extract characters representing month
macro_uncert_data$Date <- paste0(year, month)
year <- substr(return20yr[,1], 1, 4)  # Extract characters representing year
month <- substr(return20yr[,1], 5, 6)  # Extract characters representing month
return20yr[,1] <- paste0(year, month)


#Filter Date ranges to 01/1993 to 12/2023
FF <- FF[FF[,1] >= 199301 & FF[,1] <= 202312,]
FF$mkt <- FF$RF+FF$Mkt.RF
indus <- indus[indus[,1] >= 199301 & indus[,1] <= 202312,]
financial_uncert_data <- financial_uncert_data[financial_uncert_data[,1] >= 199301 & financial_uncert_data[,1] <= 202312,]
macro_uncert_data <- macro_uncert_data[macro_uncert_data[,1] >= 199301 & macro_uncert_data[,1] <= 202312,]
return20yr <- return20yr[return20yr[,1] >= 199301 & return20yr[,1] <= 202312,]
head(FF)
#Build Economic Factors
#SHOULD WE SUBTRACT RF FROM ALL FACTORS???????
short_rate <- FF[, c(1, 5)]   #Inflation proxy
colnames(short_rate) <- c("Date", "shortRate")
term_struct <- return20yr
term_struct[,2] <- term_struct[,2] - FF[,5]     #20yr - MM
colnames(term_struct) <- c("Date", "termStruct")
rmrf <- FF[, c(1, 2)]   #Market excess return
mkt <- FF[, c(6)]   #Market return
macro_uncert <- macro_uncert_data[,1:2]
colnames(macro_uncert) <- c("Date", "macroUncer")
financial_uncert <- financial_uncert_data[,1:2]
colnames(financial_uncert) <- c("Date", "financialUncer")



#Industry Excess returns
xrindus <- as.matrix(indus)
xrindus[,2:50] <- xrindus[,2:50] - FF[,5]

n_factors = 5
n_months <- nrow(xrindus)
period_lag <- 59 #1/2010 -> 12/2014
t <-n_months-period_lag-1 #Minus one because we need that last Rt+1
industry_rets <- xrindus[,2:50]


second_pass_gammas <- list()
for (j in 1:49){
  temp_df <- list()
  for (i in 1:t) {
    Rt_i <-industry_rets[i:(i+period_lag),j]
    mkt_ret_t <- rmrf$Mkt.RF[i:(i+period_lag)]
    rf_t <- short_rate$shortRate[i:(i+period_lag)]
    yield_slope_t <- term_struct$termStruct[i:(i+period_lag)]
    macro_uncert_t <- macro_uncert$macroUncer[i:(i+period_lag)]
    financial_uncert_t <- financial_uncert$financialUncer[i:(i+period_lag)]
    
    reg_i <- lm(Rt_i ~ mkt_ret_t + rf_t + yield_slope_t + macro_uncert_t + financial_uncert_t)
    Bt_i <- as.vector(reg_i$coefficients)[2:(n_factors+1)]
    
    #Rt+1 is for the forward return
    Rt_ij_1 <- industry_rets[(i+period_lag)+1,j]
    temp_df[[i]] <- c(Rt_ij_1, Bt_i)
  }
  
  second_pass_matrix_i <- matrix(unlist(temp_df), ncol = 6, byrow = TRUE)
  second_pass_reg <- lm(second_pass_matrix_i[, 1] ~ second_pass_matrix_i[, -1])
  second_pass_gammas[[j]] <- as.vector(second_pass_reg$coefficients)
  
}
# each index is industry gammas; second_pass_gammas[1] -> Agric Gammas
second_pass_gammas





