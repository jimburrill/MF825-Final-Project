#WORKING IN PERCENTAGES and EXCESS RETURNS
#Part 1: Build required data Matrices---------------------------
#Read in Data for the USA
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
indus <- indus[indus[,1] >= 199301 & indus[,1] <= 202312,]
financial_uncert_data <- financial_uncert_data[financial_uncert_data[,1] >= 199301 & financial_uncert_data[,1] <= 202312,]
macro_uncert_data <- macro_uncert_data[macro_uncert_data[,1] >= 199301 & macro_uncert_data[,1] <= 202312,]
return20yr <- return20yr[return20yr[,1] >= 199301 & return20yr[,1] <= 202312,]

#Build Economic Factors
short_rate <- FF[, c(1, 5)]   #Inflation proxy
colnames(short_rate) <- c("Date", "shortRate")
term_struct <- return20yr
term_struct[,2] <- term_struct[,2] - FF[,5]     #20yr - MM
colnames(term_struct) <- c("Date", "termStruct")
rmrf <- FF[, c(1, 2)]   #Market excess return, DO WE NEED TO ADD RF BACK TO HERE????
macro_uncert <- macro_uncert_data[,1:2]
colnames(macro_uncert) <- c("Date", "macroUncer")
financial_uncert <- financial_uncert_data[,1:2]
colnames(financial_uncert) <- c("Date", "financialUncer")

#Industry Excess returns
xrindus <- as.matrix(indus)
xrindus[,2:50] <- xrindus[,2:50] - FF[,5]


#Part 2: Compute Factor Betas---------------------------
numindus <- 49
nyear <- 2023-1993 +1		

#Generate a list to store models, and a matrix to store industry betas over time
models <- list()  #List to store models
betas_indus <- matrix(NA, nrow = 5*(nyear-4), ncol = numindus+1)
alphas_indus <- matrix(NA, nrow = (nyear-4), ncol = numindus+1)
colnames <- colnames(xrindus[,])
colnames(betas_indus) <- colnames
colnames(alphas_indus) <- colnames
years <- rep(seq(1997, 1997 + floor((nrow(betas_indus) - 1) / 5)), each = 5)
betas_indus[,1] <- years
years <- rep(seq(1997, 1997 + floor((nrow(alphas_indus) - 1))), each = 1)
alphas_indus[,1] <- years

#Use 5 years of data to estimate the betas, estimate the betas every year
for (i in 1:(nyear-4)) {		
  current_yr <- 1996 + i
  first_yr <- current_yr - 4
  first_date <- (first_yr*100) + 01
  last_date <- (current_yr*100) + 12
  
  #Build 4 factor model for current time frame
  zret <- xrindus[xrindus[,1]>=first_date & xrindus[,1]<=last_date,2:50]
  xrm  <- as.matrix(rmrf[rmrf[,1]>=first_date & rmrf[,1]<=last_date,2])
  srate <- as.matrix(short_rate[short_rate[,1]>=first_date & short_rate[,1]<=last_date,2])
  tstruct <- as.matrix(term_struct[term_struct[,1]>=first_date & term_struct[,1]<=last_date,2])
  macuncert <- as.matrix(macro_uncert[macro_uncert[,1]>=first_date & macro_uncert[,1]<=last_date,2])
  finuncert <- as.matrix(financial_uncert[financial_uncert[,1]>=first_date & financial_uncert[,1]<=last_date,2])
  
  #Calculate row indices to store results in the Beta matrix
  start_row <- (i - 1) * 5 + 1
  end_row <- start_row + 5 - 1
  
  #Run regression for each industry, store the beta coefficients
  for (j in 1:numindus){
    models[[j]] <- lm(zret[,j] ~ xrm + srate + tstruct + macuncert + finuncert)
    betas_indus[start_row:end_row, j+1] <- coef(models[[j]])[2:6]  #Exclude intercept (alpha)
    alphas_indus[i, j+1] <- coef(models[[j]])[1:1]
  }
  
  #Print summary of each model if needed for diagnostics
  #for (k in 1:numindus){
  #  print(summary(models[[k]]))
  #}
}

#Part 3: Run the 2nd Pass Regression to Extract Risk Premiums---------------------------
library(sandwich)
numindus <- 49

#Run the 2nd pass regression with industry returns of 1998 to 2023
#Run the 2nd pass regression with Betas from 1997 to 2022
#Extract relevant data periods
xrindus_2pass <- xrindus[xrindus[,1] >= 199801 & xrindus[,1] <= 202312,]
betas_2pass <- betas_indus[betas_indus[,1] >= 1997 & betas_indus[,1] <= 2022,]

#Build the matrix to store risk premias for each factor
n_years <- length(betas_2pass[,1])/5
years <- seq(1998, 2023)
frp_rmrf <- matrix(NA, nrow = n_years, ncol = 5)
colnames(frp_rmrf) <- c('Date', 'Gamma1', 'HAC Var', 'GbyP', 'Precision')
frp_rmrf[,1] <- years
frp_srate <- matrix(NA, nrow = n_years, ncol = 5)
colnames(frp_srate) <- c('Date', 'Gamma2', 'HAC Var', 'GbyP', 'Precision')
frp_srate[,1] <- years
frp_tstruct <- matrix(NA, nrow = n_years, ncol = 5)
colnames(frp_tstruct) <- c('Date', 'Gamma3', 'HAC Var', 'GbyP', 'Precision')
frp_tstruct[,1] <- years
frp_macro <- matrix(NA, nrow = n_years, ncol = 5)
colnames(frp_macro) <- c('Date', 'Gamma4', 'HAC Var', 'GbyP', 'Precision')
frp_macro[,1] <- years
frp_fin <- matrix(NA, nrow = n_years, ncol = 5)
colnames(frp_fin) <- c('Date', 'Gamma5', 'HAC Var', 'GbyP', 'Precision')
frp_fin[,1] <- years
frp_alphas <- matrix(NA, nrow = n_years, ncol = 5)
colnames(frp_alphas) <- c('Date', 'Gamma0', 'HAC Var', 'GbyP', 'Precision')
frp_alphas[,1] <- years

#Build the matrix to store yearly betas for each factor for regression
col_names <- colnames(xrindus_2pass)
betas_rmrf <- matrix(NA, nrow = 49, ncol = 1)
colnames(betas_rmrf) <-  c("rmrf")
#betas_rmrf[,1] <- col_names[2:50]
betas_srate <- matrix(NA, nrow = 49, ncol = 1)
colnames(betas_srate) <- c("srate")
#betas_srate[,1] <- col_names[2:50]
betas_tsruct <- matrix(NA, nrow = 49, ncol = 1)
colnames(betas_tsruct) <- c("tstruct")
#betas_tsruct[,1] <- col_names[2:50]
betas_macro <- matrix(NA, nrow = 49, ncol = 1)
colnames(betas_macro) <- c("macro")
#betas_macro[,1] <- col_names[2:50]
betas_fin <- matrix(NA, nrow = 49, ncol = 1)
colnames(betas_fin) <- c("fin")
#betas_fin[,1] <- col_names[2:50]

#Run the 2nd pass regression of industry returns of t on the factor betas for t-1
for (i in 1:n_years) {
  first_date <- years[i]*100 + 01
  last_date <- years[i]*100 + 12
  
  xrindus_yr <- xrindus_2pass[xrindus_2pass[,1] >= first_date & xrindus_2pass[,1] <= last_date,]
  betas_yr <- betas_2pass[betas_2pass[,1] == (years[i]-1),]
  
  #Build individual factor beta matrices with each industry for the current year
  betas_rmrf[,1] <- betas_yr[1,2:50]
  betas_srate[,1] <- betas_yr[2,2:50]
  betas_tsruct[,1] <- betas_yr[3,2:50]
  betas_macro[,1] <- betas_yr[4,2:50]
  betas_fin[,1] <- betas_yr[5,2:50]
  
  #Build mean return of each industry matrix for the year
  average_indus_returns <- colMeans(xrindus_yr[,2:50])
  avg_returns_mat <- matrix(average_indus_returns, nrow = length(average_indus_returns), ncol = 1)
  
  #Run the 2nd pass regression for each year, store the beta coefficients
  models[[i]] <- lm(avg_returns_mat[,1] ~ betas_rmrf[,1] + betas_srate[,1] + betas_tsruct[,1] + betas_macro[,1] + betas_fin[,1])
  frp_alphas[i, 2] <- coef(models[[i]])[1]
  frp_rmrf[i, 2] <- coef(models[[i]])[2]
  frp_srate[i, 2] <- coef(models[[i]])[3]
  frp_tstruct[i, 2] <- coef(models[[i]])[4]
  frp_macro[i, 2] <- coef(models[[i]])[5]
  frp_fin[i, 2] <- coef(models[[i]])[6]
  
  #Get HAC standard errors^2 for the coefficients
  hac_se <- vcovHAC(models[[i]], type = "HC")
  hac_se_factors <- diag(hac_se)
  frp_alphas[i, 3] <- hac_se_factors[1]
  frp_rmrf[i, 3] <- hac_se_factors[2]
  frp_srate[i, 3] <- hac_se_factors[3]
  frp_tstruct[i, 3] <- hac_se_factors[4]
  frp_macro[i, 3] <- hac_se_factors[5]
  frp_fin[i, 3] <- hac_se_factors[6]
}

#Gamma Aggregation Numerator (Gamma*Precision = Gamma/Var)
frp_alphas[, 4] <- frp_alphas[, 2] / frp_alphas[, 3]
frp_rmrf[, 4] <- frp_rmrf[, 2] / frp_rmrf[, 3]
frp_srate[, 4] <- frp_srate[, 2] / frp_srate[, 3]
frp_tstruct[, 4] <- frp_tstruct[, 2] / frp_tstruct[, 3]
frp_macro[, 4] <- frp_macro[, 2] / frp_macro[, 3]
frp_fin[, 4] <- frp_fin[, 2] / frp_fin[, 3]

#Gamma Aggregation Denomenator (Precision = 1/Var)
frp_alphas[, 5] <- 1 / frp_alphas[, 3]
frp_rmrf[, 5] <- 1 / frp_rmrf[, 3]
frp_srate[, 5] <- 1 / frp_srate[, 3]
frp_tstruct[, 5] <- 1 / frp_tstruct[, 3]
frp_macro[, 5] <- 1 / frp_macro[, 3]
frp_fin[, 5] <- 1 / frp_fin[, 3]

#Aggregate the gammas
gamma0 <- sum(frp_alphas[, 4]) / sum(frp_alphas[, 5])
gamma1 <- sum(frp_rmrf[, 4]) / sum(frp_rmrf[, 5])
gamma2 <- sum(frp_srate[, 4]) / sum(frp_srate[, 5])
gamma3 <- sum(frp_tstruct[, 4]) / sum(frp_tstruct[, 5])
gamma4 <- sum(frp_macro[, 4]) / sum(frp_macro[, 5])
gamma5 <- sum(frp_fin[, 4]) / sum(frp_fin[, 5])



#Rolling Gammas Calculation
n_years <- length(betas_2pass[,1])/5 - 4
years_rolling <- seq(2002, 2023)
frp_rmrf_roll <- matrix(NA, nrow = n_years, ncol = 2)
colnames(frp_rmrf_roll) <- c('Date', 'Gamma1')
frp_rmrf_roll[,1] <- years_rolling
frp_srate_roll <- matrix(NA, nrow = n_years, ncol = 2)
colnames(frp_srate_roll) <- c('Date', 'Gamma1')
frp_srate_roll[,1] <- years_rolling
frp_tstruct_roll <- matrix(NA, nrow = n_years, ncol = 2)
colnames(frp_tstruct_roll) <- c('Date', 'Gamma1')
frp_tstruct_roll[,1] <- years_rolling
frp_macro_roll <- matrix(NA, nrow = n_years, ncol = 2)
colnames(frp_macro_roll) <- c('Date', 'Gamma1')
frp_macro_roll[,1] <- years_rolling
frp_fin_roll <- matrix(NA, nrow = n_years, ncol = 2)
colnames(frp_fin_roll) <- c('Date', 'Gamma1')
frp_fin_roll[,1] <- years_rolling
frp_alphas_roll <- matrix(NA, nrow = n_years, ncol = 2)
colnames(frp_alphas_roll) <- c('Date', 'Gamma1')
frp_alphas_roll[,1] <- years_rolling

for (i in 1:(n_years)){
  first_date <- years[i]
  last_date <- years[i+4]
  frp_alphas_roll[i,2] <- sum(frp_alphas[, 4][frp_alphas[, 1] >= first_date & frp_alphas[, 1] <= last_date]) / sum(frp_alphas[, 5][frp_alphas[, 1] >= first_date & frp_alphas[, 1] <= last_date])
  frp_rmrf_roll[i,2] <- sum(frp_rmrf[, 4][frp_rmrf[, 1] >= first_date & frp_rmrf[, 1] <= last_date]) / sum(frp_rmrf[, 5][frp_rmrf[, 1] >= first_date & frp_rmrf[, 1] <= last_date])
  frp_srate_roll[i,2] <- sum(frp_srate[, 4][frp_srate[, 1] >= first_date & frp_srate[, 1] <= last_date]) / sum(frp_srate[, 5][frp_srate[, 1] >= first_date & frp_srate[, 1] <= last_date])
  frp_tstruct_roll[i,2] <- sum(frp_tstruct[, 4][frp_tstruct[, 1] >= first_date & frp_tstruct[, 1] <= last_date]) / sum(frp_tstruct[, 5][frp_tstruct[, 1] >= first_date & frp_tstruct[, 1] <= last_date])
  frp_macro_roll[i,2] <- sum(frp_macro[, 4][frp_macro[, 1] >= first_date & frp_macro[, 1] <= last_date]) / sum(frp_macro[, 5][frp_macro[, 1] >= first_date & frp_macro[, 1] <= last_date])
  frp_fin_roll[i,2] <- sum(frp_fin[, 4][frp_fin[, 1] >= first_date & frp_fin[, 1] <= last_date]) / sum(frp_fin[, 5][frp_fin[, 1] >= first_date & frp_fin[, 1] <= last_date])
}

#Part 4: Examine Accuracy of Model---------------------------
#Read in the Industry Returns at a Yearly level
FF_model_per <- read.csv(file='FF-Yearly-USA.csv', header=TRUE)
FF_model_per <- FF_model_per[FF_model_per[,1] >= 1997 & FF_model_per[,1] <= 2023,]
indus_model_per <- read.csv(file='49-Industries-ValueWeight-USA-Yearly.csv', header=TRUE)
indus_model_per <- indus_model_per[indus_model_per[,1] >= 1997 & indus_model_per[,1] <= 2023,]
xrindus_model_per <- indus_model_per
xrindus_model_per[,2:50] <- indus_model_per[,2:50] - FF_model_per[,5]

col_names <- colnames(indus_model_per)
betas_rmrf <- matrix(NA, nrow = length(betas_indus[,1])/5, ncol = 50)
colnames(betas_rmrf) <-  col_names
betas_srate <- matrix(NA, nrow = length(betas_indus[,1])/5, ncol = 50)
colnames(betas_srate) <- col_names
betas_tsruct <- matrix(NA, nrow = length(betas_indus[,1])/5, ncol = 50)
colnames(betas_tsruct) <- col_names
betas_macro <- matrix(NA, nrow = length(betas_indus[,1])/5, ncol = 50)
colnames(betas_macro) <- col_names
betas_fin <- matrix(NA, nrow = length(betas_indus[,1])/5, ncol = 50)
colnames(betas_fin) <- col_names

betas_rmrf <- betas_indus[seq(1, length(betas_indus[,1]), by = 5), ]
betas_srate <- betas_indus[seq(2, length(betas_indus[,1]), by = 5), ]
betas_tsruct <- betas_indus[seq(3, length(betas_indus[,1]), by = 5), ]
betas_macro <- betas_indus[seq(4, length(betas_indus[,1]), by = 5), ]
betas_fin <- betas_indus[seq(5, length(betas_indus[,1]), by = 5), ]

#APT Expected Excess Returns Per Year
APTxr <- matrix(NA, nrow = length(xrindus_model_per[,1]), ncol = 50)
colnames(APTxr) <- col_names
APTxr[,1] <- xrindus_model_per[,1]
APTxr[,2:50] <- gamma1*betas_rmrf[,2:50] + gamma2*betas_srate[,2:50] + gamma3*betas_tsruct[,2:50] + gamma4*betas_macro[,2:50] + gamma5*betas_fin[,2:50]
APTxr[,2:50] <- APTxr[,2:50]*100    #Get in Percents to match FF

#Build matrixes to save R-squared values
rsquared_values <- matrix(NA, nrow = 1, ncol = 49)
colnames(rsquared_values) <-  col_names[2:50]
adjusted_rsquared <- matrix(NA, nrow = 1, ncol = 49)
colnames(adjusted_rsquared) <-  col_names[2:50]

#Plot each industry APT predicated return vs. Actual Observation
finalmodel <- list()  #List to store models
for (i in 2:50){
  plot(APTxr[,1], APTxr[,i], type = "l", col = 'black', ylim = c(-100, 130), 
       main = paste("Industy Excess Return 1997 to 2023:", col_names[i]), xlab = 'Date', ylab = 'Excess Return')
  lines(xrindus_model_per[,1], xrindus_model_per[,i], type = "l", col = 'blue')
  legend("topright", legend = c("APT", "Observed"), col = c("black", "blue"), lty = 1)
  
  #Regress Observed Excess Returns on APT Returns????
  finalmodel[[i]] <- lm(xrindus_model_per[,i] ~ APTxr[,i])
  model <- lm(xrindus_model_per[,i] ~ APTxr[,i])  #Run this to grab R-square values
  
  #Get R Squared values
  rsquared_values[,i-1] <- summary(model)$r.squared
  adjusted_rsquared[,i-1] <- summary(model)$adj.r.squared
}

#Part 5: Examine Accuracy of Model Rolling---------------------------
#Read in the Industry Returns at a Yearly level
#APT Expected Excess Returns Per Year Rolling Gammas
APTxr_roll <- matrix(NA, nrow = length(xrindus_model_per[,1])-5, ncol = 50)
colnames(APTxr_roll) <- col_names
APTxr_roll[,1] <- xrindus_model_per[6:27,1]
APTxr_roll[,2:50] <- frp_rmrf_roll[,2]*betas_rmrf[6:27,2:50] + frp_srate_roll[,2]*betas_srate[6:27,2:50] + frp_tstruct_roll[,2]*betas_tsruct[6:27,2:50] + frp_macro_roll[,2]*betas_macro[6:27,2:50] + frp_fin_roll[,2]*betas_fin[6:27,2:50]
APTxr_roll[,2:50] <- APTxr_roll[,2:50]*100    #Get in Percents to match FF

#Build matrixes to save R-squared values
rsquared_values <- matrix(NA, nrow = 1, ncol = 49)
colnames(rsquared_values) <-  col_names[2:50]
adjusted_rsquared <- matrix(NA, nrow = 1, ncol = 49)
colnames(adjusted_rsquared) <-  col_names[2:50]

#Plot each industry APT predicated return vs. Actual Observation
finalmodel <- list()  #List to store models
for (i in 2:50){
  plot(APTxr_roll[,1], APTxr_roll[,i], type = "l", col = 'black', ylim = c(-100, 130), 
       main = paste("Industy Excess Return 2002 to 2023:", col_names[i]), xlab = 'Date', ylab = 'Excess Return')
  lines(xrindus_model_per[6:27,1], xrindus_model_per[6:27,i], type = "l", col = 'blue')
  legend("topright", legend = c("APT", "Observed"), col = c("black", "blue"), lty = 1)
  
  #Regress Observed Excess Returns on APT Returns????
  finalmodel[[i]] <- lm(xrindus_model_per[6:27,i] ~ APTxr_roll[,i])
  model <- lm(xrindus_model_per[6:27,i] ~ APTxr_roll[,i])  #Run this to grab R-square values
  
  #Get R Squared values
  rsquared_values[,i-1] <- summary(model)$r.squared
  adjusted_rsquared[,i-1] <- summary(model)$adj.r.squared
}


#Part 6: Does the US Model Explain Regional Index Returns---------------------------
FF <- read.csv(file='FF-Monthly-USA.csv', header=TRUE)
FF <- FF[FF[,1] >= 200401 & FF[,1] <= 202312,]
regional_indexes <- read.csv(file='Country-Index-Returns-Monthly.csv', header=TRUE)
regional_indexes <- regional_indexes[1:240,1:7]
regional_indexes[,2:7] <- regional_indexes[,2:7]*100

#Compute Excess Returns
xrregional <- matrix(NA, nrow = length(regional_indexes[,1]), ncol = 7)
xrregional <- regional_indexes
xrregional[,2:7] <- regional_indexes[,2:7] - FF[,5]

numregions <- 6
num_years <- 20

#Generate a list to store models, and a matrix to store industry betas over time
regionmodels <- list()  #List to store models
colnames <- colnames(xrregional[,])
years <- seq(from = 2008, to = 2023, by = 1)
betas_rmrf <- matrix(NA, nrow = num_years-4, ncol = numregions+1)
colnames(betas_rmrf) <-  colnames
betas_rmrf[,1] <- years
betas_srate <- matrix(NA, nrow = num_years-4, ncol = numregions+1)
colnames(betas_srate) <- colnames
betas_srate[,1] <- years
betas_tsruct <- matrix(NA, nrow = num_years-4, ncol = numregions+1)
colnames(betas_tsruct) <- colnames
betas_tsruct[,1] <- years
betas_macro <- matrix(NA, nrow = num_years-4, ncol = numregions+1)
colnames(betas_macro) <- colnames
betas_macro[,1] <- years
betas_fin <- matrix(NA, nrow = num_years-4, ncol = numregions+1)
colnames(betas_fin) <- colnames
betas_fin[,1] <- years


#Use 5 years of data to estimate the betas, estimate the betas every year
for (i in 1:(num_years-4)) {		
  current_yr <- 2007 + i
  first_yr <- current_yr - 4
  first_date <- (first_yr*100) + 01
  last_date <- (current_yr*100) + 12
  
  #Build 4 factor model for current time frame
  zret <- xrregional[xrregional[,1]>=first_date & xrregional[,1]<=last_date,2:7]
  xrm  <- as.matrix(rmrf[rmrf[,1]>=first_date & rmrf[,1]<=last_date,2])
  srate <- as.matrix(short_rate[short_rate[,1]>=first_date & short_rate[,1]<=last_date,2])
  tstruct <- as.matrix(term_struct[term_struct[,1]>=first_date & term_struct[,1]<=last_date,2])
  macuncert <- as.matrix(macro_uncert[macro_uncert[,1]>=first_date & macro_uncert[,1]<=last_date,2])
  finuncert <- as.matrix(financial_uncert[financial_uncert[,1]>=first_date & financial_uncert[,1]<=last_date,2])

  #Run regression for each industry, store the beta coefficients
  for (j in 1:numregions){
    regionmodels[[j]] <- lm(zret[,j] ~ xrm + srate + tstruct + macuncert + finuncert)
    betas_rmrf[i,j+1] <- coef(regionmodels[[j]])[2:2] 
    betas_srate[i,j+1] <- coef(regionmodels[[j]])[3:3] 
    betas_tsruct[i,j+1] <- coef(regionmodels[[j]])[4:4] 
    betas_macro[i,j+1] <- coef(regionmodels[[j]])[5:5] 
    betas_fin[i,j+1] <- coef(regionmodels[[j]])[6:6] 
  }
}

#REDO THIS LOGIG, DONT THINK CORRECT
#Aggregate regional returns to yearly frequency
xrregion_yr <- matrix(NA, nrow = num_years-4, ncol = numregions+1)
colnames(xrregion_yr) <- colnames
xrregion_yr[,1] <- seq(from = 2008, to = 2023, by = 1)
for (i in 1:num_years-4){
  current_yr <- 2007 + i
  first_date <- (current_yr*100) + 01
  last_date <- (current_yr*100) + 12
  
  temp_returns <- xrregional[xrregional[,1]>=first_date & xrregional[,1]<=last_date,2:7]
  annual_returns <- sapply(temp_returns, function(x) prod(1 + (x/100)) - 1)
  xrregion_yr[i,2] <- annual_returns[1]*100
  xrregion_yr[i,3] <- annual_returns[2]*100
  xrregion_yr[i,4] <- annual_returns[3]*100
  xrregion_yr[i,5] <- annual_returns[4]*100
  xrregion_yr[i,6] <- annual_returns[5]*100
  xrregion_yr[i,7] <- annual_returns[6]*100
}

#APT Expected Excess Returns Per Year
APTxr <- matrix(NA, nrow = length(betas_srate[,1]), ncol = 7)
colnames(APTxr) <- colnames
APTxr[,1] <- betas_srate[,1]
APTxr[,2:7] <- gamma1*betas_rmrf[,2:7] + gamma2*betas_srate[,2:7] + gamma3*betas_tsruct[,2:7] + gamma4*betas_macro[,2:7] + gamma5*betas_fin[,2:7]
APTxr[,2:7] <- APTxr[,2:7]*100    #Get in Percents to match FF

#Build matrixes to save R-squared values
rsquared_values <- matrix(NA, nrow = 1, ncol = 6)
colnames(rsquared_values) <-  colnames[2:7]
adjusted_rsquared <- matrix(NA, nrow = 1, ncol = 6)
colnames(adjusted_rsquared) <-  colnames[2:7]

#Plot each industry APT predicated return vs. Actual Observation
finalmodel <- list()  #List to store models
for (i in 2:7){
  plot(APTxr[,1], APTxr[,i], type = "l", col = 'black', ylim = c(-100, 100), 
       main = paste("Region Excess Return 1997 to 2023:", colnames[i]), xlab = 'Date', ylab = 'Excess Return')
  lines(xrregion_yr[,1], xrregion_yr[,i], type = "l", col = 'blue')
  legend("topright", legend = c("APT", "Observed"), col = c("black", "blue"), lty = 1)
  
  #Regress Observed Excess Returns on APT Returns????
  finalmodel[[i]] <- lm(xrregion_yr[,i] ~ APTxr[,i])
  model <- lm(xrregion_yr[,i] ~ APTxr[,i])  #Run this to grab R-square values
  
  #Get R Squared values
  rsquared_values[,i-1] <- summary(model)$r.squared
  adjusted_rsquared[,i-1] <- summary(model)$adj.r.squared
}


# Set up the plotting layout
par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))  # 5 rows, 5 columns, adjust margin
# Iterate over each pair of observed and predicted returns
for (i in 2:7) {
  y_max <- max(c(APTxr[,i], xrregion_yr[,i])) * 1.1
  y_min <- min(c(APTxr[,i], xrregion_yr[,i])) * 1.1
  # Create the plot
  plot(xrregion_yr[,1], xrregion_yr[,i], type = "l", col = "blue", 
       main = paste(colnames[i]), xlab = "Date", ylab = "Excess Returns", ylim = c(y_min, y_max))
  
  # Add a reference line
  lines(APTxr[,1], APTxr[,i], col = "black")
}
# Reset plotting parameters
par(mfrow = c(1, 1))

#Part 7: Does the US Model Explain Regional Index Returns Rolling Gamma---------------------------
#APT Expected Excess Returns Per Year
APTxr_roll <- matrix(NA, nrow = length(betas_srate[,1]), ncol = 7)
colnames(APTxr_roll) <- colnames
APTxr_roll[,1] <- betas_srate[,1]
APTxr_roll[,2:7] <- frp_rmrf_roll[7:22,2]*betas_rmrf[,2:7] + frp_srate_roll[7:22,2]*betas_srate[,2:7] + frp_tstruct_roll[7:22,2]*betas_tsruct[,2:7] + frp_macro_roll[7:22,2]*betas_macro[,2:7] + frp_fin_roll[7:22,2]*betas_fin[,2:7]
APTxr_roll[,2:7] <- APTxr_roll[,2:7]*100    #Get in Percents to match FF

#Build matrixes to save R-squared values
rsquared_values <- matrix(NA, nrow = 1, ncol = 6)
colnames(rsquared_values) <-  colnames[2:7]
adjusted_rsquared <- matrix(NA, nrow = 1, ncol = 6)
colnames(adjusted_rsquared) <-  colnames[2:7]

#Plot each industry APT predicated return vs. Actual Observation
finalmodel <- list()  #List to store models
for (i in 2:7){
  plot(APTxr_roll[,1], APTxr_roll[,i], type = "l", col = 'black', ylim = c(-100, 100), 
       main = paste("Region Excess Return 1997 to 2023:", colnames[i]), xlab = 'Date', ylab = 'Excess Return')
  lines(xrregion_yr[,1], xrregion_yr[,i], type = "l", col = 'blue')
  legend("topright", legend = c("APT", "Observed"), col = c("black", "blue"), lty = 1)
  
  #Regress Observed Excess Returns on APT Returns????
  finalmodel[[i]] <- lm(xrregion_yr[,i] ~ APTxr_roll[,i])
  model <- lm(xrregion_yr[,i] ~ APTxr_roll[,i])  #Run this to grab R-square values
  
  #Get R Squared values
  rsquared_values[,i-1] <- summary(model)$r.squared
  adjusted_rsquared[,i-1] <- summary(model)$adj.r.squared
}


# Set up the plotting layout
par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))  # 5 rows, 5 columns, adjust margin
# Iterate over each pair of observed and predicted returns
for (i in 2:7) {
  y_max <- max(c(APTxr_roll[,i], xrregion_yr[,i])) * 1.1
  y_min <- min(c(APTxr_roll[,i], xrregion_yr[,i])) * 1.1
  # Create the plot
  plot(xrregion_yr[,1], xrregion_yr[,i], type = "l", col = "blue", 
       main = paste(colnames[i]), xlab = "Date", ylab = "Excess Returns", ylim = c(y_min, y_max))
  
  # Add a reference line
  lines(APTxr_roll[,1], APTxr_roll[,i], col = "black")
}
# Reset plotting parameters
par(mfrow = c(1, 1))

#Part 8: Does the Model Explain Momentum Portfolio Returns---------------------------
FF <- read.csv(file='FF-Monthly-USA.csv', header=TRUE)
FF <- FF[FF[,1] >= 199301 & FF[,1] <= 202312,]
momentum_ports <- read.csv(file='25-Size-Momentum_Ports.csv', header=TRUE)
momentum_ports <- momentum_ports[momentum_ports[,1] >= 199301 & momentum_ports[,1] <= 202312,]

#Compute Excess Returns
xrmomentum <- momentum_ports
xrmomentum[,2:26] <- xrmomentum[,2:26] - FF[,5]

numports <- 25
num_years <- 2023-1993 +1

#Generate a list to store models, and a matrix to store industry betas over time
momentummodels <- list()  #List to store models
colnames <- colnames(xrmomentum[,])
years <- seq(from = 1997, to = 2023, by = 1)
betas_rmrf <- matrix(NA, nrow = num_years-4, ncol = numports+1)
colnames(betas_rmrf) <-  colnames
betas_rmrf[,1] <- years
betas_srate <- matrix(NA, nrow = num_years-4, ncol = numports+1)
colnames(betas_srate) <- colnames
betas_srate[,1] <- years
betas_tsruct <- matrix(NA, nrow = num_years-4, ncol = numports+1)
colnames(betas_tsruct) <- colnames
betas_tsruct[,1] <- years
betas_macro <- matrix(NA, nrow = num_years-4, ncol = numports+1)
colnames(betas_macro) <- colnames
betas_macro[,1] <- years
betas_fin <- matrix(NA, nrow = num_years-4, ncol = numports+1)
colnames(betas_fin) <- colnames
betas_fin[,1] <- years


#Use 5 years of data to estimate the betas, estimate the betas every year
for (i in 1:(num_years-4)) {		
  current_yr <- 1996 + i
  first_yr <- current_yr - 4
  first_date <- (first_yr*100) + 01
  last_date <- (current_yr*100) + 12
  
  #Build 4 factor model for current time frame
  zret <- xrmomentum[xrmomentum[,1]>=first_date & xrmomentum[,1]<=last_date,2:26]
  xrm  <- as.matrix(rmrf[rmrf[,1]>=first_date & rmrf[,1]<=last_date,2])
  srate <- as.matrix(short_rate[short_rate[,1]>=first_date & short_rate[,1]<=last_date,2])
  tstruct <- as.matrix(term_struct[term_struct[,1]>=first_date & term_struct[,1]<=last_date,2])
  macuncert <- as.matrix(macro_uncert[macro_uncert[,1]>=first_date & macro_uncert[,1]<=last_date,2])
  finuncert <- as.matrix(financial_uncert[financial_uncert[,1]>=first_date & financial_uncert[,1]<=last_date,2])
  
  #Run regression for each industry, store the beta coefficients
  for (j in 1:numports){
    momentummodels[[j]] <- lm(zret[,j] ~ xrm + srate + tstruct + macuncert + finuncert)
    betas_rmrf[i,j+1] <- coef(momentummodels[[j]])[2:2] 
    betas_srate[i,j+1] <- coef(momentummodels[[j]])[3:3] 
    betas_tsruct[i,j+1] <- coef(momentummodels[[j]])[4:4] 
    betas_macro[i,j+1] <- coef(momentummodels[[j]])[5:5] 
    betas_fin[i,j+1] <- coef(momentummodels[[j]])[6:6] 
  }
}


#Momentum portfolio returns at yearly frequency
xrmomentum_yr <- read.csv(file='25-Size-Momentum_Ports_Annual.csv', header=TRUE)
xrmomentum_yr <- xrmomentum_yr[xrmomentum_yr[,1] >= 1997 & xrmomentum_yr[,1] <= 2023,]
xrmomentum_yr[,2:26] <- xrmomentum_yr[,2:26] - FF_model_per[,5]

#APT Expected Excess Returns Per Year
APTxr <- matrix(NA, nrow = length(betas_srate[,1]), ncol = 26)
colnames(APTxr) <- colnames
APTxr[,1] <- betas_srate[,1]
APTxr[,2:26] <- gamma1*betas_rmrf[,2:26] + gamma2*betas_srate[,2:26] + gamma3*betas_tsruct[,2:26] + gamma4*betas_macro[,2:26] + gamma5*betas_fin[,2:26]
APTxr[,2:26] <- APTxr[,2:26]*100    #Get in Percents to match FF

#Build matrixes to save R-squared values
rsquared_values <- matrix(NA, nrow = 1, ncol = 25)
colnames(rsquared_values) <-  colnames[2:26]
adjusted_rsquared <- matrix(NA, nrow = 1, ncol = 25)
colnames(adjusted_rsquared) <-  colnames[2:26]

#Plot each industry APT predicated return vs. Actual Observation
finalmodel <- list()  #List to store models
for (i in 2:26){
  plot(APTxr[,1], APTxr[,i], type = "l", col = 'black', ylim = c(-100, 120), 
       main = paste("Momentum Port Excess Return 1997 to 2023:", colnames[i]), xlab = 'Date', ylab = 'Excess Return')
  lines(xrmomentum_yr[,1], xrmomentum_yr[,i], type = "l", col = 'blue')
  legend("topright", legend = c("APT", "Observed"), col = c("black", "blue"), lty = 1)
  
  #Regress Observed Excess Returns on APT Returns????
  finalmodel[[i]] <- lm(xrmomentum_yr[,i] ~ APTxr[,i])
  model <- lm(xrmomentum_yr[,i] ~ APTxr[,i])  #Run this to grab R-square values
  
  #Get R Squared values
  rsquared_values[,i-1] <- summary(model)$r.squared
  adjusted_rsquared[,i-1] <- summary(model)$adj.r.squared
}



# Set up the plotting layout
par(mfrow = c(5, 5), mar = c(3, 3, 1, 1))  # 5 rows, 5 columns, adjust margin

# Iterate over each pair of observed and predicted returns
for (i in 2:26) {
  # Create the plot
  plot(xrmomentum_yr[,1], xrmomentum_yr[,i], type = "l", col = "blue", 
       main = paste(colnames[i]), xlab = "Date", ylab = "Excess Returns")
  
  # Add a reference line
  lines(APTxr[,1], APTxr[,i], col = "black")
}

# Reset plotting parameters
par(mfrow = c(1, 1))



#Part 9: Does the Model Explain Momentum Portfolio Returns Rolling Gamma---------------------------
#APT Expected Excess Returns Per Year
APTxr_roll <- matrix(NA, nrow = length(betas_srate[6:27,1]), ncol = 26)
colnames(APTxr_roll) <- colnames
APTxr_roll[,1] <- betas_srate[6:27,1]
APTxr_roll[,2:26] <- frp_rmrf_roll[,2]*betas_rmrf[6:27,2:26] + frp_srate_roll[,2]*betas_srate[6:27,2:26] + frp_tstruct_roll[,2]*betas_tsruct[6:27,2:26] + frp_macro_roll[,2]*betas_macro[6:27,2:26] + frp_fin_roll[,2]*betas_fin[6:27,2:26]
APTxr_roll[,2:26] <- APTxr_roll[,2:26]*100    #Get in Percents to match FF

#Build matrixes to save R-squared values
rsquared_values <- matrix(NA, nrow = 1, ncol = 25)
colnames(rsquared_values) <-  colnames[2:26]
adjusted_rsquared <- matrix(NA, nrow = 1, ncol = 25)
colnames(adjusted_rsquared) <-  colnames[2:26]

#Plot each industry APT predicated return vs. Actual Observation
finalmodel <- list()  #List to store models
for (i in 2:26){
  #Regress Observed Excess Returns on APT Returns????
  finalmodel[[i]] <- lm(xrmomentum_yr[6:27,i] ~ APTxr_roll[,i])
  model <- lm(xrmomentum_yr[6:27,i] ~ APTxr_roll[,i])  #Run this to grab R-square values
  
  #Get R Squared values
  rsquared_values[,i-1] <- summary(model)$r.squared
  adjusted_rsquared[,i-1] <- summary(model)$adj.r.squared
}

# Set up the plotting layout
par(mfrow = c(5, 5), mar = c(3, 3, 1, 1))  # 5 rows, 5 columns, adjust margin

# Iterate over each pair of observed and predicted returns
for (i in 2:26) {
  # Create the plot
  plot(xrmomentum_yr[6:27,1], xrmomentum_yr[6:27,i], type = "l", col = "blue", 
       main = paste(colnames[i]), xlab = "Date", ylab = "Excess Returns")
  
  # Add a reference line
  lines(APTxr[6:27,1], APTxr[6:27,i], col = "black")
  lines(APTxr_roll[,1], APTxr_roll[,i], col = "red")
}

# Reset plotting parameters
par(mfrow = c(1, 1))


#Calculate asset performance using 5 year rolling gamma for Regional and Momentum