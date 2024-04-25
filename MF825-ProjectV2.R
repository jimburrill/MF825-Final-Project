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
#SHOULD WE SUBTRACT RF FROM ALL FACTORS???????
short_rate <- FF[, c(1, 5)]   #Inflation proxy
colnames(short_rate) <- c("Date", "shortRate")
term_struct <- return20yr
term_struct[,2] <- term_struct[,2] - FF[,5]     #20yr - MM
colnames(term_struct) <- c("Date", "termStruct")
rmrf <- FF[, c(1, 2)]   #Market excess return
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
numindus <- 49

#Run the 2nd pass regression with industry returns of 1998 to 2023
#Run the 2nd pass regression with Betas from 1997 to 2022
#Extract relevant data periods
xrindus_2pass <- xrindus[xrindus[,1] >= 199801 & xrindus[,1] <= 202312,]
betas_2pass <- betas_indus[betas_indus[,1] >= 1997 & betas_indus[,1] <= 2022,]

#Build the matrix to store risk premias for each factor
n_years <- length(betas_2pass[,1])/5
years <- seq(1998, 2023)
frp_rmrf <- matrix(NA, nrow = n_years, ncol = 50)
colnames(frp_rmrf) <- colnames(xrindus_2pass)
frp_rmrf[,1] <- years
frp_srate <- matrix(NA, nrow = n_years, ncol = 50)
colnames(frp_srate) <- colnames(xrindus_2pass)
frp_srate[,1] <- years
frp_tstruct <- matrix(NA, nrow = n_years, ncol = 50)
colnames(frp_tstruct) <- colnames(xrindus_2pass)
frp_tstruct[,1] <- years
frp_macro <- matrix(NA, nrow = n_years, ncol = 50)
colnames(frp_macro) <- colnames(xrindus_2pass)
frp_macro[,1] <- years
frp_fin <- matrix(NA, nrow = n_years, ncol = 50)
colnames(frp_fin) <- colnames(xrindus_2pass)
frp_fin[,1] <- years
frp_alphas <- matrix(NA, nrow = n_years, ncol = 50)
colnames(frp_alphas) <- colnames(xrindus_2pass)
frp_alphas[,1] <- years

#Build the matrix to store yearly betas (12 times) for each factor for regression
col_names <- colnames(xrindus_2pass)
betas_rmrf <- matrix(NA, nrow = 12, ncol = length(xrindus_2pass[1,])-1)
colnames(betas_rmrf) <- col_names[2:50]
betas_srate <- matrix(NA, nrow = 12, ncol = length(xrindus_2pass[1,])-1)
colnames(betas_srate) <- col_names[2:50]
betas_tsruct <- matrix(NA, nrow = 12, ncol = length(xrindus_2pass[1,])-1)
colnames(betas_tsruct) <- col_names[2:50]
betas_macro <- matrix(NA, nrow = 12, ncol = length(xrindus_2pass[1,])-1)
colnames(betas_macro) <- col_names[2:50]
betas_fin <- matrix(NA, nrow = 12, ncol = length(xrindus_2pass[1,])-1)
colnames(betas_fin) <- col_names[2:50]

#Run the 2nd pass regression of industry returns of t on the factor betas for t-1
for (i in 1:n_years) {
  first_date <- years[i]*100 + 01
  last_date <- years[i]*100 + 12
  
  xrindus_yr <- xrindus_2pass[xrindus_2pass[,1] >= first_date & xrindus_2pass[,1] <= last_date,]
  betas_yr <- betas_2pass[betas_2pass[,1] == (years[i]-1),]

  for (j in 1:12){
    betas_rmrf[j,1:49] <- betas_yr[1,2:50]
    betas_srate[j,1:49] <- betas_yr[2,2:50]
    betas_tsruct[j,1:49] <- betas_yr[3,2:50]
    betas_macro[j,1:49] <- betas_yr[4,2:50]
    betas_fin[j,1:49] <- betas_yr[5,2:50]
  }

  #Run the 2nd pass regression for each industry, store the beta coefficients
  for (k in 1:numindus){
    models[[k]] <- lm(xrindus_yr[,k+1] ~ betas_rmrf[,k] + betas_srate[,k] + betas_tsruct[,k] + betas_macro[,k] + betas_fin[,k])
    frp_alphas[i, k+1] <- coef(models[[k]])[1]
    frp_rmrf[i, k+1] <- coef(models[[k]])[2]
    frp_srate[i, k+1] <- coef(models[[k]])[3]
    frp_tstruct[i, k+1] <- coef(models[[k]])[4]
    frp_macro[i, k+1] <- coef(models[[k]])[5]
    frp_fin[i, k+1] <- coef(models[[k]])[6]
  }
}