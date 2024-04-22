#WORKING IN PERCENTAGES and EXCESS RETURNS
#Part 1: Build required data Matrices---------------------------
#Read in Data for the USA
FF <- read.csv(file='FF-Monthly-USA.csv', header=TRUE)
indus <- read.csv(file='49-Industries-ValueWeight-USA-Monthly.csv', header=TRUE)
yield10yr <- read.csv(file='10yr-Yield-USA-Monthly.csv', header=TRUE)
corp_yield <- read.csv(file='Corporate-InvestGradeIndex-Yield-Daily-Temp.csv', header=TRUE)

#Format Date Columns to be uniform
corp_yield$Date <- sub("-", "", corp_yield$Date)
year <- substr(corp_yield$Date, 3, 6)  # Extract characters representing year
month <- substr(corp_yield$Date, 1, 2)  # Extract characters representing month
corp_yield$Date <- paste0(year, month)
year_month <- substr(yield10yr$DATE, 1, 7)
yield10yr$DATE <- gsub("-", "", year_month)

#Filter Date ranges to 01/1996 to 12/2023 (Adjust idealy to 1993 once get bond yield data)
corp_yield <- corp_yield[corp_yield[,1] >= 199701 & corp_yield[,1] <= 202312,]
yield10yr <- yield10yr[yield10yr[,1] >= 199701 & yield10yr[,1] <= 202312,]
indus <- indus[indus[,1] >= 199701 & indus[,1] <= 202312,]
FF <- FF[FF[,1] >= 199701 & FF[,1] <= 202312,]

#Build Economic Factors
short_rate <- FF[, c(1, 5)]
colnames(short_rate) <- c("Date", "shortRate")
term_struct <- yield10yr
term_struct[,2] <- term_struct[,2] - FF[,5]
colnames(term_struct) <- c("Date", "termStruct")
default_prem <- corp_yield
default_prem[,2] <- default_prem[,2] - yield10yr[,2]
colnames(default_prem) <- c("Date", "defaultPrem")
rmrf <- FF[, c(1, 2)]

#Industry Excess returns
xrindus <- as.matrix(indus)
xrindus[,2:50] <- xrindus[,2:50] - FF[,5]


#Part 2: Compute Factor Betas---------------------------
numindus <- 49
nyear <- 2023-1997 +1		

#Generate a list to store models, and a matrix to store industry betas over time
models <- list()  #List to store models
betas_indus <- matrix(NA, nrow = 4*(nyear-4), ncol = numindus+1)
alphas_indus <- matrix(NA, nrow = (nyear-4), ncol = numindus+1)
colnames <- colnames(xrindus[,])
colnames(betas_indus) <- colnames
colnames(alphas_indus) <- colnames
years <- rep(seq(2001, 2001 + floor((nrow(betas_indus) - 1) / 4)), each = 4)
betas_indus[,1] <- years
years <- rep(seq(2001, 2001 + floor((nrow(alphas_indus) - 1))), each = 1)
alphas_indus[,1] <- years

#Use 5 years of data to estimate the betas, estimate the betas every year
for (i in 1:(nyear-4)) {		
  current_yr <- 2000 + i
  first_yr <- current_yr - 4
  first_date <- (first_yr*100) + 01
  last_date <- (current_yr*100) + 12
  
  #Build 4 factor model for current time frame
  zret <- xrindus[xrindus[,1]>=first_date & xrindus[,1]<=last_date,2:50]
  xrm  <- as.matrix(rmrf[rmrf[,1]>=first_date & rmrf[,1]<=last_date,2])
  srate <- as.matrix(short_rate[short_rate[,1]>=first_date & short_rate[,1]<=last_date,2])
  tstruct <- as.matrix(term_struct[term_struct[,1]>=first_date & term_struct[,1]<=last_date,2])
  defprem <- as.matrix(default_prem[default_prem[,1]>=first_date & default_prem[,1]<=last_date,2])
  
  #Calculate row indices to store results in the Beta matrix
  start_row <- (i - 1) * 4 + 1
  end_row <- start_row + 4 - 1
  
  #Run regression for each industry, store the beta coefficients
  for (j in 1:numindus){
    models[[j]] <- lm(zret[,j] ~ xrm + srate + tstruct + defprem)
    betas_indus[start_row:end_row, j+1] <- coef(models[[j]])[2:(5)]  #Exclude intercept (alpha)
    alphas_indus[i, j+1] <- coef(models[[j]])[1:1]
  }
  
  #Print summary of each model if needed for diagnostics
  #for (k in 1:numindus){
  #  print(summary(models[[k]]))
  #}
}



#Part 3: Construct Mimicking Portfolios, Extract Risk Premiums---------------------------
xrindus_mimick <- as.matrix(xrindus[xrindus[,1] >= 200101 & xrindus[,1] <= 202312,])
mimick_port_rtns <- matrix(NA, nrow = length(xrindus_mimick[,1]), ncol = 5)
colnames(mimick_port_rtns) <- c('X', 'RmRf', 'Srate', 'Tstruct', 'Dprem')
mimick_port_rtns[,1] <- xrindus_mimick[,1]

#Calculate mimicking port where weight of Brmrf=1, Bsrate=Btstruct=Bdprem=0
Betas_rmrf_yr <- betas_indus[seq(1, nrow(betas_indus), by = 4), ]
Betas_rmrf_mon <- matrix(rep(Betas_rmrf_yr, each = 12), nrow = nrow(Betas_rmrf_yr) * 12)
temp <- xrindus_mimick[,2:50] * Betas_rmrf_mon[,2:50]
mimick_port_rtns[,2] <- rowSums(temp) 
mimick_port_rtns[,2] <- mimick_port_rtns[,2] / sum(abs(mimick_port_rtns[,2]))   #Normalize

#Calculate mimicking port where weight of Bsrate=1, Brmrf=Btstruct=Bdprem=0
Betas_srate_yr <- betas_indus[seq(2, nrow(betas_indus), by = 4), ]
Betas_srate_mon <- matrix(rep(Betas_srate_yr, each = 12), nrow = nrow(Betas_srate_yr) * 12)
temp <- xrindus_mimick[,2:50] * Betas_srate_mon[,2:50]
mimick_port_rtns[,3] <- rowSums(temp) 
mimick_port_rtns[,3] <- mimick_port_rtns[,3] / sum(abs(mimick_port_rtns[,3]))   #Normalize

#Calculate mimicking port where weight of Btstruct=1, Brmrf=Bsrate=Bdprem=0
Betas_tstruct_yr <- betas_indus[seq(3, nrow(betas_indus), by = 4), ]
Betas_tstruct_mon <- matrix(rep(Betas_tstruct_yr, each = 12), nrow = nrow(Betas_tstruct_yr) * 12)
temp <- xrindus_mimick[,2:50] * Betas_tstruct_mon[,2:50]
mimick_port_rtns[,4] <- rowSums(temp) 
mimick_port_rtns[,4] <- mimick_port_rtns[,4] / sum(abs(mimick_port_rtns[,4]))   #Normalize

#Calculate mimicking port where weight of Bdprem=1, Brmrf=Bsrate=Btstruct=0
Betas_dprem_yr <- betas_indus[seq(4, nrow(betas_indus), by = 4), ]
Betas_dprem_mon <- matrix(rep(Betas_dprem_yr, each = 12), nrow = nrow(Betas_dprem_yr) * 12)
temp <- xrindus_mimick[,2:50] * Betas_dprem_mon[,2:50]
mimick_port_rtns[,5] <- rowSums(temp) 
mimick_port_rtns[,5] <- mimick_port_rtns[,5] / sum(abs(mimick_port_rtns[,5]))   #Normalize

#Calculate factor risk premia on yearly basis
factor_premia <- matrix(NA, nrow = length(Betas_dprem_yr[,1]), ncol = 5)
colnames(factor_premia) <- c('X', 'RmRf', 'Srate', 'Tstruct', 'Dprem')
factor_premia[,1] <- Betas_dprem_yr[,1]

#Extract the year from the monthly dates
years <- substr(as.character(mimick_port_rtns[, 1]), 1, 4)
#Aggregate the monthly premiums into yearly premiums
#ENSURE DOING THe AGGREGATION PROPERLY 
factor_premia[,2] <- tapply(mimick_port_rtns[, 2], INDEX = years, FUN = function(x) prod(1 + sum(x)/100) - 1)*100
factor_premia[,3] <- tapply(mimick_port_rtns[, 3], INDEX = years, FUN = function(x) prod(1 + sum(x)/100) - 1)*100
factor_premia[,4] <- tapply(mimick_port_rtns[, 4], INDEX = years, FUN = function(x) prod(1 + sum(x)/100) - 1)*100
factor_premia[,5] <- tapply(mimick_port_rtns[, 5], INDEX = years, FUN = function(x) prod(1 + sum(x)/100) - 1)*100


#Part 4: Examine APT vs. Actual Returns---------------------------
#Build required matrix for APT
APTxrtn <- matrix(NA, nrow = length(factor_premia[,1]), ncol = (numindus+1))
APTxrtn[,1] <- factor_premia[,1]
colnames <- colnames(xrindus[,])
colnames(APTxrtn) <- colnames

#Create an index vector for the rows of corrsponding factor betas
rmrf_rows <- seq(from = 1, to = length(betas_indus[,1]), by = 4)
srate_rows <- seq(from = 2, to = length(betas_indus[,1]), by = 4)
tstruct_rows <- seq(from = 3, to = length(betas_indus[,1]), by = 4)
dprem_rows <- seq(from = 4, to = length(betas_indus[,1]), by = 4)

#Calculate APT excess returns
APTxrtn[,2:50] <- alphas_indus[,2:50] + 
                  factor_premia[,2]*betas_indus[rmrf_rows,2:50] + 
                  factor_premia[,3]*betas_indus[srate_rows,2:50] + 
                  factor_premia[,4]*betas_indus[tstruct_rows,2:50] + 
                  factor_premia[,5]*betas_indus[dprem_rows,2:50]

#Extract the year from the monthly dates
filtered_xrindus <- xrindus[xrindus[,1] >= 200101,]
xrindus_yrly <- matrix(NA, nrow = length(APTxrtn[,1]), ncol = (numindus+1))
xrindus_yrly[,1] <- APTxrtn[,1]
colnames <- colnames(filtered_xrindus[,])
colnames(xrindus_yrly) <- colnames
years <- substr(as.character(filtered_xrindus[, 1]), 1, 4)
#Aggregate the monthly premiums into yearly premiums
#ENSURE DOING THE AGGREGATION PROPERLY 
for (i in 2:50){
  xrindus_yrly[,i] <- tapply(filtered_xrindus[,i], INDEX = years, FUN = function(x) prod(1 + x/100) - 1)
}
xrindus_yrly[,2:50] <- xrindus_yrly[,2:50]*100

#Create an empty matrix to store regression coefficients
betas <- matrix(NA, nrow = 49, ncol = 1)
alphas <- matrix(NA, nrow = 49, ncol = 1)
residuals_matrix <- matrix(NA, nrow = 49, ncol = 49)
for (i in 2:50) {
  model <- lm(xrindus_yrly[, i] ~ APTxrtn[, i])
  betas[i-1, 1] <- coef(model)[2]
  alphas[i-1, 1] <- coef(model)[1]
  residuals_matrix[, i-1] <- residuals(model)
}


#Maybe a good R Shiny app
#Plot Observed Industry returns
for (i in 2:50){
  plot(xrindus_yrly[,1], xrindus_yrly[,i], type = "l", col = 'Black', xlab = "Time", ylab = "Excess Returns", main = "Industry Returns Over Time")
  lines(APTxrtn[,1], APTxrtn[,i], col = 'Blue')
}