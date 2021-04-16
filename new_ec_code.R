# Analysis for changed Energy Consumption variable -------------------------

library(plm)
library(PSTR)
library(tibble)
library(timeSeries)
library(tidyverse)

# new_ec_data <- read.csv("../further analysis/PSTRVEC redo/new_ec_data.csv")
new_ec_data <- read.csv("../further analysis/PSTRVEC redo/returnedEC.csv")

colNa_ <- c(colnames(new_ec_data[, -4]), "lnec2")
new_ec_data <- data.frame(new_ec_data[, -4],
                          # sqrt(new_ec_data[, 3]),
                          log(sqrt(new_ec_data[, 4])))
colnames(new_ec_data) <- colNa_
head(new_ec_data)



# Residual of GDP model ---------------------------------------------------

linearFixedEffect_GDP_model <- plm(formula = lnrgdp ~ lnec2, data = new_ec_data, 
                         model = "within", 
                         effect = "time"
                         # effect = "twoways"
                         )
summary(linearFixedEffect_GDP_model)
residual_GDP_model <- residuals(linearFixedEffect_GDP_model)

linTestData_GDP <- as_tibble(as.data.frame(cbind(new_ec_data, residual_GDP_model)))

# eViews_data <- data.frame(new_ec_data, residual_GDP_model)
# colnames(eViews_data) <- c(colnames(new_ec_data), "residual")
# write.csv(eViews_data, file = "the_newEC2_Eviews.csv")

# Residual of Energy model ------------------------------------------------

linearFixedEffect_Energy_model <- plm(formula = lnec2 ~ lnrgdp, data = new_ec_data, 
                         model = "within", 
                         # effect = "twoways"
                         effect = "time"
                         )
summary(linearFixedEffect_Energy_model)
residual_Energy_model <- residuals(linearFixedEffect_Energy_model)

linTestData_Energy <- as_tibble(as.data.frame(cbind(new_ec_data, residual_Energy_model)))



# Linearity test for GDP model --------------------------------------------

PSTR_class <- NewPSTR(data = linTestData_GDP, dep = 3, indep = c(3, 4, 5),
                      iT = length(dataForRATS$Years), tvars = c(3, 4, 5), im = 3)
PSTR_class = LinTest(use = PSTR_class)
print(PSTR_class, "tests")



# Linearity test for Energy model -----------------------------------------

PSTR_class2 <- NewPSTR(data = linTestData_GDP, dep = 4, indep = c(3, 4, 5),
                      iT = length(dataForRATS$Years), tvars = c(3, 4, 5), im = 3)
PSTR_class2 = LinTest(use = PSTR_class2)
print(PSTR_class2, "tests")


row_num <- 38


years_col <- new_ec_data[1:row_num, 2]
gdp_mat <- matrix(new_ec_data$lnrgdp, nrow = 38)
ec_mat <- matrix(new_ec_data$lnec2, nrow = 38)
matrixOfresid_Energy <- matrix(residual_Energy_model, nrow = 38)
matrixOfresid_gdp <- matrix(residual_GDP_model, nrow = 38)


df_changedEC <- data.frame(years_col, gdp_mat[, 1], ec_mat[, 1], matrixOfresid_Energy[, 1],
                           gdp_mat[, 2], ec_mat[, 2], matrixOfresid_Energy[, 2],
                           gdp_mat[, 3], ec_mat[, 3], matrixOfresid_Energy[, 3],
                           gdp_mat[, 4], ec_mat[, 4], matrixOfresid_Energy[, 4],
                           gdp_mat[, 5], ec_mat[, 5], matrixOfresid_Energy[, 5],
                           gdp_mat[, 6], ec_mat[, 6], matrixOfresid_Energy[, 6]
                           )
colnames(df_changedEC) <- colnames(dataForRATS)
head(df_changedEC)

# write.csv(df_changedEC, file = "theNewEC.csv")


# Nonlinear tests ---------------------------------------------------------

#### Data Preparation ####

df_individualsNumber = 6
df_timeLength = row_num

df_gdp_matrix <- df_changedEC %>% 
  select(LnRgdp_Algeria, LnRgdp_Angola, LnRgdp_Egypt, LnRgdp_Gabon, 
         LnRgdp_Nigeria, LnRgdp_Rcongo) %>% 
  as.matrix()

df_energy_matrix <- df_changedEC %>% 
  select(LnEnergy_Algeria, LnEnergy_Angola, LnEnergy_Egypt, 
         LnEnergy_Gabon, LnEnergy_Nigeria, LnEnery_RCongo) %>% 
  as.matrix()

#### Definition of constant parameters ####

# Defining Tau of T ------------------------------------------------
df_tauOfT = vector(length = df_timeLength)
for (i in 1:df_timeLength) {
  df_tauOfT[i] = 1
}
View(df_tauOfT)

# Identity Matrix
df_identityMatrix = diag(x = 1, nrow = df_timeLength, ncol = df_timeLength)
View(df_identityMatrix)


# mSubsriptT calculations ----------------------------------------------------
df_mSubsriptT = df_identityMatrix - (df_tauOfT %*% solve(crossprod(df_tauOfT)) %*% 
                                       t(df_tauOfT))
View(df_mSubsriptT)


####tstat for Dependent GDP --------------------------------------------------

library(timeSeries)

# Matrix of Lagged Ys
df_matrixOfLaggedYs = matrix(nrow = df_timeLength, ncol = df_individualsNumber)
for (i in 1:df_individualsNumber) {
  useSeries = as.timeSeries(df_gdp_matrix[, i])
  #Lag of the series
  new_data = lag(x = useSeries, k = 1)
  new_data[is.na(new_data)] <- 0
  df_matrixOfLaggedYs[, i] = as.vector(new_data)
}
View(df_matrixOfLaggedYs)

# Cubes of the matrix of Lagged Ys -------------------------------------------
df_matrixOfCubesOfLaggedYs = (df_matrixOfLaggedYs)^3
View(df_matrixOfCubesOfLaggedYs)



# Defining the sigma squares -----------------------------------------------
df_vectorOfSigmaSquares = vector(length = 6)
for (i in 1:df_individualsNumber) {
  ectToUse = as.vector(df_residual_matrix[, i])
  denominator = df_timeLength - 1
  df_vectorOfSigmaSquares[i] = (t(ectToUse) %*% df_mSubsriptT %*% ectToUse) / denominator
}
View(df_vectorOfSigmaSquares)
df_vectorOfSigmas = sqrt(df_vectorOfSigmaSquares)


# Vector of the t-statistics ------------------------------------------------

df_vectorOfTStats = vector(length = df_individualsNumber)
for (i in 1:df_individualsNumber) {
  Mu = df_residual_matrix[, i]
  yCube = df_matrixOfCubesOfLaggedYs[, i]
  sigmaOfI = df_vectorOfSigmas[i]
  laggedY = df_matrixOfLaggedYs[, i]
  
  yCube[is.na(yCube)] <- 0
  laggedY[is.na(laggedY)] <- 0
  
  denominatorWithoutSigma = sqrt((t(laggedY) %*% df_mSubsriptT %*% laggedY)^3)
  denominatorWithSigma = sigmaOfI * denominatorWithoutSigma
  
  numeratorOfTStat = t(Mu) %*% df_mSubsriptT %*% yCube
  
  df_vectorOfTStats[i] = numeratorOfTStat / denominatorWithSigma
}

df_vectorOfTStats

df_tStat = sum(df_vectorOfTStats) / df_individualsNumber
df_tStat


#-------- Evaluating the p-value -------------------------------------------
pValue_gdp = 2 * pt(-abs(df_tStat), df = 
                      (df_individualsNumber * df_timeLength) - df_individualsNumber)
pValue_gdp


####tstat for Dependent Energy ------------------------------------------

# Matrix of Lagged Energy Consumption -----------------------------
df_matrixOfLaggedEnergy = matrix(nrow = df_timeLength, ncol = df_individualsNumber)
for (i in 1:df_individualsNumber) {
  useSeries = as.timeSeries(df_energy_matrix[, i])
  laggedEntity = as.vector(lag(x = useSeries, k = 1))
  laggedEntity[is.na(laggedEntity)] <- 0
  df_matrixOfLaggedEnergy[, i] = laggedEntity
}
View(df_matrixOfLaggedEnergy)

# Cubes of the matrix of Lagged Ys
df_matrixOfCubesOfLaggedEnergy = (df_matrixOfLaggedEnergy)^3
View(df_matrixOfCubesOfLaggedEnergy)


# Defining the sigma squares
df_vectorOfSigmaSquaresOfEnergy = vector(length = df_individualsNumber)
for (i in 1:df_individualsNumber) {
  ectToUse = as.vector(df_residual_matrix[, i])
  denominator = df_timeLength - 1
  df_vectorOfSigmaSquaresOfEnergy[i] = (t(ectToUse) %*% 
                                          df_mSubsriptT %*% ectToUse) / denominator
}
View(df_vectorOfSigmaSquaresOfEnergy)
df_vectorOfSigmasOfEnergy = sqrt(df_vectorOfSigmaSquaresOfEnergy)


# Vector of the t-statistics

df_vectorOfTStatsOfEnergyTest = vector(length = df_individualsNumber)
for (i in 1:df_individualsNumber) {
  Mu = df_residual_matrix[, i]
  yCube = df_matrixOfCubesOfLaggedEnergy[, i]
  sigmaOfI = df_vectorOfSigmasOfEnergy[i]
  laggedY = df_matrixOfLaggedEnergy[, i]
  
  yCube[is.na(yCube)] <- 0
  laggedY[is.na(laggedY)] <- 0
  
  denominatorWithoutSigma = sqrt((t(laggedY) %*% df_mSubsriptT %*% laggedY)^3)
  denominatorWithSigma = sigmaOfI * denominatorWithoutSigma
  
  numeratorOfTStat = t(Mu) %*% df_mSubsriptT %*% yCube
  
  df_vectorOfTStatsOfEnergyTest[i] = numeratorOfTStat / denominatorWithSigma
}

df_vectorOfTStatsOfEnergyTest

df_tStatofEnergyTest = sum(df_vectorOfTStatsOfEnergyTest) / df_individualsNumber
df_tStatofEnergyTest

#-------- Evaluating the p-value -------------------------------------------
pValue_energy = 2 * pt(-abs(df_tStatofEnergyTest), df = 
                         (df_individualsNumber * df_timeLength) - df_individualsNumber)
pValue_energy


#----------------------------------------------------------------------------------------
# Nonlinear Cointegration test ----------------------------------------------------------


# Dependent GDP ------------------------------------------------------------------
rho = vector(length = df_individualsNumber)
for (i in 1:df_individualsNumber) {
  ect = matrixOfresid_gdp[, i]
  muCap = as.timeSeries(ect)
  laggedMuCap = lag(muCap, k = 1)
  laggedMuCap[is.na(laggedMuCap)] <- 0
  
  numerator_A = as.vector(laggedMuCap) * as.vector(muCap)
  denominator_A = as.vector(laggedMuCap^2)
  rho[i] = (sum(numerator_A[-1])) / sum(denominator_A[-1])
}
rho

cointegrationStat = sum(rho) / df_individualsNumber
cointegrationStat

pValue_GDP.coin = 2 * pt(-abs(cointegrationStat), df = 
                           (df_individualsNumber * df_timeLength) - df_individualsNumber)
pValue_GDP.coin


# Dependent Energy ----------------------------------------------------------------
enrRho = vector(length = df_individualsNumber)
for (i in 1:df_individualsNumber) {
  ect = matrixOfresid_Energy[, i]
  muCap = as.timeSeries(ect)
  laggedMuCap = lag(muCap, k = 1)
  laggedMuCap[is.na(laggedMuCap)] <- 0
  
  numerator_A = as.vector(laggedMuCap) * as.vector(muCap)
  denominator_A = as.vector(laggedMuCap^2)
  enrRho[i] = (sum(numerator_A[-1])) / sum(denominator_A[-1])
}
enrRho

enr_cointegrationStat = sum(enrRho) / df_individualsNumber
enr_cointegrationStat

pValue_Ener.coin = 2 * pt(-abs(enr_cointegrationStat),
                          df = (df_individualsNumber * df_timeLength) - df_individualsNumber)
pValue_Ener.coin
