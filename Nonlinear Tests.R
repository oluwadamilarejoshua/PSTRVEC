library(timeSeries)

View(rawData)
str(rawData)
attach(rawData)

#### Data Preparation ####

df_individualsNumber = 6
df_timeLength = length(LnRgdp_Egypt)

df_gdp_matrix = as.matrix(cbind(LnRgdp_Algeria, LnRgdp_Angola, LnRgdp_Egypt, LnRgdp_Gabon, 
                                LnRgdp_Nigeria, LnRgdp_Rcongo))
df_energy_matrix = as.matrix(cbind(LnEnergy_Algeria, LnEnergy_Angola, LnEnergy_Egypt, 
                                   LnEnergy_Gabon, LnEnergy_Nigeria, LnEnery_RCongo))

df_residual_matrix = as.matrix(residualFrame)

head(matrixOfResidualsOfIndependentEnergy)
matrixOfResidualsOfDependentEnergy = matrixOfResidualsOfIndependentEnergy
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


####tstat for independent GDP --------------------------------------------------

library(timeSeries)

# Matrix of Lagged Ys
df_matrixOfLaggedYs = matrix(nrow = df_timeLength, ncol = df_individualsNumber)
for (i in 1:df_individualsNumber) {
  useSeries = as.timeSeries(df_gdp_matrix[, i])
  new_data = lag(x = useSeries, k = 1)
  new_data[is.na(new_data)] <- 0
  df_matrixOfLaggedYs[, i] = as.vector(new_data)
}
View(df_matrixOfLaggedYs)

# Cubes of the matrix of Lagged Ys -------------------------------------------
df_matrixOfCubesOfLaggedYs = (df_matrixOfLaggedYs)^3
View(df_matrixOfCubesOfLaggedYs)



# Defining the sigma squares -----------------------------------------------
df_vectorOfSigmaSquares = vector(length = 14)
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
pValue = 2 * pt(-abs(df_tStat), df = (df_individualsNumber * df_timeLength) - df_individualsNumber)
pValue


####tstat for independent Energy ------------------------------------------

# Matrix of Lagged misallignment ------------------------------------------
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




#----------------------------------------------------------------------------------------
# Nonlinear Cointegration test ----------------------------------------------------------


    # Dependent GDP ------------------------------------------------------------------
rho = vector(length = df_individualsNumber)
for (i in 1:df_individualsNumber) {
  ect = df_residual_matrix[, i]
  muCap = as.timeSeries(ect)
  laggedMuCap = lag(muCap, k = 1)
  laggedMuCap[is.na(laggedMuCap)] <- 0
  
  numerator_A = as.vector(laggedMuCap * muCap)
  denominator_A = as.vector(laggedMuCap^2)
  rho[i] = (sum(numerator_A[-1])) / sum(denominator_A[-1])
}
rho

cointegrationStat = sum(rho) / df_individualsNumber
cointegrationStat


    # Dependent Energy ----------------------------------------------------------------
enrRho = vector(length = df_individualsNumber)
for (i in 1:df_individualsNumber) {
  ect = matrixOfResidualsOfDependentEnergy[, i]
  muCap = as.timeSeries(ect)
  laggedMuCap = lag(muCap, k = 1)
  laggedMuCap[is.na(laggedMuCap)] <- 0
  
  numerator_A = as.vector(laggedMuCap * muCap)
  denominator_A = as.vector(laggedMuCap^2)
  enrRho[i] = (sum(numerator_A[-1])) / sum(denominator_A[-1])
}
enrRho

enr_cointegrationStat = sum(enrRho) / df_individualsNumber
enr_cointegrationStat

#mean(enrRho) / (var(enrRho) / sqrt(df_individualsNumber))
