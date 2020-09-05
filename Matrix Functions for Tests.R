library(timeSeries)

## Getting the data into matrices
View(originalData)
attach(originalData)

#### Data Preparation ####

gdpMatrix = matrix(data = c(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14),
                   nrow = length(y1), ncol = 14)
gdpColumns = c("y1", "y2", "y3", "y4", "y5", "y6", "y7", "y8", "y9", "y10", "y11", "y12",
               "y13", "y14")
colnames(gdpMatrix) = gdpColumns
misalMatrix = matrix(data = c(misal1, misal2, misal3, misal4, misal5, misal6, misal7, misal8,
                              misal9, misal10, misal11, misal12, misal13, misal14),
                     nrow = length(y1), ncol = 14)
misalColumns = c("misal1", "misal2", "misal3", "misal4", "misal5", "misal6", "misal7",
                 "misal8", "misal9", "misal10", "misal11", "misal12", "misal13", "misal14")
colnames(misalMatrix) = misalColumns
errorMatrix = matrix(data = c(ect1, ect2, ect3, ect4, ect5, ect6, ect7, ect8, ect9, ect10,
                              ect11, ect12, ect13, ect14), nrow = length(y1), ncol = 14)
errorColumns = c("ect1", "ect2", "ect3", "ect4", "ect5", "ect6", "ect7", "ect8", "ect9",
                 "ect10", "ect11", "ect12", "ect13", "ect14")
colnames(errorMatrix) = errorColumns

library(timeSeries)

#### Definition of constant parameters ####

timeLength = 45
numberOfIndividuals = 14

# Defining Tau of T
tauOfT = vector(length = timeLength)
for (i in 1:timeLength) {
  tauOfT[i] = 1
}
View(tauOfT)

# Identity Matrix
identityMatrix = diag(x = 1, nrow = timeLength, ncol = timeLength)
View(identityMatrix)


# mSubsriptT calculations
mSubsriptT = identityMatrix - (tauOfT %*% solve(crossprod(tauOfT)) %*% t(tauOfT))
View(mSubsriptT)

####tstat for independent GDP ####

# Matrix of Lagged Ys
matrixOfLaggedYs = matrix(nrow = timeLength, ncol = numberOfIndividuals)
for (i in 1:numberOfIndividuals) {
  useSeries = as.timeSeries(gdpMatrix[, i])
  new_data = lag(x = useSeries, k = 1)
  new_data[is.na(new_data)] <- 0
  matrixOfLaggedYs[, i] = as.vector(new_data)
}
View(matrixOfLaggedYs)

# Cubes of the matrix of Lagged Ys
matrixOfCubesOfLaggedYs = (matrixOfLaggedYs)^3
View(matrixOfCubesOfLaggedYs)


# Defining the sigma squares
vectorOfSigmaSquares = vector(length = 14)
for (i in 1:numberOfIndividuals) {
  ectToUse = as.vector(errorMatrix[, i])
  denominator = timeLength - 1
  vectorOfSigmaSquares[i] = (t(ectToUse) %*% mSubsriptT %*% ectToUse) / denominator
}
View(vectorOfSigmaSquares)
vectorOfSigmas = sqrt(vectorOfSigmaSquares)


# Vector of the t-statistics

vectorOfTStats = vector(length = 14)
for (i in 1:numberOfIndividuals) {
  Mu = errorMatrix[, i]
  yCube = matrixOfCubesOfLaggedYs[, i]
  sigmaOfI = vectorOfSigmas[i]
  laggedY = matrixOfLaggedYs[, i]
  
 # yCube[is.na(yCube)] <- 0
  # laggedY[is.na(laggedY)] <- 0
  
  denominatorWithoutSigma = sqrt((t(laggedY) %*% mSubsriptT %*% laggedY)^3)
  denominatorWithSigma = sigmaOfI * denominatorWithoutSigma
  
  numeratorOfTStat = t(Mu) %*% mSubsriptT %*% yCube
  
  vectorOfTStats[i] = numeratorOfTStat / denominatorWithSigma
}

vectorOfTStats

tStat = sum(vectorOfTStats) / numberOfIndividuals
tStat


### W subscript zbar statistic
expectedValueOfVectorOfTStats = -1.626
varianceOfVectorOfTStats = 0.727
zbar = 
  (sqrt(numberOfIndividuals) * (mean(vectorOfTStats) - expectedValueOfVectorOfTStats)) /
  sqrt(varianceOfVectorOfTStats)
zbar


####tstat for independent Misallignment ####

matrixOfMisalResiduals = matrix(data = residualForMisallignment, nrow = 45, ncol = 14,
                                byrow = F)
matrixOfMisalResiduals

# Matrix of Lagged misallignment
matrixOfLaggedMisallignment = matrix(nrow = timeLength, ncol = numberOfIndividuals)
for (i in 1:numberOfIndividuals) {
  useSeries = as.timeSeries(misalMatrix[, i])
  laggedEntity = as.vector(lag(x = useSeries, k = 1))
  laggedEntity[is.na(laggedEntity)] <- 0
  matrixOfLaggedMisallignment[, i] = laggedEntity
}
View(matrixOfLaggedMisallignment)

# Cubes of the matrix of Lagged Ys
matrixOfCubesOfLaggedMisallignment = (matrixOfLaggedMisallignment)^3
View(matrixOfCubesOfLaggedMisallignment)


# Defining the sigma squares
vectorOfSigmaSquaresOfMisal = vector(length = 14)
for (i in 1:numberOfIndividuals) {
  ectToUse = as.vector(misalMatrix[, i])
  denominator = timeLength - 1
  vectorOfSigmaSquaresOfMisal[i] = (t(ectToUse) %*% mSubsriptT %*% ectToUse) / denominator
}
View(vectorOfSigmaSquaresOfMisal)
vectorOfSigmasForMisal = sqrt(vectorOfSigmaSquaresOfMisal)


# Vector of the t-statistics

vectorOfTStatsOfMisalTest = vector(length = 14)
for (i in 1:numberOfIndividuals) {
  Mu = matrixOfMisalResiduals[, i]
  yCube = matrixOfCubesOfLaggedMisallignment[, i]
  sigmaOfI = vectorOfSigmasForMisal[i]
  laggedY = matrixOfLaggedMisallignment[, i]
  
  yCube[is.na(yCube)] <- 0
  laggedY[is.na(laggedY)] <- 0
  
  denominatorWithoutSigma = sqrt((t(laggedY) %*% mSubsriptT %*% laggedY)^3)
  denominatorWithSigma = sigmaOfI * denominatorWithoutSigma
  
  numeratorOfTStat = t(Mu) %*% mSubsriptT %*% yCube
  
  vectorOfTStatsOfMisalTest[i] = numeratorOfTStat / denominatorWithSigma
}

vectorOfTStatsOfMisalTest

tStatofMisalTest = sum(vectorOfTStatsOfMisalTest) / numberOfIndividuals
tStatofMisalTest

((sqrt(numberOfIndividuals) * tStat) / mean(vectorOfSigmas))
#### Test for Linearity ####

# Matrix of Second Lags of Ys
matrixOfSecondLagsOfYs = matrix(nrow = timeLength, ncol = numberOfIndividuals)
for (i in 1:numberOfIndividuals) {
  useSeries = as.timeSeries(gdpMatrix[, i])
  new_data = lag(x = useSeries, k = 2)
  new_data[is.na(new_data)] <- 0
  matrixOfSecondLagsOfYs[, i] = as.vector(new_data)
}
View(matrixOfSecondLagsOfYs)

wu = as.data.frame(matrixOfSecondLagsOfYs)
str(wu)
theLaggedYs = as.vector(c(wu$V1, wu$V2, wu$V3, wu$V4, wu$V5, wu$V6, wu$V7, wu$V8, 
                wu$V9, wu$V10, wu$V11, wu$V12, wu$V13, wu$V14))
str(theLaggedYs)

# Matrix of Second Lags of Misals
matrixOfSecondLagsOfMisallignment = matrix(nrow = timeLength, ncol = numberOfIndividuals)
for (i in 1:numberOfIndividuals) {
  useSeries = as.timeSeries(misalMatrix[, i])
  laggedEntity = as.vector(lag(x = useSeries, k = 2))
  laggedEntity[is.na(laggedEntity)] <- 0
  matrixOfSecondLagsOfMisallignment[, i] = laggedEntity
}
View(matrixOfSecondLagsOfMisallignment)

wo = as.data.frame(matrixOfSecondLagsOfMisallignment)
theLaggedMisal = as.vector(c(wo$V1, wo$V2, wo$V3, wo$V4, wo$V5, wo$V6, wo$V7, wo$V8, 
                             wo$V9, wo$V10, wo$V11, wo$V12, wo$V13, wo$V14))
str(theLaggedMisal)
theLaggedError = as.timeSeries(newCleanData$Residual)
theLaggedError = as.vector(lag(x = theLaggedError, k =  2))
theLaggedError[is.na(theLaggedError)] <- 0

# Combination of data to use
str(newCleanData)
cbind(newCleanData[, -c(1, 2)], theLaggedYs, theLaggedMisal, theLaggedError) -> newPanelData
str(newPanelData)

# Loading the needed packages
library(tibble)
library(PSTR)

# Making data suitable for PSTR package and creating the class PSTR
newPanelData = tibble(newPanelData)
linPstr = NewPSTR(data = newPanelData, 
                  dep = 1, indep = c(2, 3, 4, 5, 6), iT = 45, tvars = c(4, 5, 6), im = 3)
print(linPstr)
linPstr = LinTest(use = linPstr)
print(linPstr, "tests")



#### Lagged Matrix Functions ####

# Lags for y1
matrixOfLaggedY1 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y1 = as.timeSeries(y1)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y1, k = myIndex)
  matrixOfLaggedY1[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY1 = cbind(as.vector(lag(y1, k = timeLength)), matrixOfLaggedY1)
#####matrixOfLaggedDeltaY1 = cbind(matrixOfLaggedY1, y1)
View(matrixOfLaggedY1)


# Lags for y2
matrixOfLaggedY2 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y2 = as.timeSeries(y2)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y2, k = myIndex)
  matrixOfLaggedY2[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY2 = cbind(as.vector(lag(y2, k = timeLength)), matrixOfLaggedY2)
View(matrixOfLaggedY2)



# Lags for y3
matrixOfLaggedY3 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y3 = as.timeSeries(y3)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y3, k = myIndex)
  matrixOfLaggedY3[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY3 = cbind(as.vector(lag(y3, k = timeLength)), matrixOfLaggedY3)
View(matrixOfLaggedY3)


# Lags for y4
matrixOfLaggedY4 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y4 = as.timeSeries(y4)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y4, k = myIndex)
  matrixOfLaggedY4[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY4 = cbind(as.vector(lag(y4, k = timeLength)), matrixOfLaggedY4)
View(matrixOfLaggedY4)


# Lags for y5
matrixOfLaggedY5 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y5 = as.timeSeries(y5)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y5, k = myIndex)
  matrixOfLaggedY5[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY5 = cbind(as.vector(lag(y5, k = timeLength)), matrixOfLaggedY5)
View(matrixOfLaggedY5)



# Lags for y6
matrixOfLaggedY6 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y6 = as.timeSeries(y6)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y6, k = myIndex)
  matrixOfLaggedY6[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY6 = cbind(as.vector(lag(y6, k = timeLength)), matrixOfLaggedY6)
View(matrixOfLaggedY6)



# Lags for y7
matrixOfLaggedY7 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y7 = as.timeSeries(y7)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y7, k = myIndex)
  matrixOfLaggedY7[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY7 = cbind(as.vector(lag(y7, k = timeLength)), matrixOfLaggedY7)
View(matrixOfLaggedY7)


# Lags for y8
matrixOfLaggedY8 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y8 = as.timeSeries(y8)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y8, k = myIndex)
  matrixOfLaggedY8[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY8 = cbind(as.vector(lag(y8, k = timeLength)), matrixOfLaggedY8)
View(matrixOfLaggedY8)



# Lags for y9
matrixOfLaggedY9 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y9 = as.timeSeries(y9)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y9, k = myIndex)
  matrixOfLaggedY9[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY9 = cbind(as.vector(lag(y9, k = timeLength)), matrixOfLaggedY9)
View(matrixOfLaggedY9)



# Lags for y10
matrixOfLaggedY10 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y10 = as.timeSeries(y10)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y10, k = myIndex)
  matrixOfLaggedY10[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY10 = cbind(as.vector(lag(y10, k = timeLength)), matrixOfLaggedY10)
View(matrixOfLaggedY10)



# Lags for y11
matrixOfLaggedY11 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y11 = as.timeSeries(y11)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y11, k = myIndex)
  matrixOfLaggedY11[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY11 = cbind(as.vector(lag(y11, k = timeLength)), matrixOfLaggedY11)
View(matrixOfLaggedY11)



# Lags for y12
matrixOfLaggedY12 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y12 = as.timeSeries(y12)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y12, k = myIndex)
  matrixOfLaggedY12[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY14 = cbind(as.vector(lag(y12, k = timeLength)), matrixOfLaggedY12)
View(matrixOfLaggedY12)



# Lags for y13
matrixOfLaggedY13 = matrix(nrow = timeLength, ncol = timeLength - 1)

for (myIndex in 1:timeLength) {
  y13 = as.timeSeries(y13)
  columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y13, k = myIndex)
  matrixOfLaggedY13[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY14 = cbind(as.vector(lag(y13, k = timeLength)), matrixOfLaggedY13)
View(matrixOfLaggedY13)



# Lags for y14
matrixOfLaggedY14 = matrix(nrow = timeLength, ncol = timeLength - 1)
y14 = as.timeSeries(y14)
for (myIndex in 1:timeLength) {
    columnToUse = timeLength - myIndex
  laggedVariable = lag(x = y14, k = myIndex)
  matrixOfLaggedY14[,columnToUse] = as.vector(laggedVariable)
  
}
matrixOfLaggedY14 = cbind(as.vector(lag(y14, k = timeLength)), matrixOfLaggedY14)
View(matrixOfLaggedY14)












# Nonlinear Cointegration test -----------------------------------------
    # For Dependent GDP

repliRho = vector(length = numberOfIndividuals)
for (i in 1:numberOfIndividuals) {
  ect = errorMatrix[, i]
  muCap = as.timeSeries(ect)
  laggedMuCap = lag(muCap, k = 1)
  laggedMuCap[is.na(laggedMuCap)] <- 0
  
  numerator_A = as.vector(laggedMuCap * muCap)
  denominator_A = as.vector(laggedMuCap^2)
  repliRho[i] = (sum(numerator_A[-1])) / sum(denominator_A[-1])
}
repliRho

repliCointegrationStat = sum(repliRho) / numberOfIndividuals
repliCointegrationStat

