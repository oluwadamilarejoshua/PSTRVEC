library(plm)
library(timeSeries)

#-- Data importation and preparation ------------------------------------------------------

dataForTests = read.csv("newNonLinearTest_data.csv")
head(dataForTests)

df_individuals = 6
df_time = 38

#---- Generating the residuals ----------------------------------------------------------

df_modelForWesterlund = plm(lnec ~ lnperinc + es + open, effect = "time", 
                         model = "within", data = dataForTests)
df_residualForWesterlund = residuals(df_modelForWesterlund)
summary(df_modelForWesterlund)
head(df_residualForWesterlund)

df_matrixOfWesterlundResidual = matrix(data = df_residualForWesterlund, nrow = df_time)
head(df_matrixOfWesterlundResidual)


#---- E subscript i,t computations -----------------------------------------------

capitalE_subIT = function(dataToUse){
  require(timeSeries)
  useData = as.timeSeries(dataToUse)
  lagged_useData = lag(useData, k = 1)
  lagged_useData[is.na(lagged_useData)] <- 0
  differencedSeries = as.vector(useData) - as.vector(lagged_useData)
  nameOfMatrixToCreate = as.matrix(rbind(t(as.timeSeries(lagged_useData)),
                                         t(as.timeSeries(differencedSeries))))
  return(nameOfMatrixToCreate)
}

capitalE_individual_1 = capitalE_subIT(df_matrixOfWesterlundResidual[, 1])
dim(capitalE_individual_1)

capitalE_individual_2 = capitalE_subIT(df_matrixOfWesterlundResidual[, 2])
dim(capitalE_individual_2)

capitalE_individual_3 = capitalE_subIT(df_matrixOfWesterlundResidual[, 3])
dim(capitalE_individual_3)

capitalE_individual_4 = capitalE_subIT(df_matrixOfWesterlundResidual[, 4])
dim(capitalE_individual_4)

capitalE_individual_5 = capitalE_subIT(df_matrixOfWesterlundResidual[, 5])
dim(capitalE_individual_5)

capitalE_individual_6 = capitalE_subIT(df_matrixOfWesterlundResidual[, 6])
dim(capitalE_individual_6)


#---- sigma squared and sigma computations -------------------------------------------


df_sigma_squared = vector()
df_sigma = vector()
for (error in 1:df_individuals) {
  variableInUse = df_matrixOfWesterlundResidual[, error]
  df_sigma_squared[error] = mean(variableInUse ^ 2)
  df_sigma[error] = sqrt(mean(variableInUse ^ 2))
}
df_sigma_squared
df_sigma

df_populationSigma = sqrt(mean(df_sigma_squared))


#---- The E matrices ---------------------------------------------------------------

capitalEsubI = function(matrixInUse){
  reducedMatrix = matrixInUse %*% t(matrixInUse)
  return(as.matrix(reducedMatrix))
}

EsubI_individual_1 = capitalEsubI(capitalE_individual_1)
dim(EsubI_individual_1)

EsubI_individual_2 = capitalEsubI(capitalE_individual_2)
dim(EsubI_individual_2)

EsubI_individual_3 = capitalEsubI(capitalE_individual_3)
dim(EsubI_individual_3)

EsubI_individual_4 = capitalEsubI(capitalE_individual_4)
dim(EsubI_individual_4)

EsubI_individual_5 = capitalEsubI(capitalE_individual_5)
dim(EsubI_individual_5)

EsubI_individual_6 = capitalEsubI(capitalE_individual_6)
dim(EsubI_individual_6)


#---- The four (4) test statistics ----------------------------------------------------

df_listOfE_matrices = list(EsubI_individual_1, EsubI_individual_2, EsubI_individual_3, 
                           EsubI_individual_4, EsubI_individual_5, EsubI_individual_6)
df_vectorOf_Ei11 = vector(length = df_individuals)
df_vectorOf_Ei12 = vector(length = df_individuals)

for (eMatrices in 1:length(df_listOfE_matrices)) {
    df_vectorOf_Ei11[eMatrices] = df_listOfE_matrices[[eMatrices]][1, 1]
  df_vectorOf_Ei12[eMatrices] = df_listOfE_matrices[[eMatrices]][1, 2]
}
df_vectorOf_Ei11
df_vectorOf_Ei12

    # EP subscript gamma

EP_gamma = sum(df_vectorOf_Ei12) / sum(df_vectorOf_Ei11)
EP_gamma
p_value_EP_gamma = pnorm(-abs(EP_gamma))
p_value_EP_gamma

    # EP subscript t

EP_t = sum(df_vectorOf_Ei12) / (df_populationSigma * sqrt(sum(df_vectorOf_Ei11)))
EP_t
p_value_EP_t = pnorm(-abs(EP_t))
p_value_EP_t

    # EG subscript gamma

EG_gamma = sum(df_vectorOf_Ei12 / df_vectorOf_Ei11)
EG_gamma
p_value_EG_gamma = pnorm(-abs(EG_gamma))
p_value_EG_gamma

    # EG subscript t

EG_t = sum(df_vectorOf_Ei12 / (df_sigma * sqrt(df_vectorOf_Ei11)))
EG_t
p_value_EG_t = pnorm(-abs(EG_t))
p_value_EG_t
