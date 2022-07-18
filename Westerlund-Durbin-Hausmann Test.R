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


# Defining the E ----------------------------------------------------------


errorMatFun <- function(u) {
  e_it_cap <- df_matrixOfWesterlundResidual[, u]
  e_it_lag <- lag(e_it_cap, k = 1)
  e_it_lag[is.na(e_it_lag)] <- 0
  def_e_it <- e_it_cap - e_it_lag
  
  return(t(as.matrix(cbind(e_it_cap, e_it_lag, def_e_it))))
}

e_1 <- errorMatFun(1)
e_2 <- errorMatFun(2)
e_3 <- errorMatFun(3)
e_4 <- errorMatFun(4)
e_5 <- errorMatFun(5)
e_6 <- errorMatFun(6)


theEIFun <- function(E_i){
  return(E_i %*% t(E_i))
}

E_1 <- theEIFun(e_1)
E_2 <- theEIFun(e_2)
E_3 <- theEIFun(e_3)
E_4 <- theEIFun(e_4)
E_5 <- theEIFun(e_5)
E_6 <- theEIFun(e_6)


