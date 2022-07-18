library(plm)
library(timeSeries)

#-- Data importation and preparation ------------------------------------------------------
dataForTests = read.csv("newNonLinearTest_data.csv")
head(dataForTests)

df_individuals = 6
df_time = 38


#-- The tau and the identity matrix--------------------------------------------------------
df_tau = vector(length = df_time)
for (i in 1:df_time) {
  df_tau[i] = 1
}
head(df_tau)

df_theIdentityMatrix = diag(x = 1, nrow = df_time, ncol = df_time)
head(df_theIdentityMatrix)


#-- The Observed y -----------------------------------------------------------------------
observed_Y = dataForTests[, 3]
head(observed_Y)

df_matrixOfObservedY = matrix(data = observed_Y, nrow = df_time)
head(df_matrixOfObservedY)

df_matrixOfLaggedObservedY = matrix(nrow = df_time, ncol = df_individuals)
for (i in 1:df_individuals) {
  yToUse = as.timeSeries(df_matrixOfObservedY[, i])
  df_laggedY = lag(x = yToUse, k = 1)
  df_laggedY[is.na(df_laggedY)] <- 0
  
  df_matrixOfLaggedObservedY[, i] = as.vector(df_laggedY)
}
head(df_matrixOfLaggedObservedY)

df_meanOfCrossSectionalObservedY = vector(length = df_time)
for (i in 1:df_time) {
  df_meanOfCrossSectionalObservedY[i] = mean(df_matrixOfLaggedObservedY[i,])
}
head(df_meanOfCrossSectionalObservedY)



#--  The Intercept model and it's residuals -------------------------------------------------
df_model_intercept = plm(lnec ~ lnperinc + es + open, effect = "individual", 
                         model = "within", data = dataForTests)
df_residual_intercept = residuals(df_model_intercept)
summary(df_model_intercept)
head(df_residual_intercept)

df_matrixOfInterceptResidual = matrix(data = df_residual_intercept, nrow = df_time)
head(df_matrixOfInterceptResidual)

        # Mean of cross sectional residuals of the intercept only
df_intercept_meanOfCrossSectionalResiduals = vector(length = df_time)
for (i in 1:df_time) {
  df_intercept_meanOfCrossSectionalResiduals[i] = mean(df_matrixOfInterceptResidual[i,])
}
head(df_intercept_meanOfCrossSectionalResiduals)



#-- The Intercept & Trend model and it's residuals -----------------------------------------
df_model_inteceptAndTrend = plm(lnec ~ lnperinc + es + open, effect = "twoways",
                                model = "within", data = dataForTests)
df_residual_inteceptAndTrend = residuals(df_model_inteceptAndTrend)
summary(df_residual_inteceptAndTrend)
head(df_residual_inteceptAndTrend)

df_matrixOfInterceptAndTrendResidual = matrix(data = df_residual_inteceptAndTrend,
                                              nrow = df_time)
head(df_matrixOfInterceptAndTrendResidual)

        # Mean of cross sectional residuals of the intercept and trend
df_interceptAndTrend_meanOfCrossSectionalResiduals = vector(length = df_time)
for (i in 1:df_time) {
  df_interceptAndTrend_meanOfCrossSectionalResiduals[i] = 
    mean(df_matrixOfInterceptAndTrendResidual[i,])
}
head(df_interceptAndTrend_meanOfCrossSectionalResiduals)


#-- The w_bar -----------------------------------------------------------------------------
          # Intercept
df_intercept_wBar = as.matrix(cbind(df_tau, df_intercept_meanOfCrossSectionalResiduals,
                     df_meanOfCrossSectionalObservedY))
head(df_intercept_wBar)

          # Intecept & Trend
df_interceptAndTrend_wBar = as.matrix(cbind(df_tau, df_interceptAndTrend_meanOfCrossSectionalResiduals,
  df_meanOfCrossSectionalObservedY))
head(df_interceptAndTrend_wBar)


#-- The G -----------------------------------------------------------------------------
          # Intercept

df_intercept_matrixOfG_firstIndividual = as.matrix(cbind(df_intercept_wBar,
                                                         df_matrixOfLaggedObservedY[, 1]))
df_intercept_matrixOfG_secondIndividual = as.matrix(cbind(df_intercept_wBar,
                                                         df_matrixOfLaggedObservedY[, 2]))
df_intercept_matrixOfG_thirdIndividual = as.matrix(cbind(df_intercept_wBar,
                                                         df_matrixOfLaggedObservedY[, 3]))
df_intercept_matrixOfG_fourthIndividual = as.matrix(cbind(df_intercept_wBar,
                                                         df_matrixOfLaggedObservedY[, 4]))
df_intercept_matrixOfG_fifthIndividual = as.matrix(cbind(df_intercept_wBar,
                                                         df_matrixOfLaggedObservedY[, 5]))
df_intercept_matrixOfG_sixthIndividual = as.matrix(cbind(df_intercept_wBar,
                                                         df_matrixOfLaggedObservedY[, 6]))


          # Intecept & Trend 

df_interceptAndTrend_matrixOfG_firstIndividual = as.matrix(cbind(df_interceptAndTrend_wBar,
                                                         df_matrixOfLaggedObservedY[, 1]))
df_interceptAndTrend_matrixOfG_secondIndividual = as.matrix(cbind(df_interceptAndTrend_wBar,
                                                          df_matrixOfLaggedObservedY[, 2]))
df_interceptAndTrend_matrixOfG_thirdIndividual = as.matrix(cbind(df_interceptAndTrend_wBar,
                                                         df_matrixOfLaggedObservedY[, 3]))
df_interceptAndTrend_matrixOfG_fourthIndividual = as.matrix(cbind(df_interceptAndTrend_wBar,
                                                          df_matrixOfLaggedObservedY[, 4]))
df_interceptAndTrend_matrixOfG_fifthIndividual = as.matrix(cbind(df_interceptAndTrend_wBar,
                                                         df_matrixOfLaggedObservedY[, 5]))
df_interceptAndTrend_matrixOfG_sixthIndividual = as.matrix(cbind(df_interceptAndTrend_wBar,
                                                         df_matrixOfLaggedObservedY[, 6]))

#-- The M sub(i, w) --------------------------------------------------------------------
            # The intercept
df_intercept_mSubiw_firstIndividual = df_theIdentityMatrix - (df_intercept_matrixOfG_firstIndividual %*% 
                                                                solve(crossprod(
                                                                  df_intercept_matrixOfG_firstIndividual)) %*% 
                                                                t(df_intercept_matrixOfG_firstIndividual))

df_intercept_mSubiw_secondIndividual = df_theIdentityMatrix - (df_intercept_matrixOfG_secondIndividual %*% 
                                                                 solve(crossprod(
                                                                   df_intercept_matrixOfG_secondIndividual)) %*% 
                                                                 t(df_intercept_matrixOfG_secondIndividual))

df_intercept_mSubiw_thirdIndividual = df_theIdentityMatrix - (df_intercept_matrixOfG_thirdIndividual %*%
                                                                solve(crossprod(
                                                                  df_intercept_matrixOfG_thirdIndividual)) %*% 
                                                                t(df_intercept_matrixOfG_thirdIndividual))

df_intercept_mSubiw_fourthIndividual = df_theIdentityMatrix - (df_intercept_matrixOfG_fourthIndividual %*%
                                                                solve(crossprod(
                                                                  df_intercept_matrixOfG_fourthIndividual)) %*% 
                                                                t(df_intercept_matrixOfG_fourthIndividual))

df_intercept_mSubiw_fifthIndividual = df_theIdentityMatrix - (df_intercept_matrixOfG_fifthIndividual %*%
                                                                 solve(crossprod(
                                                                   df_intercept_matrixOfG_fifthIndividual)) %*% 
                                                                 t(df_intercept_matrixOfG_fifthIndividual))

df_intercept_mSubiw_sixthIndividual = df_theIdentityMatrix - (df_intercept_matrixOfG_sixthIndividual %*%
                                                                solve(crossprod(
                                                                  df_intercept_matrixOfG_sixthIndividual)) %*% 
                                                                t(df_intercept_matrixOfG_sixthIndividual))

          # Intecept and Trend
df_interceptAndTrend_mSubiw_firstIndividual = df_theIdentityMatrix - (df_interceptAndTrend_matrixOfG_firstIndividual %*% 
                                                                solve(crossprod(
                                                                  df_interceptAndTrend_matrixOfG_firstIndividual)) %*% 
                                                                t(df_interceptAndTrend_matrixOfG_firstIndividual))

df_interceptAndTrend_mSubiw_secondIndividual = df_theIdentityMatrix - (df_interceptAndTrend_matrixOfG_secondIndividual %*% 
                                                                        solve(crossprod(
                                                                          df_interceptAndTrend_matrixOfG_secondIndividual)) %*% 
                                                                        t(df_interceptAndTrend_matrixOfG_secondIndividual))

df_interceptAndTrend_mSubiw_thirdIndividual = df_theIdentityMatrix - (df_interceptAndTrend_matrixOfG_thirdIndividual %*% 
                                                                         solve(crossprod(
                                                                           df_interceptAndTrend_matrixOfG_thirdIndividual)) %*% 
                                                                         t(df_interceptAndTrend_matrixOfG_thirdIndividual))

df_interceptAndTrend_mSubiw_fourthIndividual = df_theIdentityMatrix - (df_interceptAndTrend_matrixOfG_fourthIndividual %*% 
                                                                        solve(crossprod(
                                                                          df_interceptAndTrend_matrixOfG_fourthIndividual)) %*% 
                                                                        t(df_interceptAndTrend_matrixOfG_fourthIndividual))

df_interceptAndTrend_mSubiw_fifthIndividual = df_theIdentityMatrix - (df_interceptAndTrend_matrixOfG_fifthIndividual %*% 
                                                                         solve(crossprod(
                                                                           df_interceptAndTrend_matrixOfG_fifthIndividual)) %*% 
                                                                         t(df_interceptAndTrend_matrixOfG_fifthIndividual))

df_interceptAndTrend_mSubiw_sixthIndividual = df_theIdentityMatrix - (df_interceptAndTrend_matrixOfG_sixthIndividual %*% 
                                                                        solve(crossprod(
                                                                          df_interceptAndTrend_matrixOfG_sixthIndividual)) %*% 
                                                                        t(df_interceptAndTrend_matrixOfG_sixthIndividual))

#-- The m_bar subscript w -------------------------------------------------------------------------------------
        # Intercept
df_intercept_mBarSubW = df_theIdentityMatrix - (df_intercept_wBar %*% solve(crossprod(df_intercept_wBar)) %*% 
                                                  t(df_intercept_wBar))

        # Intercept and Trend
df_interceptAndTrend_mBarSubW = df_theIdentityMatrix - (df_interceptAndTrend_wBar %*% solve(crossprod(df_interceptAndTrend_wBar)) %*% 
                                                          t(df_interceptAndTrend_wBar))

#-- The sigma computation -------------------------------------------------------------------------------------
df_degreeOfFreedom = df_time - 4

        # Intercept
dim(t(df_matrixOfInterceptResidual[, 1]) %*% df_intercept_mBarSubW %*% df_matrixOfInterceptResidual[, 1])

df_intercept_sigma = vector()
for (res in 1:df_individuals) {
  residualToUse = df_matrixOfInterceptResidual[, res]
  sigmaSquared = (t(residualToUse) %*% df_intercept_mBarSubW %*% residualToUse) / df_degreeOfFreedom
  df_intercept_sigma[res] = sqrt(sigmaSquared)
}
df_intercept_sigma

        # Inercept and Trend
df_interceptAndTrend_sigma = vector()
for (res in 1:df_individuals) {
  residualToUse = df_matrixOfInterceptAndTrendResidual[, res]
  sigmaSquared = (t(residualToUse) %*% df_intercept_mBarSubW %*% residualToUse) / df_degreeOfFreedom
  df_interceptAndTrend_sigma[res] = sqrt(sigmaSquared)
}
df_interceptAndTrend_sigma

#-- The t-stats and CIPS stat computations --------------------------------------------------------------------
        # Intercept
df_intercept_vectorOfTStats = vector()
for (individual in 1:df_individuals) {
  residualToUse = df_matrixOfInterceptResidual[, individual]
  laggedYToUse = df_matrixOfLaggedObservedY[, individual]
  sigmaToUse = df_intercept_sigma[individual]
  
  df_denominator = sigmaToUse / sqrt(t(laggedYToUse) %*% df_intercept_mBarSubW %*% laggedYToUse)
  df_numerator = t(residualToUse) %*% df_intercept_mBarSubW %*% laggedYToUse
  
  df_intercept_vectorOfTStats[individual] = df_numerator / df_denominator
}
df_intercept_vectorOfTStats

df_intercept_CIPS = mean(df_intercept_vectorOfTStats)
df_intercept_CIPS

df_intercept_pValue = pt(-abs(df_intercept_CIPS), df = 222)
df_intercept_pValue

        # Intercept and Trend
df_interceptAndTrend_vectorOfTStats = vector()
for (individual in 1:df_individuals) {
  residualToUse = df_matrixOfInterceptAndTrendResidual[, individual]
  laggedYToUse = df_matrixOfLaggedObservedY[, individual]
  sigmaToUse = df_interceptAndTrend_sigma[individual]
  
  df_denominator = sigmaToUse / sqrt(t(laggedYToUse) %*% df_intercept_mBarSubW %*% laggedYToUse)
  df_numerator = t(residualToUse) %*% df_intercept_mBarSubW %*% laggedYToUse
  
  df_interceptAndTrend_vectorOfTStats[individual] = df_numerator / df_denominator
}
df_interceptAndTrend_vectorOfTStats

df_interceptAndTrend_CIPS = mean(df_interceptAndTrend_vectorOfTStats)
df_interceptAndTrend_CIPS

df_interceptAndTrend_pValue = pt(-abs(df_interceptAndTrend_CIPS), df = 222)
df_interceptAndTrend_pValue
