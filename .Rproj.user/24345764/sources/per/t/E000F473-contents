#----------------- Getting data ready for Clustering analysis---------------------------
View(newClean)
clusteringData = newClean[, -c(1, 2)]
head(clusteringData)

clusteringData$lnrgdp = as.numeric(clusteringData$lnrgdp)
clusteringData$lnenr = as.numeric(clusteringData$lnenr)
clusteringData$resid = as.numeric(clusteringData$resid)

#---------------- Clusterring Analysis -------------------------------------------------
clusteringResults = kmeans(clusteringData$lnrgdp, centers = 2)
clusteringResults

length(clusteringResults$cluster) == length(clusteringData$lnrgdp)

str(clusteringResults$cluster)
classifierFromClustering = clusteringResults$cluster

classifier = vector(length = length(classifierFromClustering))
for (i in 1:length(classifierFromClustering)) {
  if(classifierFromClustering[i] == 1){
    classifier[i] = 0
  } else {
    classifier[i] = 1
  }
}
classifier

logitData = as.data.frame(cbind(classifier, clusteringData[, -c(1,3)]))
head(logitData)

#---------------- Logistic Function ---------------------------------------------------
logitModel = glm(classifier ~ lnenr, data = logitData, family = "binomial")
plot(logitModel)

plot(classifier, logitData$lnrgdp)


vectorOfMeans = as.vector(logitModel$fitted.values)
mean(vectorOfMeans)

logitModel$coefficients

#--------------- Energy Consumption transformation -----------------------------------
transformableEnergy = read.csv("theDataForWinRATS.csv")
head(transformableEnergy)
header = list(2, 5, 8, 11, 14, 17)
for (i in header) {
  transformableEnergy[, i] <- log(abs(transformableEnergy[, i])^(1/2))
}

head(transformableEnergy)
write.csv(file = "justTransformed.csv", x = transformableEnergy)

#--------------- Analysis of the transformed data ------------------------------------
newlyTransformed = read.csv("structuredForR.csv")
head(newlyTransformed)
str(newlyTransformed)


library(plm)
linearModelOfTransformed = 
  plm(logGDP ~ trans_ENR, data = newlyTransformed, model = "within", effect = "time")
summary(linearModelOfTransformed)
newBaseResidual = residuals(linearModelOfTransformed)
head(newBaseResidual)

#--------------- Data preparation for RATS -------------------------------------------
newRealData = read.csv("justTransformed.csv", header = T)
head(newRealData)
summary(newRealData)
newRealData = newRealData[, -c(4, 7, 10, 13, 16, 19)]
str(newRealData)

# structurng the estimated residual into a matrix and coercing it into a data frame
newResidualGDP = 
  as.data.frame(matrix(data = newBaseResidual, nrow = length(newRealData[, 1]),
                       ncol = (length(newBaseResidual) / length(newRealData[, 1]))))
head(newResidualGDP)
residualNames = 
  c("res_Algeria", "res_Angola", "res_Egypt", "res_Gabon", "res_Nigeria", "res_Congo")
colnames(newResidualGDP, do.NULL = F)
colnames(newResidualGDP) <- residualNames

res_Algeria = newResidualGDP[, 1]
res_Angola = newResidualGDP[, 2]
res_Egypt = newResidualGDP[, 3]
res_Gabon = newResidualGDP[, 4]
res_Nigeria = newResidualGDP[, 5]
res_Congo = newResidualGDP[, 6]

dataOfAlgeria = cbind(newRealData[, c(2, 3)], res_Algeria)
dataOfAngola = cbind(newRealData[, c(4, 5)], res_Angola)
dataOfEgypt = cbind(newRealData[, c(6, 7)], res_Egypt)
dataOfGabon = cbind(newRealData[, c(8, 9)], res_Gabon)
dataOfNigeria = cbind(newRealData[, c(10, 11)], res_Nigeria)
dataOfCongo = cbind(newRealData[, c(12, 13)], res_Congo)
Years = rawData[, 1]

transformedDataForRATS = 
  cbind(Years, dataOfAlgeria, dataOfAngola, dataOfEgypt, dataOfGabon, dataOfNigeria,
        dataOfCongo)

str(transformedDataForRATS)

write.csv(transformedDataForRATS, file = "transformedDataForRATS.csv")

#----------------- Clustering and Loogistic functions -----------------------------
head(newlyTransformed)
head(newBaseResidual)

gammaC_Data = as.data.frame(cbind(newlyTransformed[, -c(1, 2)], newBaseResidual))
str(gammaC_Data)
head(gammaC_Data)


transCluster = kmeans(x = newlyTransformed[, c(3, 4)], centers = 2)
summary(transCluster)
transCluster$centers
classifiedByClusters = transCluster$cluster

df_classifier = vector(length = length(classifiedByClusters))
for (i in 1:length(classifiedByClusters)) {
  if(classifiedByClusters[i] == 1){
    df_classifier[i] = 0
  } else {
    df_classifier[i] = 1
  }
}
df_classifier

df_logitData = cbind(df_classifier, newlyTransformed[, c(3, 4)])
str(df_logitData)

df_logitData$df_classifier = as.factor(df_logitData$df_classifier)

df_logitModel = glm(df_classifier ~ ., data = df_logitData, family = "binomial")
df_logitModel$R





