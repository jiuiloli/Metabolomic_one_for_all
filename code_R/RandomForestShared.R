

library(ggplot2)
library(lattice)
library(caret)
library(C50)
library(kernlab)
library(mlbench)
library(randomForest)
library(caretEnsemble)
library(MASS)
library(klaR)
library(nnet)
library(readxl)
library(HDclassif)
library(kknn)
library(gbm)

rm(list=ls())

name.split = function(x){
as.numeric(unlist(strsplit(x,"X"))[2])  
}

dat = data.frame(Class=mydata$Class, mydata[,9:ncol(mydata)])

set.seed(17)
# Stratified sampling
TrainingDataIndex <- createDataPartition(dat$Class, p=0.65, list = FALSE)
# Create Training Data 
trainingData <- dat[TrainingDataIndex,]
testData <- dat[-TrainingDataIndex,]
TrainingParameters <- trainControl(method = "repeatedcv", number = 10, repeats=10)

table(trainingData$Class)
table(testData$Class)

# training the model 
tgrid <- expand.grid(
  .mtry = ceiling(sqrt(2000)),
  .splitrule = "gini",
  .min.node.size = 4)

m_rf <- train(Class~., data=trainingData, method="ranger", trControl=TrainingParameters, preProc=c("center","scale"),
              tuneGrid = tgrid, num.trees = 100, importance = "permutation")

print(m_rf)

# validation using test data
m_rfPredictions <-predict(m_rf, testData)
(cmm_rf <-confusionMatrix(m_rfPredictions, testData$Class))







