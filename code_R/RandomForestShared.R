

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



set.seed(17)
# Stratified sampling
TrainingDataIndex <- createDataPartition(data$rn, p=0.7, list = FALSE)
# Create Training Data 
trainingData <- data[TrainingDataIndex,]
testData <- data[-TrainingDataIndex,]
TrainingParameters <- trainControl(method = "repeatedcv", number = 10, repeats=3)


# training the model 


m_rf <- train(rn~., data=trainingData, method="ranger", trControl=TrainingParameters, preProc=c("center","scale"),
               num.trees = 500, importance = "permutation")

print(m_rf)

# validation using test data
m_rfPredictions <-predict(m_rf, testData)
(cmm_rf <-confusionMatrix((table(m_rfPredictions, testData$rn)) ))








