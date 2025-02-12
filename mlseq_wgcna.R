#loading the package
BiocManager::install("MLSeq")
BiocManager::install("S4Vectors")
library("MLSeq")
library("S4Vectors")
library(DESeq2)
library(MLmetrics)

#setting path and loading file
setwd("path/to/my/data.csv")

#loading complete data
mydata <- read.csv("path/to/my/data.csv", header=TRUE)



#making gene ids as row names
rownames(mydata) <- mydata$GeneID
mydata<-mydata[,-1]
View(mydata)


# Read in rownames from file- sybsetting wrt to wgcna driver genes
rownames_subset <- readLines("wgcna.txt")

# Subset of data frame based on rownames read in from file and same columns
subset_df <- mydata[rownames(mydata) %in% 
                  rownames_subset, 
                c("1",
                  "2",
                  "3",
                  "4" ,
                  "5",
                  "6",
                  "7",
                  "8", 
                  "9",
                  "10")]



#CREATING CLASS DATA
#creating a variable class which will contain our labels
class_df <- data.frame(condition = factor(rep(c("M","N"), c(5, 5))))


##SPLITTING THE DATA
set.seed(2128)
n_df <- ncol(subset_df)
p_df <- nrow(subset_df)


#here we select 30% columns as test data 
nTest_df <- ceiling(n_df*0.3)

#we create the same 30% sample split in class data 
#to maintain uniformity
ind_df <- sample(n_df, nTest_df, FALSE)

#TRAIN SET
rownames_df <- rownames(subset_df)
data.train_df <- as.matrix(subset_df[ ,-ind_df] + 1)
rownames(data.train_df) <- rownames_df


#TEST SET
data.test_df <- as.matrix(subset_df[ ,ind_df] + 1)
rownames(data.test_df) <- rownames_df



#CLASSIFICATION MODELS
#converting class- train to a dataframe
classtr_df <- data.frame(condition = class_df[-ind_df, ], row.names = row.names(class_df)[-ind_df])

#converting class- test to a dataframe
classts_df <- data.frame(condition = class_df[ind_df, ], row.names = row.names(class_df)[ind_df])




#creating input training set for ML using DESeq object
data.trainS4_sub <- DESeqDataSetFromMatrix(countData = data.train_df, 
                                       colData = classtr_df, 
                                       formula(~1))

#creating input testing set for ML using DESeq object
data.testS4_sub = DESeqDataSetFromMatrix(countData = data.test_df, 
                                     colData = classts_df,
                                     design = formula(~1))

set.seed(2128)

#Normalization and transformation

#1.CARRET BASED 

#rf-classifier
#normalization-deseq
#transformation-vst
rf_sub <- classify(data = data.trainS4_sub, 
               method = "rf", 
               preProcessing = "deseq-vst", 
               ref = "N", 
               control = trainControl(method = "repeatedcv", 
                                      number = 5, 
                                      repeats = 2, 
                                      classProbs = TRUE))

rf_sub
MLSeq::plot(rf_sub)
trained(rf_sub)

#Predictions
pred.rf_sub <- predict(rf_sub, 
                       data.testS4_sub)
pred.rf_sub

pred.rf_sub <- relevel(pred.rf_sub, 
                       ref = "N")
actual.rf_sub <- relevel(classts_df$condition,
                         ref = "N")
tbl.rf_sub <- table(Predicted = pred.rf_sub, 
                    Actual = actual.rf_sub)
confusionMatrix(tbl.rf_sub,
                positive = "N")


rf2_sub <- classify(data = data.trainS4_sub, 
                method = "rf", 
                preProcessing = "deseq-vst", 
                ref = "N", 
                control = trainControl(method = "repeatedcv", 
                                       number = 5, 
                                       repeats = 4, 
                                       classProbs = TRUE))

rf2_sub
MLSeq::plot(rf2)

#Predictions
pred.rf2_sub <- predict(rf2_sub, data.testS4_sub)
pred.rf2_sub

pred.rf2_sub <- relevel(pred.rf2_sub, ref = "N")
actual.rf2_sub <- relevel(classts_df$condition, ref = "N")
tbl.rf2_sub <- table(Predicted = pred.rf2_sub, Actual = actual.rf2_sub)
confusionMatrix(tbl.rf2_sub, positive = "N")


#rf-classifier
#normalization-deseq
#transformation-rlog
rf3_sub <- classify(data = data.trainS4_sub, 
                method = "rf", 
                preProcessing = "deseq-rlog", 
                ref = "N", 
                control = trainControl(method = "repeatedcv", 
                                       number = 5,
                                       repeats = 2, 
                                       classProbs = TRUE))

rf3_sub
MLSeq::plot(rf3_sub)

#Predictions
pred.rf3_sub <- predict(rf3_sub, data.testS4_sub)
pred.rf3_sub

pred.rf3_sub <- relevel(pred.rf3_sub, ref = "N")
actual.rf3_sub <- relevel(classts_df$condition, ref = "N")
tbl.rf3_sub <- table(Predicted = pred.rf3_sub, Actual = actual.rf3_sub)
confusionMatrix(tbl.rf3_sub, positive = "N")




show(fit_sub)
MLSeq::plot(fit_sub)

#Predictions
pred.fit_sub <- predict(fit_sub, data.testS4_sub)
pred.fit_sub

pred.fit_sub <- relevel(pred.fit_sub, ref = "N")
actual.fit_sub <- relevel(classts_df$condition, ref = "N")
tbl.fit_sub <- table(Predicted = pred.fit_sub, Actual = actual.fit_sub)
confusionMatrix(tbl.fit_sub, positive = "N")

set.seed(1235)
fit2_sub <- classify(data = data.trainS4_sub,
                 method = "svmRadial",
                 preProcessing = "deseq-rlog",
                 ref = "N",
                 control = trainControl(method = "repeatedcv", 
                                        number = 5,
                                        repeats = 4, 
                                        classProbs = TRUE))

show(fit2_sub)
MLSeq::plot(fit2_sub)

#Predictions
pred.fit2_sub <- predict(fit2_sub, data.testS4_sub)
pred.fit2_sub

pred.fit2_sub <- relevel(pred.fit2_sub, ref = "N")
actual.fit2_sub <- relevel(classts_df$condition, ref = "N")
tbl.fit2_sub <- table(Predicted = pred.fit2_sub, Actual = actual.fit2_sub)
confusionMatrix(tbl.fit2_sub, positive = "N")


#MODEL OPTIMIZATION
# Support vector machines with radial basis function kernel
#classifier-svm
#normalization-deseq
#transformation-vst
set.seed(1236)
fit.svm_sub <- classify(data = data.trainS4_sub,
                    method = "svmRadial",
                    preProcessing = "deseq-vst", 
                    ref = "N", 
                    tuneLength = 10,
                    control = trainControl(method = "repeatedcv", 
                                           number = 5,
                                           repeats = 10, 
                                           classProbs = TRUE))
show(fit.svm_sub)
trained(fit.svm_sub)
plot(fit.svm_sub)

#Predictions
pred.fit.svm_sub <- predict(fit.svm_sub, data.testS4_sub)
pred.fit.svm_sub

pred.fit.svm_sub <- relevel(pred.fit.svm_sub, ref = "N")
actual.fit.svm_sub <- relevel(classts_df$condition, ref = "N")
tbl.fit.svm_sub <- table(Predicted = pred.fit.svm_sub, Actual = actual.fit.svm_sub)
confusionMatrix(tbl.fit.svm_sub, positive = "N")



# Define control list
ctrl.svm <- trainControl(method = "repeatedcv", 
                         number = 5, repeats = 1)



fit.svm2 <- classify(data = data.trainS4_sub,
                     method = "svmRadial",
                     preProcessing = "deseq-vst", 
                     ref = "N", tuneLength = 10,
                     control = ctrl.svm)
show(fit.svm2)
trained(fit.svm2)

#Predicted class labels
pred.svm <- predict(fit.svm2, data.testS4)
pred.svm


pred.svm <- relevel(pred.svm, ref = "N")
actual <- relevel(classts$condition, ref = "N")
tbl <- table(Predicted = pred.svm, Actual = actual)
confusionMatrix(tbl, positive = "N")


# Define control lists.- COMPARING MODELS PERFORMANCE
ctrl.continuous <- trainControl(method = "repeatedcv", 
                                number = 5, repeats = 10)




fit.svm3 <- classify(data = data.trainS4_sub, 
                     method = "svmRadial",
                     preProcessing = "deseq-vst",
                     ref = "N", 
                     tuneLength = 10,
                     control = ctrl.continuous)


# Define control lists.- COMPARING MODELS PERFORMANCE
ctrl.continuous2_sub <- trainControl(method = "repeatedcv",
                                 number = 5, repeats = 5)


set.seed(1234)
fit.NSC_sub <- classify(data = data.trainS4_sub, 
                    method = "pam",
                    preProcessing = "deseq-vst", 
                    ref = "N", 
                    tuneLength = 10,
                    control = ctrl.continuous2_sub)
show(fit.NSC_sub)
trained(fit.NSC_sub)
MLSeq::plot(fit.NSC_sub)

#Predictions
pred.svm <- predict(fit.svm, data.testS4)
pred.NSC_sub <- predict(fit.NSC_sub, data.testS4_sub)
pred.NSC_sub

pred.NSC_sub <- relevel(pred.NSC_sub, ref = "N")
actual.nsc_sub <- relevel(classts_df$condition, ref = "N")
tbl.nsc_sub <- table(Predicted = pred.NSC_sub, Actual = actual.nsc_sub)
confusionMatrix(tbl.nsc_sub, positive = "N")


selectedGenes(fit.NSC_sub)
