#loading the package
BiocManager::install("MLSeq")
BiocManager::install("S4Vectors")
library("MLSeq")
library("S4Vectors")
library(DESeq2)
library(MLmetrics)

#setting path and loading file
setwd("C:/Users/dua_n/Desktop/Dissertation")

#loading complete data
mydata <- read.csv("C:/Users/dua_n/Desktop/Dissertation/count_data.csv", header=TRUE)


#data manipulation- PREPARING THE RAW COUNTS DATA 
colnames(mydata)[1]<-"GeneID"
colnames(mydata)[2]<-"GSM4313065"
colnames(mydata)[3]<-"GSM4313067"
colnames(mydata)[4]<-"GSM4313069"
colnames(mydata)[5]<-"GSM4313071"
colnames(mydata)[6]<-"GSM4313073"
colnames(mydata)[7]<-"GSM4313075"
colnames(mydata)[8]<-"GSM4313077"
colnames(mydata)[9]<-"GSM4313079"
colnames(mydata)[10]<-"GSM4313081"
colnames(mydata)[11]<-"GSM4313082"
View(mydata)


#making gene ids as row names
rownames(mydata) <- mydata$GeneID
mydata<-mydata[,-1]
View(mydata)

#just checking the loaded data
head(mydata[ ,1:10]) 


# Read in rownames from file- sybsetting wrt to wgcna driver genes
rownames_subset <- readLines("subset_wgcna.txt")

# Subset of data frame based on rownames read in from file and same columns
subset_df <- mydata[rownames(mydata) %in% 
                  rownames_subset, 
                c("GSM4313065",
                  "GSM4313067",
                  "GSM4313069",
                  "GSM4313071" ,
                  "GSM4313073",
                  "GSM4313075",
                  "GSM4313077",
                  "GSM4313079", 
                  "GSM4313081",
                  "GSM4313082")]

# Output the subset data frame
View(subset_df)


#CREATING CLASS DATA
#creating a variable class which will contain our labels
class_df <- data.frame(condition = factor(rep(c("M","N"), c(5, 5))))
rownames(class_df) <- c("GSM4313065", "GSM4313067", "GSM4313069","GSM4313071","GSM4313073","GSM4313075","GSM4313077","GSM4313079","GSM4313081","GSM4313082")  # set row names
View(class_df)


##SPLITTING THE DATA
set.seed(2128)
n_df <- ncol(subset_df)
p_df <- nrow(subset_df)


#here we select 30% columns as test data in oral cancer
nTest_df <- ceiling(n_df*0.3)

#we create the same 30% sample split in class data 
#to maintain uniformity
ind_df <- sample(n_df, nTest_df, FALSE)

#TRAIN SET
# extract the row names
rownames_df <- rownames(subset_df)

#putting oral cancer dataset as matrix and adding 1 
#to values so that there are no 0 values
data.train_df <- as.matrix(subset_df[ ,-ind_df] + 1)

# assign back the row names since in previous step rownames
#get lost in the dataframe.
rownames(data.train_df) <- rownames_df


#TEST SET
#putting oral cancer dataset as matrix and adding 1 
#to values so that there are no 0 values
data.test_df <- as.matrix(subset_df[ ,ind_df] + 1)

# assign back the row names since in previous step rownames
#get lost in the dataframe.
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

#1.CARRET BASED thus using

#pre-processing and not normalize. 
#deseq-vst: Normalization is applied with deseq 
#median ratio method. Variance stabiling 
#transformation is applied to the normalized data.

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


# 2. Discrete classifiers:
# Poisson Linear Discriminant Analysis
#pmodel <- classify(data = data.trainS4_sub, method = "PLDA", ref = "N",
#                  class.labels = "condition",normalize = "deseq",
#                 control = discreteControl(number = 5, repeats = 2,
#                                          tuneLength = 10, parallel = TRUE))
#pmodel

set.seed(2128)
# 3. Voom based Nearest Shrunken Centroids.
#this is a voom based classifier thus using normalize
#fit <- classify(data = data.trainS4_sub,
#               method = "voomNSC",
#              normalize = "deseq", 
#             ref = "N",
#            control = voomControl(tuneLength = 20))


# Support Vector Machines with Radial Kernel
#classifier-svm
#normalization-deseq
#transformation-rlog
fit_sub <- classify(data = data.trainS4_sub,
                method = "svmRadial",
                preProcessing = "deseq-rlog",
                ref = "N",
                control = trainControl(method = "repeatedcv", 
                                       number = 5,
                                       repeats = 3, 
                                       classProbs = TRUE))

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


#The final values used for the model were sigma
#= 8.992966e-06 and C = 0.25.
## An object of class "MLSeq"
## Model Description: Support Vector Machines with 
#Radial Basis Function Kernel (svmRadial)


# Define control list
ctrl.svm <- trainControl(method = "repeatedcv", 
                         number = 5, repeats = 1)


#trying to optimize by changing hyperparameters
#Support vector machines with radial basis function kernel
#svm, deseq, vst
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
#the accuracy and specificity is very poor of this model.


# Define control lists.- COMPARING MODELS PERFORMANCE
ctrl.continuous <- trainControl(method = "repeatedcv", 
                                number = 5, repeats = 10)
# Continuous classifiers, SVM and NSC


#svm,deseq,vst
fit.svm3 <- classify(data = data.trainS4_sub, 
                     method = "svmRadial",
                     preProcessing = "deseq-vst",
                     ref = "N", 
                     tuneLength = 10,
                     control = ctrl.continuous)


# Define control lists.- COMPARING MODELS PERFORMANCE
ctrl.continuous2_sub <- trainControl(method = "repeatedcv",
                                 number = 5, repeats = 5)

#pam,deseq,vst
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
