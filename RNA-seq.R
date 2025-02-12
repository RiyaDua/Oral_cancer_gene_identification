#RNA-seq analysis on HTP expression profiling data on GSE145272
#DESeq2 workflow

#load libraries
library(dplyr)
library(tidyr)
library(cli)
library(tidyverse)
library(GEOquery)
library(DESeq2)
library(ggplot2)



#read in sample information
gse<- getGEO(GEO='GSE145272', GSEMatrix=TRUE)
metadata<-pData(phenoData(gse[[1]]))
head(metadata)
View(metadata)



#exploring gene expression data
counts_data<-read.csv("path/to/my/dataset.csv")
head(counts_data)

head(metadata_renames)

#making genes as row names in counts data so that all columns and rows match in the next step
rownames(counts_data)<-counts_data$GeneID
head(counts_data)
View(counts_data)


#ordering the columns of counts data to metadata's row names
rownames(metadata_renames)<-colnames(counts_data2)

all(colnames(counts_data2)%in%rownames(metadata_renames))
all(colnames(counts_data2)==rownames(metadata_renames))

#creating DESeq2 objects and round() since it convts to integers
dds<-DESeqDataSetFromMatrix(countData = round(counts_data2),
                            colData = metadata_renames,
                            design = ~tumor_stage)
dds


#pre-filtering: removing rows with low gene counts(lower than 10)
keep<-rowSums(counts(dds))>=5
keep

dds<-dds[keep,]
dds

dds$tumor_stage<-relevel(dds$tumor_stage, ref="Non-Metastatic")

#running DESEQ2
dds<-DESeq(dds)
res <-results(dds)
res

summary(res)
plotMA(res)


#contrast
resultsNames(dds)

#experimenting with alpha value
res0.01<-results(dds, alpha=0.01)
res0.01
summary(res0.01)
plotMA(res0.01)


res0.9<-results(dds, alpha=0.9)
res0.9
summary(res0.9)
plotMA(res0.9)


