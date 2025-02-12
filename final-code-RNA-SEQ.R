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

#load("GSE145272.RData")
#count <- fc$counts
#head(count)
#write.csv(as.data.frame(count), file="count_data.csv")

#read in sample information
gse<- getGEO(GEO='GSE145272', GSEMatrix=TRUE)
metadata<-pData(phenoData(gse[[1]]))
head(metadata)
View(metadata)

#modifying metadata according to our needs- selecting required columns
metadata.modified<-metadata %>% 
  select(1,10,12)

#renaming columns
colnames(metadata.modified)[2] <- "tumor_stage"
colnames(metadata.modified)[3]<-"tissue"
head(metadata.modified)

#renaming cells and removing extra texts
metadata_renames<-mutate(metadata.modified,tumor_stage = gsub("tumor stage:","",tumor_stage), tissue=gsub("tissue:","",tissue))
colnames(metadata_renames)[1]<-"title"
head(metadata_renames)
View(metadata_renames)

#exploring gene expression data
counts_data<-read.csv("C:/Users/dua_n/Desktop/Dissertation/count_data.csv")
head(counts_data)
colnames(counts_data)[1]<-"GeneID"
colnames(counts_data)[2]<-"TM1"
colnames(counts_data)[3]<-"TM2"
colnames(counts_data)[4]<-"T3M"
colnames(counts_data)[5]<-"TC"
colnames(counts_data)[6]<-"LLP"
colnames(counts_data)[7]<-"B4M"
colnames(counts_data)[8]<-"B7M"
colnames(counts_data)[9]<-"B17M"
colnames(counts_data)[10]<-"T8M"
colnames(counts_data)[11]<-"T11M"

View(counts_data)


#rearranging columns in metadata renames to match in the next step
rownames(metadata_renames) <- metadata_renames$title
metadata_renames<- metadata_renames[,-1]

#changing t1n1 and all to metatstatic and non-metastatic
metadata_renames[1,1]<-"Metastatic"
metadata_renames[2,1]<-"Metastatic"
metadata_renames[3,1]<-"Metastatic"
metadata_renames[4,1]<-"Metastatic"
metadata_renames[5,1]<-"Metastatic"
metadata_renames[6,1]<-"Non-Metastatic"
metadata_renames[7,1]<-"Non-Metastatic"
metadata_renames[8,1]<-"Non-Metastatic"
metadata_renames[9,1]<-"Non-Metastatic"
metadata_renames[10,1]<-"Non-Metastatic"
head(metadata_renames)

#making genes as row names in counts data so that all column and rows match in the next step
rownames(counts_data)<-counts_data$GeneID
head(counts_data)
View(counts_data)


#deleting existing genes column 
counts_data2<- counts_data[,-1]
head(counts_data2)
View(counts_data2)

#ordering the columns of counts data to metadata's row names
#ordered_counts<- counts_data2[,c(5,6,3,7,8,1,9,10,4,2)]
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

#res0.1<-as.data.frame(res)
#write.csv(res0.1,'C:/Users/dua_n/Desktop/Dissertation/RNA-seq analysis/res0.1.csv')

#contrast
resultsNames(dds)

#changing alpha value
res0.01<-results(dds, alpha=0.01)
res0.01
summary(res0.01)
plotMA(res0.01)
#res001<-as.data.frame(res0.01)
#write.csv(res001,'C:/Users/dua_n/Desktop/Dissertation/RNA-seq analysis/res0.01.csv')


res0.05<-results(dds, alpha=0.05)
res0.05
summary(res0.05)
plotMA(res0.05)
res005<-as.data.frame(res0.05)
write.csv(res005,'C:/Users/dua_n/Desktop/Dissertation/RNA-seq analysis/res0.05.csv')



res0.09<-results(dds, alpha=0.09)
res0.09
summary(res0.09)
plotMA(res0.09)
res009<-as.data.frame(res0.09)
write.csv(res009,'C:/Users/dua_n/Desktop/Dissertation/RNA-seq analysis/res0.09.csv')


res0.5<-results(dds, alpha=0.5)
res0.5
summary(res0.5)
plotMA(res0.5)
res05<-as.data.frame(res0.5)
write.csv(res05,'C:/Users/dua_n/Desktop/Dissertation/RNA-seq analysis/res0.5.csv')


#final
res0.9<-results(dds, alpha=0.9)
res0.9
summary(res0.9)
plotMA(res0.9)

#volcanoplot
