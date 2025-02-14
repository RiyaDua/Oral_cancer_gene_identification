#WGCNA
#LOADING LIBRARIES
#BiocManager::install("WGCNA")
#install.packages("remotes")
#remotes::install_github("kevinblighe/CorLevelPlot")
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)


#loading the dataset
counts_data<-read.csv('path/to/my/data.csv')
View(counts_data)

#getting metadata
phenoData<-read.csv('path/to/my/data.csv')
head(phenoData)


#making gene ids as row names
rownames(phenoData) <- phenoData$geo_accession
phenoData<-phenoData[,-4]
View(phenoData)

#making geo accession as row names in phenodata
rownames(counts_data) <- counts_data$GeneID
counts_data<-counts_data[,-1]

#identifying outliers
gsg<-goodSamplesGenes(t(counts_data))
summary(gsg)

gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

data<-counts_data[gsg$goodGenes == TRUE,]


#detecting outliers using hierarchial clustering- method 1
htree<-hclust(dist(t(counts_data)),method= 'average')
plot(htree)


#principle component analysis
pca<-prcomp(t(counts_data))
pca_data<-pca$x
View(pca_data)
pca_data<- as.data.frame(pca_data)

#calculating variance
pca_var<-pca$sdev^2
pca_var_percent<-round(pca_var/sum(pca_var)*100 , digits=2)
pca_var_percent

#creating a plot between PC1 and PC2 and labeling them as their percentage
ggplot(pca_data, aes(PC1, PC2))+
  geom_point()+
  geom_text(label= row.names(pca_data))+
  labs(x=paste0('PC1',pca_var_percent[1],'%'),
       y=paste0('PC2',pca_var_percent[2],'%'))

#exclude the samples with outliers
samples_to_be_excluded<- c('1','2')
data_subset<-counts_data[,!(colnames(counts_data) %in% samples_to_be_excluded)]
View(data_subset)

#excluding outliers from phenodata as well
colData<-phenoData%>%
  filter(!row.names(.)%in%samples_to_be_excluded)

#making sure that col names and row names in pheno and counts is same
all(colnames(data_subset)%in%rownames(colData))
all(colnames(data_subset)==rownames(colData))


#normalization using DESeq2-since WGCNA requires normalized data


#creating DESeq2 dataset
dds<-DESeqDataSetFromMatrix(countData = data_subset,
                            colData = colData,
                            design = ~ 1) #not specifying design
dds


dds75<-dds[rowSums(counts(dds) >=15) >=6,]
nrow(dds75)

#performing variance stabilization
dds_norm<-vst(dds75)


#normalization- it gives normalized values
norm_counts<-assay(dds_norm)%>%
  t()
View(norm_counts)


#NETWORK CONSTRUCTION
#choosing soft threshold values.
power<-c(c(1:10),seq(from=12, to=50, by=2))
power


#calling network topology analysis function
sft<-pickSoftThreshold(norm_counts,
                  powerVector = power,
                  networkType = 'signed',
                  verbose=5)

sft_data<-sft$fitIndices
View(sft_data)

#visualize to pick a power- soft threshold value
a1<-ggplot(sft_data, aes(Power,SFT.R.sq, label= Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept= 0.8, color='red')+
  labs(x='Power',y='scale free topology fit, signed R^2')+
  theme_classic()

#plot for mean connectivity
a2<-ggplot(sft_data, aes(Power,mean.k., label=Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  labs(x='Power',y='mean connectivity')+
  theme_classic()

grid.arrange(a1, a2, nrow=2)



norm_counts[]<-sapply(norm_counts, as.numeric)
soft_power<-28
temp_cor<-cor
cor<-WGCNA::cor


#memory estimate wrt block size
bwnet<- blockwiseModules(norm_counts,
                 maxBlockSize = 22000,
                 TOMType = 'signed',
                 power=soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3)
cor<-temp_cor


#module eigengene information
module_eigengenes<-bwnet$MEs
head(module_eigengenes)
View(module_eigengenes)


#to get number of genes in same module
modules<-table(bwnet$colors)
View(modules)


# Using updated color names in plotDendroAndColors() function
plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(bwnet$unmergedColors, bwnet$colors),
                    c('unmerged', 'merged'),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)


traits <- colData%>%
  mutate(tumor_stage_bin = ifelse(grepl('Metastatic', tumor_stage), 1, 0))%>%
  select(4)


# Defining numbers of genes and samples
nSamples <- nrow(norm_counts)
nGenes <- ncol(norm_counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)



# visualizing module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')




CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[65],
             y = names(heatmap.data)[1:64],
             col = c("blue1", "skyblue", "white", "pink", "red"))




# Identifying driver genes
module.membership.measure <- cor(module_eigengenes, norm_counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals %>% 
  as.data.frame() %>% 
  head(25)

module.membership.measure.pvals[1:10,1:10]

gene.signf.corr <- cor(norm_counts, traits$tumor_stage_bin, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)
