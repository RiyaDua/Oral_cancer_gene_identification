#WGCNA
#LOAD LIBRARIES
#BiocManager::install("WGCNA")
#install.packages("remotes")
#remotes::install_github("kevinblighe/CorLevelPlot")
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)


#load the dataset
counts_data<-read.csv('C:/Users/dua_n/Desktop/Dissertation/count_data.csv')
View(counts_data)
#getting metadata
phenoData<-read.csv('C:/Users/dua_n/Desktop/Dissertation/metadata.csv')
head(phenoData)

#MUTATING COL NAMES OF COUNTS DATA
colnames(counts_data)[1]<-"GeneID"
colnames(counts_data)[2]<-"GSM4313065"
colnames(counts_data)[3]<-"GSM4313067"
colnames(counts_data)[4]<-"GSM4313069"
colnames(counts_data)[5]<-"GSM4313071"
colnames(counts_data)[6]<-"GSM4313073"
colnames(counts_data)[7]<-"GSM4313075"
colnames(counts_data)[8]<-"GSM4313077"
colnames(counts_data)[9]<-"GSM4313079"
colnames(counts_data)[10]<-"GSM4313081"
colnames(counts_data)[11]<-"GSM4313082"
View(counts_data)

#making gene ids as row names
rownames(phenoData) <- phenoData$geo_accession
phenoData<-phenoData[,-4]
View(phenoData)

#making geo accession as row names in phenodata
rownames(counts_data) <- counts_data$GeneID
counts_data<-counts_data[,-1]
#quality control-outlier detection - this function needs data as row as samples and columns as genes 
#thus we have to transpose data using function t

gsg<-goodSamplesGenes(t(counts_data))
summary(gsg)

#if allok says true then there are no outliers
#if it says false- we have outlier genes

gsg$allOK

#check how many genes are outliers
table(gsg$goodGenes)
#9145 genes are outliers as per result
#they need to be excluded from the analysis

#check the same for sample
table(gsg$goodSamples)
#true- thus all samples are good to go

#exclude outlier genes from the data
data<-counts_data[gsg$goodGenes == TRUE,]


#detecting outliers using hierarchial clustering- method 1
htree<-hclust(dist(t(counts_data)),method= 'average')
plot(htree)


#principle component analysis- method 2 to detect outliers
pca<-prcomp(t(counts_data))
pca_data<-pca$x
View(pca_data)
pca_data<- as.data.frame(pca_data)
#this has principle component computed for each sample


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
samples_to_be_excluded<- c('GSM4313073','GSM4313082')
data_subset<-counts_data[,!(colnames(counts_data) %in% samples_to_be_excluded)]
View(data_subset)

#excluding outliers from phenodata as well
colData<-phenoData%>%
  filter(!row.names(.)%in%samples_to_be_excluded)

#making sure that col names and row names in pheno and counts is same
all(colnames(data_subset)%in%rownames(colData))
all(colnames(data_subset)==rownames(colData))


#normalization using DESeq2-since WGCNA requires normalized data
#in case of fpkm data we just need to log transform the data

#creating DESeq2 dataset
dds<-DESeqDataSetFromMatrix(countData = data_subset,
                            colData = colData,
                            design = ~ 1) #not specifying design
dds

#REMOVING GENES WITH LOW COUNTS- gene having 15 counts or more in 75% of samples
#i.e. genes with 15 or more counts in 6 samples will be left
#21952
dds75<-dds[rowSums(counts(dds) >=15) >=6,]
nrow(dds75)

#performing variance stabilization
dds_norm<-vst(dds75)


#normalization- it gives normalized values
#these values now need to be reversed like samples as rows and genes as columns
#thus we transpose it
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


#sof-threshold=28
norm_counts[]<-sapply(norm_counts, as.numeric)
soft_power<-28
temp_cor<-cor
cor<-WGCNA::cor

# unsigned -> nodes with positive & negative correlation are treated equally 
# signed -> nodes with negative correlation are considered *unconnected*, treated as zero
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


#write.csv(modules,'C:/Users/dua_n/Desktop/Dissertation/WGCNA/modules.csv')
#modules_data<-read.csv("C:/Users/dua_n/Desktop/Dissertation/WGCNA/modules.csv")
#View(modules_data)

#grey modules= genes that doesnt fall into any module are kept in grey module


#to get module colors dendrogram
#mergedColors = labels2colors(bwnet$colors)
#plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
 #                   "Module colors",
  #                  dendroLabels = FALSE, hang = 0.03,
   #                 addGuide = TRUE, guideHang = 0.05)



#plot dendrogram and colors before and merging underneath
#not working- LABELS NAMES TAKEN BY THE FUNCTION HAD SOME ERRORS.
#THUS, CHANGING EACH LABEL NAME.
bwnet$colors <- gsub("1$", "", bwnet$colors)
bwnet$unmergedColors <- gsub("1$", "", bwnet$unmergedColors)
bwnet$unmergedColors[bwnet$unmergedColors == 'magenta41'] <- 'magenta3'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'magenta41'] <- 'magenta3'  # Replace with a valid color name or code

# Check bwnet object for color 'purple.1' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'purple.1'] <- 'purple'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'purple.1'] <- 'purple'  # Replace with a valid color name or code

# Check bwnet object for color 'orangered1.1' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'orangered41'] <- 'orangered1'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'orangered41'] <- 'orangered1'  # Replace with a valid color name or code

# Check bwnet object for color 'lavenderblush3.1' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'lavenderblush11'] <- 'lavenderblush1'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'lavenderblush11'] <- 'lavenderblush1'  # Replace with a valid color name or code

# Check bwnet object for color 'darkmagenta.1' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'darkmagenta.1'] <- 'darkmagenta'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'darkmagenta.1'] <- 'darkmagenta'  # Replace with a valid color name or code

# Check bwnet object for color 'chocolate4.1' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'chocolate4.1'] <- 'chocolate'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'chocolate'] <- 'chocolate'  # Replace with a valid color name or code

# Check bwnet object for color 'chocolate4.1' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'yellow3.1'] <- 'yellow'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'yellow3.1'] <- 'yellow'  # Replace with a valid color name or code

# Check bwnet object for color 'chocolate4.1' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'midnightblue.1'] <- 'midnightblue'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'midnightblue.1'] <- 'midnightblue'  # Replace with a valid color name or code

# Check bwnet object for color 'darkolivegreen11' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'darkolivegreen21'] <- 'darkolivegreen2'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'darkolivegreen21'] <- 'darkolivegreen2'  # Replace with a valid color name or code

# Check bwnet object for color 'darkolivegreen11' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'grey60'] <- 'grey'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'grey60'] <- 'grey'  # Replace with a valid color name or code

# Check bwnet object for color 'greenyellow1' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'greenyellow1'] <- 'greenyellow'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'greenyellow1'] <- 'greenyellow'  # Replace with a valid color name or code

# Check bwnet object for color 'tan41' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'tan41'] <- 'tan'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'tan41'] <- 'tan'  # Replace with a valid color name or code

# Check bwnet object for color 'tan41' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'thistle31'] <- 'thistle3'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'thistle31'] <- 'thistle3'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'mediumpurple11'] <- 'mediumpurple1'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'mediumpurple11'] <- 'midiumpurple1'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'palevioletred21'] <- 'palevioletred'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'palevioletred21'] <- 'palevioletred'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'blueviolet1'] <- 'blueviolet'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'blueviolet1'] <- 'blueviolet'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'lightslateblue1'] <- 'lightslateblue'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'lightslateblue1'] <- 'lightslateblue1'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'lightslateblue1'] <- 'lightslateblue'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'lightslateblue1'] <- 'lightslateblue1'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'lightsteelblue11'] <- 'lightsteelblue'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'lightsteelblue11'] <- 'lightsteelblue'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'skyblue11'] <- 'skyblue1'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'skyblue11'] <- 'skyblue1'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'lightblue31'] <- 'lightblue3'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'lightblue31'] <- 'lightblue3'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'palevioletred31'] <- 'palevioletred3'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'palevioletred31'] <- 'palevioletred3'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'darkolivegreen41'] <- 'darkolivegreen3'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'darkolivegreen41'] <- 'darkolivegreen3'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'lightpink31'] <- 'lightpink2'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'lightpink31'] <- 'lightpink2'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'lightpink2'] <- 'lightpink2'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'lightpink21'] <- 'lightpink2'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'firebrick31'] <- 'firebrick3'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'firebrick31'] <- 'firebrick3'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'coral11'] <- 'coral1'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'coral11'] <- 'coral1'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'green41'] <- 'green3'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'green41'] <- 'green3'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'sienna41'] <- 'sienna1'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'sienna41'] <- 'sienna1'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'honeydew11'] <- 'honeydew1'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'honeydew11'] <- 'honeydew1'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'darkviolet1'] <- 'darkviolet'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'darkviolet1'] <- 'darkviolet'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'antiquewhite4'] <- 'antiquewhite3'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'antiquewhite4'] <- 'antiquewhite3'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'yellowgreen1'] <- 'yellowgreen'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'yellowgreen1'] <- 'yellowgreen'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'floralwhite1'] <- 'floralwhite'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'floralwhite1'] <- 'floralwhite'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'darkslateblue1'] <- 'darkslateblue'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'darkslateblue1'] <- 'darkslateblue'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'lightskyblue41'] <- 'lightskyblue'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'lightskyblue41'] <- 'lightskyblue'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'pink31'] <- 'pink3'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'pink31'] <- 'pink3'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'darkgreen1'] <- 'darkgreen'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'darkgreen1'] <- 'darkgreen'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'blue21'] <- 'blue2'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'blue21'] <- 'blue2'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'yellow21'] <- 'yellow2'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'yellow21'] <- 'yellow2'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'deeppink11'] <- 'deeppink1'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'deeppink11'] <- 'deeppink1'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'violet1'] <- 'violet'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'violet1'] <- 'violet'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'lightcyan11'] <- 'lightcyan'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'lightcyan11'] <- 'lightcyan'  # Replace with a valid color name or code

# Check bwnet object for color 'mediumpurple31' and replace it with a valid color
bwnet$unmergedColors[bwnet$unmergedColors == 'salmon4'] <- 'salmon'  # Replace with a valid color name or code
bwnet$colors[bwnet$colors == 'salmon4'] <- 'salmon'  # Replace with a valid color name or code



# Remove '.' from color names in bwnet object
#bwnet$unmergedColors <- gsub("\\.", "", bwnet$unmergedColors)
#bwnet$colors <- gsub("\\.", "", bwnet$colors)

# Convert numeric colors to character in bwnet object
#bwnet$unmergedColors <- as.character(bwnet$unmergedColors)
#bwnet$colors <- as.character(bwnet$colors)

# Use updated color names in plotDendroAndColors() function
plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(bwnet$unmergedColors, bwnet$colors),
                    c('unmerged', 'merged'),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

# 6A. Relate modules to traits --------------------------------------------------
# module trait associations


colData[5,2]<-"Non-metastatic"
colData[6,2]<-"Non-metastatic"
colData[7,2]<-"Non-metastatic"
colData[8,2]<-"Non-metastatic"

# create traits file - binarize categorical variables
traits <- colData%>%
  mutate(tumor_stage_bin = ifelse(grepl('Metastatic', tumor_stage), 1, 0))%>%
  select(4)


# Define numbers of genes and samples
nSamples <- nrow(norm_counts)
nGenes <- ncol(norm_counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)



# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')




CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[65],
             y = names(heatmap.data)[1:64],
             col = c("blue1", "skyblue", "white", "pink", "red"))


#not working
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()

# 6B. Intramodular analysis: Identifying driver genes ---------------



# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm_counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals %>% 
  as.data.frame() %>% 
  head(25)

module.membership.measure.pvals[1:10,1:10]


# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm_counts, traits$tumor_stage_bin, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)


# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.