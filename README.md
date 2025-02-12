# 🧬 Identification of Oral Cancer Candidate Genes using Machine Learning Approaches

This repository contains reusable scripts and workflow details from my Master's thesis project, where we developed a **machine learning-based model** to identify candidate genes for **Oral Cancer** using **RNA-seq data** and network analysis.  
*(Full study is under review for publication.)*  

## 🔍 Research Overview
- **Objective**: Develop a **machine learning model** to classify oral cancer genes based on **RNA-seq differential expression data**.
- **Data Source**: Publicly available dataset *(Preprocessed data not included due to publication restrictions).*
- **Modeling Approach**: Feature selection using **DESeq2**, network-based filtering via **WGCNA**, and classification using **MLseq**.

## 🏆 Key Methods & Tools Used
### 🔹 **Upstream RNA-Seq Processing (Linux)**
- Quality Control: `FastQC`
- Alignment & Preprocessing: `STAR`, `HPC`

### 🔹 **Downstream Analysis & Machine Learning (R)**
- **Differential Expression Analysis**: `DESeq2`, `ggplot2`, `tidyverse`
- **Feature Selection & Normalization**: `DESeq-VST`, `WGCNA`
- **Network Analysis**: `Cytoscape`
- **Machine Learning Classification**: `MLseq`
  - **Best Performing Model**: SVM (`Support Vector Machine`)
  - **Final Model Accuracy**: **99%** *(Under validation, pending publication.)*

