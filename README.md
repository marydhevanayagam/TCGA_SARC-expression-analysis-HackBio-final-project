# **Transcriptomic Profiling of Older Age Sarcoma Patients using TCGA RNA-seq data**

## 1. **Project Overview**

In this project, RNA-seq data for sarcoma (SARC) from The Cancer Genome Atlas (TCGA) was analyzed to characterize the transcriptomic profile of older age (OA: ≥65 years) compared to younger age (YA: 18-65 years) sarcoma patients. The significant differentially expressed genes, transcription factors, hub genes, and pathways that characterise the OA sarcoma patients were identified. 

## 2. **GitHub repository folders**

* Code \- contains the R script of the code for data preprocessing, differential gene expression analysis (DGEA), and functional enrichment analysis (FEA).  
* Data \- contains the data obtained after pre-processing, differential expression analysis, and network analysis.  
* Images \- contains the results generated from DGEA, FEA, and protein-protein interaction (PPI) network.  
* Report \- contains the project report.

## 3. **Requirements**

   The following R libraries were used

* TCGAbiolinks   
* SummarizedExperiment  
* dplyr  
* gplots   
* ggplot2    
* EnhancedVolcano   
* org.Hs.eg.db  
* reshape2  
  
These can be installed by running:  
install.packages(c(“dplyr”, “ggplot2”, “gplots”,“reshape2”))

To install Bioconductor packages:
* install.packages("BiocManager")  
* BiocManager::install("TCGAbiolinks")   
* BiocManager::install("SummarizedExperiment")  
* BiocManager::install("org.Hs.eg.db")  
* BiocManager::install("EnhancedVolcano")

## 4. **Methodology \- Code**

### **4.1.  Data collection and preprocessing** 

* RNA-seq and clinical data of sarcoma patient samples ("**TCGA-SARC**" project) were acquired from TCGA using the **TCGAbiolinks** package functions in R.   
* A query was prepared to retrieve "**Gene Expression Quantification**" data from the "**Transcriptome Profiling**" data category of the "**TCGA-SARC**" project, for the tissue types **"Primary Tumor", "Recurrent Tumor", "Metastatic"** using **GDCquery()** function. 
* Using **GDCdownload()** and **GDCprepare()** functions of the TCGAbiolinks package, the sample sets were downloaded and prepared for analysis.  
* From the retrieved data, metadata was obtained with the sample "**barcode**" and "**age\_at\_diagnosis**" fields.  
* Samples with missing or inaccurate age at diagnosis values were excluded.  
* The values of the "**age\_at\_diagnosis**" subgroup were converted from days to years, and the samples were divided into two age groups based on the threshold of 23,741.25 days for 65 years and 6574.5 for 18 years: **OA** (Older age: ≥65 years) and **YA** (Younger age: 18-65 years). 
* The unstranded dataset was selected for analysis.  
* The **TCGAanalyze\_Normalization()** function was used to normalize the gene expression data by gene length and read depth.  
* The **TCGAanalyze\_Filtering()** function was used to eliminate low-expression genes from the normalized data with the cut-off set at the first quantile (0.25).  
* The final filtered data was analyzed downstream for differential gene expression analysis (DGEA), functional enrichment analysis (FEA), transcription factor enrichment analysis (TFEA), gene-specific survival analysis, and network analysis.
  
### **4.2.  Differential gene expression analysis** 

* The **TCGAanalyze\_DEA()** function was used to compare the gene expression levels between the age groups ≥65 and 18-65 using the **edgeR** pipeline.   
* Genes were categorized as significantly upregulated and significantly downregulated based on a log fold-change (logFC) threshold of \>1.5 or \<(-1.5) and a false discovery rate (FDR) cut-off of \<0.005.
* The **mapIds()** function in the **org.Hs.eg.db** package was used to convert the Ensembl IDs of the DGEA results to gene symbols.    
* A volcano plot was generated for the differentially expressed genes (DEGs) using the **ggplot()** function and an enhanced volcano plot using the **EnhancedVolcano()** function.   
* A heatmap with diverging color palette was generated using the **heatmap.2()** function to visualize the significant DEGs. 
* The top 10 significantly up-regulated and down-regulated genes based on the logFC values were obtained, and the gene information was validated using the [GeneCards](https://www.genecards.org/) database.

### **4.3.  Functional enrichment analysis and pathway analysis**

* The **TCGAanalyze\_EAcomplete()** function was used to conduct functional enrichment analysis and pathway analysis for the genes with significant upregulation and downregulation.  
* The results of enrichment analysis were presented as bar plots using the **TCGAvisualize\_EAbarplot()** function to highlight the top 5 most enriched terms based on fold enrichment and FDR values for biological processes (BP), cellular components (CC), molecular functions (MF) and pathways.  

### **4.4. Network analysis**

* The protein-protein interaction (PPI) network of 538 out of 733 differentially expressed genes was visualized in [Cytoscape version 3.10.2](https://cytoscape.org/download.html) using the [STRING](https://string-db.org/) database.
* The network was subjected to ‘Network Analysis’ to obtain the network parameters.
* Additionally, the top 10 gene hubs were identified using Cytohubba plug-in from the network based on the closeness, degree, and MCC algorithms.
* The results from these three algorithms were analyzed using a Venn diagram to select the 10 top-ranking genes
