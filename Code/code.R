library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(gplots)
library(ggplot2)
library(reshape2)
library(biomaRt)

gdcprojects<- getGDCprojects()

getProjectSummary("TCGA-SARC") 

# Prepare the query
SARCq<- GDCquery(project = "TCGA-SARC",                    
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification",
                 sample.type = c("Metastatic", "Primary Tumor", "Recurrent Tumor")) 

# Download and prepare the dataset
GDCdownload(SARCq)
sarc.data<-GDCprepare(SARCq) 


# Create a metadata & obtain the sub-groups
sarcMetadata <- data.frame("barcode" = sarc.data$barcode,
                           "age_at_diagnosis" = sarc.data$age_at_diagnosis)

# Check for missing values and remove na's
sum(is.na(sarcMetadata$age_at_diagnosis))
sarcMetadata_na_omitted <- na.omit(sarcMetadata)
sum(is.na(sarcMetadata_na_omitted$age_at_diagnosis))
sum(sarcMetadata_na_omitted$age_at_diagnosis<0)

# Convert age_at_diagnosis from days to years (conversion factor: 1 year = 365.25 days, therefore  23,741.25 as the cutoff for 65 years and 6574.5 for 18 years)
sarcMetadata_na_omitted$age_group <- ifelse(sarcMetadata_na_omitted$age_at_diagnosis<23741.25 & 
                                            sarcMetadata_na_omitted$age_at_diagnosis>=6574.5, "YA", "OA")

# Group the samples based on age
YA_metadata<- sarcMetadata_na_omitted[sarcMetadata_na_omitted$age_group == 'YA',] #grouping samples with age between 18-65
OA_metadata<- sarcMetadata_na_omitted[sarcMetadata_na_omitted$age_group == 'OA',] #grouping samples with age >=65

# Select the unstranded dataset
sarc.raw.data <- assays(sarc.data)
dim(sarc.raw.data$unstranded)

# Obtain the unstranded data of the selected IDH groups
selectedBarcodes <- c(sample(YA_metadata$barcode), sample(OA_metadata$barcode))

write.csv(selectedBarcodes, file ="selected_barcodes.csv", row.names =TRUE)

selectedData <- sarc.raw.data$unstranded[,c(selectedBarcodes)]
dim(selectedData)

# Data normalization 
normData <- TCGAanalyze_Normalization(tabDF = selectedData, geneInfo = geneInfoHT, method = "geneLength")
dim(normData)

# Data filtering
filtData <- TCGAanalyze_Filtering(tabDF = normData, method = "quantile", qnt.cut = 0.25)
dim(filtData)

write.csv(filtData, file ="filtered_data.csv", row.names =TRUE)


# -----------------Differential expression analysis (DEA)-----------------------

YA_metadata <- sarcMetadata_na_omitted[sarcMetadata_na_omitted$age_group == 'YA',]
OA_metadata <- sarcMetadata_na_omitted[sarcMetadata_na_omitted$age_group == 'OA',]

selectResults<-TCGAanalyze_DEA(mat1 = filtData[, c(selectedBarcodes)[1:154]], #1st group is 18-65 (YA)
                               mat2 = filtData[, c(selectedBarcodes)[155:262]], #2nd group is >=65 (OA)
                               Cond1type = "YA", #Defining condition 1
                               Cond2type = "OA", #Defining condition 2
                               pipeline = "edgeR")

#Differential expression levels for the different conditions adds the average values for each group
selectResults.levels <- TCGAanalyze_LevelTab(selectResults,"YA" ,"OA", 
                                             filtData[,c(selectedBarcodes)[1:154]],
                                             filtData[,c(selectedBarcodes)[155:262]])

dim(selectResults)
dim(selectResults.levels)

# Set the logFC and p-value filter
selectResults.levels$diff_exp <- "No"
selectResults.levels$diff_exp[selectResults.levels$logFC > 1.5 & selectResults.levels$FDR <0.005] <-"UP"
selectResults.levels$diff_exp[selectResults.levels$logFC < (-1.5) & selectResults.levels$FDR <0.005] <-"DOWN"

table(selectResults.levels$diff_exp)
# Obtain the gene names as gene symbols
converted_gene_names<- mapIds(org.Hs.eg.db, 
                              keys = selectResults.levels$mRNA, 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL", 
                              multiVals = "first")

# Merge the conversion results back to the original dataframe
selectResults.levels$gene<-converted_gene_names
# Assign ensemble IDs to genes without gene names (lncRNA's etc.)
selectResults.levels$gene <- ifelse(is.na(selectResults.levels$gene), selectResults.levels$mRNA, selectResults.levels$gene)
sum(is.na(selectResults.levels$gene))

write.csv(selectResults.levels, file ="All_DGEA_results.csv", row.names = TRUE)

# Generate a volcano plot
ggplot(data = selectResults.levels, aes(x = logFC, y = (- log10(FDR)), col = diff_exp)) +  #basic Volcano plot
  geom_vline(xintercept = c(-1.5, 1.5), col = "blue", linetype = 'dashed') + #setting the threshold to 1
  geom_hline(yintercept = -log10(0.005), col = "red", linetype = 'dashed') + #setting the significance to 0.005
  geom_point() +  #make a continuous plot
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colors of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels we want to overwrite the categories from the dataframe (UP, DOWN, No)
  labs(color = 'Gene condition', x= expression("log"[2]*"FoldChange"), y = expression("-log"[10]*"p-adj value"))+
  coord_cartesian(ylim = c(0, 40), xlim = c(-10, 10)) + # to set the limits of the axis
  ggtitle("Volcano plot")


# Obtain the list of significant DEGs 
DE_results <- selectResults.levels[selectResults.levels$diff_exp == "UP" | selectResults.levels$diff_exp == "DOWN",]

write.csv(DE_results, file ="Sig_DGEA_results.csv", row.names =TRUE)

# Obtain the upregulated and downregulated genes
upreg.genes <- rownames(subset(selectResults.levels[selectResults.levels$diff_exp =='UP',]))
dnreg.genes <- rownames(subset(selectResults.levels[selectResults.levels$diff_exp =='DOWN',]))

write.csv(upreg.genes, "upregulated_genes.csv", row.names = TRUE)
write.csv(dnreg.genes, "downregulated_genes.csv", row.names = TRUE)

upreg_gene_symbols <- selectResults.levels$gene[selectResults.levels$diff_exp =='UP']
dnreg_gene_symbols <- selectResults.levels$gene[selectResults.levels$diff_exp =='DOWN']

write.csv(upreg_gene_symbols, "upregulated_gene_symbols.csv", row.names = TRUE)
write.csv(dnreg_gene_symbols, "downregulated_gene_symbols.csv", row.names = TRUE)


# Heatmap visualisation
dim(DE_results)
heat.data <- filtData[rownames(DE_results),] # selecting the genes that are significantly differentiated from the filtered data

dim(heat.data)
write.csv(heat.data, file ="heat_data.csv", row.names =TRUE)

# Color based on the age groups of the samples - column colors
cancer.type<-c(rep("YA",154), rep("OA",108))

ccodes<-c()
for(i in cancer.type)
{
  if(i == "OA")
    ccodes <- c(ccodes,"red")
  else
    ccodes <- c(ccodes, "blue")
}
ccodes


# Plot the Heatmap
par(oma = c(1,1,1,1)) #Setting outer margins
par(mar = c(1,1,1,1)) #setting inner plot margins
par(cex.main = 0.75) #size of the title
png(file= "heatmap.png", width =1500, height =1500, res=150)
heatmap.2(as.matrix(heat.data),
          col = hcl.colors(100, palette = "Blue-Red 3"), # Diverging palette
          Colv = F,                         # Cluster columns
          Rowv = F,                         # Cluster rows
          dendrogram = "none",              # No cluster both rows and columns
          trace = "none",                   # Remove trace lines
          scale = "row",                    # Standardizes rows (genes) across samples
          sepcolor = "black",               # Separate the columns
          key = TRUE,                       # Show color key
          cexRow = 0.5,                     # Adjust row label size
          cexCol = 0.5,                     # Adjust column label size
          margins = c(9, 7),                # Adjust margins to fit labels
          main = "Heatmap", #Title
          xlab = "Samples",                 # x-axis label
          ylab = "Genes",                   # y-axis label
          key.title = "Expression Level")
#ColSideColors = ccodes)  # columns are the samples and are color-coded based on the previous for loop
legend("topright", legend = c("OA group", "YA group"), fill = c("red", "blue"), title = "Column Colors", cex = 0.8)
dev.off()


# -------------------Functional enrichment analysis-----------------------------

up.EA <- TCGAanalyze_EAcomplete(TFname = "Upregulated", upreg_gene_symbols) # produces result based on BP, CC, MF and Pathways(P)
dn.EA <- TCGAanalyze_EAcomplete(TFname = "Downregulated", dnreg_gene_symbols)

# Visualization
EAbarplot_upreg_genes <- (TCGAvisualize_EAbarplot(tf = rownames(up.EA$ResBP), #Row names
                                                  GOBPTab = up.EA$ResBP, #results for BP
                                                  GOMFTab = up.EA$ResMF, #results for MF
                                                  GOCCTab = up.EA$ResCC, #results for CC
                                                  PathTab = up.EA$ResPat, #results for Pathway
                                                  nRGTab = upreg_gene_symbols, #number of genes in the list
                                                  nBar = 5, #max number of bars is 5 but can be increased to 10
                                                  text.size = 2, # 2 
                                                  fig.width = 30, # size of figure
                                                  fig.height = 15) #generates a pdf in the working directory
)


EAbarplot_downreg_genes <- (TCGAvisualize_EAbarplot(tf = rownames(dn.EA$ResBP),
                                                    GOBPTab = dn.EA$ResBP, 
                                                    GOMFTab = dn.EA$ResMF, 
                                                    GOCCTab = dn.EA$ResCC, 
                                                    PathTab = dn.EA$ResPat, 
                                                    nRGTab = dnreg_gene_symbols, 
                                                    nBar = 5, 
                                                    text.size = 2, 
                                                    fig.width = 30,
                                                    fig.height = 15))

