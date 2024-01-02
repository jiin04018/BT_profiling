####### All Data Integration

library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)

library(Seurat)
library(data.table)
library(tools)

library(gplots)

library(RColorBrewer)
library(yarrr)
library(colorspace)
library(ggpubr)

#### Load the MNG dataset
MNG1.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE183655_MNG_sc/samples/GSM5567093/")
MNG2.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE183655_MNG_sc/samples/GSM5567094/")
MNG3.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE183655_MNG_sc/samples/GSM5567095/")
MNG4.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE183655_MNG_sc/samples/GSM5567096/")
MNG5.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE183655_MNG_sc/samples/GSM5567098/")
MNG6.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE183655_MNG_sc/samples/GSM5567101/")


#### Create Seurat object
MNG1 <- CreateSeuratObject(counts = MNG1.data, project = "gse183655", min.cells = 3, min.features = 200)
MNG2 <- CreateSeuratObject(counts = MNG2.data, project = "gse183655", min.cells = 3, min.features = 200)
MNG3 <- CreateSeuratObject(counts = MNG3.data, project = "gse183655", min.cells = 3, min.features = 200)
MNG4 <- CreateSeuratObject(counts = MNG4.data, project = "gse183655", min.cells = 3, min.features = 200)
MNG5 <- CreateSeuratObject(counts = MNG5.data, project = "gse183655", min.cells = 3, min.features = 200)
MNG6 <- CreateSeuratObject(counts = MNG6.data, project = "gse183655", min.cells = 3, min.features = 200)



#### Add disease to metadata
disease = rep("Meningioma", ncol(MNG1))
MNG1@meta.data$disease = disease

disease = rep("Meningioma", ncol(MNG2))
MNG2@meta.data$disease = disease

disease = rep("Meningioma", ncol(MNG3))
MNG3@meta.data$disease = disease

disease = rep("Meningioma", ncol(MNG4))
MNG4@meta.data$disease = disease

disease = rep("Meningioma", ncol(MNG5))
MNG5@meta.data$disease = disease

disease = rep("Meningioma", ncol(MNG6))
MNG6@meta.data$disease = disease

MNG3[["percent.mt"]] <- PercentageFeatureSet(MNG3, pattern = "^MT-")
VlnPlot(MNG3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)




########################################################################################################################

####### Load PA datasets
mat = fread("/home/jiin/BT_immune/scPortal_PA/rawdata.txt")

meta = read.table("/home/jiin/BT_immune/scPortal_PA/metadata3_20181216.txt", header=T, sep="\t", as.is=T, row.names=1)
meta=meta[-1,]
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)

##PA = CreateSeuratObject(counts = mat, project = "nowakowski", meta.data=meta)

PA = CreateSeuratObject(counts = mat, project = "scPortal_PA", meta.data=meta)


##### Add metadata
disease = rep("Pilocytic Astrocytoma", ncol(PA))
PA@meta.data$disease = disease

PA[["percent.mt"]] <- PercentageFeatureSet(PA, pattern = "^MT-")
VlnPlot(PA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


########################################################################################################################

###### EPN datasets

mat = fread("/data1/CGL_data1/jiin/BT_data/gse125969_EPN/GSE125969_count_matrix.tsv")

meta = read.table("/data1/CGL_data1/jiin/BT_data/gse125969_EPN/GSE125969_cell_metadata.tsv", header=T, sep="\t", as.is=T, row.names=1)

genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)



EPN = CreateSeuratObject(counts = mat, project = "GSE125969_EPN", meta.data=meta)


#### Add disease
disease = rep("Ependymoma", ncol(EPN))
EPN@meta.data$disease = disease


EPN[["percent.mt"]] <- PercentageFeatureSet(EPN, pattern = "^MT-")
VlnPlot(EPN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

########################################################################################################################



###### MB datasets

mat = fread("/data1/CGL_data1/jiin/BT_data/gse155446_MB/GSE155446_human_raw_counts.csv")

meta = read.csv("/data1/CGL_data1/jiin/BT_data/gse155446_MB/GSE155446_human_cell_metadata.csv", header = T, as.is = T, row.names = 1)

genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)

##tiss = CreateSeuratObject(counts = mat, project = "nowakowski", meta.data=meta)

MB = CreateSeuratObject(counts = mat, project = "gse155446", meta.data=meta)


test = subset(MB, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)


### Add disease
disease = rep("Medulloblastoma", ncol(MB))
MB@meta.data$disease = disease

MB[["percent.mt"]] <- PercentageFeatureSet(MB, pattern = "^MT-")
VlnPlot(MB, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

########################################################################################################################



#### Load the GBM dataset *****
GBM1.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE162631_GBM_sc/gse162631_GBM/R1_N/")
GBM2.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE162631_GBM_sc/gse162631_GBM/R1_T/")
GBM3.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE162631_GBM_sc/gse162631_GBM/R2_N/")
GBM4.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE162631_GBM_sc/gse162631_GBM/R2_T/")
GBM5.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE162631_GBM_sc/gse162631_GBM/R3_N/")
GBM6.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE162631_GBM_sc/gse162631_GBM/R3_T/")
GBM7.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE162631_GBM_sc/gse162631_GBM/R4_N/")
GBM8.data <- Read10X(data.dir = "/data1/CGL_data1/BT_immune/GSE162631_GBM_sc/gse162631_GBM/R4_T/")




#### Create Seurat object
GBM1 <- CreateSeuratObject(counts = GBM1.data, project = "gse162631", min.cells = 3, min.features = 200)
GBM2 <- CreateSeuratObject(counts = GBM2.data, project = "gse162631", min.cells = 3, min.features = 200)
GBM3 <- CreateSeuratObject(counts = GBM3.data, project = "gse162631", min.cells = 3, min.features = 200)
GBM4 <- CreateSeuratObject(counts = GBM4.data, project = "gse162631", min.cells = 3, min.features = 200)
GBM5 <- CreateSeuratObject(counts = GBM5.data, project = "gse162631", min.cells = 3, min.features = 200)
GBM6 <- CreateSeuratObject(counts = GBM6.data, project = "gse162631", min.cells = 3, min.features = 200)
GBM7 <- CreateSeuratObject(counts = GBM7.data, project = "gse162631", min.cells = 3, min.features = 200)
GBM8 <- CreateSeuratObject(counts = GBM8.data, project = "gse162631", min.cells = 3, min.features = 200)


#### Add disease to metadata
disease = rep("Glioblastoma", ncol(GBM1))
GBM1@meta.data$disease = disease

disease = rep("Glioblastoma", ncol(GBM2))
GBM2@meta.data$disease = disease

disease = rep("Glioblastoma", ncol(GBM3))
GBM3@meta.data$disease = disease

disease = rep("Glioblastoma", ncol(GBM4))
GBM4@meta.data$disease = disease

disease = rep("Glioblastoma", ncol(GBM5))
GBM5@meta.data$disease = disease

disease = rep("Glioblastoma", ncol(GBM6))
GBM6@meta.data$disease = disease

disease = rep("Glioblastoma", ncol(GBM7))
GBM7@meta.data$disease = disease

disease = rep("Glioblastoma", ncol(GBM8))
GBM8@meta.data$disease = disease

GBM4[["percent.mt"]] <- PercentageFeatureSet(GBM4, pattern = "^MT-")
VlnPlot(GBM4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



########################################################################################################################


########## Load LGG datasets
mat = fread("/home/jiin/BT_immune/gse89567_LGG_sc/GSE89567_IDH_A_processed_data.txt")

meta = read.table("/home/jiin/BT_immune/gse89567_LGG_sc/IDH_A_cell_type_assignment_portal_v2.txt", header=T, sep="\t", as.is=T, row.names=1)
meta = data.frame(meta[-1,])


genes = mat[,1][[1]]
genes = gsub("'", "", genes)
mat = data.frame(mat[,-1], row.names=genes)

##PA = CreateSeuratObject(counts = mat, project = "nowakowski", meta.data=meta)

LGG = CreateSeuratObject(counts = mat, project = "LGG", meta.data=meta)


#### Add disease to metadata
disease = rep("Low Grade Glioma", ncol(LGG))
LGG@meta.data$disease = disease


LGG[["percent.mt"]] <- PercentageFeatureSet(LGG, pattern = "^MT-")
VlnPlot(LGG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)











rm(MNG1.data,MNG2.data,MNG3.data,MNG4.data,MNG5.data,MNG6.data)
rm(GBM1.data,GBM2.data,GBM3.data,GBM4.data,GBM5.data,GBM6.data,GBM7.data,GBM8.data)
### List
disease.list <- list(MNG1, MNG2, MNG3, MNG4, MNG5, MNG6, PA, EPN, MB, 
                     GBM1, GBM2, GBM3, GBM4, GBM5, GBM6,GBM7, GBM8, LGG)


for (i in disease.list){
  print(i)
}

#### QC and selecting cells for further analysis
disease.list <- lapply(X = disease.list, FUN = function(x) {
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 30)
})



# Normalize and identify variable features for each dataset independently
disease.list <- lapply(X = disease.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = disease.list)


####### Perform integration
immune.anchors <- FindIntegrationAnchors(object.list = disease.list, anchor.features = features) 




# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)



######################################### Integration workflow
# Run the standard workflow for visualization and clustering
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
ElbowPlot(immune.combined, ndims = 30)
#DefaultAssay(immune.combined) <- "RNA"

########################### ***** Elbowplot : quantitative approach
library(ggplot2)
library(ggrepel)
# Determine percent of variation associated with each PC
pct <- immune.combined[["pca"]]@stdev / sum(immune.combined[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()






immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:23)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:23)

rm(disease.list,EPN,GBM1,GBM2,GBM3,GBM4,GBM5,GBM6,GBM7,GBM8,i,
   LGG,mat,MB,meta,MNG1,MNG2,MNG3,MNG4,MNG5,MNG6,PA,disease,genes)

DimPlot(immune.combined, reduction = "umap", label = F,raster = F) +
  font("legend.title", size = 25) +
  font("legend.text", size = 25) + 
  guides(color = guide_legend(override.aes = list(size=10), ncol=1)) + 
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 25),
        title = element_text(size = 20))
  

immune.combined@meta.data$disease.abb <- ifelse(immune.combined@meta.data$disease == "Ependymoma","EPN",
                                                ifelse(immune.combined@meta.data$disease == "Glioblastoma","GBM",
                                                       ifelse(immune.combined@meta.data$disease == "Low Grade Glioma","LGG",
                                                              ifelse(immune.combined@meta.data$disease == "Medulloblastoma","MED",
                                                                     ifelse(immune.combined@meta.data$disease == "Meningioma","MNG","PA")))))
  
DimPlot(immune.combined, reduction = "umap", label = F, split.by = "disease.abb",raster = F,cols = "Paired")+
  guides(color = guide_legend(override.aes = list(size=7), ncol=1) ) +
  labs(title = "") +
  # font("xlab",size = 18) + 
  # font("ylab",size = 25) + 
  # font("xy.text", size = 18) + 
  font("legend.title", size = 21) +
  font("legend.text", size = 21) +
  font("subtitle", size = 40) + 
  theme(plot.subtitle = element_text(size = 20),
        axis.title = element_text(size=20))




immune.combined@meta.data <- immune.combined@meta.data[,-c(6:26)]



# find markers for every cluster compared to all remaining cells, report only the positive
# ones
DefaultAssay(immune.combined) <- "RNA"
BT.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25)
## logfc.threshold = 0.25
dim(BT.markers)
table(BT.markers$cluster)


BT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> Top.BT.markers


DotPlot(immune.combined, features = c("CX3CR1", "P2RY12", 
                                      "APOC1", "CD163", "CD68",
                                      "TOP2A","HDAC2","STMN2",
                                      "CXCR2","IL1R2","FPR2",
                                      "CD3D", "CD3E", "GZMK", 
                                      "VWF","CLDN5","CD34","PECAM1",
                                      "RGS5", "NOTCH3", "TAGLN",
                                      "S100B", "GFAP",
                                     "HLA-DQA1", "HLA-DPB1"
                                      ), 
        col.min = 0.5, dot.scale = 10,cols = c("white","#DF536B")) + RotatedAxis() +
  # font("xlab",size = 18) +
  # font("ylab",size = 25) +
  font("xy.text", size = 20) +
  font("legend.title", size = 25) +
  font("legend.text", size = 25) +
  labs(x = "", y = "") + theme(legend.key.size = unit(3,"line")) + 
  theme(axis.text.y = element_text(size = 30),
        axis.text.x = element_text(size = 20,angle = 90,vjust = .5))

####### Assigning cell type identity to clusters
new.cluster.ids <- c("Microglia", "Macrophage", "Neoplastic cell", "Neoplastic cell", "Microglia",
                     "Neoplastic cell",  "Neoplastic cell", "Neutrophil", "Neoplastic cell", "T cell",
                     "Neoplastic cell", "Macrophage", "Neoplastic cell", "Endothelial cell", "Neoplastic cell",
                     "Mural cell", "Neoplastic cell", "Astrocyte", "Neoplastic cell", "Dendritic cell",
                     "Macrophage", "Neoplastic cell", "Macrophage", "Neoplastic cell", "T cell")

names(new.cluster.ids) <- levels(immune.combined)
immune.combined = RenameIdents(immune.combined, new.cluster.ids)
DimPlot(immune.combined, reduction = 'umap', label = F,raster = F)
immune.combined[["cell.type.annotation"]] = Idents(object = immune.combined)

DimPlot(immune.combined, reduction = "umap", label = F,raster = F,cols = "Paired") +
  font("legend.title", size = 25) +
  font("legend.text", size = 25) + 
  guides(color = guide_legend(override.aes = list(size=10), ncol=1)) + 
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 25),
        title = element_text(size = 20))

dev.off()


Idents(immune.combined) = "disease"
DimPlot(immune.combined)
MNG = subset(immune.combined, idents = "Meningioma")
PA = subset(immune.combined, idents = "Pilocytic Astrocytoma")
EPN = subset(immune.combined, idents = "Ependymoma")
MB = subset(immune.combined, idents = "Medulloblastoma")
GBM = subset(immune.combined, idents = "Glioblastoma")
LGG = subset(immune.combined, idents = "Low Grade Glioma")

Idents(MNG) = "cell.type.annotation"
Idents(PA) = "cell.type.annotation"
Idents(EPN) = "cell.type.annotation"
Idents(MB) = "cell.type.annotation"
Idents(GBM) = "cell.type.annotation"
Idents(LGG) = "cell.type.annotation"


## Stacked bar plot by sample

disease = c("MNG","PA","EPN","MED","GBM","LGG")

dataset_disease = data.frame()

# for (sam in disease){
#   # proportion
  celltype = c("Microglia","Macrophage","Neutrophil","T cell","Endothelial cell","Mural cell",
               "Astrocyte","Dendritic cell")
  
  
  samp = c(rep("BT", length(celltype)))
  
  value = c(ncol(subset(x = immune.combined, idents = "Microglia", subset = disease.abb == "LGG")), 
            ncol(subset(x = immune.combined, idents = "Macrophage", subset = disease.abb == "LGG")),
            ncol(subset(x = immune.combined, idents = "Neutrophil", subset = disease.abb == "LGG")),
            ncol(subset(x = immune.combined, idents = "T cell", subset = disease.abb == "LGG")),
            ncol(subset(x = immune.combined, idents = "Endothelial cell", subset = disease.abb == "LGG")),
            ncol(subset(x = immune.combined, idents = "Mural cell", subset = disease.abb == "LGG")),
            ncol(subset(x = immune.combined, idents = "Astrocyte", subset = disease.abb == "LGG")),
            ncol(subset(x = immune.combined, idents = "Dendritic cell", subset = disease.abb == "LGG")))
  
  
  sam <- data.frame(samp,celltype,value)
  
  sam$percent = sprintf("%1.2f", sam$value/sum(sam$value))
  
  #data_sam$percent = paste0(sprintf("%1.2f", 100 * data$value/sum(data$value)),"%")
  
  sam$Celltype = factor(sam$celltype, levels = celltype)
  dataset_disease = rbind(dataset_disease,sam)

# }
dataset_disease$celltype = dataset_disease$Celltype  
dataset_disease$celltype = factor(dataset_disease$celltype,levels = c("Microglia","Macrophage","Neutrophil",
                                                                      "T cell","Endothelial cell","Mural cell",
                                                                      "Astrocyte","Dendritic cell"))  
########## Plotting
library(ggplot2)
library(dplyr)
ggplot(dataset_disease, aes(fill=celltype, y=value, x=samp)) + 
  geom_bar(position="fill", stat = "identity") + scale_y_continuous(labels = scales::percent) + RotatedAxis()+
  labs(fill = 'Cell type') + ylab('Percent') +  
  #scale_fill_brewer(palette = "Paired") + 
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6")) + 
  scale_x_discrete(limits=c("MNG","PA","EPN","MED","GBM","LGG")) +
  theme_classic() +
  theme(axis.text = element_text(size = 17),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.title.x = element_blank())
