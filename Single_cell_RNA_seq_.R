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

MNG1[["percent.mt"]] <- PercentageFeatureSet(MNG1, pattern = "^MT-")


disease = rep("MNG1", ncol(MNG1))
MNG1@meta.data$sample = disease

disease = rep("MNG2", ncol(MNG2))
MNG2@meta.data$sample = disease

disease = rep("MNG3", ncol(MNG3))
MNG3@meta.data$sample = disease

disease = rep("MNG4", ncol(MNG4))
MNG4@meta.data$sample = disease

disease = rep("MNG5", ncol(MNG5))
MNG5@meta.data$sample = disease

disease = rep("MNG6", ncol(MNG6))
MNG6@meta.data$sample = disease

rm(MNG1.data,MNG2.data,MNG3.data,MNG4.data,MNG5.data,MNG6.data)
########################################################################################################################

####### Load PA datasets
mat = fread("/data1/CGL_data1/jiin/BT_immune/scPortal_PA/rawdata.txt")

meta = read.table("/data1/CGL_data1/jiin/BT_immune/scPortal_PA/metadata3_20181216.txt", header=T, sep="\t", as.is=T, row.names=1)

meta=meta[-1,]
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)

##PA = CreateSeuratObject(counts = mat, project = "nowakowski", meta.data=meta)

PA = CreateSeuratObject(counts = mat, project = "scPortal_PA", meta.data=meta)
PA@meta.data

##### Add metadata
disease = rep("Pilocytic Astrocytoma", ncol(PA))
PA@meta.data$disease = disease

PA@meta.data$sample = PA@meta.data$orig.ident

PA[["percent.mt"]] <- PercentageFeatureSet(PA, pattern = "^MT-")
VlnPlot(PA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


########################################################################################################################

###### EPN datasets

mat = fread("/data1/CGL_data1/jiin/BT_data/gse125969_EPN/GSE125969_count_matrix.tsv")

meta = read.table("/data1/CGL_data1/jiin/BT_data/gse125969_EPN/GSE125969_cell_metadata.tsv", header=T, sep="\t", as.is=T, row.names=1)

genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)

unique(substr(rownames(EPN@meta.data),1,11))

EPN = CreateSeuratObject(counts = mat, project = "GSE125969_EPN", meta.data=meta)

EPN@meta.data
#### Add disease
disease = rep("Ependymoma", ncol(EPN))
EPN@meta.data$disease = disease

EPN@meta.data$sample = substr(rownames(EPN@meta.data),1,11)
unique(EPN@meta.data$sample)

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

MB@meta.data
unique(MB@meta.data$orig.ident)
MB@meta.data$sample = MB@meta.data$orig.ident
MB@meta.data$sample[which(MB@meta.data$sample == "X966.2")] = "X966"
MB@meta.data
unique(MB@meta.data$sample)
MB@meta.data$sample = as.character(MB@meta.data$sample)
### Add disease
disease = rep("Medulloblastoma", ncol(MB))
MB@meta.data$disease = disease

MB[["percent.mt"]] <- PercentageFeatureSet(MB, pattern = "^MT-")
VlnPlot(MB, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

##################################################################################################

GBM1 = fread("/data1/CGL_data1/jiin/BT_immune/GSE103224_RAW_GBM/GSM2758472_PJ017.filtered.matrix.txt")
GBM1 = GBM1[,-1]
GBM1 <- GBM1[!duplicated(GBM1[,1]),]
GBM1 = data.frame(GBM1[,-1], row.names=GBM1[,1][[1]])
GBM1 = CreateSeuratObject(counts = GBM1, project = "GBM")
GBM1@meta.data$disease = rep("Glioblastoma",ncol(GBM1))
GBM1@meta.data$sample = rep("GSM2759472",ncol(GBM1))

GBM2 = fread("/data1/CGL_data1/jiin/BT_immune/GSE103224_RAW_GBM/GSM2758473_PJ018.filtered.matrix.txt")
GBM2 = GBM2[,-1]
GBM2 <- GBM2[!duplicated(GBM2[,1]),]
GBM2 = data.frame(GBM2[,-1], row.names=GBM2[,1][[1]])
GBM2 = CreateSeuratObject(counts = GBM2, project = "GBM")
GBM2@meta.data$disease = rep("Glioblastoma",ncol(GBM2))
GBM2@meta.data$sample = rep("GSM2759473",ncol(GBM2))

GBM3 = fread("/data1/CGL_data1/jiin/BT_immune/GSE103224_RAW_GBM/GSM2758474_PJ025.filtered.matrix.txt")
GBM3 = GBM3[,-1]
GBM3 <- GBM3[!duplicated(GBM3[,1]),]
GBM3 = data.frame(GBM3[,-1], row.names=GBM3[,1][[1]])
GBM3 = CreateSeuratObject(counts = GBM3, project = "GBM")
GBM3@meta.data$disease = rep("Glioblastoma",ncol(GBM3))
GBM3@meta.data$sample = rep("GSM2759474",ncol(GBM3))

GBM4 = fread("/data1/CGL_data1/jiin/BT_immune/GSE103224_RAW_GBM/GSM2758476_PJ032.filtered.matrix.txt")
GBM4 = GBM4[,-1]
GBM4 <- GBM4[!duplicated(GBM4[,1]),]
GBM4 = data.frame(GBM4[,-1], row.names=GBM4[,1][[1]])
GBM4 = CreateSeuratObject(counts = GBM4, project = "GBM")
GBM4@meta.data$disease = rep("Glioblastoma",ncol(GBM4))
GBM4@meta.data$sample = rep("GSM2759476",ncol(GBM4))

GBM5 = fread("/data1/CGL_data1/jiin/BT_immune/GSE103224_RAW_GBM/GSM2758477_PJ035.filtered.matrix.txt")
GBM5 = GBM5[,-1]
GBM5 <- GBM5[!duplicated(GBM5[,1]),]
GBM5 = data.frame(GBM5[,-1], row.names=GBM5[,1][[1]])
GBM5 = CreateSeuratObject(counts = GBM5, project = "GBM")
GBM5@meta.data$disease = rep("Glioblastoma",ncol(GBM5))
GBM5@meta.data$sample = rep("GSM2759477",ncol(GBM5))

GBM6 = fread("/data1/CGL_data1/jiin/BT_immune/GSE103224_RAW_GBM/GSM2940098_PJ048.filtered.matrix.txt")
GBM6 = GBM6[,-1]
GBM6 <- GBM6[!duplicated(GBM6[,1]),]
GBM6 = data.frame(GBM6[,-1], row.names=GBM6[,1][[1]])
GBM6 = CreateSeuratObject(counts = GBM6, project = "GBM")
GBM6@meta.data$disease = rep("Glioblastoma",ncol(GBM6))
GBM6@meta.data$sample = rep("GSM2940098",ncol(GBM6))

########## Load LGG datasets
mat = fread("/data/CGL_data1/jiin/BT_immune/gse89567_LGG_sc/GSE89567_IDH_A_processed_data.txt")

meta = read.table("/data/CGL_data1/jiin/BT_immune/gse89567_LGG_sc/IDH_A_cell_type_assignment_portal_v2.txt", header=T, sep="\t", as.is=T, row.names=1)
meta = data.frame(meta[-1,])


genes = mat[,1][[1]]
genes = gsub("'", "", genes)
mat = data.frame(mat[,-1], row.names=genes)

##PA = CreateSeuratObject(counts = mat, project = "nowakowski", meta.data=meta)

LGG = CreateSeuratObject(counts = mat, project = "LGG", meta.data=meta)
LGG@meta.data
LGG@meta.data$sample = toupper(substr(rownames(LGG@meta.data),1,6))
LGG@meta.data$sample = ifelse(LGG@meta.data$sample == "MGH42_","MGH42",
                              ifelse(LGG@meta.data$sample == "MGH61_","MGH61",
                                     ifelse(LGG@meta.data$sample == "MGH43_","MGH43",
                                            ifelse(LGG@meta.data$sample == "MGH45_","MGH45",
                                                   ifelse(LGG@meta.data$sample == "MGH44_","MGH44",
                                                          ifelse(LGG@meta.data$sample == "MGH56_","MGH56",
                                                                 ifelse(LGG@meta.data$sample == "MGH64_" | LGG@meta.data$sample == "MGH64.","MGH64",
                                                                        ifelse(LGG@meta.data$sample == "MGH57_" | LGG@meta.data$sample == "X57_P1","MGH57",LGG@meta.data$sample))))))))
unique(LGG@meta.data$sample)
#### Add disease to metadata
disease = rep("Low Grade Glioma", ncol(LGG))
LGG@meta.data$disease = disease


LGG[["percent.mt"]] <- PercentageFeatureSet(LGG, pattern = "^MT-")
VlnPlot(LGG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

LGG@meta.data









rm(MNG1.data,MNG2.data,MNG3.data,MNG4.data,MNG5.data,MNG6.data)
rm(GBM1.data,GBM2.data,GBM3.data,GBM4.data,GBM5.data,GBM6.data,GBM7.data,GBM8.data)
### List
disease.list <- list(MNG1, MNG2, MNG3, MNG4, MNG5, MNG6, PA, EPN, MB, 
                     GBM1, GBM2, GBM3, GBM4, GBM5, GBM6, LGG)


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
#1/23 4:28 PM


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

rm(disease.list,EPN,GBM1,GBM2,GBM3,GBM4,GBM5,GBM6,i,
   LGG,mat,MB,meta,MNG1,MNG2,MNG3,MNG4,MNG5,MNG6,PA,disease,genes)

DimPlot(immune.combined, reduction = "umap", label = F,raster = F) +
  font("legend.title", size = 25) +
  font("legend.text", size = 25) + 
  guides(color = guide_legend(override.aes = list(size=10), ncol=3)) + 
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 25),
        title = element_text(size = 20))

DimPlot(subset(immune.combined,subset = disease.abb == "GBM"), reduction = "umap", label = F,raster = F,cols = "Paired") +
  ggtitle("GBM")+
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

Idents(immune.combined) <- "disease.abb"


immune.combined@meta.data <- immune.combined@meta.data[,-c(7:27)]

immune.combined@meta.data$cell.type.annotation.neoplastic = ifelse(immune.combined@meta.data$cell.type.annotation == "Neoplastic cell",
                                                                   "Neoplastic cell","Other cell")
Idents(immune.combined) <- "cell.type.annotation.neoplastic"
DimPlot(immune.combined)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
DefaultAssay(immune.combined) <- "RNA"
BT.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25)
## logfc.threshold = 0.25
dim(BT.markers)
table(BT.markers$cluster)


BT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) -> Top.BT.markers


DotPlot(immune.combined, features = c("OLIG1","OLIG2",
                                      "TOP2A","HDAC2",
                                      "CX3CR1","P2RY12",
                                      "CD3D","CD3E","GZMK",
                                      "RGS5","NOTCH3","TAGLN",
                                      "VWF","CLDN5","CD34","PECAM1",
                                      "CD68","CD14","APOC1",
                                      "HLA-DQA1","HLA-DPB1",
                                      "GFAP","S100B",
                                      "CD79A","CD79B","MS4A1"
), 
col.min = 0.2, dot.scale = 10,cols = c("white","#DF536B")) + RotatedAxis() +
  # font("xlab",size = 18) +
  # font("ylab",size = 25) +
  font("xy.text", size = 20) +
  font("legend.title", size = 25) +
  font("legend.text", size = 25) +
  labs(x = "", y = "") + theme(legend.key.size = unit(3,"line")) + 
  theme(axis.text.y = element_text(size = 30),
        axis.text.x = element_text(size = 20,angle = 90,vjust = .5))


## featureplot
FeaturePlot(immune.combined,features = c("HDAC2"), raster = F,cols = c("lightgrey","red"),max.cutoff = 3,min.cutoff = 0.7)



####### Assigning cell type identity to clusters
new.cluster.ids <- c("OPC","Neoplastic cell","Microglia","Neoplastic cell","T cell",
                     "Neoplastic cell","Mural cell","Neoplastic cell","Neoplastic cell","Neoplastic cell",
                     "Neoplastic cell","Neoplastic cell","Neoplastic cell","Endothelial cell","Macrophage",
                     "Neoplastic cell","Neoplastic cell","Dendritic cell","Astrocyte","Macrophage",
                     "Neoplastic cell","Neoplastic cell","B cell","Astrocyte")

names(new.cluster.ids) <- levels(immune.combined)
immune.combined = RenameIdents(immune.combined, new.cluster.ids)
DimPlot(immune.combined, reduction = 'umap', label = F,raster = F)
immune.combined[["cell.type.annotation.new"]] = Idents(object = immune.combined)

DimPlot(immune.combined, reduction = "umap", label = F,raster = F,cols = "Paired") +
  font("legend.title", size = 25) +
  font("legend.text", size = 25) + 
  guides(color = guide_legend(override.aes = list(size=10), ncol=1)) + 
  theme(axis.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        title = element_text(size = 20))

dev.off()

##################################### ploting(percent)
Idents(immune.combined) = "cell.type.annotation"
celltype = c(levels(immune.combined))

samp = c(rep("ALL", length(celltype)))

value = c()


for (i in celltype){
  value = append(value, values = ncol(subset(immune.combined, idents = i)))
}


sam <- data.frame(samp,celltype,value)
sam <- sam[!(sam$celltype == "Neoplastic cell" ), ]

#sam$percent = sprintf("%1.2f", sam$value/sum(sam$value))
sam$percent <- (sam$value / sum(sam$value)) * 100

sam$Celltype = factor(sam$celltype, levels = celltype)


brewer.pal(n=20,name="Paired")

label <- paste(round(sam$percent,1),"%")

library(ggpubr)
ggbarplot(sam,"Celltype","percent",fill = "Celltype",
          palette = c("#A6CEE3", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")) + 
  labs(x = "",y = "Percent") + ylim(0,50) + 
  geom_text(aes(label = label),vjust=-0.2,size=5) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20))



## Stacked bar plot by sample
immune.combined.non <- subset(immune.combined, idents = c("OPC","Microglia","T cell","Mural cell",
                                                          "Endothelial cell","Macrophage","Dendritic cell","Astrocyte","B cell"))
disease = c(unique(immune.combined.non@meta.data$disease.abb))

dataset_disease = data.frame()

for (sam in disease){
  # proportion
  celltype = c("OPC","Microglia","T cell","Mural cell",
               "Endothelial cell","Macrophage","Dendritic cell","Astrocyte","B cell")
  
  
  samp = c(rep(sam, length(celltype)))
  
  value = c(tryCatch(ncol(subset(x = immune.combined.non, idents = "OPC", subset = disease.abb == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Microglia", subset = disease.abb == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "T cell", subset = disease.abb == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Mural cell", subset = disease.abb == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Endothelial cell", subset = disease.abb == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Macrophage", subset = disease.abb == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Dendritic cell", subset = disease.abb == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Astrocyte", subset = disease.abb == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "B cell", subset = disease.abb == sam)),error = function(e){print(0)}))
  
  
  
  sam <- data.frame(samp,celltype,value)
  
  sam$percent = sprintf("%1.2f", sam$value/sum(sam$value))
  
  #data_sam$percent = paste0(sprintf("%1.2f", 100 * data$value/sum(data$value)),"%")
  
  sam$Celltype = factor(sam$celltype, levels = celltype)
  dataset_disease = rbind(dataset_disease,sam)
  
}
dataset_disease$celltype = dataset_disease$Celltype  
dataset_disease$celltype = factor(dataset_disease$celltype,levels = c("OPC","Microglia","T cell","Mural cell",
                                                                      "Endothelial cell","Macrophage","Dendritic cell","Astrocyte","B cell"))  
########## Plotting
library(ggplot2)
library(dplyr)
ggplot(dataset_disease, aes(fill=celltype, y=value, x=samp)) + 
  geom_bar(position="fill", stat = "identity") + scale_y_continuous(labels = scales::percent) + RotatedAxis()+
  labs(fill = 'Cell type') + ylab('Percent') +  
  #scale_fill_brewer(palette = "Paired") + 
  scale_fill_manual(values = c("#A6CEE3", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")) + 
  ggtitle("")+
  scale_x_discrete(limits=c("MNG","PA","EPN","MED","GBM","LGG")) +
  theme_classic() +
  theme(title = element_text(size = 25),
        axis.text = element_text(size = 17),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.title.x = element_blank())


## Stacked bar plot (immune cells)
immune.combined.non <- subset(immune.combined, idents = c("Microglia","T cell","Macrophage","Dendritic cell","B cell"))
disease = c(unique(immune.combined.non@meta.data$disease.abb))

dataset_disease = data.frame()

for (sam in disease){
  # proportion
  celltype = c("Microglia","Macrophage","T cell","B cell","Dendritic cell")
  
  
  samp = c(rep(sam, length(celltype)))
  
  value = c(tryCatch(ncol(subset(x = immune.combined.non, idents = "Microglia", subset = disease.abb == sam)),error = function(e){print(0)}), 
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Macrophage", subset = disease.abb == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "T cell", subset = disease.abb == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "B cell", subset = disease.abb == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Dendritic cell", subset = disease.abb == sam)),error = function(e){print(0)}))
  
  
  
  sam <- data.frame(samp,celltype,value)
  
  sam$percent = sprintf("%1.2f", sam$value/sum(sam$value))
  
  #data_sam$percent = paste0(sprintf("%1.2f", 100 * data$value/sum(data$value)),"%")
  
  sam$Celltype = factor(sam$celltype, levels = celltype)
  dataset_disease = rbind(dataset_disease,sam)
  
}
dataset_disease$celltype = dataset_disease$Celltype  
dataset_disease$celltype = factor(dataset_disease$celltype,levels = c("Microglia","Macrophage","T cell","B cell","Dendritic cell"))  
########## Plotting
library(ggplot2)
library(dplyr)
ggplot(dataset_disease, aes(fill=celltype, y=value, x=samp)) + 
  geom_bar(position="fill", stat = "identity") + scale_y_continuous(labels = scales::percent) + RotatedAxis()+
  labs(fill = 'Cell type') + ylab('Percent') +  
  #scale_fill_brewer(palette = "Paired") + 
  scale_fill_manual(values = c("#B2DF8A","#FF7F00", "#FB9A99","#FFFF99","#CAB2D6")) + 
  ggtitle("")+
  scale_x_discrete(limits=c("MNG","PA","EPN","MED","GBM","LGG")) +
  theme_classic() +
  theme(title = element_text(size = 25),
        axis.text = element_text(size = 17),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.title.x = element_blank())


## disease
MNG = subset(immune.combined,subset = disease.abb == "PA")
Idents(MNG) <- "sample"
DimPlot(MNG) + ggtitle("PA") 



LGG = subset(immune.combined,subset = disease.abb == "LGG")
LGG@meta.data
LGG@meta.data$sample.new = toupper(substr(rownames(LGG@meta.data),1,6))
LGG@meta.data$sample.new = ifelse(LGG@meta.data$sample.new == "MGH42_","MGH42",
                                  ifelse(LGG@meta.data$sample.new == "MGH61_","MGH61",
                                         ifelse(LGG@meta.data$sample.new == "MGH43_","MGH43",
                                                ifelse(LGG@meta.data$sample.new == "MGH45_","MGH45",
                                                       ifelse(LGG@meta.data$sample.new == "MGH44_","MGH44",
                                                              ifelse(LGG@meta.data$sample.new == "MGH56_","MGH56",
                                                                     ifelse(LGG@meta.data$sample.new == "MGH64_" | LGG@meta.data$sample.new == "MGH64.","MGH64",
                                                                            ifelse(LGG@meta.data$sample.new == "MGH57_" | LGG@meta.data$sample.new == "X57_P1","MGH57",LGG@meta.data$sample.new))))))))
unique(LGG@meta.data$sample.new)

## Stacked bar plot by sample (Cell types)

immune.combined.non <- subset(immune.combined,idents = "Other cell")
DimPlot(immune.combined.non)

Idents(immune.combined.non) <- "cell.type.annotation"
DimPlot(immune.combined.non)


disease = c(unique(immune.combined.non@meta.data$sample))


dataset_disease = data.frame()

for (sam in disease){
  # proportion
  celltype = c("OPC","Microglia","T cell","Mural cell",
               "Endothelial cell","Macrophage","Dendritic cell","Astrocyte","B cell")
  
  
  samp = c(rep(sam, length(celltype)))
  
  value = c(tryCatch(ncol(subset(x = immune.combined.non, idents = "OPC", subset = sample == sam)),error = function(e){print(0)}), 
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Microglia", subset = sample == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "T cell", subset = sample == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Mural cell", subset = sample == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Endothelial cell", subset = sample == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Macrophage", subset = sample == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Dendritic cell", subset = sample == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "Astrocyte", subset = sample == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined.non, idents = "B cell", subset = sample == sam)),error = function(e){print(0)}))
  
  
  
  sam <- data.frame(samp,celltype,value)
  
  sam$percent = sprintf("%1.2f", sam$value/sum(sam$value))
  
  #data_sam$percent = paste0(sprintf("%1.2f", 100 * data$value/sum(data$value)),"%")
  
  sam$Celltype = factor(sam$celltype, levels = celltype)
  dataset_disease = rbind(dataset_disease,sam)
  
}
dataset_disease$celltype = dataset_disease$Celltype  
dataset_disease$celltype = factor(dataset_disease$celltype,levels = c("OPC","Microglia","T cell","Mural cell",
                                                                      "Endothelial cell","Macrophage","Dendritic cell","Astrocyte","B cell"))  


########## Plotting
library(ggplot2)
library(dplyr)
ggplot(dataset_disease, aes(fill=celltype, y=value, x=samp)) + 
  geom_bar(position="fill", stat = "identity") + scale_y_continuous(labels = scales::percent) + RotatedAxis()+
  labs(fill = 'Cell type') + ylab('Percent') +  
  #scale_fill_brewer(palette = "Paired") + 
  #scale_fill_manual(values = c("#B2DF8A","lightblue4")) + 
  #ggtitle("LGG")+
  # scale_x_discrete(limits=c("MNG","PA","EPN","MED","GBM","LGG")) +
  theme_classic() +
  theme(title = element_text(size = 25),
        axis.text = element_text(size = 17),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.title.x = element_blank())

## Boxplots
dataset_disease$disease <- c(rep("MNG",54),rep("PA",54),rep("EPN",234),rep("MED",261),rep("GBM",54),rep("LGG",90))
dataset_disease$percent <- as.numeric(dataset_disease$percent)
dataset_disease$per <- dataset_disease$percent * 100
dendritic <- dataset_disease[dataset_disease$celltype == "B cell",]
dendritic$disease = factor(dendritic$disease,levels = c("MNG","PA","EPN","MED","GBM","LGG"))

ggboxplot(dendritic, x = "disease", y = "per",
          color = "black", palette =c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33"),
          , x.text.angle = 30, xlab = "Disease", ylab = "Proportion (B cell)", outlier.size = 0.05,fill = "disease"
)+
  scale_y_continuous(label = function(x) paste0(x,"%"))+
  stat_compare_means(label.x = 1.4,label.y = 7.5,size = 4) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        legend.position = "none") 


## Stacked bar plot by sample (Neoplastic cell vs Other cells)
immune.combined@meta.data$neoplas <- ifelse(immune.combined@meta.data$cell.type.annotation.new == "Neoplastic cell","Neoplastic cell","Other cell")
Idents(immune.combined) <- "neoplas"
DimPlot(immune.combined)
EPN <- subset(immune.combined,subset = disease.abb == "PA")
DimPlot(EPN)

disease = c(unique(EPN@meta.data$sample))


dataset_disease = data.frame()

for (sam in disease){
  # proportion
  celltype = c("Neoplastic cell", "Other cell")
  
  
  samp = c(rep(sam, length(celltype)))
  
  value = c(tryCatch(ncol(subset(x = immune.combined, idents = "Neoplastic cell", subset = sample == sam)),error = function(e){print(0)}),
            tryCatch(ncol(subset(x = immune.combined, idents = "Other cell", subset = sample == sam)),error = function(e){print(0)}))
  
  
  sam <- data.frame(samp,celltype,value)
  
  sam$percent = sprintf("%1.2f", sam$value/sum(sam$value))
  
  #data_sam$percent = paste0(sprintf("%1.2f", 100 * data$value/sum(data$value)),"%")
  
  sam$Celltype = factor(sam$celltype, levels = celltype)
  dataset_disease = rbind(dataset_disease,sam)
  
}
dataset_disease$celltype = dataset_disease$Celltype  
dataset_disease$celltype = factor(dataset_disease$celltype,levels = c("Neoplastic cell", "Other cell"))  
########## Plotting
library(ggplot2)
library(dplyr)
ggplot(dataset_disease, aes(fill=celltype, y=value, x=samp)) + 
  geom_bar(position="fill", stat = "identity") + scale_y_continuous(labels = scales::percent) + RotatedAxis()+
  labs(fill = 'Cell type') + ylab('Percent') +  
  #scale_fill_brewer(palette = "Paired") + 
  scale_fill_manual(values = c("#B2DF8A","lightblue4")) + 
  ggtitle("PA")+
  # scale_x_discrete(limits=c("MNG","PA","EPN","MED","GBM","LGG")) +
  theme_classic() +
  theme(title = element_text(size = 25),
        axis.text = element_text(size = 17),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.title.x = element_blank())