## Sample information (supplementary)
library(readxl)
sample_info <- read_excel("/Users/knu_cgl3/Desktop/all/final/Supplementary_files_revise/Supplementary_table.xlsx",
                          sheet = 'TableS1')
colnames(sample_info) <- sample_info[1,]
sample_info <- sample_info[-1,]
sample_info <- sample_info[,-c(4,5,6)]
sample_info$`Accession number`[c(11,12)] <- "GSE4290"
sample_info$`Accession number`[c(14)] <- "GSE44971"
sample_info$`Accession number`[c(16,17,18,19)] <- "GSE50161"
sample_info <- sample_info[-c(12,14,19),]
sample_info$`Tumor type` <- factor(sample_info$`Tumor type`,
                                   levels = c("MNG","PA","EPN","MED","GBM","LGG"))
sample_info$Sample <- as.numeric(sample_info$Sample)

library(ggpubr)
ggbarplot(sample_info,x='Accession number',y='Sample',color = "Tumor type",
          palette = "Set1",fill = 'Tumor type') +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 20),
          axis.text.x = element_text(size = 20,angle = 45,vjust = 1,hjust = 1),
          axis.title = element_text(size = 20),
          axis.text.y = element_text(size = 20))

# Heatmap
library(ggplot2)
library(gplots)
library(pheatmap)
library(heatmap.plus)
y = read.table("/Users/knu_cgl3/Desktop/all/CIBERSORT_result/immune_expression.txt",header = T, sep = "\t", quote = "")

s = y

y = as.matrix(y)
s = as.matrix(s)

breaks = seq(-2, 2, 0.001)
col = colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

dis_info = unlist(lapply(colnames(y),function(x){
  if(grepl("MNG",x)) "MNG" 
  else if(grepl('PA',x)) "PA" 
  else if(grepl('EPN',x)) "EPN"
  else if(grepl('MED',x)) "MED" 
  else if(grepl('GBM', x)) "GBM" 
  else if(grepl('LGG', x)) "LGG"
  
}))


annotation <- data.frame(Disease = dis_info)

rownames(annotation) = colnames(y)

annotation_colors = list(Disease = c(MNG = "#E41A1C", PA = "#377EB8", EPN = "#4DAF4A", MED = "#984EA3", 
                                     GBM = "#FF7F00", LGG = "#FFFF33"))
par(mar = c(4,4,2,2))
rownames(y) = gsub("\\."," ",rownames(y))
pheatmap(y, color = col, breaks = breaks, scale = "row", cluster_cols = F, show_colnames = F,
         annotation_col = annotation, annotation_colors = annotation_colors, fontsize = 15)


# Violin & Box plots
library(ggpubr)

kk = read.table("/Users/knu_cgl3/Desktop/all/CIBERSORT_result/immune_expression.txt",header = T, sep = "\t")


B.cells.naive = kk[1,]
B.cells.memory = kk[2,]
Plasma.cells = kk[3,]
T.cells.CD8 = kk[4,]
T.cells.CD4.memory.resting = kk[5,]
T.cells.follicular.helper = kk[6,]
T.cells.regulatory.Tregs = kk[7,]
T.cells.gamma.delta = kk[8,]
NK.cells.resting = kk[9,]
NK.cells.activated = kk[10,]
Monocytes = kk[11,]
Macrophages.M0 = kk[12,]
Macrophages.M1 = kk[13,]
Macrophages.M2 = kk[14,]
Dendritic.cells.resting = kk[15,]
Dendritic.cells.activated = kk[16,]
Mast.cells.resting = kk[17,]
Mast.cells.activated = kk[18,]
Eosinophils = kk[19,]
Neutrophils = kk[20,]
#Total.Immune.Infiltration = kk[21,]


Value = c(B.cells.naive, B.cells.memory, Plasma.cells, T.cells.CD8, T.cells.CD4.memory.resting, T.cells.follicular.helper,
          T.cells.regulatory.Tregs, T.cells.gamma.delta, NK.cells.resting, NK.cells.activated, Monocytes, Macrophages.M0,
          Macrophages.M1, Macrophages.M2, Dendritic.cells.resting, Dendritic.cells.activated, Mast.cells.resting,
          Mast.cells.activated, Eosinophils, Neutrophils) #Total.Immune.Infiltration

cluster = c(rep(c("MNG"), times = 68), rep(c("PA"), times = 152), rep(c("EPN"), times = 130), rep(c("MED"), times = 160), rep(c("GBM"), times = 227), rep(c("LGG"), times = 76))
Cluster = rep(cluster, times = 20)

Cell.type = c(rep(c("B.cells.naive"), times = 813), rep(c("B.cells.memory"), times = 813), rep(c("Plasma.cells"), times = 813), rep(c("T.cells.CD8"), times = 813),
              rep(c("T.cells.CD4.memory.resting"), times = 813), rep(c("T.cells.follicular.helper"), times = 813), rep(c("T.cells.regulatory.Tregs"), times = 813), 
              rep(c("T.cells.gamma.delta"), times = 813), rep(c("NK.cells.resting"), times = 813), rep(c("NK.cells.activated"), times = 813),
              rep(c("Monocytes"), times = 813), rep(c("Macrophages.M0"), times = 813), rep(c("Macrophages.M1"), times = 813), rep(c("Macrophages.M2"), times = 813),
              rep(c("Dendritic.cells.resting"), times = 813), rep(c("Dendritic.cells.activated"), times = 813), rep(c("Mast.cells.resting"), times = 813),
              rep(c("Mast.cells.activated"), times = 813), rep(c("Eosinophils"), times = 813), rep(c("Neutrophils"), times = 813)
) #rep(c("Total.Immune.Infiltration"), times = 813)

input = cbind(Value, Cluster, Cell.type)

input1 = as.data.frame(input, stringsAsFactors = True)
input1$Value = as.numeric(input1$Value)
input1$Cluster = factor(input1$Cluster, levels = c("MNG", "PA", "EPN", "MED", "GBM", "LGG"))
input1$Cell.type = factor(input1$Cell.type, levels = c("B.cells.naive", "B.cells.memory","Plasma.cells", "T.cells.CD8", "T.cells.CD4.memory.resting",
                                                       "T.cells.follicular.helper", "T.cells.regulatory.Tregs", "T.cells.gamma.delta", "NK.cells.resting",
                                                       "NK.cells.activated", "Monocytes", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2",
                                                       "Dendritic.cells.resting", "Dendritic.cells.activated", "Mast.cells.resting",
                                                       "Mast.cells.activated", "Eosinophils", "Neutrophils")) #, "Total.Immune.Infiltration"



input1$Cell.type <- gsub("\\."," ",input1$Cell.type)

p1 <- ggboxplot(input1, x = "Cell.type", y = "Value",
                color = "black", fill = "Cluster",palette =c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "Goldenrod"),
                , x.text.angle = 45, xlab = "", ylab = "Abundance", outlier.size = 0.05,
) + 
  font("xlab",size = 10) + 
  font("ylab",size = 25) + 
  font("xy.text", size = 15) + 
  font("legend.title", size = 15) + 
  font("legend.text", size = 20) +
  theme(legend.title = element_blank(),
        legend.position = "right")
p1



p1 + stat_compare_means(aes(group = Cluster), label = "p.signif", size = 5,label.y = 0.32)

######################################
input2 = input1[
  input1$Cell.type == "Macrophages.M2" | 
    input1$Cell.type == "T.cells.CD8" | 
    input1$Cell.type == "T.cells.CD4.memory.resting",]

p1 <- ggviolin(input2, x = "Cell.type", y = "Value",
               color = "Cluster",,palette =c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "Goldenrod"),
               , x.text.angle = 30, xlab = "", ylab = "Abundance", outlier.size = 0.05,
               add = "boxplot") + 
  scale_x_discrete(labels = c("CD8 T cells","CD4 T cells","M2 macrophages")) + 
  font("xlab",size = 20) + 
  font("ylab",size = 25) + 
  font("xy.text", size = 15) + 
  font("legend.title", size = 15) + 
  font("legend.text", size = 20) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(size = 18)) + 
  ylim(0,0.3)
p1



p1 + stat_compare_means(aes(group = Cluster), label = "p.signif", size = 5,label.y = 0.32)

ggviolin(input2, x = "Cluster",y = "Value",
         color = "Cluster",fill = "Cluster",palette =c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "Goldenrod"),
         add = "boxplot", add.params = list(fill = "white"),facet.by = "Cell.type",
         ylab = "Abundance") + 
  ylim(0,0.3) + 
  #stat_compare_means(aes(group = Cluster),size = 6,label.x = 1.7,label.y = 0.29) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(size = 18,angle = 30,vjust = .9,hjust = .9),
        axis.text.y = element_text(size = 15),
        legend.position = "none",
        text = element_text(size = 23))


##########
input2 = input1[input1$Cell.type == "Macrophages.M1",]

ggviolin(input2,x = "Cluster",y = "Value", color = "Cluster", fill = "Cluster",
         palette =c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "Goldenrod"),
         xlab = "",ylab = "Macrophage.M1 signatures",outlier.size = 0.05,
         add = "boxplot",add.params = list(fill="white")) + 
  theme(legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 20)) + 
  stat_compare_means(size = 6,label.y = 0.0,label.x = 1.3) 
#ylim(0,0.35)

##########
input2 = input1[input1$Cell.type == "T.cells.CD8",]

ggviolin(input2,x = "Cluster",y = "Value", color = "black", fill = "Cluster",
         palette =c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "Goldenrod"),
         xlab = "",ylab = "T.cells.CD8 signatures",outlier.size = 0.05,
         add = "boxplot",add.params = list(fill="white")) + 
  theme(legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 20)) + 
  stat_compare_means(size = 6,label.y = 0.3,label.x = 1.3) 
#ylim(0,0.35)

##########
input2 = input1[input1$Cell.type == "T.cells.CD4.memory.resting",]

ggviolin(input2,x = "Cluster",y = "Value", color = "black", fill = "Cluster",
         palette =c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "Goldenrod"),
         xlab = "",ylab = "T.cells.CD4 signatures",outlier.size = 0.05,
         add = "boxplot",add.params = list(fill="white")) + 
  theme(legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 20)) + 
  stat_compare_means(size = 6,label.y = 0.3,label.x = 1.8) 
#ylim(0,0.35)


########## Total immune infiltration
kk = read.table("/Users/knu_cgl3/Desktop/all/CIBERSORT_result/immune_cells.txt",header = T, sep = "\t")


total_immune = c(kk[,27])

cluster = c(rep(c("MNG"), times = 68), rep(c("PA"), times = 152), rep(c("EPN"), times = 130), rep(c("MED"), times = 160), rep(c("GBM"), times = 227), rep(c("LGG"), times = 76))

immune_input = cbind(total_immune, cluster)

immune_input = as.data.frame(immune_input)

immune_input$total_immune = as.numeric(immune_input$total_immune)
immune_input$cluster = factor(immune_input$cluster, levels = c("MNG", "PA", "EPN", "MED", "GBM", "LGG"))

p1 <- ggboxplot(immune_input, x = "cluster", y = "total_immune",
                color = "black", palette =c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "Goldenrod"),
                , x.text.angle = 0, xlab = "", ylab = "Total immune infiltration", outlier.size = 0.1, 
                fill = "cluster") + 
  font("xlab",size = 13) + 
  font("ylab",size = 25) + 
  font("xy.text", size = 17) + 
  font("legend.title", size = 21) + 
  font("legend.text", size = 21) + 
  theme(legend.title = element_blank(),
        legend.position = "none")
p1



p1 + stat_compare_means(size = 6, label.y = 0.97,label.x = 1.5)


p1 <- ggviolin(immune_input, x = "cluster", y = "total_immune",
               color = "cluster", palette =c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "Goldenrod"),
               , x.text.angle = 0, xlab = "", ylab = "Total immune infiltration", outlier.size = 0.1, 
               fill = "cluster",add = "boxplot",add.params = list(fill="white")) + 
  font("xlab",size = 13) + 
  font("ylab",size = 25) + 
  font("xy.text", size = 17) + 
  font("legend.title", size = 21) + 
  font("legend.text", size = 21) + 
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylim(0.75,0.95)
p1



p1 + stat_compare_means(size = 6, label.y = 0.76,label.x = 1.2)









