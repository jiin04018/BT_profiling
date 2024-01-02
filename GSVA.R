gene_matrix = read.table("/Users/knu_cgl3/Desktop/all/mixturefile.txt",header = T, sep = "\t")

DC = gene_matrix[gene_matrix$gene_symbols == "HLA-DQA1" | gene_matrix$gene_symbols == "HLA-DPB1",]

T_cell = gene_matrix[gene_matrix$gene_symbols == "CD3D" | gene_matrix$gene_symbols == "CD3E" |
                       gene_matrix$gene_symbols == "GZMK",]

Mural_cell = gene_matrix[gene_matrix$gene_symbols == "RGS5" | gene_matrix$gene_symbols == "NOTCH3" |
                           gene_matrix$gene_symbols == "TAGLN",]

Endothelial_cell = gene_matrix[gene_matrix$gene_symbols == "VWF" | gene_matrix$gene_symbols == "CLDN5" |
                                 gene_matrix$gene_symbols == "CD34" | gene_matrix$gene_symbols == "PECAM1",]

rownames(DC) = DC[,1]
DC = DC[,-1]

rownames(T_cell) = T_cell[,1]
T_cell = T_cell[,-1]

rownames(Mural_cell) = Mural_cell[,1]
Mural_cell = Mural_cell[,-1]

rownames(Endothelial_cell) = Endothelial_cell[,1]
Endothelial_cell = Endothelial_cell[,-1]

Total = rbind(DC, T_cell, Mural_cell, Endothelial_cell)

######################################## GSVA
library(GSVA)
library(ggpubr)

gs <- as.list(sample(10:100, size=4, replace=TRUE))
gs[[1]] = rownames(Total)[1:2]
gs[[2]] = rownames(Total)[3:5]
gs[[3]] = rownames(Total)[6:8]
gs[[4]] = rownames(Total)[9:12]

names(gs) <- c("DC.markers", "T.cells.markers", "Mural.cells.markers", "Endothelial.cells.markers")

Total = as.matrix(Total)
gsva.es <- gsva(Total, gs, verbose=FALSE) #, method = "ssgsea"
dim(gsva.es)

gsva.es[1:4, 1:5]

Total = gsva.es

DC.markers = Total[1,]
T.cells.markers = Total[2,]
Mural.cells.markers = Total[3,]
Endothelial.cells.markers = Total[4,]



Value = c(DC.markers, T.cells.markers, Mural.cells.markers, Endothelial.cells.markers)


cluster = c(rep(c("MNG"), times = 68), rep(c("PA"), times = 152), rep(c("EPN"), times = 130),
            rep(c("MED"), times = 160),rep(c("GBM"), times = 227),rep(c("LGG"), times = 76))
Cluster = rep(cluster, times = 4)

Class = c(rep(c("DC.markers"), times = 813), rep(c("T.cells.markers"), times = 813), 
          rep(c("Mural.cells.markers"), times = 813), rep(c("Endothelial.cells.markers"), times = 813))

input = cbind(Value, Cluster, Class)


input1 = as.data.frame(input)
input1$Value = as.numeric(input1$Value)
input1$Cluster = factor(input1$Cluster, levels = c("MNG", "PA", "EPN", "MED", "GBM" ,"LGG"))
input1$Class = factor(input1$Class, levels = c("DC.markers", "T.cells.markers",
                                               "Mural.cells.markers", "Endothelial.cells.markers")) 

##### Mural cell
mural_input = input1[1627:2439,]

p1 <- ggboxplot(mural_input, x = "Cluster", y = "Value",
                color = "black", palette =c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "Goldenrod"),
                ,  xlab = "", ylab = "GSVA score", outlier.size = 0.05,
                title = "Mural cell markers",fill = "Cluster") + 
  ylim(-1,1.3)+
  font("xlab",size = 18) + 
  font("ylab",size = 25) + 
  font("xy.text", size = 20) + 
  font("legend.title", size = 21) + 
  font("legend.text", size = 21) +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 25,hjust = 0.5),
        legend.position = "none")
  
my_comparisons = list(c("MNG", "LGG"))
p1 + stat_compare_means(size = 6,label.y = 1.1,label.x = 1.5)


##### DC marker
DC.markers = Total[1,]

Value_DC = c(DC.markers)

Cluster_DC = c(rep(c("Supratentorial tumors"), times = 68), rep(c("Infratentorial tumors"), times = 152), rep(c("Infratentorial tumors"), times = 130),
            rep(c("Infratentorial tumors"), times = 160),rep(c("Supratentorial tumors"), times = 227),rep(c("Supratentorial tumors"), times = 76))


Class_DC = c(rep(c("DC.markers"), times = 813))

input_DC = cbind(Value_DC, Cluster_DC, Class_DC)

# library(dplyr)
# input_DC = arrange(input_DC,Cluster_DC)

input_DC = as.data.frame(input_DC)
input_DC$Value_DC = as.numeric(input_DC$Value)
input_DC$Cluster_DC = factor(input_DC$Cluster, levels = c("Supratentorial tumors", "Infratentorial tumors" ))
input_DC$Class_DC = factor(input_DC$Class, levels = c("DC.markers")) 


p1 <- ggboxplot(input_DC, x = "Cluster_DC", y = "Value_DC",
                color = "black", palette =c("darkseagreen4", "deep sky blue"),
                title = "Dendritic cell markers",  xlab = "", ylab = "GSVA score", outlier.size = 0.05,
                fill = "Cluster_DC") + 
  ylim(-1,1.3)+
  scale_x_discrete(labels = c("ST","IT"))+
  font("xlab",size = 20) + 
  font("ylab",size = 25) + 
  font("xy.text", size = 20) + 
  font("legend.title", size = 21) + 
  font("legend.text", size = 21) +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 25,hjust = 0.5),
        legend.position = "none")

my_comparisons = list(c("Supratentorial tumors", "Infratentorial tumors"))
p1 + stat_compare_means(method = "t.test", size = 6,
                        label.y = 1.1,comparisons = my_comparisons)


#axis.text.x = element_text(size = 20, angle = 30,vjust=0.6)


