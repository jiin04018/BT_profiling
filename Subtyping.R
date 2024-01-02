IPS = read.table("/Users/knu_cgl3/Desktop/all/IPS/IPS_genes.txt",header = T, sep = "\t")

index = IPS[,1]

gene_matrix = read.table("/Users/knu_cgl3/Desktop/all/mixturefile.txt",header = T, sep = "\t")

rownames(gene_matrix) = gene_matrix[,1]
gene_matrix = gene_matrix[,-1]


gene_matrix = as.data.frame(gene_matrix[index, ])

gene_matrix = na.omit(gene_matrix)




library(ConsensusClusterPlus)
d = gene_matrix
d = as.matrix(d)
d = sweep(d, 1, apply(d, 1, median, na.rm=T))
results = ConsensusClusterPlus(d, maxK = 6, reps = 1000, title = "/Users/knu_cgl3/Desktop/all/IPS", pItem = 0.8, pFeature = 1,
                               clusterAlg = "km", distance = "euclidean", plot = "png")



samples = results[[2]][["consensusClass"]]

k = gene_matrix
#colnames(k) = samples

k[155,] = colnames(k)

# k = t(k)

rownames(k)[155] = "Sample"

k[156,] = samples
rownames(k)[156] = "Cluster"


k = t(k)

setwd("/Users/knu_cgl3/Desktop/all/IPS")

write.table(k, file = "IPS_Gene_based_consensus_sample.txt", quote = F, sep = "\t")
# k = as.matrix(k)
# k = as.data.frame(k)
k = t(k)
k = as.data.frame(k)
library(dplyr)
clus1 = k %>% filter(Cluster == "1")
clus2 = k %>% filter(Cluster == "2")


fin = rbind(clus1, clus2)
fin = t(fin)
colnames(fin) = fin[156,]
fin = fin[1:20,]
setwd("/Users/knu_cgl3/Desktop/all/IPS")

write.table(fin, file = "IPS_Gene_based_Consensus_cluster.txt", quote = F, sep = "\t")



##### Consensus clustering result
consen = read.table("/Users/knu_cgl3/Desktop/all/IPS/IPS_Gene_based_Consensus_cluster.txt",header = T, sep = "\t", quote = "")

consen = as.matrix(consen)
dis = consen
colnames(dis) = dis[155,]
dis = as.matrix(dis)
class(consen) = "numeric"
class(dis) = "numeric"
consen = consen[1:154,]
dis = dis[1:154,]




breaks = seq(-2, 2, 0.001)
col = colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

cluster_info = unlist(lapply(colnames(consen),function(x){
  if(grepl("X1",x)) "ImmuneCluster.1"
  else if(grepl('X2',x)) "ImmuneCluster.2"
  
  
}))


dis_info = unlist(lapply(colnames(dis),function(x){
  if(grepl("MNG",x)) "MNG" 
  else if(grepl('PA',x)) "PA" 
  else if(grepl('EPN',x)) "EPN"
  else if(grepl('MED',x)) "MED" 
  else if(grepl('GBM', x)) "GBM" 
  else if(grepl('LGG', x)) "LGG"
  
}))




annotation <- data.frame(TMEcluster = cluster_info, Disease = dis_info)

rownames(annotation) = colnames(consen)
library(RColorBrewer)
cccc = brewer.pal(8, "Accent")
annotation_colors = list(TMEcluster = c(ImmuneCluster.1="Coral", ImmuneCluster.2="Cyan"),
                         Disease = c(MNG = "#E41A1C", PA = "#377EB8", EPN = "#4DAF4A", MED = "#984EA3", 
                                     GBM = "#FF7F00", LGG = "#FFFF33"),
                         Class = c(MHC = "#7FC97F", CP = "#BEAED4",
                                       Act.CD4 = "#FDC086", Act.CD8 = "#FFFF99",
                                       Tem.CD4 = "#386CB0", Tem.CD8 = "#F0027F",
                                       MDSC = "#BF5B17", Treg = "#666666"))


library(pheatmap)

par(mar = c(4,4,2,2))
pheatmap(consen, color = col, breaks = breaks, scale = "row", cluster_cols = F, show_colnames = F,
         annotation_col = annotation, annotation_colors = annotation_colors, show_rownames = F)

dev.off()

ll = read.table("/Users/knu_cgl3/Desktop/all/IPS/IPS_genes.txt",header = T, sep = "\t", quote = "")

a = ll[1:113,]
b = ll[116:162,]
ll = rbind(a,b)

index = rownames(consen)
rownames(ll) = ll[,1]

ll = as.data.frame(ll[index, ])
ll = as.data.frame(ll[,2:3])
rname = rownames(ll)
# ll = as.data.frame(ll[,-1])
# rownames(ll) = rname
# colnames(ll) = "Cell.type"

rownames(ll)[107:108] = rownames(consen)[107:108]


ll[107,1] = "Tem CD8"
ll[108,1] = "Tem CD8"
ll[107,2] = "EC"
ll[108,2] = "EC"

MHC = as.data.frame(ll[1:9,2])
CP = as.data.frame(ll[10:19,2])
others = as.data.frame(ll[20:154, 1])
colnames(MHC) = "Cell.type"
colnames(CP) = "Cell.type"
colnames(others) = "Cell.type"
annot_row = rbind(MHC, CP, others)
rownames(annot_row) = rownames(ll)
colnames(annot_row) = "Class"

annot_row[20:42,1] = "Act.CD4"
annot_row[43:68,1] = "Act.CD8"
annot_row[69:91,1] = "Tem.CD4"
annot_row[92:115,1] = "Tem.CD8"

rownames(annot_row)[107] = rownames(consen)[107]
rownames(annot_row)[108] = rownames(consen)[108]


library(ComplexHeatmap)
annotation$Disease = factor(annotation$Disease,levels = c("MNG","PA","EPN","MED","GBM","LGG"))
annotation$TMEcluster = factor(annotation$TMEcluster,levels = c("ImmuneCluster.1","ImmuneCluster.2"))
pheatmap(consen, color = col, breaks = breaks, scale = "row", cluster_cols = F, show_colnames = F,
         annotation_col = annotation, annotation_colors = annotation_colors, show_rownames = F,
         annotation_row = annot_row, cluster_rows = T,
         column_split = annotation$TMEcluster,name = " ",
         fontsize = 20
         )


dev.off()

######### cluster ratio (disease) 
## "Coral" "Cyan"
library(ggpubr)
kk = read.table("/Users/knu_cgl3/Desktop/all/IPS/IPS_Gene_based_Consensus_sample.txt",header = T, sep = "\t", quote = "")
kk = kk[155:156,]

MNG = kk[,1:68]
PA = kk[,69:220]
EPN = kk[,221:350]
MED = kk[,351:510]
GBM = kk[,511:737]
LGG = kk[,738:813]

MNG = as.data.frame(t(MNG))
PA = as.data.frame(t(PA))
EPN = as.data.frame(t(EPN))
MED = as.data.frame(t(MED))
GBM = as.data.frame(t(GBM))
LGG = as.data.frame(t(LGG))

Cluster = c("ImmuneCluster.1", "ImmuneCluster.2")

value1 = c()
value2 = c()

for (i in (MNG$Cluster)){
  if (i == "1"){
    value1 = append(value1, values = i)
  }
  else if (i == "2"){
    value2 = append(value2, values = i)
  }
}

value = c(length(value1),length(value2))

MNG_input = as.data.frame(cbind(Cluster, value))
MNG_input$value = as.numeric(MNG_input$value)
MNG_input$percent = (MNG_input$value / sum(MNG_input$value)) * 100
#MNG_input = MNG_input[1:3,]
MNG_input$disease = "MNG"

################################################ PA

Cluster = c("ImmuneCluster.1", "ImmuneCluster.2")

value1 = c()
value2 = c()

for (i in (PA$Cluster)){
  if (i == "1"){
    value1 = append(value1, values = i)
  }
  else if (i == "2"){
    value2 = append(value2, values = i)
  }
}

value = c(length(value1),length(value2))

PA_input = as.data.frame(cbind(Cluster, value))
PA_input$value = as.numeric(PA_input$value)
PA_input$percent = (PA_input$value / sum(PA_input$value)) * 100
#PA_input = PA_input[1:3,]
PA_input$disease = "PA"
 
################################################ EPN

Cluster = c("ImmuneCluster.1", "ImmuneCluster.2")

value1 = c()
value2 = c()

for (i in (EPN$Cluster)){
  if (i == "1"){
    value1 = append(value1, values = i)
  }
  else if (i == "2"){
    value2 = append(value2, values = i)
  }
}

value = c(length(value1),length(value2))

EPN_input = as.data.frame(cbind(Cluster, value))
EPN_input$value = as.numeric(EPN_input$value)
EPN_input$percent = (EPN_input$value / sum(EPN_input$value)) * 100
#EPN_input = EPN_input[1:3,]
EPN_input$disease = "EPN"

################################################ MED

Cluster = c("ImmuneCluster.1", "ImmuneCluster.2")

value1 = c()
value2 = c()

for (i in (MED$Cluster)){
  if (i == "1"){
    value1 = append(value1, values = i)
  }
  else if (i == "2"){
    value2 = append(value2, values = i)
  }
}

value = c(length(value1),length(value2))

MED_input = as.data.frame(cbind(Cluster, value))
MED_input$value = as.numeric(MED_input$value)
MED_input$percent = (MED_input$value / sum(MED_input$value)) * 100
#MED_input = MED_input[1:3,]
MED_input$disease = "MED"

################################################ GBM

Cluster = c("ImmuneCluster.1", "ImmuneCluster.2")

value1 = c()
value2 = c()

for (i in (GBM$Cluster)){
  if (i == "1"){
    value1 = append(value1, values = i)
  }
  else if (i == "2"){
    value2 = append(value2, values = i)
  }
}

value = c(length(value1),length(value2))

GBM_input = as.data.frame(cbind(Cluster, value))
GBM_input$value = as.numeric(GBM_input$value)
GBM_input$percent = (GBM_input$value / sum(GBM_input$value)) * 100
#GBM_input = GBM_input[1:3,]
GBM_input$disease = "GBM"

################################################ LGG

Cluster = c("ImmuneCluster.1", "ImmuneCluster.2")

value1 = c()
value2 = c()

for (i in (LGG$Cluster)){
  if (i == "1"){
    value1 = append(value1, values = i)
  }
  else if (i == "2"){
    value2 = append(value2, values = i)
  }
}

value = c(length(value1),length(value2))

LGG_input = as.data.frame(cbind(Cluster, value))
LGG_input$value = as.numeric(LGG_input$value)
LGG_input$percent = (LGG_input$value / sum(LGG_input$value)) * 100
#LGG_input = LGG_input[1:3,]
LGG_input$disease = "LGG"

######## ggbarplot
df2 = rbind(MNG_input, PA_input, EPN_input, MED_input, GBM_input, LGG_input)

library(ggrepel)
labs <- paste0(" (", sprintf("%.2f",df2$percent), "%)")
ggbarplot(df2, "disease", "percent", fill = "Cluster", palette = c("Coral", "Cyan"),
           ylab = "Percent (%)", xlab = "Disease") +
  font("xlab",size = 18) + 
  font("ylab",size = 25) + 
  font("xy.text", size = 18) + 
  font("legend.title", size = 21) + 
  font("legend.text", size = 21) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 23))



######################################## GSVA
library(GSVA)
kk = read.table("/Users/knu_cgl3/Desktop/all/IPS/IPS_Gene_based_Consensus_cluster.txt",header = T, sep = "\t", quote = "")
kk = kk[1:154,]
for (i in 1:813){
  kk[,i] = as.numeric(kk[,i])
}

gs <- as.list(sample(10:100, size=8, replace=TRUE))
gs[[1]] = rownames(kk)[1:9]
gs[[2]] = rownames(kk)[10:19]
gs[[3]] = rownames(kk)[20:42]
gs[[4]] = rownames(kk)[43:68]
gs[[5]] = rownames(kk)[69:91]
gs[[6]] = rownames(kk)[92:115]
gs[[7]] = rownames(kk)[116:135]
gs[[8]] = rownames(kk)[136:154]
names(gs) <- c("MHC", "CP", "Act.CD4", "Act.CD8", "Tem.CD4", "Tem.CD8", "MDSC", "Treg")

kk = as.matrix(kk)
gsva.es <- gsva(kk, gs, verbose=FALSE,method = "gsva") #, method = "ssgsea"
dim(gsva.es)

gsva.es[1:5, 1:5]

kk = gsva.es

MHC = kk[1,]
CP = kk[2,]
Act.CD4 = kk[3,]
Act.CD8 = kk[4,]
Tem.CD4 = kk[5,]
Tem.CD8 = kk[6,]
MDSC = kk[7,]
Treg = kk[8,]


Value = c(MHC, CP, Act.CD4, Act.CD8, Tem.CD4, Tem.CD8, MDSC, Treg)


cluster = c(rep(c("ImmuneCluster.1"), times = 421), rep(c("ImmuneCluster.2"), times = 392))
Cluster = rep(cluster, times = 8)

Class = c(rep(c("MHC"), times = 813), rep(c("CP"), times = 813), rep(c("Act.CD4"), times = 813),
          rep(c("Act.CD8"), times = 813), rep(c("Tem.CD4"), times = 813), rep(c("Tem.CD8"), times = 813),
          rep(c("MDSC"), times = 813), rep(c("Treg"), times = 813))

input = cbind(Value, Cluster, Class)


input1 = as.data.frame(input)
input1$Value = as.numeric(input1$Value)
input1$Cluster = factor(input1$Cluster, levels = c("ImmuneCluster.1", "ImmuneCluster.2"))
input1$Class = factor(input1$Class, levels = c("MHC", "CP", "Act.CD4", "Act.CD8",
                                               "Tem.CD4", "Tem.CD8", "MDSC", "Treg")) 


ggviolin(input1, x = "Class", y = "Value",
          color = "Cluster", palette =c("Coral", "Cyan"),
          , x.text.angle = 30, xlab = "Class", ylab = "GSVA score", outlier.size = 0.05,fill = "Cluster",
         add = "boxplot",add.params = list(fill="white")) + 
  font("xlab",size = 18) + 
  font("ylab",size = 25) + 
  font("xy.text", size = 18) + 
  font("legend.title", size = 21) + 
  font("legend.text", size = 21) + 
  ylim(-1,1) + stat_compare_means(aes(group = Cluster), label = "p.signif",
                                  size = 7,label.y = 0.97) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 25),
        legend.title = element_blank())


####### DEG
mix = read.table("/Users/knu_cgl3/Desktop/all/mixturefile.txt",header = T, sep = "\t", quote = "")
clus = read.table("/Users/knu_cgl3/Desktop/all/IPS/IPS_Gene_based_consensus_sample.txt",header = T, sep = "\t", quote = "")

rownames(mix) <- mix$gene_symbols
mix = mix[,-1]

grp = unlist(lapply(clus[156,],function(x){
  if(grepl(1,x)) "ImmuneCluster.1"
  else if(grepl(2,x)) "ImmuneCluster.2"
}))

colnames(mix) = grp

library(limma)
design <- model.matrix(~0 + grp)
colnames(design)

colnames(design) <- c("test","ctrl")

fit <- lmFit(mix,design)
cont <- makeContrasts(test-ctrl,levels=design)
fit.cont <- contrasts.fit(fit,cont)
fit.cont <- eBayes(fit.cont)
res <- topTable(fit.cont,number=Inf)
head(res)

res$logFC = 10^res$logFC
res$logFC = log2(res$logFC)


res$logP = -log(res$P.Value)
res$label = "no difference"
res$label[res$logP > -log10(0.05) & res$logFC > 0.5] <- "Up in ImmuneCluster.1"
res$label[res$logP > -log10(0.05) & res$logFC < -0.5] <- "Up in ImmuneCluster.2"

up_Immune1= res[res$label == "Up in ImmuneCluster.1",]
list_Immune1_DEG = as.character(rownames(up_Immune1))

up_Immune2= res[res$label == "Up in ImmuneCluster.2",]
list_Immune2_DEG = as.character(rownames(up_Immune2))

res$Gene = rownames(res)

res$genelabels <- ifelse(res$logP > -log10(0.05) & res$logFC > 1 |
                           res$logP > -log10(0.05) & res$logFC < -1 , T, F)
library(ggplot2)
library(ggrepel)

ggplot(data = res, aes(x = logFC, y = logP, col=label)) + geom_point() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "red") +
  scale_color_manual(breaks = c("Up in ImmuneCluster.1","nodiff","Up in ImmuneCluster.2"),
                     values = c("red","grey","blue")) + 
  theme_classic() +
  scale_x_continuous(limits = c(-2.5,2)) + 
  geom_text_repel(aes(x = logFC, y = logP, col=label), label = ifelse(res$genelabels, res$Gene, "")) + 
  labs( x = "log2FoldChange", y = "-logP") + theme(axis.title = element_text(size = 50, color = "black"),
                                                   axis.text = element_text(size = 18, color = "black"),
                                                   legend.title = element_blank(),
                                                   legend.text = element_text(size = 15),
                                                   legend.position = "none")

#### GSEA
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

matching = mapIds(org.Hs.eg.db, c(rownames(res)), 'ENTREZID', 'SYMBOL')
d = as.data.frame(cbind(matching, 10^(res$logFC)))
geneList = as.numeric(d[,2])
names(geneList) = as.character(d[,1])
geneList = sort(geneList, decreasing = T)

#ggo = gseGO(geneList, OrgDb = "org.Hs.eg.db", ont = "BP")
library(msigdbr)
#m_df = msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
m_df = msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::select(gs_name, entrez_gene)
ggo1 = GSEA(geneList, TERM2GENE = m_df)

gsea_res = ggo1@result

gseaplot2(ggo1,geneSetID = 280,title = ggo1@result$Description[280],base_size = 15) +
  theme(text = element_text(size = 15))

ggo2 = gseWP(geneList,organism = "Homo sapiens")
ggo_res2 = as.data.frame(ggo2@result)
ggo_res2$logp <- -log10(ggo_res2$p.adjust)

gseaplot2(ggo2,geneSetID = 102,title = ggo2@result$Description[102],base_size = 15) +
  theme(text = element_text(size = 15))


