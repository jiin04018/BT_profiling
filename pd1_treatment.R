############# GBM PD-1 treatment
gbm_pd1 <- read.table("/data1/NGS_process/RNASeq/tpmcalculator/SRP155030/PRJNA482620_GBM_TPM_lg2.gct",
                      sep = "\t",header = T,skip = 2,fill = T)
IPS = read.table("/data1/CGL_data1/jiin/BT_data/IPS_surv_.txt")
IPS = IPS[c(1:154),]
idx = rownames(IPS)
gbm_pd1 <- gbm_pd1[gbm_pd1$NAME %in% idx,]
rownames(gbm_pd1) <- gbm_pd1$NAME
gbm_pd1 <- gbm_pd1[,-c(1,2)]

library(ConsensusClusterPlus)
d = gbm_pd1
d = as.matrix(d)
d = sweep(d, 1, apply(d, 1, median, na.rm=T))
results = ConsensusClusterPlus(d, maxK = 6, reps = 1000, title = "/data1/CGL_data1/jiin/BT_data/SRP155030_PD1/consensus/", pItem = 0.8, pFeature = 1,
                               clusterAlg = "km", distance = "euclidean", plot = "png")



samples = results[[2]][["consensusClass"]]
gbm_pd1 = t(gbm_pd1)
gbm_pd1 = as.data.frame(gbm_pd1)
gbm_pd1$cluster = samples

##### mechanism
library(limma)
gbm_pd1_a <- read.table("/data1/NGS_process/RNASeq/tpmcalculator/SRP155030/PRJNA482620_GBM_TPM_lg2.gct",
                        sep = "\t",header = T,skip = 2,fill = T)
gbm_pd1_a <- gbm_pd1_a[c(5:nrow(gbm_pd1_a)),]
rownames(gbm_pd1_a) <- gbm_pd1_a[,1]
t = gbm_pd1_a[-which(duplicated(gbm_pd1_a$NAME)),]
rownames(t) = t[,1]
gbm_pd1_a = t
rm(t)
gbm_pd1_a <- gbm_pd1_a[,-c(1,2)]

grp = unlist(lapply(samples,function(x){
  if(grepl(1,x)) "ImmuneCluster.1"
  else if(grepl(2,x)) "ImmuneCluster.2"
}))

# gbm_pd1_a <- t(gbm_pd1)
# gbm_pd1_a <- as.data.frame(gbm_pd1_a)

design = model.matrix(~0 + grp)
colnames(design)

colnames(design) <- c("test","ctrl")

fit = lmFit(gbm_pd1_a,design)
cont = makeContrasts(test-ctrl,levels = design)
fit.cont <- contrasts.fit(fit,cont)
fit.cont <- eBayes(fit.cont)
res <- topTable(fit.cont,number=Inf)
head(res)

res$logFC = 10^res$logFC
res$logFC = log2(res$logFC)

res$logP = -log(res$P.Value)
res$label = "no difference"
res$label[res$logP > -log10(0.05) & res$logFC > 2] <- "Up in ImmuneCluster.1"
res$label[res$logP > -log10(0.05) & res$logFC < -2] <- "Up in ImmuneCluster.2"

up_Immune1= res[res$label == "Up in ImmuneCluster.1",]
list_Immune1_DEG = as.character(rownames(up_Immune1))

up_Immune2= res[res$label == "Up in ImmuneCluster.2",]
list_Immune2_DEG = as.character(rownames(up_Immune2))

res$Gene = rownames(res)

res$genelabels <- ifelse(res$logP > -log10(0.05) & res$logFC > 2 |
                           res$logP > -log10(0.05) & res$logFC < -2 , T, F)
library(ggplot2)
library(ggrepel)

ggplot(data = res, aes(x = logFC, y = logP, col=label)) + geom_point() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "red") +
  scale_color_manual(breaks = c("Up in ImmuneCluster.1","nodiff","Up in ImmuneCluster.2"),
                     values = c("red","grey","blue")) + 
  theme_classic() +
  #scale_x_continuous(limits = c(-2.5,2)) + 
  #  geom_text_repel(aes(x = logFC, y = logP, col=label), label = ifelse(res$genelabels, res$Gene, "")) + 
  labs( x = "log2FoldChange", y = "-logP") + theme(axis.title = element_text(size = 18, color = "black"),
                                                   axis.text = element_text(size = 16, color = "black"),
                                                   legend.title = element_blank(),
                                                   legend.text = element_text(size = 15),
                                                   legend.position = c(0.25,0.8))


#### GSEA
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

matching = mapIds(org.Hs.eg.db, c(rownames(res)), 'ENTREZID', 'SYMBOL')
d = as.data.frame(cbind(matching,(res$logFC)))
geneList = as.numeric(d[,2])
names(geneList) = as.character(d[,1])
geneList = sort(geneList, decreasing = T)

#ggo = gseGO(geneList, OrgDb = "org.Hs.eg.db", ont = "BP")
library(msigdbr)
#m_df = msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
m_df = msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::select(gs_name, entrez_gene)
ggo1 = GSEA(geneList, TERM2GENE = m_df)
ggo1_res = ggo1@result

gseaplot2(ggo1,geneSetID = 128,title = ggo1@result$Description[128],base_size = 20) +
  theme(text = element_text(size = 20))

ggo1_wp = gseWP(geneList,organism = "Homo sapiens")
ggo1_wp_res = ggo1_wp@result

gseaplot2(ggo1_wp,geneSetID = 186,title = ggo1_wp@result$Description[186],base_size = 15) +
  theme(text = element_text(size = 15))

##### immune cluster 1 (pre vs post)
design = model.matrix(~0 + immune_cluster1_treat)
colnames(design)

colnames(design) <- c("test","ctrl")

fit = lmFit(immune_cluster1,design)
cont = makeContrasts(test-ctrl,levels = design)
fit.cont <- contrasts.fit(fit,cont)
fit.cont <- eBayes(fit.cont)
res <- topTable(fit.cont,number=Inf)
head(res)

res$logFC = 10^res$logFC
res$logFC = log2(res$logFC)

res$logP = -log(res$P.Value)
res$label = "no difference"
res$label[res$logP > -log10(0.05) & res$logFC > 2] <- "Up in Post"
res$label[res$logP > -log10(0.05) & res$logFC < -2] <- "Up in Pre"

up_post_immune1= res[res$label == "Up in Post",]
list_post_DEG_immune1 = as.character(rownames(up_post_immune1))

up_pre_immune1 = res[res$label == "Up in Pre",]
list_pre_DEG_immune1 = as.character(rownames(up_pre_immune1))

res$Gene = rownames(res)

res$genelabels <- ifelse(res$logP > -log10(0.05) & res$logFC > 2 |
                           res$logP > -log10(0.05) & res$logFC < -2 , T, F)
library(ggplot2)
library(ggrepel)

ggplot(data = res, aes(x = logFC, y = logP, col=label)) + geom_point() + 
  # geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  # geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  scale_color_manual(breaks = c("Up in Post","no difference","Up in Pre"),
                     values = c("#E64B35FF","grey","#00A087FF")) + 
  theme_classic() +
  #scale_x_continuous(limits = c(-2.5,2)) + 
  #  geom_text_repel(aes(x = logFC, y = logP, col=label), label = ifelse(res$genelabels, res$Gene, "")) + 
  labs(title = "ImmuneCluster1", x = "Log2(Fold Change)", y = "-Log10(adj.P.Val)") + 
  theme(axis.title = element_text(size = 18, color = "black"),
        plot.title = element_text(size = 20,hjust = .5),
        axis.text = element_text(size = 20, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "none")

rm(up_post_immune1,up_pre_immune1)



##### immune cluster 2 (pre vs post)
design = model.matrix(~0 + immune_cluster2_treat)
colnames(design)

colnames(design) <- c("test","ctrl")

fit = lmFit(immune_cluster2,design)
cont = makeContrasts(test-ctrl,levels = design)
fit.cont <- contrasts.fit(fit,cont)
fit.cont <- eBayes(fit.cont)
res <- topTable(fit.cont,number=Inf)
head(res)

res$logFC = 10^res$logFC
res$logFC = log2(res$logFC)

res$logP = -log(res$P.Value)
res$label = "no difference"
res$label[res$logP > -log10(0.05) & res$logFC > 2] <- "Up in Post"
res$label[res$logP > -log10(0.05) & res$logFC < -2] <- "Up in Pre"

up_post_immune2= res[res$label == "Up in Post",]
list_post_DEG_immune2 = as.character(rownames(up_post_immune2))

up_pre_immune2 = res[res$label == "Up in Pre",]
list_pre_DEG_immune2 = as.character(rownames(up_pre_immune2))

res$Gene = rownames(res)

res$genelabels <- ifelse(res$logP > -log10(0.05) & res$logFC > 2 |
                           res$logP > -log10(0.05) & res$logFC < -2 , T, F)
library(ggplot2)
library(ggrepel)

ggplot(data = res, aes(x = logFC, y = logP, col=label)) + geom_point() + 
  # geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  # geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  scale_color_manual(breaks = c("Up in Post","no difference","Up in Pre"),
                     values = c("#E64B35FF","grey","#00A087FF")) + 
  theme_classic() +
  #scale_x_continuous(limits = c(-2.5,2)) + 
  #  geom_text_repel(aes(x = logFC, y = logP, col=label), label = ifelse(res$genelabels, res$Gene, "")) + 
  labs(title = "ImmuneCluster2", x = "Log2(Fold Change)", y = "-Log10(adj.P.Val)") + 
  theme(axis.title = element_text(size = 18, color = "black"),
        plot.title = element_text(size = 20,hjust = .5),
        axis.text = element_text(size = 20, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "none")

rm(up_post_immune2,up_pre_immune2)

## GO
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

matching = mapIds(org.Hs.eg.db, list_post_DEG_immune2, 'ENTREZID', 'SYMBOL')


ggo = enrichGO(matching, OrgDb = "org.Hs.eg.db", ont = "BP", readable = T)


barplot(ggo, showCategory = 20, title = "", 
        font.size = 13) + theme(plot.title = element_text(size = 25),
                                legend.title = element_text(size = 20),
                                legend.text =element_text(size = 20),
                                legend.key.size = unit(5, "line"))

library(ggpubr)
ggo_res <- ggo@result[,c(1:6)]
ggo_res$logFDR <- -log10(ggo_res$p.adjust)
ggo_res <- ggo_res[c(1:15),]

ggbarplot(ggo_res, x = "Description", y = "logFDR",
          fill = "cyan", rotate = T, sort.val = "asc",
          color = "cyan", xlab = "", ylab = "-Log10(adj.P.Val)",title = "ImmuneCluster2 (Post vs Pre)") +
  geom_text(aes(x = Description, y = logFDR, label = Description),size = 8,
            xpd = TRUE,cex=0.9,adj=0,y=-log10(0.05)-1.3,colour="black") +
  theme(legend.position = "none",
        title=element_text(size = 25),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_blank(),
        # axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.ticks.y = element_line(size=0)) +
  geom_hline(yintercept = c(-log10(0.05)), linetype = "dashed", color = "black")


###### PD-1 therapy datasets
therapy = readxl::read_excel("/Users/knu_cgl3/Desktop/all/PD-1/GSE121810.xlsx")
therapy = as.data.frame(therapy)
therapy = therapy[-23876,]
therapy = therapy[-23874,]

rownames(therapy) = therapy[,1]
therapy = therapy[,-1]

cpm = apply(therapy, 2, function(x) (x/sum(x))*1000000)
log.cpm = log(cpm + 1)
log.cpm = as.data.frame(log.cpm)

IPS = read.table("/Users/knu_cgl3/Desktop/all/IPS/IPS_genes.txt",header = T, sep = "\t")

index = IPS[,1]


gene_matrix = as.data.frame(log.cpm[index, ])

gene_matrix = na.omit(gene_matrix)




library(ConsensusClusterPlus)
d = gene_matrix
d = as.matrix(d)
d = sweep(d, 1, apply(d, 1, median, na.rm=T))
results = ConsensusClusterPlus(d, maxK = 6, reps = 1000, title = "/Users/knu_cgl3/Desktop/all/PD-1", pItem = 0.8, pFeature = 1,
                               clusterAlg = "km", distance = "euclidean", plot = "png")

samples = results[[2]][["consensusClass"]]
gene_matrix[162,] = samples
rownames(gene_matrix)[162] = "Cluster"
gene_matrix[162,] = as.character(gene_matrix[162,])


setwd("/Users/knu_cgl3/Desktop/all/PD-1")

write.table(gene_matrix, file = "IPS_Gene_based_consensus_PD-1.txt", quote = F, sep = "\t")

new_clus = ifelse(gene_matrix[162,] == "1", "ImmuneCluster.1", "ImmuneCluster.2")
gene_matrix[162,] = new_clus

cls = as.data.frame(new_clus)

write.table(cls, file = "CLuster_cls.txt", quote = F, sep = "\t")

write.table(log.cpm, file = "GSEA_input.txt", quote = F, sep = "\t")

description = rep("NA", times = 37034)
Hugo_Symbol = rownames(log.cpm)

final = cbind(Hugo_Symbol, description, log.cpm)

write.table(final, file = "GSEA_input.txt", quote = F, sep = "\t", row.names = F)

#### GSEA
kk = read.table("/Users/knu_cgl3/Desktop/all/PD-1/GSEA_input.txt",header = T, sep = "\t")
cluster = read.table("/Users/knu_cgl3/Desktop/all/PD-1/cluster_information.txt",header = F,sep = "\t")
rownames(kk) <- kk[,1]
kk = kk[,-c(1,2)]
colnames(kk) <- cluster

library(limma)

grp = unlist(lapply(colnames(kk),function(x){
  if(grepl("ImmuneCluster.1",x)) "ImmuneCluster.1"
  else if(grepl('ImmuneCluster.2',x)) "ImmuneCluster.2"
  
  
}))

library(limma)
design <- model.matrix(~0 + grp)
colnames(design)

colnames(design) <- c("test","ctrl")

fit <- lmFit(kk,design)
cont <- makeContrasts(test-ctrl,levels=design)
fit.cont <- contrasts.fit(fit,cont)
fit.cont <- eBayes(fit.cont)
res <- topTable(fit.cont,number=Inf)
head(res)

res$logFC = 10^res$logFC
res$logFC = log2(res$logFC)


res$logP = -log(res$P.Value)
res$label = "no difference"
res$label[res$logP > -log10(0.05) & res$logFC > 1] <- "Up in ImmuneCluster.1"
res$label[res$logP > -log10(0.05) & res$logFC < -1] <- "Up in ImmuneCluster.2"

up_Immune1= res[res$label == "Up in ImmuneCluster.1",]
list_Immune1_DEG = as.character(rownames(up_Immune1))

up_Immune2= res[res$label == "Up in ImmuneCluster.2",]
list_Immune2_DEG = as.character(rownames(up_Immune2))

res$Gene = rownames(res)

res$genelabels <- ifelse(res$logP > -log10(0.05) & res$logFC > 1.5 |
                           res$logP > -log10(0.05) & res$logFC < -1 , T, F)
library(ggplot2)
library(ggrepel)

ggplot(data = res, aes(x = logFC, y = logP, col=label)) + geom_point() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  scale_color_manual(breaks = c("Up in ImmuneCluster.1","nodiff","Up in ImmuneCluster.2"),
                     values = c("red","grey","blue")) + 
  theme_classic() +
  #scale_x_continuous(limits = c(-2.5,2)) + 
  #geom_text_repel(aes(x = logFC, y = logP, col=label), label = ifelse(res$genelabels, res$Gene, "")) + 
  labs( x = "log2FoldChange", y = "-logP") + theme(axis.title = element_text(size = 18, color = "black"),
                                                   axis.text = element_text(size = 16, color = "black"),
                                                   legend.title = element_blank(),
                                                   legend.text = element_text(size = 15),
                                                   legend.position = c(0.25,0.8))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

matching = mapIds(org.Hs.eg.db, c(rownames(res)), 'ENTREZID', 'SYMBOL')
d = as.data.frame(cbind(matching, 2^(res$logFC)))
geneList = as.numeric(d[,2])
names(geneList) = as.character(d[,1])
geneList = sort(geneList, decreasing = T)

library(msigdbr)
#m_df = msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
m_df = msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
ggo1 = GSEA(geneList, TERM2GENE = m_df)

gseaplot2(ggo1,geneSetID = 32,title = ggo1@result$Description[32],base_size = 15) +
  theme(text = element_text(size = 15))
