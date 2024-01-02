library(GEOquery)
library(limma)
library(umap)
library(affy)
library(hgu133a.db)
library(hgu133acdf)

gse16581 = list.files("/Users/knu_cgl3/Desktop/all/rawdata/GSE16581_RAW", full.names = T)
gse5675 = list.files("/Users/knu_cgl3/Desktop/all/rawdata/GSE5675_RAW", full.names = T)
gse73066 = list.files("/Users/knu_cgl3/Desktop/all/rawdata/GSE73066_RAW", full.names = T)
gse16155 = list.files("/Users/knu_cgl3/Desktop/all/rawdata/GSE16155_RAW", full.names = T)
gse50385 = list.files("/Users/knu_cgl3/Desktop/all/rawdata/GSE50385_RAW", full.names = T)
gse10327 = list.files("/Users/knu_cgl3/Desktop/all/rawdata/GSE10327_RAW", full.names = T)
gse37418 = list.files("/Users/knu_cgl3/Desktop/all/rawdata/GSE37418_RAW", full.names = T)
gse36245 = list.files("/Users/knu_cgl3/Desktop/all/rawdata/GSE36245_RAW", full.names = T)
gse53733 = list.files("/Users/knu_cgl3/Desktop/all/rawdata/GSE53733_RAW", full.names = T)
gse4290 = list.files("/Users/knu_cgl3/Desktop/all/rawdata/GSE4290_RAW", full.names = T)
gse44971 = list.files("/Users/knu_cgl3/Desktop/all/rawdata/GSE44971_RAW", full.names = T)
gse50161 = list.files("/Users/knu_cgl3/Desktop/all/rawdata/GSE50161_RAW", full.names = T)


## preprocessing and normalization of Affymetrix expression data
library(affy)
affy.data = ReadAffy(filenames = c(gse16581, gse5675, gse73066, gse16155, gse50385, gse10327,
                                   gse37418, gse36245, gse53733, gse4290, gse44971, gse50161))
data.rma.norm = rma(affy.data)
rma.nologs = exprs(data.rma.norm)

colnames(rma.nologs) = c(paste("gse16581", 1:68, sep = "_"), paste("gse5675", 1:41, sep = "_"), paste("gse73066", 1:47, sep = "_"),
                         paste("gse16155", 1:19, sep = "_"), paste("gse50385", 1:65, sep = "_"), 
                         paste("gse10327", 1:62, sep = "_"), paste("gse37418", 1:76, sep = "_"),
                         paste("gse36245", 1:46, sep = "_"), paste("gse53733", 1:70, sep = "-"),
                         paste("gse4290", 1:180, sep = "_"), paste("gse44971", 1:58, sep = "_"),
                         paste("gse50161", 1:130, sep = "_"))

rma = log(rma.nologs, 2)
head(rma)

setwd("/Users/knu_cgl3/Desktop/all")
write.table(rma, file = "rma_normalization.txt", quote = F, sep = "\t")


################# phenoData
phenodata = pData(affy.data)


colnames(phenodata) = c(paste("gse16581", 1:68, sep = "_"), paste("gse5675", 1:41, sep = "_"), paste("gse73066", 1:47, sep = "_"),
                         paste("gse16155", 1:19, sep = "_"), paste("gse50385", 1:65, sep = "_"), 
                         paste("gse10327", 1:62, sep = "_"), paste("gse37418", 1:76, sep = "_"),
                         paste("gse36245", 1:46, sep = "_"), paste("gse53733", 1:70, sep = "-"),
                         paste("gse4290", 1:180, sep = "_"), paste("gse44971", 1:58, sep = "_"),
                         paste("gse50161", 1:130, sep = "_"))


setwd("/Users/knu_cgl3/Desktop/all/clinical")
write.table(phenodata, file = "phenoData.txt", quote = F, sep = "\t")

