stem_index<-read.table("/md01/nieyg/project/metabolism/stemness_index.txt")
rownames(stem_index)<-stem_index[,1]
o<-stem_index[,2]
names(o)<-stem_index[,1]
stem_index<-o

tumor_index<-read.table("/md01/nieyg/project/metabolism/tumor_index.txt")
rownames(tumor_index)<-tumor_index[,1]
o<-tumor_index[,2]
names(o)<-tumor_index[,1]
tumor_index<-o

counts<-read.table("/md01/nieyg/project/metabolism/data/cross-species/human/4_featureCounts/counts-addgenesymbol.txt")
## Reduce the signature to the common set of genes
overlapgene<-intersect(names(stem_index),rownames(counts))
stem_index <- stem_index[ overlapgene ]
data<-counts[overlapgene,]
####### Score via Spearman correlation
s <- apply( data, 2, function(z) {cor( z, stem_index, method = "sp", use = "complete.obs" )} )
## Scale the scores to be between 0 and 1
s <- s - min(s)
s <- s / max(s)
stemness<-s
    WS_rep1     WS_rep2     WT_rep1     WT_rep2 
0.003573656 0.000000000 1.000000000 0.815273449 

overlapgene<-intersect(names(tumor_index),rownames(counts))
tumor_index <- tumor_index[ overlapgene ]
data<-counts[overlapgene,]
####### Score via Spearman correlation
s <- apply( data, 2, function(z) {cor( z, tumor_index, method = "sp", use = "complete.obs" )} )
## Scale the scores to be between 0 and 1
s <- s - min(s)
s <- s / max(s)

tumorness<-s

data<-rbind(stemness,tumorness)
library(pheatmap)

pdf("hMSC-stem-tumor-index-withoutFeatureSelecton.pdf")
#count<-na.omit(count)
pheatmap(data,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)
dev.off()

###top gene####


colData <- data.frame(row.names=colnames(counts),
                      condition=condition)
colData
colData$condition <- relevel(colData$condition, ref="WT")

countData <- counts[apply(counts, 1, sum) > 1 , ] 
countData<-na.omit(countData)
dds<-DESeqDataSetFromMatrix(countData,colData, 
                            formula(~condition)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)


stem_index<-stem_index[order(stem_index)]

gene_nega<-names(stem_index[1:100])
gene_posi<-names(stem_index[12848:12948])

gene_nega<-gene_nega[which(gene_nega%in% rownames(normalized_counts))]

count=t(scale(t(normalized_counts[gene_nega,]),scale = T,center = T))
head(count)
pdf("hMSC-stem-tumor-nagetive_gene-withoutFeatureSelecton.pdf",width=6,height=12)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)
dev.off()

gene_posi<-gene_posi[which(gene_posi%in% rownames(normalized_counts))]


count=t(scale(t(counts[gene_posi,]),scale = T,center = T))
head(count)
pdf("hMSC-stem-tumor-positive_gene-withoutFeatureSelecton.pdf",width=6,height=12)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)
dev.off()



data<-read.table("GSE95755_MultiCellularRNAseq_EdgeR_CPM.txt",row.names=4,sep="\t")

colnames(data)<-data[1,]
data<-data[-1,]
data<-data[,-1:-9]
data<-data[,grep("Myo",colnames(data),value=T)]

counts<-data
rownames(counts)<-toupper(rownames(counts))
overlapgene<-intersect(names(stem_index),rownames(counts))
stem_index <- stem_index[ overlapgene ]
data1<-counts[overlapgene,]
 data<-as.matrix(data1)

data<-apply(data,2,as.numeric)
rownames(data)<-rownames(data1)
####### Score via Spearman correlation
s <- apply( data, 2, function(z) {cor( z, stem_index, method = "sp", use = "complete.obs" )} )
## Scale the scores to be between 0 and 1
s <- s - min(s)
s <- s / max(s)
stemness<-s
ShamP1_Myo_1  ShamP1_Myo_2  ShamP1_Myo_3  ShamP1_Myo_4 ShamP56_Myo_1 
   0.76848713    0.20165803    0.00000000    0.03307210    1.00000000 
ShamP56_Myo_2 ShamP56_Myo_3 ShamP56_Myo_4    MIP1_Myo_1    MIP1_Myo_2 
   0.78541617    0.61293460    0.91433259    0.09944631    0.16729145 
   MIP1_Myo_3    MIP1_Myo_4   MIP56_Myo_1   MIP56_Myo_2   MIP56_Myo_3 
   0.54564411    0.55411759    0.35147182    0.67539369    0.49841178 
  MIP56_Myo_4 
   0.27891458 

overlapgene<-intersect(names(tumor_index),rownames(counts))
tumor_index <- tumor_index[ overlapgene ]
data1<-counts[overlapgene,]
data<-as.matrix(data1)
data<-apply(data,2,as.numeric)
rownames(data)<-rownames(data1)
####### Score via Spearman correlation
s <- apply( data, 2, function(z) {cor( z, tumor_index, method = "sp", use = "complete.obs" )} )
## Scale the scores to be between 0 and 1
s <- s - min(s)
s <- s / max(s)

tumorness<-s
ShamP1_Myo_1  ShamP1_Myo_2  ShamP1_Myo_3  ShamP1_Myo_4 ShamP56_Myo_1 
   0.07341136    0.24157249    0.05534595    0.06504870    0.39263136 
ShamP56_Myo_2 ShamP56_Myo_3 ShamP56_Myo_4    MIP1_Myo_1    MIP1_Myo_2 
   0.85890377    0.89283294    0.62992350    0.58933013    0.41703088 
   MIP1_Myo_3    MIP1_Myo_4   MIP56_Myo_1   MIP56_Myo_2   MIP56_Myo_3 
   0.40278599    0.00000000    0.85769460    0.74442458    1.00000000 
  MIP56_Myo_4 
   0.76158747 


data<-rbind(stemness,tumorness)
library(pheatmap)

pdf("mm10-stem-tumor-index-withoutFeatureSelecton.pdf",width=20,height=5)
#count<-na.omit(count)
pheatmap(data,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)
dev.off()



stem_index2<-read.table("/md01/nieyg/project/metabolism/stemness_index_withFS.txt")
#rownames(stem_index2)<-stem_index2[,1]
o<-stem_index2[,1]
names(o)<-rownames(stem_index2)
stem_index2<-o

