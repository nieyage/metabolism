













 
                     ifelse(gene.df$ALIAS %in% cc_gene,'cc_gene',"None")))
                     ifelse(gene.df$ALIAS %in% negative_info$gene,'negative_gene',
                .hj = "BP", ###BP,MF,CC
                gene = negative_gene.df$ENTREZID,
                gene = positive_gene.df$ENTREZID,
                keyType = "kegg",
                keyType = "kegg",
                keyType = 'ENTREZID',
                keyType = 'ENTREZID',
                ont = "BP", ###BP,MF,CC
                organism  = 'hsa',
                organism  = 'hsa',
                OrgDb = org.Hs.eg.db);
                OrgDb = org.Hs.eg.db,
                OrgDb = org.Hs.eg.db,
                pAdjustMethod  = "BH",
                pAdjustMethod  = "BH",
                pAdjustMethod = "BH",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.5,
                pvalueCutoff  = 0.5,
                pvalueCutoff = 0.5,
                pvalueCutoff = 0.5,
                qvalueCutoff  = 0.5);
                qvalueCutoff  = 0.5);
                qvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
                readable = TRUE)
                toType = c("ENSEMBL","SYMBOL", "ENTREZID","GENENAME"),
         #annotation_colors = ann_colors,
         #annotation_colors = ann_colors,
         annotation_col = annotation_col, 
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_row = annotation_row,
         breaks=bk,
         breaks=bk,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         legend_breaks=seq(-2,2,1),
         show_rownames=F,show_colnames=F)
         show_rownames=F,show_colnames=F)
    ## Drop empty duds
    ## Make sure there was at least one mapping
    ID
    ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )
    if( length(j) > 0 ) ID <- ID[-j,]
    if( nrow(ID) < 1 ) top( "No IDs mapped successfully" )
    j <- which( ID[,2] == "" )
    stopifnot( all( ID[,1] %in% v ) )
  # Gene-Concept Network
  # Gene-Concept Network
  #Reactome pathway
  #Reactome pathway
  #WikiPathways 
  #WikiPathways 
  dev.off()
  dev.off()
  ego <- enrichGO(negative_gene.df$ENTREZID,
  ego <- enrichGO(positive_gene.df$ENTREZID,
  ekk<-setReadable(kegg,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  ekk<-setReadable(kegg,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  kegg <- enrichKEGG(
  kegg <- enrichKEGG(
  p_ego <- barplot(ego, showCategory=20)
  p_ego <- barplot(ego, showCategory=20)
  p_ego_go <- goplot(ego)
  p_ego_go <- goplot(ego)
  p_kegg <- barplot(kegg, showCategory=20)
  p_kegg <- barplot(kegg, showCategory=20)
  p_kegg_cnet <- cnetplot(ekk, foldChange=rowSums_cc_cor, cex_label_category=0.8)
  p_kegg_cnet <- cnetplot(ekk, foldChange=rowSums_cc_cor, cex_label_category=0.8)
  p_Reactome<-barplot(Reactome, showCategory=20)
  p_Reactome<-barplot(Reactome, showCategory=20)
  p_WP<-barplot(WP, showCategory=20)
  p_WP<-barplot(WP, showCategory=20)
  pdf("negative_gene-functional_annotation.pdf"),width=12,height=8)
  pdf("positive_gene-functional_annotation.pdf"),width=12,height=8)
  print(p_ego);
  print(p_ego);
  print(p_ego_go);
  print(p_ego_go);
  print(p_kegg);
  print(p_kegg);
  print(p_kegg_cnet);
  print(p_kegg_cnet);
  print(p_Reactome);
  print(p_Reactome);
  print(p_WP);
  print(p_WP);·
  Reactome <- enrichPathway(gene=negative_gene.df$ENTREZID, pvalueCutoff = 0.5, readable=TRUE)
  Reactome <- enrichPathway(gene=positive_gene.df$ENTREZID, pvalueCutoff = 0.5, readable=TRUE)
  WP2<-setReadable(WP,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  WP2<-setReadable(WP,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  WP<-enrichWP(negative_gene.df$ENTREZID,pvalueCutoff  = 0.5, organism = "Homo sapiens") 
  WP<-enrichWP(positive_gene.df$ENTREZID,pvalueCutoff  = 0.5, organism = "Homo sapiens") 
  write.csv(ego,"negative_gene-GO-BP.csv"))
  write.csv(ego,"positive_gene-GO-BP.csv")
  write.csv(ekk,"negative_gene-kegg.csv"))
  write.csv(ekk,"positive_gene-kegg.csv"))
  write.table(Reactome,"negative_gene-Reactome.csv"))
  write.table(Reactome,"positive_gene-Reactome.csv")
  write.table(WP2,"negative_gene-WikiPathways.csv"))
  write.table(WP2,"positive_gene-WikiPathways.csv"))
# Fix the missing labels by hand
# ID trans 
# Load PCBC RNAseq data to crate stemness index
# Maps ENSEMBL IDs to HUGO
# negative gene upsetR 
# positive gene upsetR 
# Retrieve metadata
# Retrieve the labels from the metadata
# Stem cell data 
# the top gene results from data 
# Use srcType = "ensembl_gene_id" for Ensembl IDs
# Use srcType = "entrezgene" for Entrez IDs
# validate the top100 gene 
# validation in human he1234567890-art data 
# validation in other data 
# Verify in independent data
## Drop the splice form ID from the gene names
';LKJJHGFDSAzxvbnm,./'
annotation_col = data.frame(Sample_type = type)
annotation_row = data.frame(GeneClass = factor(genedf_X$gene_type))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
cc_gene<- gene.df[which(gene.df$gene_type=="cc_gene"),]$SYMBOL
cc_gene<-c(s.genes,g2m.genes)
Convert('/data/R02/nieyg/project/metabolism/data/human_heart/global_raw.h5ad', "h5seurat",overwrite = TRUE,assay = "RNA")
count<-na.omit(count)
count=t(scale(t(X[validation_gene,]),scale = T,center = T))
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dittoHeatmap(scRNA_usful, genes%in% rownames(scRNA_usful), annot.by = c("cell_states", "age_group"), order.by = "cell_states")
dittoHeatmap(scRNA_usful, genes%in% rownames(scRNA_usful), annot.by = c("cell_states", "age_group"), order.by = "cell_states")

ensembl <- biomaRt::useMart( "ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl" )
g2m.genes<-cc.genes$g2m.genes
gene.df <- bitr(validation_gene, fromType = c("ALIAS"),
gene.df$gene_type<-  ifelse(gene.df$ALIAS %in% positive_info$gene , "positive_gene",
gene.df<- read.csv("/data/R02/nieyg/project/metabolism/proliferation/top_genesets_df.csv")
genedf_X <- gene.df[which(gene.df$SYMBOL%in% rownames(X)),]
genedf_X <- genedf_X[!duplicated(genedf_X$SYMBOL),]
genes2hugo <- function( v, srcType = "ensembl_gene_id" )
genes<-c(negative_gene,cc_gene)
genes<-c(positive_gene,cc_gene)
groupInfo <- split(negative_info$gene, negative_info$BP_term)
groupInfo <- split(positive_info$gene, positive_info$BP_term)
head(V)
head(y)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
Idents(scRNA)<- scRNA$cell_type
library(biomaRt)
library(clusterProfiler) 
library(dittoSeq)
library(dplyr)
library(gelnet)
library(org.Hs.eg.db) 
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(ReactomePA)
library(Seurat)
library(Seurat)
library(SeuratDisk)
library(STRINGdb) 
library(synapser)
library(tidyverse) 
library(UpSetR)
names(y) <- colnames(X)
negative_gene.df <- gene.df[gene.df$gene_type=="negative_gene",]
negative_gene<- gene.df[which(gene.df$gene_type=="negative_gene"),]$SYMBOL
negative_info <- top_genesets[top_genesets$gene_type=="negative_gene",]
pdf("./negative_gene_human_heart_heatmap.pdf",width=12,height=8)
pdf("./positive_gene_human_heart_heatmap.pdf",width=12,height=8)
pdf("./top_gene_negative_upsetR.pdf")
pdf("./top_gene_positive_upsetR.pdf")
pdf("./top_genesets_validation_in_PCBC_datasets.pdf",width=10,height=10)
pheatmap(count,cluster_cols = T,cluster_rows = F,
pheatmap(count,cluster_cols = T,cluster_rows = T,
positive_gene.df <- gene.df[gene.df$gene_type=="positive_gene",]
positive_gene<- gene.df[which(gene.df$gene_type=="positive_gene"),]$SYMBOL
positive_info <- top_genesets[top_genesets$gene_type=="positive_gene",]
rownames(annotation_col) = factor(colnames(count))
rownames(annotation_row) = factor(genedf_X$SYMBOL)
rownames(X) <- v
rownames(X) <- V[,2]
rownames(Y1)<-rownames(Y)
s.genes<-cc.genes$s.genes
scRNA <- LoadH5Seurat("/data/R02/nieyg/project/metabolism/data/human_heart/global_raw.h5seurat")
scRNA_usful <- subset(scRNA,idents = levels(scRNA)[1:12])
synLogin("fraya","nyg789654")
synMeta <- synTableQuery( "SELECT UID, Diffname_short FROM syn3156503" )
synRNA <- synGet( "syn2701943", downloadLocation = "./data/PCBC" )
top_genesets <- read.csv("related_BP_top_gene_data.csv",header=T)
type <- factor(y[colnames(count)])
upset(fromList(groupInfo), nintersects = 30, mb.ratio = c(0.5, 0.5),order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
upset(fromList(groupInfo), nintersects = 30, mb.ratio = c(0.5, 0.5),order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
V <- genes2hugo( rownames(X) )
v <- strsplit( rownames(X), "\\." ) %>% lapply( "[[", 1 ) %>% unlist()
validation_gene <- rownames(annotation_row)
validation_gene = c(unique(positive_info$gene),unique(negative_info$gene),unique(cc_gene))
write.csv(gene.df,"top_genesets_df.csv")
X <- read.delim( synRNA$path ) %>%tibble::column_to_rownames( "tracking_id" ) %>% as.matrix
X <- X[V[,1],]
Y <- read.delim(synMeta$filepath,sep = ",") %>%  mutate( UID = gsub("-", ".", UID) ) %>%  tibble::column_to_rownames( "UID" )
y <- Y1[colnames(X),]
Y1<-data.frame(Diffname_short=Y$Diffname_short)
y["SC11.014BEB.133.5.6.11"] <- "EB"
y["SC12.039ECTO.420.436.92.16"] <- "ECTO"
{   ## Retrieve the EMSEMBL -> HUGO mapping
}




