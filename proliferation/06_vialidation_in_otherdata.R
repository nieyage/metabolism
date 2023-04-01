# validation in other data 
# the top gene results from data 
top_genesets <- read.csv("related_BP_top_gene_data.csv",header=T)
# positive gene upsetR 
library(UpSetR)
positive_info <- top_genesets[top_genesets$gene_type=="positive_gene",]
groupInfo <- split(positive_info$gene, positive_info$BP_term)
pdf("./top_gene_positive_upsetR.pdf")
upset(fromList(groupInfo), nintersects = 30, mb.ratio = c(0.5, 0.5),order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
dev.off()

# negative gene upsetR 
negative_info <- top_genesets[top_genesets$gene_type=="negative_gene",]
groupInfo <- split(negative_info$gene, negative_info$BP_term)
pdf("./top_gene_negative_upsetR.pdf")
upset(fromList(groupInfo), nintersects = 30, mb.ratio = c(0.5, 0.5),order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
dev.off()

# ID trans 
library(Seurat)
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
cc_gene<-c(s.genes,g2m.genes)

validation_gene = c(unique(positive_info$gene),unique(negative_info$gene),unique(cc_gene))

library(ReactomePA)
library(STRINGdb) 
library(tidyverse) 
library(clusterProfiler) 
library(org.Hs.eg.db) 
gene.df <- bitr(validation_gene, fromType = c("ALIAS"),
                toType = c("ENSEMBL","SYMBOL", "ENTREZID","GENENAME"),
                OrgDb = org.Hs.eg.db);
gene.df$gene_type<-  ifelse(gene.df$ALIAS %in% positive_info$gene , "positive_gene",
                     ifelse(gene.df$ALIAS %in% negative_info$gene,'negative_gene',
                     ifelse(gene.df$ALIAS %in% cc_gene,'cc_gene',"None")))

# Stem cell data 
library(gelnet)
library(dplyr)
library(biomaRt)
library(synapser)
synLogin("fraya","nyg789654")
# Maps ENSEMBL IDs to HUGO
# Use srcType = "ensembl_gene_id" for Ensembl IDs
# Use srcType = "entrezgene" for Entrez IDs
ensembl <- biomaRt::useMart( "ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl" )
httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes2hugo <- function( v, srcType = "ensembl_gene_id" )
{   ## Retrieve the EMSEMBL -> HUGO mapping
    ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )
    ## Make sure there was at least one mapping
    if( nrow(ID) < 1 ) top( "No IDs mapped successfully" )
    ## Drop empty duds
    j <- which( ID[,2] == "" )
    if( length(j) > 0 ) ID <- ID[-j,]
    stopifnot( all( ID[,1] %in% v ) )
    ID
}

# Load PCBC RNAseq data to crate stemness index
synRNA <- synGet( "syn2701943", downloadLocation = "./data/PCBC" )
X <- read.delim( synRNA$path ) %>%tibble::column_to_rownames( "tracking_id" ) %>% as.matrix
# Retrieve metadata
synMeta <- synTableQuery( "SELECT UID, Diffname_short FROM syn3156503" )
Y <- read.delim(synMeta$filepath,sep = ",") %>%  mutate( UID = gsub("-", ".", UID) ) %>%  tibble::column_to_rownames( "UID" )
Y1<-data.frame(Diffname_short=Y$Diffname_short)
rownames(Y1)<-rownames(Y)
# Retrieve the labels from the metadata
y <- Y1[colnames(X),]
names(y) <- colnames(X)
# Fix the missing labels by hand
y["SC11.014BEB.133.5.6.11"] <- "EB"
y["SC12.039ECTO.420.436.92.16"] <- "ECTO"
## Drop the splice form ID from the gene names
v <- strsplit( rownames(X), "\\." ) %>% lapply( "[[", 1 ) %>% unlist()
rownames(X) <- v
head(y)
V <- genes2hugo( rownames(X) )
head(V)
X <- X[V[,1],]
rownames(X) <- V[,2]

# validate the top100 gene 
# Verify in independent data
library(pheatmap)
library(RColorBrewer)
type <- factor(y[colnames(count)])
annotation_col = data.frame(Sample_type = type)
rownames(annotation_col) = factor(colnames(count))
genedf_X <- gene.df[which(gene.df$SYMBOL%in% rownames(X)),]
genedf_X <- genedf_X[!duplicated(genedf_X$SYMBOL),]
annotation_row = data.frame(GeneClass = factor(genedf_X$gene_type))
rownames(annotation_row) = factor(genedf_X$SYMBOL)
validation_gene <- rownames(annotation_row)
count=t(scale(t(X[validation_gene,]),scale = T,center = T))
count<-na.omit(count)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf("./top_genesets_validation_in_PCBC_datasets.pdf",width=10,height=10)
pheatmap(count,cluster_cols = T,cluster_rows = F,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         #annotation_colors = ann_colors,
         show_rownames=F,show_colnames=F)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         #annotation_colors = ann_colors,
         show_rownames=F,show_colnames=F)
dev.off()

write.csv(gene.df,"top_genesets_df.csv")
positive_gene.df <- gene.df[gene.df$gene_type=="positive_gene",]
negative_gene.df <- gene.df[gene.df$gene_type=="negative_gene",]

  ego <- enrichGO(positive_gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.05,
                readable = TRUE)
  p_ego <- barplot(ego, showCategory=20)
  p_ego_go <- goplot(ego)
  write.csv(ego,"positive_gene-GO-BP.csv")
  kegg <- enrichKEGG(
                gene = positive_gene.df$ENTREZID,
                keyType = "kegg",
                organism  = 'hsa',
                pvalueCutoff  = 0.5,
                pAdjustMethod  = "BH",
                qvalueCutoff  = 0.5);
  p_kegg <- barplot(kegg, showCategory=20)
  ekk<-setReadable(kegg,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  write.csv(ekk,"positive_gene-kegg.csv"))
  # Gene-Concept Network
  p_kegg_cnet <- cnetplot(ekk, foldChange=rowSums_cc_cor, cex_label_category=0.8)
  #WikiPathways 
  WP<-enrichWP(positive_gene.df$ENTREZID,pvalueCutoff  = 0.5, organism = "Homo sapiens") 
  p_WP<-barplot(WP, showCategory=20)
  WP2<-setReadable(WP,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  write.table(WP2,"positive_gene-WikiPathways.csv"))
  #Reactome pathway
  Reactome <- enrichPathway(gene=positive_gene.df$ENTREZID, pvalueCutoff = 0.5, readable=TRUE)
  p_Reactome<-barplot(Reactome, showCategory=20)
  write.table(Reactome,"positive_gene-Reactome.csv")
  pdf("positive_gene-functional_annotation.pdf"),width=12,height=8)
  print(p_ego);
  print(p_ego_go);
  print(p_kegg);
  print(p_kegg_cnet);
  print(p_WP);
  print(p_Reactome);
  dev.off()
 

  ego <- enrichGO(negative_gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.05,
                readable = TRUE)
  p_ego <- barplot(ego, showCategory=20)
  p_ego_go <- goplot(ego)
  write.csv(ego,"negative_gene-GO-BP.csv"))
  kegg <- enrichKEGG(
                gene = negative_gene.df$ENTREZID,
                keyType = "kegg",
                organism  = 'hsa',
                pvalueCutoff  = 0.5,
                pAdjustMethod  = "BH",
                qvalueCutoff  = 0.5);
  p_kegg <- barplot(kegg, showCategory=20)
  ekk<-setReadable(kegg,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  write.csv(ekk,"negative_gene-kegg.csv"))
  # Gene-Concept Network
  p_kegg_cnet <- cnetplot(ekk, foldChange=rowSums_cc_cor, cex_label_category=0.8)
  #WikiPathways 
  WP<-enrichWP(negative_gene.df$ENTREZID,pvalueCutoff  = 0.5, organism = "Homo sapiens") 
  p_WP<-barplot(WP, showCategory=20)
  WP2<-setReadable(WP,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  write.table(WP2,"negative_gene-WikiPathways.csv"))
  #Reactome pathway
  Reactome <- enrichPathway(gene=negative_gene.df$ENTREZID, pvalueCutoff = 0.5, readable=TRUE)
  p_Reactome<-barplot(Reactome, showCategory=20)
  write.table(Reactome,"negative_gene-Reactome.csv"))
  pdf("negative_gene-functional_annotation.pdf"),width=12,height=8)
  print(p_ego);
  print(p_ego_go);
  print(p_kegg);
  print(p_kegg_cnet);
  print(p_WP);
  print(p_Reactome);
  dev.off()


# validation in human heart data 
library(SeuratDisk)
library(patchwork)
library(dittoSeq)
Convert('/data/R02/nieyg/project/metabolism/data/human_heart/global_raw.h5ad', "h5seurat",
        overwrite = TRUE,assay = "RNA")
scRNA <- LoadH5Seurat("/data/R02/nieyg/project/metabolism/data/human_heart/global_raw.h5seurat")
scRNA

dittoHeatmap(sce, genes,
    annot.by = c("label", "donor"),
    order.by = "donor")