library(ggplot2)
library(RColorBrewer)
validation_data<-read.csv("./version3/validation_data.csv",row.names=1)
cor_data<- read.csv("./version3/highCV_gene_correalation.csv",row.names=1)
# gene related to cell cycle in seurat 
library(Seurat)
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
cc_gene<-c(s.genes,g2m.genes)

# gene related to metabolism/metabolic enzyme/TF in kegg 
data <- read.csv("~/project/metabolism/Enzymes/KEGG_pathway_ko_uniq.txt",sep="\t")
metabolism_data <-data[data$level1_pathway_id=="Metabolism",]
metabolism_gene <- unique(metabolism_data$ko)
metabolism_enzyme <- unique(metabolism_data[metabolism_data$ec!="-",]$ko)
tf_data <- data[data$level3_pathway_id=="Transcription factors",]
tf_gene <- unique(tf_data$ko)

# gene type barplot 
gene_type<-as.data.frame(table(metabolism_data$level2_pathway_id))
pdf("./version3/metabolism_gene_type_barplot.pdf",width=10,height=6)
ggplot(gene_type, aes(x=as.factor(Var1), y=Freq ,fill=as.factor(Var1))) + 
  geom_bar( stat = "identity", width = 0.6) +
  xlab("Classification of metabolic genes") +
  ylab("# of genes") + 
  scale_fill_brewer(palette = "Set3")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5));
dev.off()

gene_type2<-data.frame(type=c("Transcription factors","Metabolic genes","Enzymes"),
  gene_number=c(1351,1670,1484))
pdf("./version3/metabolism_gene_type2_barplot.pdf",width=5,height=6)
ggplot(gene_type2, aes(x=as.factor(type), y=gene_number ,fill=as.factor(type))) + 
  geom_bar( stat = "identity", width = 0.6) +
  xlab("Classification of genes") +
  ylab("# of genes") + 
  scale_fill_brewer(palette = "Set2")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5));
dev.off()

# all metabolism gene in kegg 
cc_gene_hCV<-cc_gene[cc_gene %in% colnames(cor_data)];
cc_gene_cor_data <- cor_data[,cc_gene_hCV]
cc_gene_cor_data <- cc_gene_cor_data[-which(rownames(cc_gene_cor_data)%in%cc_gene),]
rowSums_cc_cor <- rowSums(cc_gene_cor_data);

pdf("./version3/TCGA-PanCancer-pearson-density.pdf",width=7,height=5)
  Freq<-as.data.frame(rowSums_cc_cor)
  ggplot(Freq, aes(x = rowSums_cc_cor)) +
  geom_density()+
  geom_vline(aes(xintercept=5.729450 ),linetype="dashed",col="#FB8072")+
  geom_vline(aes(xintercept=-6.908023),linetype="dashed",col="#8DD3C7")+
  scale_x_continuous(breaks = c(-10,-6.9,-5,0,5,5.7,10))+
  theme_classic()
  dev.off()

metabolism_gene_cor <- na.omit(rowSums_cc_cor[metabolism_gene])
metabolism_gene_cor  <-metabolism_gene_cor[order(metabolism_gene_cor,decreasing=TRUE)]
top100_cc_gene_cor_positive<-names(metabolism_gene_cor)[1:100]
metabolism_gene_cor  <-metabolism_gene_cor[order(metabolism_gene_cor,decreasing=FALSE)]
top100_cc_gene_cor_negative<-names(metabolism_gene_cor)[1:100]
metabolism_gene_cor_positive <- top100_cc_gene_cor_positive 
metabolism_gene_cor_negative <- top100_cc_gene_cor_negative 

pdf("./version3/metabolism_gene_cor-pearson-density.pdf",width=7,height=5)
  Freq<-as.data.frame(metabolism_gene_cor)
  ggplot(Freq, aes(x = metabolism_gene_cor)) +
  geom_density()+
  geom_vline(aes(xintercept=2.180336 ),linetype="dashed",col="#FB8072")+
  geom_vline(aes(xintercept=-5.010242),linetype="dashed",col="#8DD3C7")+
  #scale_x_continuous(breaks = c(-10,-6.9,-5,0,5,5.7,10))+
  theme_classic()
  dev.off()
#Ranking plot ;

library(ggrepel)
data<-as.data.frame(metabolism_gene_cor);
data$ranking<- 1:nrow(data);
nega<-nrow(data)-100;
data$cutoff = ifelse(data$ranking <= 100, "positive",ifelse(data$ranking >nega ,"negative","normal"))
data$label = ""
data$label[1:5] = rownames(data)[1:5];
p<-nrow(data)-5
data$label[p:nrow(data)] = rownames(data)[p:nrow(data)]

pdf("./version3/metabolism_gene_dotplot_featureselcetion_cor_result.pdf",width=7,height=5)
ggplot(data=data, aes(x=ranking, y=metabolism_gene_cor, color=cutoff)) +
  geom_point(alpha=1,size=0.5) + 
  geom_hline(aes(yintercept=2.180336 ),linetype="dashed",col="#FB8072") +
  geom_hline(aes(yintercept=-5.010242),linetype="dashed",col="#8DD3C7") +
  theme_classic(base_size=15)+  
  scale_colour_manual(values = c("#8DD3C7",'grey',"#FB8072")) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+ 
  theme(plot.title = element_text(size=15,hjust = 0.5))+
   geom_text_repel(
    data =data[data$label!="",],
    aes(label = label),
    size = 3,
    max.overlaps= 15,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
dev.off()


metabolism_enzyme_cor <- na.omit(rowSums_cc_cor[metabolism_enzyme])
metabolism_enzyme_cor  <-metabolism_enzyme_cor[order(metabolism_enzyme_cor,decreasing=TRUE)]
top100_cc_gene_cor_positive<-names(metabolism_enzyme_cor)[1:100]
metabolism_enzyme_cor  <-metabolism_enzyme_cor[order(metabolism_enzyme_cor,decreasing=FALSE)]
top100_cc_gene_cor_negative<-names(metabolism_enzyme_cor)[1:100]
metabolism_enzyme_cor_positive <- top100_cc_gene_cor_positive 
metabolism_enzyme_cor_negative <- top100_cc_gene_cor_negative 

tf_gene_cor <- na.omit(rowSums_cc_cor[tf_gene])
tf_gene_cor  <-tf_gene_cor[order(tf_gene_cor,decreasing=TRUE)]
top100_cc_gene_cor_positive<-names(tf_gene_cor)[1:100]
tf_gene_cor  <-tf_gene_cor[order(tf_gene_cor,decreasing=FALSE)]
top100_cc_gene_cor_negative<-names(tf_gene_cor)[1:100]
tf_gene_cor_positive <- top100_cc_gene_cor_positive 
tf_gene_cor_negative <- top100_cc_gene_cor_negative 

pdf("./version3/TF_cor-pearson-density.pdf",width=7,height=5)
  Freq<-as.data.frame(tf_gene_cor)
  ggplot(Freq, aes(x = tf_gene_cor)) +
  geom_density()+
  geom_vline(aes(xintercept=3.334937 ),linetype="dashed",col="#FB8072")+
  geom_vline(aes(xintercept=-3.319746),linetype="dashed",col="#8DD3C7")+
  #scale_x_continuous(breaks = c(-10,-6.9,-5,0,5,5.7,10))+
  theme_classic()
  dev.off()
data<-as.data.frame(tf_gene_cor);
data$ranking<- 1:nrow(data);
nega<-nrow(data)-100;
data$cutoff = ifelse(data$ranking <= 100, "positive",ifelse(data$ranking >nega ,"negative","normal"))
data$label = ""
data$label[1:5] = rownames(data)[1:5];
p<-nrow(data)-5
data$label[p:nrow(data)] = rownames(data)[p:nrow(data)]

pdf("./version3/TF_dotplot_featureselcetion_cor_result.pdf",width=7,height=5)
ggplot(data=data, aes(x=ranking, y=tf_gene_cor, color=cutoff)) +
  geom_point(alpha=1,size=0.5) + 
  geom_hline(aes(yintercept=3.334937  ),linetype="dashed",col="#FB8072") +
  geom_hline(aes(yintercept=-3.319746),linetype="dashed",col="#8DD3C7") +
  theme_classic(base_size=15)+  
  scale_colour_manual(values = c("#8DD3C7",'grey',"#FB8072")) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+ 
  theme(plot.title = element_text(size=15,hjust = 0.5))+
   geom_text_repel(
    data =data[data$label!="",],
    aes(label = label),
    size = 3,
    max.overlaps= 15,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
dev.off()
#the upsetR among 6 genesets
library(UpSetR)
listInput <- list(
        metabolism_gene_posi= metabolism_gene_cor_positive, 
        metabolism_enzyme_posi = metabolism_enzyme_cor_positive, 
        tf_gene_cor_posi = tf_gene_cor_positive,
        metabolism_gene_nega= metabolism_gene_cor_negative, 
        metabolism_enzyme_nega = metabolism_enzyme_cor_negative, 
        tf_gene_cor_nega = tf_gene_cor_negative)
pdf("./version3/6genesets_upsetR.pdf")
upset(fromList(listInput), nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
dev.off()

positive_gene<-c(metabolism_gene_cor_positive,metabolism_enzyme_cor_positive,tf_gene_cor_positive)
negative_gene<-c(metabolism_gene_cor_negative,metabolism_enzyme_cor_negative,tf_gene_cor_negative)
positive_gene<-unique(positive_gene)
negative_gene<-unique(negative_gene)
# validate the top100 gene 
# Verify in independent data
Tcga_Anno<- read.csv("./version3/Tcga_Anno_info.csv",row.names=1)
rownames(Tcga_Anno)<-Tcga_Anno[,1]
validation_gene=c(positive_gene,negative_gene,cc_gene_hCV)
count=t(scale(t(validation_data[validation_gene,]),scale = T,center = T))

colnames(validation_data)<-gsub("\\.","-",colnames(validation_data))
type <- factor(Tcga_Anno[match(colnames(validation_data),Tcga_Anno$barcode),2])
annotation_col = data.frame(
    Sample_type = type
)
rownames(annotation_col) = factor(colnames(count))
annotation_row = data.frame(
  GeneClass = factor(c(
    rep("positive_gene",length(positive_gene)),
    rep("negative_gene",length(negative_gene)),
    rep("cc_gene_hCV",length(cc_gene_hCV))
    ))
)
rownames(annotation_row) = factor(validation_gene)
 ann_colors= list(
  Sample_type = c("tumor"=brewer.pal(7, "Set3")[7],"non_tumor"=brewer.pal(7, "Set3")[5]),
  GeneClass=c("positive_gene"="#FB8072",
    "negative_gene"="#8DD3C7",
    "cc_gene_hCV"=brewer.pal(7, "Set3")[2]))

count<-na.omit(count)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,.by=0.01))
pdf("./version3/merge_gene_cor_heatmap_tumor_test_set_data.pdf",width=10,height=7)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         #cellwidth = 10, cellheight = 10,
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         show_rownames=F,show_colnames=F)
dev.off()

# the GO and KEGG 
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

### positive_gene
gene.df <- bitr(positive_gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

pdf("./version3/positive-gene-GO.pdf")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
goplot(ego)
write.csv(ego,"./version3/positive_gene-BP.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"./version3/cc_gene_cor_positive-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"./version3/cc_gene_cor_positive-CC.csv")
dev.off()

######KEGG##########
pdf("./version3/cc_gene_cor_positive-KEGG.pdf")
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ekk<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")

write.table(ekk,"./version3/cc_gene_cor_positive-KEGG.csv")
dev.off()

edox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
# Gene-Concept Network
p1 <- cnetplot(edox, foldChange=rowSums_cc_cor, cex_label_category=0.8)
#p3 <- cnetplot(edox, foldChange=rowSums_cc_cor, circular = TRUE, colorEdge = TRUE) 
pdf("./version3/Gene-Concept-Network.pdf")
p1
#p3
dev.off()

#library(enrichplot)
#
#edox2 <- pairwise_termsim(edox)
#p1 <- treeplot(edox2)
#p2 <- treeplot(edox2, hclust_method = "average")
#pdf("./test.pdf")
#p1
#p2
#dev.off()

pdf("./version3/cc_gene_cor_positive-mkk-KEGG.pdf",width=10,height=6)
mkk <- enrichMKEGG(gene = gene.df$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mkk)             
barplot(mkk, showCategory=20)
mkk_r<-setReadable(mkk,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
write.table(mkk_r,"./version3/cc_gene_cor_positive-mkk_r-KEGG.csv")
dev.off()

#browseKEGG(kk, 'hsa04110')
#library("pathview")
#hsa04110 <- pathview(gene.data  = geneList,
#                     pathway.id = "hsa04110",
#                     species    = "hsa",
#                     limit      = list(gene=max(abs(geneList)), cpd=1))

#WikiPathways
pdf("./version3/positive-WikiPathways.pdf")
WP<-enrichWP(gene.df$ENTREZID, organism = "Homo sapiens") 
WP
barplot(WP, showCategory=20)
WP2<-setReadable(WP,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")

write.table(WP2,"./version3/cc_gene_cor_positive-WikiPathways.csv")
dev.off()

#Reactome pathway
library(ReactomePA)
pdf("./version3/positive-ReactomePathways.pdf")
Reactome <- enrichPathway(gene=gene.df$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
Reactome
barplot(Reactome, showCategory=20)
write.table(Reactome,"./version3/cc_gene_cor_positive-ReactomePathways.csv")
dev.off()


#String 
library(STRINGdb) 
library(tidyverse) 
library(clusterProfiler) 
library(org.Hs.eg.db) 
library(igraph) 
library(ggraph)

string_db <- STRINGdb$new( version="11.5", species=9606)   
STRINGdb$methods() 
 
gene <-positive_gene
gene <- gene %>% bitr(fromType = "SYMBOL",                       
 toType = "ENTREZID",
 OrgDb = "org.Hs.eg.db",
  drop = T)  

data_mapped <- gene %>% string_db$map(my_data_frame_id_col_names = "ENTREZID",     
                                   removeUnmappedRows = TRUE)  
string_db$plot_network( data_mapped$STRING_id )

data_links <- data_mapped$STRING_id %>% string_db$get_interactions() 
# 转换stringID为Symbol，只取前两列和最后一列 
links <- data_links %>% 
mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"])
 %>% mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%    
 dplyr::select(from, to , last_col()) %>%   
 dplyr::rename(weight = combined_score)
 # 节点数据
 nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct() 
# 创建网络图 # 根据links和nodes创建 
net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
 # 添加一些参数信息用于后续绘图 # V和E是igraph包的函数，分别用于修改网络图的节点（nodes）和连线(links) 
igraph::V(net)$deg <- igraph::degree(net) 
# 每个节点连接的节点数 
igraph::V(net)$size <- igraph::degree(net)/5 
# igraph::E(net)$width <- igraph::E(net)$weight/10 
##########################################################  
##########################################################  
ggraph(net,layout = "linear", circular = TRUE)+ 
geom_edge_arc(aes(edge_width=width),
 color = "lightblue", show.legend = F)+ 
geom_node_point(aes(size=size), color="orange", alpha=0.7)+ 
geom_node_text(aes(filter=deg>1, label=name), size = 3, repel = F)+ 
scale_edge_width(range = c(0.2,1))+ scale_size_continuous(range = c(1,10) )+ 
guides(size=F)+ theme_graph() 


# negative genes functional enrichment analysis
gene.df <- bitr(negative_gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

pdf("./version3/negative-gene-GO.pdf")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
goplot(ego)
write.csv(ego,"./version3/negative_gene-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"./version3/cc_gene_cor_negative-MF.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"./version3/cc_gene_cor_negative-CC.csv")
dev.off()

######KEGG##########
pdf("./version3/cc_gene_cor_negative-KEGG.pdf")
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ekk<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
write.table(ekk,"./version3/cc_gene_cor_negative-KEGG.csv")
dev.off()

edox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
# Gene-Concept Network
p1 <- cnetplot(edox, foldChange=rowSums_cc_cor, cex_label_category=0.8)
pdf("./version3/Gene-Concept-Network.pdf")
p1
dev.off()

pdf("./version3/cc_gene_cor_negative-mkk-KEGG.pdf",width=10,height=6)
mkk <- enrichMKEGG(gene = gene.df$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mkk)             
barplot(mkk, showCategory=20)
mkk_r<-setReadable(mkk,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
write.table(mkk_r,"./version3/cc_gene_cor_negative-mkk_r-KEGG.csv")
dev.off()

#WikiPathways
pdf("./version3/negative-WikiPathways.pdf",height=7,width=10)
WP<-enrichWP(gene.df$ENTREZID, organism = "Homo sapiens") 
WP
barplot(WP, showCategory=20)
WP2<-setReadable(WP,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
write.table(WP2,"./version3/cc_gene_cor_negative-WikiPathways.csv")
dev.off()

#Reactome pathway
library(ReactomePA)
pdf("./version3/negative-ReactomePathways.pdf")
Reactome <- enrichPathway(gene=gene.df$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
Reactome
barplot(Reactome, showCategory=20)
write.table(Reactome,"./version3/cc_gene_cor_negative-ReactomePathways.csv")
dev.off()

# the go plot 
gene.df <- bitr(positive_gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
pdf("./version3/positive-gene-GOplot.pdf",height=10,width=15)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
#ego
#barplot(ego, showCategory=20)
goplot(ego)
dev.off()
gene.df <- bitr(negative_gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
pdf("./version3/nagetive-gene-GOplot.pdf",height=10,width=15)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
goplot(ego)
dev.off()

write.csv(positive_gene,"./version3/positive_gene_205.csv",row.names=F)
write.csv(negative_gene,"./version3/negative_gene_201.csv",row.names=F)

library("pathview")
data<-data.frame(gene=names(rowSums_cc_cor),value=rowSums_cc_cor)
write.csv(data,"./version3/rowSums_cc_cor.csv")
