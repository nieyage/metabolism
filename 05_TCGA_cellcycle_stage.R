
#Step1: load TCGA tumor and normal data 

library(gelnet)
library(dplyr)
library(biomaRt)
library(synapser)
synLogin("fraya","nyg789654")

# Load TCGA pancancer RNAseq data to crate tumor index

s <- synGet( "syn4976369", downloadLocation = "./data/pancan" )

# Auxiliary function: Reduces HUGO|POSITION gene IDs to just HUGO
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )

X <- read.delim( s$path, as.is=TRUE, check.names=FALSE ) %>%    ## Read the raw values
     filter( !grepl( "\\?", gene_id ) ) %>%     ## Drop genes with no mapping to HUGO
     mutate( gene_id = f( gene_id ) ) #%>%       ## Clip gene ids to HUGO
     #filter( gene_id %in% names(w) )            ## Reduce to the signature's gene set
## SLC35E2 has multiple entries with the same HUGO id
  ## Keep the first entry only
  j <- grep( "SLC35E2", X[,1] )
  if( length(j) > 1 )
    X <- X[-j[-1],]
  
  ## Convert to a matrix
  rownames(X) <- NULL
  X <- X %>% tibble::column_to_rownames( "gene_id" ) %>% as.matrix()

## get pancancer metadata 
library(TCGAbiolinks)
pandata<-PanCancerAtlas_subtypes()
pandata<-as.data.frame(pandata)
pandata$sample<-NULL
for (i in 1:length(pandata$pan.samplesID)){
	pandata$sample[i]<-paste(unlist(strsplit( pandata$pan.samplesID[i], "-" ))[1],
		unlist(strsplit( pandata$pan.samplesID[i], "-" ))[2],
		unlist(strsplit( pandata$pan.samplesID[i], "-" ))[3],sep="-")
}

barcode<-colnames(X)
non_tumor<-grep("^TCGA-..-....-1[0-9].-...-....-..",barcode,value=T)
    tumor<-grep("^TCGA-..-....-0[0-9].-...-....-..",barcode,value=T)
Tcga_Anno<-data.frame(barcode,type="tumor")
for (i in 1:length(barcode)){
  if (Tcga_Anno[i,]$barcode %in% non_tumor){Tcga_Anno[i,]$type="non_tumor"}
}


#Step2: divide the training set and test set;
set.seed(1234)
# remove the NA data
 X <- na.omit(X)
 tumor_barcode <- Tcga_Anno[Tcga_Anno$type=="tumor",1];
 non_tumor_barcode <- Tcga_Anno[Tcga_Anno$type=="non_tumor",1];
 training_set_barcode<-sample(tumor_barcode,as.integer(0.75*length(tumor_barcode)))
 test_set_barcode<-setdiff(tumor_barcode,training_set_barcode)


 training_set_data <- X[,training_set_barcode]
 test_set_data <- X[,test_set_barcode]
 non_tumor_data <- X[,non_tumor_barcode]
 validation_data <- cbind(test_set_data,non_tumor_data)
#Step3: Calculate the CV and filter genes

cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}
gene_CV<-apply(training_set_data, 1, cal_cv)

# density plot 
data<-as.data.frame(gene_CV);
range(data$gene_CV)
# [1]  0.2260352 51.7948408
library(ggplot2)
pdf("./version2/gene_cv.pdf",width=10,height=4)
ggplot(data, aes(x=gene_CV)) + xlab("gene_CV")+
              geom_density(alpha=.25) + theme_classic() 
dev.off()
data$gene<-rownames(data)

#select gene_cv>1 as cutoff
gene_hCV<-data[data$gene_CV>1,2]

#the high CV gene correlation matrix 
training_set_data<-training_set_data[which(rownames(training_set_data)%in%gene_hCV),]
cor_data<-cor(t(training_set_data));

write.csv(cor_data,"./version3/highCV_gene_correalation.csv")
write.csv(validation_data,"./version3/validation_data.csv")
write.csv(training_set_data,"./version3/training_set_data.csv")
write.csv(test_set_data,"./version3/test_set_data.csv")
write.csv(Tcga_Anno,"./version3/Tcga_Anno_info.csv")

gene_meta<-rownames(training_set_data)[which(rownames(training_set_data)%in%gene_related_Meta$`metabolic process`)]
gene_hCV_meta<-intersect(gene_hCV,gene_meta)

#Step4: Feature Selection 
# gene related to cell cycle 
library(Seurat)
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
cc_gene<-c(s.genes,g2m.genes)
gene_hCV_meta<-unique(c(gene_hCV_meta,cc_gene))

training_set_data<-training_set_data[which(rownames(training_set_data)%in%gene_hCV_meta),]
#test_set_data<-test_set_data[which(rownames(training_set_data)%in%gene_hCV_meta),]

# Correlation Matrix
cor_data<-cor(t(training_set_data));
cc_gene_hCV<-cc_gene[cc_gene %in% gene_hCV];
cc_gene_cor_data <- cor_data[,cc_gene_hCV]
cc_gene_cor_data <- cc_gene_cor_data[-which(rownames(cc_gene_cor_data)%in%cc_gene),]

rowSums_cc_cor <- rowSums(cc_gene_cor_data);

pdf("./version2/TCGA-PanCancer-pearson-density.pdf",width=7,height=5)
  Freq<-as.data.frame(rowSums_cc_cor)
  ggplot(Freq, aes(x = rowSums_cc_cor)) +
  geom_density()+
  geom_vline(aes(xintercept=5.729450 ),linetype="dashed",col="#FB8072")+
  geom_vline(aes(xintercept=-6.908023),linetype="dashed",col="#8DD3C7")+
  scale_x_continuous(breaks = c(-10,-6.9,-5,0,5,5.7,10))+
  theme_classic()
  dev.off()

#Seletc a cutoff:
range(rowSums_cc_cor)
[1] -11.97110  22.94465

all_gene<-names(rowSums_cc_cor)
data<-data.frame()
cut<- seq(-11,22,length=1000)
for( cutoff in cut ){
  if(cutoff>0){
    gene<-names(rowSums_cc_cor[rowSums_cc_cor>cutoff]);
    random<-sample(all_gene,length(gene))
    gene_cor<-sum(rowSums_cc_cor[rowSums_cc_cor>cutoff])/length(gene)
    random_cor<-sum(rowSums_cc_cor[which(names(rowSums_cc_cor)%in%random)])/length(gene)  
  }
  if(cutoff<0){
    gene<-names(rowSums_cc_cor[rowSums_cc_cor<cutoff]);
    random<-sample(all_gene,length(gene))
    gene_cor<-sum(rowSums_cc_cor[rowSums_cc_cor<cutoff])/length(gene)
    random_cor<-sum(rowSums_cc_cor[which(names(rowSums_cc_cor)%in%random)])/length(gene)  
  };
  data_subset<-data.frame(cutoff,gene_cor,random_cor)
  data<-rbind(data,data_subset)
}

#data trans 
library(reshape2)  
colnames(data)<-c("cutoff","related_gene","random_gene")
data_long_m<-melt(data, id.vars = c("cutoff"), 
                  measure.vars = c('related_gene','random_gene'),
                  variable.name = c('group'),
                  value.name = 'mean_cor_sum')

pdf("./version2/select_cutoff.pdf",width=20,height=10)
ggplot(data_long_m, aes(x=cutoff, y=mean_cor_sum, colour=group,group=group)) + 
  geom_line(size=2)+theme_classic()+
  geom_point(size=1)
dev.off()

rowSums_cc_cor  <-rowSums_cc_cor[order(rowSums_cc_cor,decreasing=TRUE)]
top300_cc_gene_cor_positive<-names(rowSums_cc_cor)[1:300]

rowSums_cc_cor  <-rowSums_cc_cor[order(rowSums_cc_cor,decreasing=F)]
top300_cc_gene_cor_negative<-names(rowSums_cc_cor)[1:300]

# cell cycle gene 
gene_related_cellcycle<-getGO("GO:0007049")
cc_gene_cor_positive <- setdiff(top300_cc_gene_cor_positive,gene_related_cellcycle$`cell cycle`)
cc_gene_cor_negative <- setdiff(top300_cc_gene_cor_negative,gene_related_cellcycle$`cell cycle`)

#Ranking plot ;

library(ggrepel)

data<-as.data.frame(rowSums_cc_cor);
data$ranking<- 1:nrow(data);
nega<-nrow(data)-200;
data$cutoff = ifelse(data$ranking <= 200, "positive",ifelse(data$ranking >nega ,"negative","normal"))
data$label = ""
data$label[1:5] = rownames(data)[1:5];
data$label[3540:nrow(data)] = rownames(data)[3540:nrow(data)]

data$label[which(rownames(data)=="FOXM1")]="FOXM1";
data$label[which(rownames(data)=="TEAD4")]="TEAD4";
data$label[which(rownames(data)=="MYC")]="MYC";
data$label[which(rownames(data)=="NANOG")]="NANOG";
data$label[which(rownames(data)=="POU5F1")]="OCT4";


pdf("./version2/dotplot_featureselcetion_cor_result.pdf",width=7,height=5)
ggplot(data=data, aes(x=ranking, y=rowSums_cc_cor, color=cutoff)) +
  geom_point(alpha=1,size=0.5) + 
  geom_hline(aes(yintercept=5.729450 ),linetype="dashed",col="#FB8072") +
  geom_hline(aes(yintercept=-6.908023),linetype="dashed",col="#8DD3C7") +
  theme_classic(base_size=15)+  
  scale_colour_manual(values = c("#8DD3C7",'grey',"#FB8072")) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+ 
  theme(plot.title = element_text(size=15,hjust = 0.5))+
   geom_text_repel(
    #data = subset(data, data$ranking <= 5 | data$ranking >= 3540),
    data =data[data$label!="",],
    aes(label = label),
    size = 3,
    max.overlaps= 15,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )

dev.off()


#top gene correalation 

gene<-data.frame(type=c(rep("positive",length(cc_gene_cor_positive)),
  rep("negative",length(cc_gene_cor_negative))),
  gene=c(cc_gene_cor_positive,cc_gene_cor_negative))
library('biomaRt')
library("curl")

mart <- useDataset("hsapiens_gene_ensembl",useMart("ensembl"))
hs_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',
                                 'ensembl_transcript_id',"description"),
                    filters = 'external_gene_name', 
                    values= gene$gene, mart = mart)

gene$Description<-hs_symbols[match(gene$gene,hs_symbols$external_gene_name),4]

write.csv(gene,"./version2/gene_correalation.csv")

#GO and KEGG 


library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

### cc_gene_cor_positive
gene.df <- bitr(cc_gene_cor_positive, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

pdf("./version2/cc_gene_cor_positive-GO.pdf")
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
write.csv(ego,"./version2/cc_gene_cor_positive-BP.csv")

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
write.csv(ego,"./version2/cc_gene_cor_positive-MF.csv")

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
write.csv(ego,"./version2/cc_gene_cor_positive-CC.csv")
dev.off()

######KEGG##########
pdf("./version2/cc_gene_cor_positive-KEGG.pdf")
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

write.table(ekk,"./version2/cc_gene_cor_positive-KEGG.csv")
dev.off()
pdf("./version2/cc_gene_cor_positive-mkk-KEGG.pdf")

mkk <- enrichMKEGG(gene = gene.df$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mkk)             
barplot(mkk, showCategory=20)
mkk_r<-setReadable(mkk,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
write.table(mkk_r,"./version2/cc_gene_cor_positive--mkk_r-KEGG.csv")
dev.off()

#browseKEGG(kk, 'hsa04110')

#library("pathview")
#hsa04110 <- pathview(gene.data  = geneList,
#                     pathway.id = "hsa04110",
#                     species    = "hsa",
#                     limit      = list(gene=max(abs(geneList)), cpd=1))

#WikiPathways
pdf("positive-WikiPathways.pdf")
WP<-enrichWP(gene.df$ENTREZID, organism = "Homo sapiens") 
WP
barplot(WP, showCategory=20)
WP2<-setReadable(WP,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")

write.table(WP2,"./version2/cc_gene_cor_positive-WikiPathways.csv")
dev.off()

#Reactome pathway
library(ReactomePA)
pdf("./version2/positive-ReactomePathways.pdf")
Reactome <- enrichPathway(gene=gene.df$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
Reactome
barplot(Reactome, showCategory=20)
write.table(Reactome,"./version2/cc_gene_cor_positive-ReactomePathways.csv")
dev.off()

library(AnnotationHub)
library(MeSHDbi)
ah <- AnnotationHub(localHub=TRUE)
hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
file_hsa <- hsa[[1]]
db <- MeSHDbi::MeSHDb(file_hsa)
library(meshes)
data(geneList, package="DOSE")
de <- names(geneList)[1:100]
x <- enrichMeSH(de, MeSHDb = db, database='gendoo', category = 'C')

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

p1 <- cnetplot(edox, node_label="category", 
        cex_label_category = 1.2) 
p2 <- cnetplot(edox, node_label="gene", 
        cex_label_gene = 0.8) 
p3 <- cnetplot(edox, node_label="all") 
p4 <- cnetplot(edox, node_label="none", 
        color_category='firebrick', 
        color_gene='steelblue') 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

edo <- pairwise_termsim(edo)
p1 <- emapplot(edo)
p2 <- emapplot(edo, cex_category=1.5)
p3 <- emapplot(edo, layout="kk")
p4 <- emapplot(edo, cex_category=1.5,layout="kk") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

#String 
library(STRINGdb) 
library(tidyverse) 
library(clusterProfiler) 
library(org.Hs.eg.db) 
library(igraph) 
library(ggraph)

string_db <- STRINGdb$new( version="11.5", species=9606)   
STRINGdb$methods() 
 
gene <-cc_gene_cor_positive
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


### negative top 200 
gene.df <- bitr(cc_gene_cor_negative, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

pdf("./version2/cc_gene_cor_negative-GO.pdf")
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
write.csv(ego,"./version2/cc_gene_cor_negative-BP.csv")

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
write.csv(ego,"./version2/cc_gene_cor_negative-MF.csv")

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
write.csv(ego,"./version2/cc_gene_cor_negative-CC.csv")
dev.off()

# KEGG##########
pdf("./version2/cc_gene_cor_negative-KEGG.pdf")
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.1,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.1)
ego
barplot(ego, showCategory=20)
ekk<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")

write.table(ekk,"./version2/cc_gene_cor_negative-KEGG.csv")
dev.off()

pdf("./version2/cc_gene_cor-KEGG.pdf")
gene.df <- bitr(gene$gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

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

write.table(ekk,"./version2/cc_gene_cor-KEGG.csv")
dev.off()

goplot(ego)
# Verify in independent data

library(RColorBrewer)
rownames(Tcga_Anno)<-Tcga_Anno[,1]
validation_gene=c(cc_gene_cor_positive,cc_gene_cor_negative,cc_gene_hCV)
count=t(scale(t(validation_data[validation_gene,]),scale = T,center = T))

type <- factor(Tcga_Anno[match(colnames(validation_data),Tcga_Anno$barcode),2])
annotation_col = data.frame(
    Sample_type = type
)
rownames(annotation_col) = factor(colnames(count))
annotation_row = data.frame(
  GeneClass = factor(c(
    rep("cc_gene_cor_positive",length(cc_gene_cor_positive)),
    rep("cc_gene_cor_negative",length(cc_gene_cor_negative)),
    rep("cc_gene_hCV",length(cc_gene_hCV))
    ))
)
rownames(annotation_row) = factor(validation_gene)
 ann_colors= list(
  Sample_type = c("tumor"=brewer.pal(7, "Set3")[7],"non_tumor"=brewer.pal(7, "Set3")[5]),
  GeneClass=c("cc_gene_cor_positive"="#FB8072",
    "cc_gene_cor_negative"="#8DD3C7",
    "cc_gene_hCV"=brewer.pal(7, "Set3")[2]))

count<-na.omit(count)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf("./version2/top300-filter_cc_gene_cor_heatmap_tumor_test_set_data.pdf",width=10,height=7)
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



count=t(scale(t(validation_data[cc_gene_cor_negative,]),scale = T,center = T))

type <- factor(Tcga_Anno[match(colnames(validation_data),Tcga_Anno$barcode),2])
annotation_col = data.frame(
    Sample_type = type
)
rownames(annotation_col) = factor(colnames(count))
annotation_row = data.frame(
  GeneClass = factor(rep("cc_gene_cor_negative",length(cc_gene_cor_negative)))
)
rownames(annotation_row) = factor(cc_gene_cor_negative)

ann_colors= list(
  timepoint = c("tumor"=brewer.pal(7, "Set2")[1],"non_tumor"=brewer.pal(7, "Set2")[2]),
  GeneClass=c("cc_gene_cor_negative"=brewer.pal(7, "Set2")[7]))

count<-na.omit(count)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf("./version2/cc_gene_cor_negative_heatmap_tumor_test_set_data.pdf")
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

# CM differentiation 

CM_diff<-read.table("/data/R02/nieyg/project/metabolism/data/RNAseq/CMdifferentiation/Total.txt",header=TRUE)
#trans ID to genesymbol#
library('biomaRt')
library("curl")
library(ggrepel)
library(dplyr)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))# mmusculus_gene_ensembl
my_ensembl_gene_id<-CM_diff$GID
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
  filters = 'ensembl_gene_id', values= my_ensembl_gene_id, mart = mart)
CM_diff$symbol<-mms_symbols[match(CM_diff$GID,mms_symbols$ensembl_gene_id),2]
CM_diff<-na.omit(CM_diff)
CM_diff<-CM_diff[!duplicated(CM_diff$symbol),]
CM_diff<-CM_diff[,c(2:13,15)]
rownames(CM_diff)<-CM_diff[,13]
CM_diff<-CM_diff[,-13]




validation_gene=c(cc_gene_cor_positive,cc_gene_cor_negative,cc_gene_hCV)
type<-c(rep("positive",200),rep("negative",200),rep("cc_gene",40))
top_gene<-data.frame(gene=validation_gene,type);
write.table(top_gene,"./version2/top_gene.txt",row.names=F)
annotation_row = data.frame(
  GeneClass = factor(c(
    rep("cc_gene_cor_positive",length(cc_gene_cor_positive)),
    rep("cc_gene_cor_negative",length(cc_gene_cor_negative)),
    rep("cc_gene_hCV",length(cc_gene_hCV))
    ))
)
rownames(annotation_row) = factor(validation_gene)

count=t(scale(t(CM_diff[validation_gene,]),scale = T,center = T))
head(count)
timepoint <- factor(c(rep("Day0",2),rep("Day2",2),rep("Day5",2),rep("Day7",2),
                      rep("Day15",2),rep("Day80",2)))
annotation_col = data.frame(
    timepoint = timepoint
)
rownames(annotation_col) = factor(colnames(count))
library(RColorBrewer)

ann_colors= list(
  timepoint = c("Day0"=brewer.pal(7, "Set2")[1],"Day2"=brewer.pal(7, "Set2")[2],"Day5"=brewer.pal(7, "Set2")[3],
    "Day7"=brewer.pal(7, "Set2")[4],"Day15"=brewer.pal(7, "Set2")[5],"Day80"=brewer.pal(7, "Set2")[6]),
  GeneClass=c("cc_gene_cor_positive"="#FB8072",
    "cc_gene_cor_negative"="#8DD3C7",
    "cc_gene_hCV"=brewer.pal(7, "Set3")[2]))

count<-na.omit(count)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf("./version2/cc_gene_heatmap_CM_diff.pdf")
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         #cellwidth = 10, cellheight = 10,
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         show_rownames=T,show_colnames=T)
dev.off()

count=t(scale(t(CM_diff[cc_gene_cor_negative,]),scale = T,center = T))
head(count)
timepoint <- factor(c(rep("Day0",2),rep("Day2",2),rep("Day5",2),rep("Day7",2),
                      rep("Day15",2),rep("Day80",2)))
annotation_col = data.frame(
    timepoint = timepoint
)
rownames(annotation_col) = factor(colnames(count))
annotation_row = data.frame(
  GeneClass = factor(rep("cc_gene_cor_negative",length(cc_gene_cor_negative)))
)
rownames(annotation_row) = factor(cc_gene_cor_negative)
library(RColorBrewer)

ann_colors= list(
  timepoint = c("Day0"=brewer.pal(7, "Set2")[1],"Day2"=brewer.pal(7, "Set2")[2],"Day5"=brewer.pal(7, "Set2")[3],
    "Day7"=brewer.pal(7, "Set2")[4],"Day15"=brewer.pal(7, "Set2")[5],"Day80"=brewer.pal(7, "Set2")[6]),
  GeneClass=c("cc_gene_cor_negative"=brewer.pal(7, "Set2")[7]))

count<-na.omit(count)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf("./version2/cc_gene_cor_negative_heatmap_CM_diff.pdf")
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         #cellwidth = 10, cellheight = 10,
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         show_rownames=T,show_colnames=T)
dev.off()



