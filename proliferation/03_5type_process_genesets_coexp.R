library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(cowplot)
library(clusterProfiler)
library(org.Hs.eg.db)
validation_data<-read.csv("./proliferation/validation_data.csv",row.names=1)
cor_data<- read.csv("./proliferation/highCV_gene_correalation.csv",row.names=1)
# gene related to cell cycle in seurat 
library(Seurat)
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
cc_gene<-c(s.genes,g2m.genes)

# five type genesets
# gene related to metabolism/TF/immune/signal in kegg 
data <- read.csv("~/project/metabolism/Enzymes/KEGG_pathway_ko_uniq.txt",sep="\t")
metabolism_data <-data[data$level1_pathway_id=="Metabolism",]
metabolism_gene <- unique(metabolism_data$ko)
tf_data <- data[data$level3_pathway_id=="Transcription factors",]
tf_gene <- unique(tf_data$ko)
signal_gene <- unique(data[data$level2_pathway_id=="Signal transduction",]$ko)
immune_gene <- unique(data[data$level2_pathway_id=="Immune system",]$ko)

# gene related to Epifactor in Epifactor database (https://epifactors.autosome.org/) 
Epidata<- read.table("./proliferation/Epi_HGNC_symbol.txt",header=T)
Epi_gene <- Epidata$HGNC_symbol
genesets_list<- list(metabolism_gene,tf_gene,signal_gene,immune_gene,Epi_gene)
# gene related to cell adhesion in GO
adhesion_gene <- getGO("GO:0007155")
adhesion_gene <- unique(adhesion_gene[[1]])

gene_type<-data.frame(type=factor(c("Signal transduction","Epigenetic factors","Transcription factors",
	"Metabolic genes","Cell adhesion","Immune system")),
  gene_number=c(length(signal_gene),length(Epi_gene),length(tf_gene),length(metabolism_gene),length(adhesion_gene),length(immune_gene)))

gene_type$type <-factor(gene_type$type,levels=c("Signal transduction","Epigenetic factors","Transcription factors",
	"Metabolic genes","Cell adhesion","Immune system"))
pdf("./proliferation/5_genetype_barplot.pdf",width=7,height=5)
ggplot(gene_type, aes(x=as.factor(type), y=gene_number ,fill=as.factor(type))) + 
  geom_bar( stat = "identity", width = 0.6) +
  xlab("Classification of genes") +
  ylab("# of genes") + 
  scale_fill_brewer(palette = "Set3")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5));
dev.off()


Tcga_Anno<- read.csv("./version3/Tcga_Anno_info.csv",row.names=1)
rownames(Tcga_Anno)<-Tcga_Anno[,1]
gene_sets_list<- list(signal_gene,Epi_gene,tf_gene,metabolism_gene,adhesion_gene,immune_gene)
genetype<-c("Signal transduction","Epigenetic factors","Transcription factors","Metabolic genes","Cell adhesion","Immune system")
BP_gene_data<-data.frame()
gene_related_data<-data.frame()
pdf("./proliferation/proliferation_related_BP_results.pdf",width=16,height=8)
for (i in 1:length(genetype)){
  # get genesets by GO
  BP<- genetype[i]
  print(BP)
  gene_related <- gene_sets_list[[i]]

  cc_gene_hCV<-cc_gene[cc_gene %in% colnames(cor_data)];
  cc_gene_cor_data <- cor_data[,cc_gene_hCV]
  cc_gene_cor_data <- cc_gene_cor_data[-which(rownames(cc_gene_cor_data)%in%cc_gene),]
  rowSums_cc_cor <- rowSums(cc_gene_cor_data);
  BP_gene_cor <- na.omit(rowSums_cc_cor[which(gene_related%in%names(rowSums_cc_cor))])
  # get the positive genesets and negativa genesets
  if(length(BP_gene_cor)>100){
  BP_gene_cor <- BP_gene_cor[order(BP_gene_cor,decreasing=TRUE)]
  if(BP_gene_cor[100]<0){top100_cc_gene_cor_positive<-names(which(BP_gene_cor>0))   
  }else{top100_cc_gene_cor_positive<-names(BP_gene_cor)[1:100]}
  BP_gene_cor  <-BP_gene_cor[order(BP_gene_cor,decreasing=FALSE)]
  if(BP_gene_cor[100]>0){top100_cc_gene_cor_negative<-names(which(BP_gene_cor<0))   
  }else{top100_cc_gene_cor_negative<-names(BP_gene_cor)[1:100]}
  }
  if(length(BP_gene_cor)<100){
    top100_cc_gene_cor_positive<-names(which(BP_gene_cor>0));
    top100_cc_gene_cor_negative<-names(which(BP_gene_cor<0))
  }
  min_positive_cor_value <- BP_gene_cor[top100_cc_gene_cor_positive[length(top100_cc_gene_cor_positive)]]
  max_negative_cor_value <- BP_gene_cor[top100_cc_gene_cor_negative[length(top100_cc_gene_cor_negative)]]
  Freq<-as.data.frame(rowSums_cc_cor)
  #plot the rowSums_cc_cor density
  p1<-ggplot(Freq, aes(x = rowSums_cc_cor)) + geom_density()+
  geom_vline(aes(xintercept=min_positive_cor_value ),linetype="dashed",col="#FB8072")+
  geom_vline(aes(xintercept=max_negative_cor_value ),linetype="dashed",col="#8DD3C7")+
  #scale_x_continuous(breaks = c(-10,-6.9,-5,0,5,5.7,10))+
  theme_classic()
  #the Ranking plot
  BP_gene_cordata<-as.data.frame(BP_gene_cor);
  BP_gene_cordata$ranking<- 1:nrow(BP_gene_cordata);
  BP_gene_cordata$type = ifelse( rownames(BP_gene_cordata)%in% top100_cc_gene_cor_positive, "positive",ifelse(rownames(BP_gene_cordata)%in% top100_cc_gene_cor_negative ,"negative","normal"))
  BP_gene_cordata$label = ""
  BP_gene_cordata$label[1:5] = rownames(BP_gene_cordata)[1:5];
  p<-nrow(BP_gene_cordata)-4
  BP_gene_cordata$label[p:nrow(BP_gene_cordata)] = rownames(BP_gene_cordata)[p:nrow(BP_gene_cordata)]
  p2<-ggplot(data=BP_gene_cordata, aes(x=ranking, y=BP_gene_cor, color=type)) +
      geom_point(alpha=1,size=0.5) + 
      geom_hline(aes(yintercept=min_positive_cor_value),linetype="dashed",col="#FB8072") +
      geom_hline(aes(yintercept=max_negative_cor_value),linetype="dashed",col="#8DD3C7") +
      theme_classic(base_size=15)+  
      scale_colour_manual(values = c("#8DD3C7",'grey',"#FB8072")) + 
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"))+ 
      theme(plot.title = element_text(size=15,hjust = 0.5))+
       geom_text_repel(
        data =BP_gene_cordata[BP_gene_cordata$label!="",],
        aes(label = label),
        size = 3,
        max.overlaps= 15,
        box.padding = unit(0.5, "lines"),
        point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
  # validate the top100 gene 
  # Verify in independent data
  validation_gene=c(top100_cc_gene_cor_positive,top100_cc_gene_cor_negative,cc_gene_hCV)
  count=t(scale(t(validation_data[validation_gene,]),scale = T,center = T))
  colnames(validation_data)<-gsub("\\.","-",colnames(validation_data))
  type <- factor(Tcga_Anno[match(colnames(validation_data),Tcga_Anno$barcode),2])
  annotation_col = data.frame(Sample_type = type)
  rownames(annotation_col) = factor(colnames(count))
  annotation_row = data.frame(
  GeneClass = factor(c(
    rep("positive_gene",length(top100_cc_gene_cor_positive)),
    rep("negative_gene",length(top100_cc_gene_cor_negative)),
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
  bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
  p3<-pheatmap(count,cluster_cols = T,cluster_rows = T,
           color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
           legend_breaks=seq(-2,2,1),
           breaks=bk,
           annotation_col = annotation_col, 
           annotation_row = annotation_row,
           annotation_colors = ann_colors,
           show_rownames=F,show_colnames=F)
  # merge the three plots 
  left<-plot_grid(p1,p2,labels = c("A","B"),ncol = 1)
  #right<-plot_grid(top_right,p4$gtable,ncol = 1,labels = c(" ","D"))
  last<-plot_grid(left, p3$gtable, labels = c('', 'C'),rel_widths = c(1, 2.5), label_size = 12, ncol = 2)
  title <- ggdraw() + 
    draw_label(
      paste("Biological process: ",BP),
      fontface = 'bold',
      x = 0,
      hjust = 0
    )
  add_title<-plot_grid(title,last,ncol = 1 , rel_heights = c(0.1,1))
  add_title2<-plot_grid(title,p3$gtable,ncol = 1 , rel_heights = c(0.1,1))
  add_title3<-plot_grid(title,p2,ncol = 1 , rel_heights = c(0.1,1))
  print(add_title)
  print(add_title2)
  print(add_title3)
  # Data managment 
  # result 
  all_length <- length(top100_cc_gene_cor_positive) + length(top100_cc_gene_cor_negative)
  BP_gene_data_subset <- data.frame(BP_term=c(rep(BP,all_length)),
    gene_type=c(rep("positive_gene",length(top100_cc_gene_cor_positive)),rep("negative_gene",length(top100_cc_gene_cor_negative))),
    gene=c(top100_cc_gene_cor_positive,top100_cc_gene_cor_negative))
  BP_gene_data <- rbind(BP_gene_data,BP_gene_data_subset)
  # related genesets
  genesets_length<-length(gene_related)
  gene_related_data_subset<-data.frame(BP_term=c(rep(BP,genesets_length)),gene=gene_related)
  gene_related_data <- rbind(gene_related_data,gene_related_data_subset)
}
dev.off()
write.csv(BP_gene_data,"./proliferation/related_BP_top_gene_data.csv")
write.csv(gene_related_data,"./proliferation/related_BP_genesets.csv")
# the functional annotation for gene results
library(ReactomePA)
library(STRINGdb) 
library(tidyverse) 
library(clusterProfiler) 
library(org.Hs.eg.db) 
library(igraph) 
library(ggraph)
head(BP_gene_data)
cc_gene_hCV<-cc_gene[cc_gene %in% colnames(cor_data)];
cc_gene_cor_data <- cor_data[,cc_gene_hCV]
cc_gene_cor_data <- cc_gene_cor_data[-which(rownames(cc_gene_cor_data)%in%cc_gene),]
rowSums_cc_cor <- rowSums(cc_gene_cor_data);

genetype<-c("Signal transduction","Epigenetic factors","Transcription factors","Metabolic genes","Cell adhesion","Immune system")


for (BP in genetype){
  print(BP)
  input_BP_data<- BP_gene_data[which(BP_gene_data$BP_term==BP),]
  positive_gene<- input_BP_data[which(input_BP_data$gene_type=="positive_gene"),]$gene;
  negative_gene<- input_BP_data[which(input_BP_data$gene_type=="negative_gene"),]$gene;
  type<-c("positive","negative")
  for (type_gene in type){
  print(type_gene)
  chr <- paste(BP,type_gene,sep="_")
  print(chr)
  if (type_gene=="positive"){input_gene=positive_gene}else{input_gene=negative_gene};

  gene.df <- bitr(input_gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db);
  ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.05,
                readable = TRUE)
  p_ego <- barplot(ego, showCategory=20)
  p_ego_go <- goplot(ego)
  write.csv(ego,paste0("./proliferation/",chr,"-GO-BP.csv"))
  kegg <- enrichKEGG(
                gene = gene.df$ENTREZID,
                keyType = "kegg",
                organism  = 'hsa',
                pvalueCutoff  = 0.5,
                pAdjustMethod  = "BH",
                qvalueCutoff  = 0.5);
  p_kegg <- barplot(kegg, showCategory=20)
  ekk<-setReadable(kegg,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  write.csv(ekk,paste0("./proliferation/",chr,"-kegg.csv"))
  # Gene-Concept Network
  p_kegg_cnet <- cnetplot(ekk, foldChange=rowSums_cc_cor, cex_label_category=0.8)
  #WikiPathways 
  WP<-enrichWP(gene.df$ENTREZID,pvalueCutoff  = 0.5, organism = "Homo sapiens") 
  p_WP<-barplot(WP, showCategory=20)
  WP2<-setReadable(WP,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  write.table(WP2,paste0("./proliferation/",chr,"-WikiPathways.csv"))
  #Reactome pathway
  Reactome <- enrichPathway(gene=gene.df$ENTREZID, pvalueCutoff = 0.5, readable=TRUE)
  p_Reactome<-barplot(Reactome, showCategory=20)
  write.table(Reactome,paste0("./proliferation/",chr,"-Reactome.csv"))
  pdf(paste0("./proliferation/",chr,"functional_annotation.pdf"),width=12,height=8)
  print(p_ego);
  print(p_ego_go);
  print(p_kegg);
  print(p_kegg_cnet);
  print(p_WP);
  print(p_Reactome);
  dev.off()
}
}
 