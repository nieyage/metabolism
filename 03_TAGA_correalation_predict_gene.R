#####TAGA data correalation predict gene related to regenaration 

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







#create a dataframe to store the correalation value 
library(Seurat)
> cc.genes
$s.genes
 [1] "MCM5"     "PCNA"     "TYMS"     "FEN1"     "MCM2"     "MCM4"    
 [7] "RRM1"     "UNG"      "GINS2"    "MCM6"     "CDCA7"    "DTL"     
[13] "PRIM1"    "UHRF1"    "MLF1IP"   "HELLS"    "RFC2"     "RPA2"    
[19] "NASP"     "RAD51AP1" "GMNN"     "WDR76"    "SLBP"     "CCNE2"   
[25] "UBR7"     "POLD3"    "MSH2"     "ATAD2"    "RAD51"    "RRM2"    
[31] "CDC45"    "CDC6"     "EXO1"     "TIPIN"    "DSCC1"    "BLM"     
[37] "CASP8AP2" "USP1"     "CLSPN"    "POLA1"    "CHAF1B"   "BRIP1"   
[43] "E2F8"    

$g2m.genes
 [1] "HMGB2"   "CDK1"    "NUSAP1"  "UBE2C"   "BIRC5"   "TPX2"    "TOP2A"  
 [8] "NDC80"   "CKS2"    "NUF2"    "CKS1B"   "MKI67"   "TMPO"    "CENPF"  
[15] "TACC3"   "FAM64A"  "SMC4"    "CCNB2"   "CKAP2L"  "CKAP2"   "AURKB"  
[22] "BUB1"    "KIF11"   "ANP32E"  "TUBB4B"  "GTSE1"   "KIF20B"  "HJURP"  
[29] "CDCA3"   "HN1"     "CDC20"   "TTK"     "CDC25C"  "KIF2C"   "RANGAP1"
[36] "NCAPD2"  "DLGAP5"  "CDCA2"   "CDCA8"   "ECT2"    "KIF23"   "HMMR"   
[43] "AURKA"   "PSRC1"   "ANLN"    "LBR"     "CKAP5"   "CENPE"   "CTCF"   
[50] "NEK2"    "G2E3"    "GAS2L3"  "CBX5"    "CENPA"  





cor<-matrix(ncol=24,nrow=20500);
colnames(cor)=cancer.type;
rownames(cor)=rownames(X)[-which(rownames(X)=="MKI67")]

#cancer.type[i]
for(i in 1:length(cancer.type)){
	data<-X[,tumor[which(tumor$cancer.type==cancer.type[i]),1]];
	Ki67<-data["MKI67",]
    data<-data[which(rownames(data)!="MKI67"),]
    s <- apply( data, 1, function(z) {cor( z, Ki67, method = "pearson" )} )
    cor[,cancer.type[i]]<-s;
}
write.csv(cor,"TCGA-PanCancer-pearson-Ki67.csv");
cor_density<-na.omit(as.numeric(cor));
cor<-as.data.frame(cor);
#trans NA to 0 
cor[is.na(cor)]<-0
pdf("TCGA-PanCancer-pearson-density.pdf",width=20,height=5)
  Freq<-as.data.frame(cor_density)
  ggplot(Freq, aes(x = cor_density)) +geom_density()
  data<-na.omit(rowSums(cor))
  Freq<-as.data.frame(data)
  ggplot(Freq, aes(x = data)) +geom_density()
  dev.off()

pdf("TCGA-PanCancer-pearson-density.pdf",width=10,height=5)
ggplot(Freq, aes(x=data)) + xlab("cor")+
              geom_density(alpha=.25) + theme_classic() 

d <- density(Freq$data)
d$x[which.min(abs(diff(d$y)))]

hist(Freq$data,prob=TRUE)
lines(d, col="red", lty=2)
#v <- optimize(approxfun(d$x,d$y),interval=c(0,1))$minimum
#abline(v=v, col="blue")
#
#Kmeans
df<-Freq
km <- kmeans(df$data,centers=3)
df$clust <- as.factor(km$cluster)
library(ggplot2)
ggplot(df, aes(x=data)) + 
  geom_histogram(aes(fill=clust,y=..count../sum(..count..)),
                 binwidth=0.5, color="grey50")+
  stat_density(geom="line", color="red")
dev.off()


cor$rowSums<-rowSums(cor[,1:24])
cor$sd<-apply(cor[,1:24],1,sd)
cor$mean<-apply(cor[,1:24],1,mean)
pdf("TCGA-PanCancer-sd-density.pdf",width=10,height=5)

d <- density(cor$sd)
d$x[which.min(abs(diff(d$y)))]

hist(cor$sd,prob=TRUE)
lines(d, col="red", lty=2)
df<-cor
km <- kmeans(df$sd,centers=2)
df$clust <- as.factor(km$cluster)
library(ggplot2)
ggplot(df, aes(x=sd)) + 
  geom_histogram(aes(fill=clust),
                 binwidth=0.01, color="grey50")+
  stat_density(geom="line", color="red")+theme_classic()
dev.off()

#Random select some gene to plot sd density tree 
random<-sample(cor$sd,2000);
pdf("TCGA-PanCancer-sd-density-random.pdf",width=20,height=5)
d <- density(random)
d$x[which.min(abs(diff(d$y)))]

hist(random,prob=TRUE)
lines(d, col="red", lty=2)
df<-data.frame(random)
km <- kmeans(df$random,centers=2)
df$clust <- as.factor(km$cluster)
library(ggplot2)
ggplot(df, aes(x=random)) + 
  geom_histogram(aes(fill=clust),
                 binwidth=0.01, color="grey50")+
  stat_density(geom="line", color="red")
dev.off()

pdf("TCGA-PanCancer-mean-density.pdf",width=10,height=5)
d <- density(cor$mean)
d$x[which.min(abs(diff(d$y)))]

hist(cor$mean,prob=TRUE)
lines(d, col="red", lty=2)
df<-cor
km <- kmeans(df$mean,centers=3)
df$clust <- as.factor(km$cluster)
library(ggplot2)
ggplot(df, aes(x=mean)) + 
  geom_histogram(aes(fill=clust),
                 binwidth=0.01, color="grey50")+
  stat_density(geom="line", color="red")+theme_classic()
dev.off()

#range(apply(cor[,1:24],1,sd))
#[1] 0.000000 0.318706


cutoff<-min(cor$mean)
all_gene<-rownames(cor)
data<-data.frame()
cut<-c(-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8)
for( cutoff in cut ){
  if(cutoff<0){
    gene<-rownames(cor[cor$mean< (cutoff)&cor$sd < 0.125,]);
    random<-sample(all_gene,length(gene))
    if(length(gene)>1){
    gene_cor<-abs(mean(cor[gene,27]));
    random_cor<-abs(mean(cor[random,27]));
    }
    if(length(gene)==1){
      gene_cor<-abs(cor[gene,27]);
    random_cor<-abs(cor[random,27]);
    }
    
  }
  if(cutoff>0){
    gene<-rownames(cor[cor$mean>cutoff&cor$sd < 0.125,])
    random<-sample(all_gene,length(gene));
    if(length(gene)>1){
    gene_cor<-abs(mean(cor[gene,27]));
    random_cor<-abs(mean(cor[random,27]));
    }
    if(length(gene)==1){
      gene_cor<-abs(cor[gene,27]);
    random_cor<-abs(cor[random,27]);
    }
  };
  data_subset<-data.frame(cutoff,gene_cor,random_cor)
  data<-rbind(data,data_subset)
}










#positive correlation
#mean>0.5,sd<0.15
p_gene<-rownames(cor[cor$mean>0.3&cor$sd<0.125,])

#negative correlation
#mean<-0.3,sd<0.15
n_gene<-rownames(cor[cor$mean< (-0.3)&cor$sd<0.125,])

library(pheatmap)

test_data<-Tcga_Anno[!(Tcga_Anno$type=="tumor"&!is.na(Tcga_Anno$cancer.type)),]
test<-X[p_gene,test_data$barcode];
count=t(scale(t(test),scale = T,center = T))
annotation_col = data.frame(
  Group = test_data$type 
)
rownames(annotation_col) = factor(test_data$barcode)
count<-na.omit(count)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf("test_data-p-gene-selectcutoff.pdf")

pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         annotation_col = annotation_col, 
         show_rownames=T,show_colnames=F);
dev.off();


test<-X[n_gene,test_data$barcode];
count=t(scale(t(test),scale = T,center = T))
annotation_col = data.frame(
  Group = test_data$type 
)
rownames(annotation_col) = factor(test_data$barcode)
count<-na.omit(count)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf("test_data-n-gene-selectcutoff.pdf")

pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         annotation_col = annotation_col, 
         show_rownames=T,show_colnames=F);
dev.off();




#create a dataframe to store the correalation value 
cor<-matrix(ncol=24,nrow=20500);
colnames(cor)=cancer.type;
rownames(cor)=rownames(X)[-which(rownames(X)=="CCNB1")]

for(i in 1:length(cancer.type)){
  data<-X[,tumor[which(tumor$cancer.type==cancer.type[i]),1]];
  CCNB1<-data["CCNB1",]
    data<-data[which(rownames(data)!="CCNB1"),]
    s <- apply( data, 1, function(z) {cor( z, CCNB1, method = "pearson" )} )
    cor[,cancer.type[i]]<-s;
}
write.csv(cor,"TCGA-PanCancer-pearson-CCNB1.csv");
cor_density<-na.omit(as.numeric(cor));
cor<-as.data.frame(cor);
cor$rowSums<-rowSums(cor[,1:24])
cor$sd<-apply(cor[,1:24],1,sd)
cor$mean<-apply(cor[,1:24],1,mean)
#trans NA to 0 
cor[is.na(cor)]<-0


pdf("TCGA-PanCancer-sd-density-CCNB1.pdf",width=10,height=5)

d <- density(cor$sd)
d$x[which.min(abs(diff(d$y)))]

hist(cor$sd,prob=TRUE)
lines(d, col="red", lty=2)
df<-cor
km <- kmeans(df$sd,centers=2)
df$clust <- as.factor(km$cluster)
library(ggplot2)
ggplot(df, aes(x=sd)) + 
  geom_histogram(aes(fill=clust),
                 binwidth=0.01, color="grey50")+
  stat_density(geom="line", color="red")+theme_classic()
dev.off()


pdf("TCGA-PanCancer-mean-density-CCNB1.pdf",width=10,height=5)
d <- density(cor$mean)
d$x[which.min(abs(diff(d$y)))]

hist(cor$mean,prob=TRUE)
lines(d, col="red", lty=2)
df<-cor
km <- kmeans(df$mean,centers=3)
df$clust <- as.factor(km$cluster)
library(ggplot2)
ggplot(df, aes(x=mean)) + 
  geom_histogram(aes(fill=clust),
                 binwidth=0.01, color="grey50")+
  stat_density(geom="line", color="red")+theme_classic()
dev.off()

#range(apply(cor[,1:24],1,sd))
#[1] 0.000000 0.318706


cutoff<-min(cor$mean)
all_gene<-rownames(cor)
data<-data.frame()
cut<-c(-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8)
for( cutoff in cut ){
  if(cutoff<0){
    gene<-rownames(cor[cor$mean< (cutoff)&cor$sd < 0.125,]);
    random<-sample(all_gene,length(gene))
    if(length(gene)>1){
    gene_cor<-abs(mean(cor[gene,27]));
    random_cor<-abs(mean(cor[random,27]));
    }
    if(length(gene)==1){
      gene_cor<-abs(cor[gene,27]);
    random_cor<-abs(cor[random,27]);
    }
    
  }
  if(cutoff>0){
    gene<-rownames(cor[cor$mean>cutoff&cor$sd < 0.125,])
    random<-sample(all_gene,length(gene));
    if(length(gene)>1){
    gene_cor<-abs(mean(cor[gene,27]));
    random_cor<-abs(mean(cor[random,27]));
    }
    if(length(gene)==1){
      gene_cor<-abs(cor[gene,27]);
    random_cor<-abs(cor[random,27]);
    }
  };
  data_subset<-data.frame(cutoff,gene_cor,random_cor)
  data<-rbind(data,data_subset)
}