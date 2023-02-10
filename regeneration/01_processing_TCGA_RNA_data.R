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
library(ggplot2)
pdf("./gene_cv.pdf",width=10,height=4)
ggplot(data, aes(x=gene_CV)) + xlab("gene_CV")+
              geom_density(alpha=.25) + theme_classic() 
dev.off()
data$gene<-rownames(data)

#select gene_cv>1 as cutoff
gene_hCV<-data[data$gene_CV>1,2]

#the high CV gene correlation matrix 
training_set_data<-training_set_data[which(rownames(training_set_data)%in%gene_hCV),]
cor_data<-cor(t(training_set_data));

write.csv(cor_data,"highCV_gene_correalation.csv")
write.csv(validation_data,"validation_data.csv")
write.csv(training_set_data,"training_set_data.csv")
write.csv(test_set_data,"test_set_data.csv")
write.csv(Tcga_Anno,"Tcga_Anno_info.csv")

