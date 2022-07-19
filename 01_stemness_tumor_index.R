library(gelnet)
library(dplyr)
library(biomaRt)
library(synapser)
synLogin("fraya","nyg789654")
# Maps ENSEMBL IDs to HUGO
# Use srcType = "ensembl_gene_id" for Ensembl IDs
# Use srcType = "entrezgene" for Entrez IDs
ensembl <- biomaRt::useMart( "ENSEMBL_MART_ENSEMBL", 
  host="www.ensembl.org", 
  dataset="hsapiens_gene_ensembl" )
httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes2hugo <- function( v, srcType = "ensembl_gene_id" )
{
    ## Retrieve the EMSEMBL -> HUGO mapping
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
X <- read.delim( synRNA$path ) %>%
tibble::column_to_rownames( "tracking_id" ) %>% as.matrix

# Retrieve metadata
synMeta <- synTableQuery( "SELECT UID, Diffname_short FROM syn3156503" )
Y <- read.delim(synMeta$filepath,sep = ",") %>%
  mutate( UID = gsub("-", ".", UID) ) %>%
  tibble::column_to_rownames( "UID" )

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


#f_m<-as.data.frame(t(X));
#f_m$class<-"Non-SC";
#scsample<-names(y[y=="SC"])
#for (i in 1:length(y)){
#  if(rownames(f_m)[i]%in% scsample){f_m$class="SC"}
#}

#if(!is.null(fnGenes)){
#  vGenes <- read.delim( fnGenes, header=FALSE ) %>% as.matrix() %>% drop()
#  VE <- genes2hugo( vGenes, "entrezgene" )
#  X <- X[intersect( rownames(X), VE[,2] ),]
#}

m <- apply( X, 1, mean )
X <- X - m
X[1:3,1:3]

j <- which( y == "SC" )
X.tr <- X[,j]
X.tr[1:3,1:3]

X.bk <- X[,-j]
X.bk[1:3,1:3]

f_m_tr<-as.data.frame(t(X.tr));
f_m_tr$class<-"SC";

f_m_bk<-as.data.frame(t(X.bk));
f_m_bk$class<-"Non-SC";

f_m<-rbind(f_m_bk,f_m_tr)
f_m$class<-factor(f_m$class)
# install.packages('Boruta')
library(Boruta)
#> Loading required package: ranger
boruta_output <- Boruta(class ~ ., data=na.omit(f_m),
 doTrace=0,maxRuns=300)
boruta_signif <- getSelectedAttributes(boruta_output, 
  withTentative = TRUE)
pdf("stemness_index_Feature_select.pdf",width=20,height=6)
plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")
dev.off();

X.tr<-X.tr[boruta_signif[-141],]

mm <- gelnet( t(X.tr), NULL, 0, 1 )
X.bk<-X.bk[boruta_signif[-141],]
## Perform leave-one-out cross-validation

auc <- c()
for(i in 1:ncol(X.tr)){
  ## Train a model on non-left-out data
  X1 <- X.tr[,-i]
  m1 <- gelnet( t(X1), NULL, 0, 1 )
  ## Score the left-out sample against the background
  s.bk <- apply( X.bk, 2, function(z) {cor( m1$w, z, method="sp" )} )
  s1 <- cor( m1$w, X.tr[,i], method="sp" )

  ## AUC = P( left-out sample is scored above the background )
  auc[i] <- sum( s1 > s.bk ) / length(s.bk)
  cat( "Current AUC: ", auc[i], "\n" )
  cat( "Average AUC: ", mean(auc), "\n" )
}

## stemness index 

stemness_index<-mm$w
stemness_gene<-boruta_signif[-141];
write.table(stemness_index,"stemness_index_withFS.txt")
# Load TCGA pancancer RNAseq data to crate tumor index

s <- synGet( "syn4976369", downloadLocation = "./data/pancan" )

# Auxiliary function: Reduces HUGO|POSITION gene IDs to just HUGO
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )

X <- read.delim( s$path, as.is=TRUE, check.names=FALSE ) %>%    ## Read the raw values
     filter( !grepl( "\\?", gene_id ) ) %>%     ## Drop genes with no mapping to HUGO
     mutate( gene_id = f( gene_id ) ) #%>%       ## Clip gene ids to HUGO
     #filter( gene_id %in% names(w) )            ## Reduce to the signature's gene set
X

  ## SLC35E2 has multiple entries with the same HUGO id
  ## Keep the first entry only
# j <- grep( "SLC35E2", X[,1] )
# if( length(j) > 1 )
#   X <- X[-j[-1],]
# 
# ## Convert to a matrix
# rownames(X) <- NULL
# X <- X %>% tibble::column_to_rownames( "gene_id" ) %>% as.matrix()
# 
# ## Reduce the signature to the common set of genes
# stopifnot( all( rownames(X) %in% names(w) ) )
# w <- w[ rownames(X) ]
# 
# ####### Score via Spearman correlation
# s <- apply( X, 2, function(z) {cor( z, w, method = "sp", use = "complete.obs" )} )
# 
# ## Scale the scores to be between 0 and 1
# s <- s - min(s)
# s <- s / max(s)
# 
# write.table(cbind(s), file = fnOut, sep = "\t", quote = FALSE, col.names = FALSE)


## get pancancer metadata 
sample<-colnames(X)
non_tumor<-grep("^TCGA-..-....-1[0-9].-...-....-..",sample,value=T)
tumor<-grep("^TCGA-..-....-0[0-9].-...-....-..",sample,value=T)
Tcga_Anno<-data.frame(sample,type="tumor")
for (i in 1:length(sample)){
  if (Tcga_Anno[i,]$sample %in% non_tumor){Tcga_Anno[i,]$type="non_tumor"}
}

X<-X[!duplicated(X$gene_id),]
rownames(X)<-X[,1]
X<-X[,-1]


m <- apply( X, 1, mean )
X <- X - m
X[1:3,1:3]

#j <- which( y == "SC" )
X.tr <- X[,tumor]
X.tr[1:3,1:3]

X.bk <- X[,non_tumor]
X.bk[1:3,1:3]

f_m_tr<-as.data.frame(t(X.tr));
f_m_tr$class<-"tumor";

f_m_bk<-as.data.frame(t(X.bk));
f_m_bk$class<-"Non-tumor";

f_m<-rbind(f_m_bk,f_m_tr)
f_m$class<-factor(f_m$class)
# install.packages('Boruta')
library(Boruta)
#> Loading required package: ranger
;
tm1<-f_m[,-which(apply(f_m,2,function(x) all(is.na(x))))]

f_m<-na.omit(tm1)
boruta_output <- Boruta(class ~ ., data=na.omit(f_m),
 doTrace=0,maxRuns=11)
boruta_signif <- getSelectedAttributes(boruta_output, 
  withTentative = TRUE)
pdf("stemness_index_Feature_select.pdf",width=20,height=6)
plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")
dev.off();

X.tr<-X.tr[boruta_signif[-141],]

mm <- gelnet( t(X.tr), NULL, 0, 1 )
X.bk<-X.bk[boruta_signif[-141],]
## Perform leave-one-out cross-validation

auc <- c()
for(i in 1:ncol(X.tr)){
  ## Train a model on non-left-out data
  X1 <- X.tr[,-i]
  m1 <- gelnet( t(X1), NULL, 0, 1 )
  ## Score the left-out sample against the background
  s.bk <- apply( X.bk, 2, function(z) {cor( m1$w, z, method="sp" )} )
  s1 <- cor( m1$w, X.tr[,i], method="sp" )

  ## AUC = P( left-out sample is scored above the background )
  auc[i] <- sum( s1 > s.bk ) / length(s.bk)
  cat( "Current AUC: ", auc[i], "\n" )
  cat( "Average AUC: ", mean(auc), "\n" )
}

## stemness index 

## tumor index 

tumor_index<-mm$w

write.table(tumor_index, file = "tumor_index.txt", sep = "\t", quote = FALSE, col.names = FALSE)
write.table(stemness_index, file = "stemness_index.txt", sep = "\t", quote = FALSE, col.names = FALSE)
