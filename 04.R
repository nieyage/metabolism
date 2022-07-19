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
