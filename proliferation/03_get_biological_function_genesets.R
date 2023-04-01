library(DOSE)
library(GOSemSim)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(GO.db)
get_GO_data <- function(OrgDb, ont, keytype) {
  GO_Env <- get_GO_Env()
  use_cached <- FALSE
  if (exists("organism", envir=GO_Env, inherits=FALSE) &&
      exists("keytype", envir=GO_Env, inherits=FALSE)) {
    org <- get("organism", envir=GO_Env)
    kt <- get("keytype", envir=GO_Env)
    if (org == DOSE:::get_organism(OrgDb) &&
        keytype == kt &&
        exists("goAnno", envir=GO_Env, inherits=FALSE)) {
      use_cached <- TRUE
    }
  }
  if (use_cached) {
    goAnno <- get("goAnno", envir=GO_Env)
  } else {
    OrgDb <- GOSemSim:::load_OrgDb(OrgDb)
    kt <- keytypes(OrgDb)
    if (! keytype %in% kt) {
      stop("keytype is not supported...")
    }    
    kk <- keys(OrgDb, keytype=keytype)
    goAnno <- suppressMessages(
      AnnotationDbi::select(OrgDb, keys=kk, keytype=keytype,
             columns=c("GOALL", "ONTOLOGYALL")))
    goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ]) 
    assign("goAnno", goAnno, envir=GO_Env)
    assign("keytype", keytype, envir=GO_Env)
    assign("organism", DOSE:::get_organism(OrgDb), envir=GO_Env)
  }
  if (ont == "ALL") {
    GO2GENE <- unique(goAnno[, c(2,1)])
  } else {
    GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
  }
  GO_DATA <- DOSE:::build_Anno(GO2GENE, get_GO2TERM_table())
  goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
  goOnt <- goOnt.df[,2]
  names(goOnt) <- goOnt.df[,1]
  assign("GO2ONT", goOnt, envir=GO_DATA)
  return(GO_DATA)
}
get_GO_Env <- function () {
  if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
  }
  get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}
get_GO2TERM_table <- function() {
  GOTERM.df <- get_GOTERM()
  GOTERM.df[, c("go_id", "Term")] %>% unique
}
get_GOTERM <- function() {
  pos <- 1
  envir <- as.environment(pos)
  if (!exists(".GOTERM_Env", envir=envir)) {
    assign(".GOTERM_Env", new.env(), envir)
  }
  GOTERM_Env <- get(".GOTERM_Env", envir = envir)
  if (exists("GOTERM.df", envir = GOTERM_Env)) {
    GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
  } else {
    GOTERM.df <- toTable(GOTERM)
    assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
  }
  return(GOTERM.df)
}
GO_DATA <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
findGO <- function(pattern, method = "key"){
    if(!exists("GO_DATA"))
        load("GO_DATA.RData")
    if(method == "key"){
        pathways = cbind(GO_DATA$PATHID2NAME[grep(pattern, GO_DATA$PATHID2NAME)])
    } else if(method == "gene"){
        pathways = cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
    }
    colnames(pathways) = "pathway"
    if(length(pathways) == 0){
        cat("No results!\n")
    } else{
        return(pathways)
    }
}
getGO <- function(ID){
    if(!exists("GO_DATA"))
        load("GO_DATA.RData")
    allNAME = names(GO_DATA$PATHID2EXTID)
    if(ID %in% allNAME){
        geneSet = GO_DATA$PATHID2EXTID[ID]
        names(geneSet) = GO_DATA$PATHID2NAME[ID]
        return(geneSet)     
    } else{
        cat("No results!\n")
    }
}

# the related biological process:
gene_related_Meta<-getGO("GO:0008152")
gene_related_cellcycle<-getGO("GO:0007049")

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