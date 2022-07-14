##' @title Regional Plot Association Results of Pig
##' @description This is a plotting tool to simultaneously display the association p-value,
##'              the LD, recombination rate and annotated gene in region of interesting in pig.
##'              The interpretation and visualization tool of association results will
##'              contribute to further identification of causal mutations.
##'
##' @details nothing
##'
##' @param geneName a character of gene SYMBOL name
##' @param GPL a character vector of GPL number
##' @param expr a GSE expression dataframe where probes would to be selected based on their probes max expression
##' @return a character of probe
##' @importFrom stringr str_detect str_split



getProbe <- function(geneName,GPL,expr){

  #geneName=gene;GPL=GSE_GPL[id,"GPL"];expr=exprSet
  anno_raw <- readLines(paste0(GPL,".soft"))
  anno_gene <- anno_raw[str_detect(anno_raw,paste(' ',geneName,' |\t',geneName,"\t",sep=""))]
  anno_gene <- str_split(anno_gene,"\t", simplify = T)
  # get matrix of specific gene
  if(length(anno_gene)==0){
    probe_name <- NA
  } else if(length(anno_gene)==1){
    probe_name <- anno_gene[,1]
  } else if(length(anno_gene)>1){
    probe_names <- anno_gene[,1]
    temp <- expr[probe_names,]
    probe_name <- rownames(temp)[which.max(rowMeans(temp))]
  }

  # get full annotation
  if(nrow(anno_gene)==0) {
    anno <- cbind(geneName,GPL,nrow(expr),"NA")
    colnames(anno) <- c("Gene","GPL","detectedExprNum","probe")
  } else{
    anno <- cbind(geneName,GPL,nrow(expr),anno_gene)
    title <- anno_raw[str_detect(anno_raw,"^ID")]
    title <- str_split(anno_gene,"\t", simplify = T)
    colnames(anno) <- c("Gene","GPL","detectedExprNum",title)
  }

  return(list(probe_name,anno))
}
