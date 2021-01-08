##' @title Regional Plot Association Results of Pig
##' @description This is a plotting tool to simultaneously display the association p-value,
##'              the LD, recombination rate and annotated gene in region of interesting in pig.
##'              The interpretation and visualization tool of association results will
##'              contribute to further identification of causal mutations.
##'
##' @details nothing
##'
##' @param GSE_GPL a character vactor comparised of GSE and GPL number.Note if GSE was sequenced on only one platform, please just give GSE number.
##' @return download GSE expression matrix and GPL annotation,
##'         and return a matrix with colanmes of "GSE" and "GPL" and row names of GSE_GPL


getProbe = function(geneName,GPL,expr){
  ##' @param geneName a character of gene SYMBOL name
  ##' @param GPL a character vector of GPL number
  ##' @param expr a GSE expression dataframe where probes would to be selected based on their probes max expression
  ##' @return a character of probe
  #geneName=gene;GPL=GSE_GPL[id,"GPL"];expr=exprSet
  annoRaw = readLines(paste0(GPL,".soft"))
  annoGene = annoRaw[str_detect(annoRaw,paste(' ',geneName,' |\t',geneName,"\t",sep=""))]
  annoGene = str_split(annoGene,"\t", simplify = T)
  # get matrix of specific gene
  if(length(annoGene)==0){
    probeName = NA
  } else if(length(annoGene)==1){
    probeName = annoGene[,1]
  } else if(length(annoGene)>1){
    probeNames = annoGene[,1]
    temp = expr[probeNames,]
    probeName = rownames(temp)[which.max(rowMeans(temp))]
  }

  # get full annotation

  if(nrow(annoGene)==0) {
    anno = cbind(geneName,GPL,nrow(expr),"NA")
    colnames(anno)=c("Gene","GPL","detectedExprNum","probe")
  } else{
    anno = cbind(geneName,GPL,nrow(expr),annoGene)
    title = annoRaw[str_detect(annoRaw,"^ID")]
    title = str_split(annoGene,"\t", simplify = T)
    colnames(anno)=c("Gene","GPL","detectedExprNum",title)
  }

  return(list(probeName,anno))
}
