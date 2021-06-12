##' @title Get expression matrix without annotation
##' @description nothing
##'
##' @details nothing
##'
##' @param GSE A character, the number of GSE.
##' @param destdir A character, the path to download GSE related files.
##' @return a character of probe
##' @importFrom GEOquery getGEO
##' @importFrom Biobase exprs
##' @importFrom utils write.table
##'
getExpr <- function(GSE,destdir="tmp"){
  if(!file.exists(destdir)) dir.create(destdir)
  # download expresstion matrix
  gse = getGEO(GSE, destdir = destdir, getGPL = FALSE)
  for(i in 1:length(gse)){
    eSet = gse[[i]]
    exprFileName = paste0(destdir,"/",GSE,"-",GPL,"-matrix.raw")
    exprSet = exprs(eSet)
    if(file.exists(exprFileName)) return(gse)
    write.table(exprSet,file = exprFileName,sep="\t",quote = TRUE)
  }
  return(gse)
}

##' @title Get phenotype data
##' @description nothing
##'
##' @details nothing
##'
##' @param GSE A character, the number of GSE.
##' @param destdir A character, the path to download GSE related files.
##' @return a character of probe
##' @importFrom GEOquery getGEO
##' @importFrom Biobase pData
##' @importFrom utils write.table
getPhe <- function(GSE,destdir="tmp"){
  if(!file.exists(destdir)) dir.create(destdir)
  # download expresstion matrix
  gse = getGEO(GSE, destdir = destdir, GSEMatrix = FALSE, getGPL = FALSE)
  for(i in 1:length(gse)){
    eSet = gse[[i]]
    pheFileName = paste0(destdir,"/",GSE,"-",GPL,"-phe.raw")
    pdata = pData(eSet)
    if(file.exists(pheFileName)) return(gse)
    write.table(pdata,file = pheFileName,sep="\t",quote = TRUE)
  }

  return(gse)
}

