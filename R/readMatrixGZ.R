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


readMatrixGZ = function(GSE){
  ##' @param GSE a character of GSE number
  ##' @return a GSE expression dataframe

  #GSE="GSE151839"
  # get the number of redundant line
  exprSet = readLines(paste(GSE,"_series_matrix.txt.gz",sep=""))
  skipNum = which(str_detect(exprSet,"ID_REF"))-1

  # read matrix
  exprSet = try(fread(paste(GSE,"_series_matrix.txt.gz",sep=""),head=T,skip=skipNum),silent=T)
  if("try-error" %in% class(exprSet)){
    exprSet= read.table(paste(GSE,"_series_matrix.txt.gz",sep=""),head=T,skip=skipNum,sep="\t",comment.char = "!",row.names = "ID_REF")
  } else {
    probes = exprSet$ID_REF
    exprSet = as.data.frame(exprSet[,-1])
    row.names(exprSet) = probes
  }

  # replace null using 0
  if(any(as.vector(exprSet)=="null",na.rm = T)){
    temp=row.names(exprSet)
    exprSet = apply(exprSet,2,function(x) {x[x=="null"]=0 ; as.numeric(x)})
    exprSet = as.data.frame(exprSet)
    row.names(exprSet) = temp
  }
  return(exprSet)
}

