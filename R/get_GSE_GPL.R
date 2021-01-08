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
##' @export
##' @importFrom stringr  str_split_fixed
##' @importFrom GEOquery getGEO

get_GSE_GPL = function(GSE_GPL){
  #GSE_GPL= row.names(GSE_GPL_ids)[1:5]
  GSE_GPL_summ = str_split_fixed(GSE_GPL,"-",n=2)
  row.names(GSE_GPL_summ) = GSE_GPL
  colnames(GSE_GPL_summ)=c("GSE","GPL")
  for(i in 1:nrow(GSE_GPL_summ)) {
    #i=1
    GSE=GSE_GPL_summ[i,"GSE"]
    print(GSE)

    # download expresstion matrix
    eSet = getGEO(GSE, destdir = ".", getGPL = T)
    # exprSet = exprs(eSet[[1]])
    # pdata = pData(eSet[[1]])

    #get GPL
    GPLs = sapply(eSet, annotation)
    if(length(GPLs)==1) GSE_GPL_summ[i,"GPL"] = GPLs
  }
  return(GSE_GPL_summ)
}

