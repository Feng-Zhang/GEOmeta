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

get_gene_expr_GSEs = function(geneName,GSE_GPL){
  ##' @param geneName a character of gene SYMBOL name
  ##' @param GSE_GPL a dataframe of GSE and GPL number
  ##' @return a dataframe of sepecific gene
  #geneName=gene;GSE_GPL_ids
  geneExpr = NULL
  if(file.exists(paste0(geneName,"/annoGene.txt"))) file.remove(paste0(geneName,"/annoGene.txt"))
  for(id in row.names(GSE_GPL)) {
    #id="GSE152991"
    print(id)

    # read matrix
    exprSet = readMatrixGZ(id)

    # get probe name of specific gene
    probe = getProbe(geneName,GSE_GPL[id,"GPL"],exprSet)
    probeName=probe[[1]]
    write.table(probe[[2]],file=paste0(geneName,"/annoGene.txt"),sep="\t",append = T,col.names=T,row.names = F)

    # get expression
    if(is.na(probeName)){
      temp=data.frame(row.names = colnames(exprSet),rep(NA,ncol(exprSet)))
    } else {
      temp = t(exprSet[probeName,])
    }

    colnames(temp)=geneName
    geneExpr = rbind(geneExpr,temp)
  }
  geneExpr$GSM = row.names(geneExpr)
  return(geneExpr)
}
