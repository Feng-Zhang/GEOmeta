##' @title Regional Plot Association Results of Pig
##' @description This is a plotting tool to simultaneously display the association p-value,
##'              the LD, recombination rate and annotated gene in region of interesting in pig.
##'              The interpretation and visualization tool of association results will
##'              contribute to further identification of causal mutations.
##'
##' @details nothing
##'
##' @param geneName a character of gene SYMBOL name
##' @param metaDat a dataframe includes GSE_GPL, GSE, GSM, Group
##' @return a dataframe of meta input data
##' @export

meta_gene_input = function(metaDat,geneName){

  #metaDat=dat;geneName=gene
  metaDat = metaDat[!is.na(metaDat[,geneName]),]
  summ = as.data.frame(table(metaDat$GSE_GPL),stringsAsFactors =F)
  metaInput=NULL
  for(i in 1:nrow(summ)){
    print(summ[i,"Var1"])
    temp = metaDat[metaDat$GSE_GPL==summ[i,"Var1"], ]

    id=summ[i,"Var1"]
    caseInd = temp[temp$Group == "case","GSM"]
    expr_case = as.numeric(temp[temp$GSM%in%caseInd,geneName])
    n_case = length(caseInd)
    mean_case = mean(expr_case,na.rm = T)
    sd_case = sd(expr_case,na.rm=T)

    controlInd = temp[temp$Group == "control","GSM"]
    expr_control =  as.numeric(temp[temp$GSM%in%controlInd,geneName])
    n_control = length(controlInd)
    mean_control = mean(expr_control,na.rm = T)
    sd_control = sd(expr_control,na.rm=T)
    tempInput = data.frame(study=id,n_case,mean_case,sd_case,n_control,mean_control,sd_control)
    metaInput = rbind(metaInput,tempInput)

  }
  return(metaInput)
}
