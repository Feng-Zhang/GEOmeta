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

meta_gene = function(metaInputData,geneName){
  ##' @param geneName a character of gene SYMBOL name
  ##' @param metaInputData a dataframe of meta input data
  ##' @return a dataframe of meta results for one gene
  #metaInputData=metaInput;geneName=gene
  metaAna = metacont(studlab=study,n_case,mean_case,sd_case,n_control,mean_control,sd_control, data=metaInputData, sm="SMD",comb.fixed=F)
  forest(metaAna)
  res= data.frame(metaAna)
  res[,c("fixed_TE","fixed_lower","fixed_upper","fixed_z","fixed_pvalue","random_TE","random_lower","random_upper","random_z","random_pvalue","I2","tao2","heter_pvalue")] = NA
  res[1,c("fixed_TE","fixed_lower","fixed_upper","fixed_z","fixed_pvalue")] = data.frame(summary(metaAna)$fixed)[c("TE", "lower", "upper","statistic","p")]
  res[1,c("random_TE","random_lower","random_upper","random_z","random_pvalue")] = data.frame(summary(metaAna)$random)[c("TE", "lower", "upper","statistic","p")]
  res[1,c("I2","tao2","heter_pvalue")] = c(metaAna$I2*100,metaAna$tau^2,metaAna$pval.Q)
  return(res)
}


