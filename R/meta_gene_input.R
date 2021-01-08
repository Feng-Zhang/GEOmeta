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

meta_gene_input = function(metaDat,geneName){
  ##' @param geneName a character of gene SYMBOL name
  ##' @param metaDat a dataframe includes GSE_GPL, GSE, GSM, Group
  ##' @return a dataframe of meta input data
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




# gene="CD38"
# GSE_group = read.csv(paste0(gene,"/GSE_group.csv"))
#
# # step1: download GSE expression and GPL annotation data
# GSE_GPL_ids = as.data.frame(table(GSE_group$GSE_GPL),stringsAsFactors =F)
# GSE_GPL_ids = get_GSE_GPL(GSE_GPL_ids[,"Var1"])
# #GSE_GPL_ids = read.csv(paste0(gene,"/GSE_ids_gpl.csv"),row.names = "GSE_GPL")
# #GSE_GPL_ids = get_GSE_GPL(row.names(GSE_GPL_ids)[1:5])
#
# # step2: get expression based on gene name
# geneExpr = get_gene_expr_GSEs(gene,GSE_GPL_ids)
#
# # step3: meta input data
# dat = merge(GSE_group,geneExpr,by.x="GSM")
# metaInput = meta_gene_input(dat,gene)
#
# # step4: meta analysis and plot
# png(paste0(gene,"_forest.png"),width = 3500, height =2000,res=300)
# res = meta_gene(metaInput)
# dev.off()
