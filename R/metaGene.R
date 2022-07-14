##' @title Meta analysis for one gene
##' @description nothing
##'
##' @details nothing
##'
##' @param pheGene A dataframe, includes GSE, GSM, Group, GeneExpr
##' @param geneName a character of gene SYMBOL name
##' @return a dataframe of meta input data
##' @export
##' @importFrom stats sd

metaGene <- function(pheGene,geneName){
  study <- NULL
  group_summ <- as.data.frame(table(pheGene$GSE,pheGene$Group))
  GSEs <- unique(pheGene$GSE)
  metaInput <- NULL
  for(GSE in GSEs){
    temp <- pheGene[pheGene$GSE==GSE, ]
    caseInd <- temp[temp$Group == "treatment","GSM"]
    expr_case <- as.numeric(temp[temp$Group == "treatment",geneName])
    n_case <- length(caseInd)
    mean_case <- mean(expr_case,na.rm = TRUE)
    sd_case <- sd(expr_case,na.rm=TRUE)

    controlInd <- temp[temp$Group == "control","GSM"]
    expr_control <- as.numeric(temp[temp$Group == "control",geneName])
    n_control <- length(controlInd)
    mean_control <- mean(expr_control,na.rm = TRUE)
    sd_control <- sd(expr_control,na.rm=TRUE)
    tempInput <- data.frame(study=GSE,n_case,mean_case,sd_case,n_control,mean_control,sd_control)
    metaInput <- rbind(metaInput,tempInput)

  }
  metaData <- metacont(studlab=study,n_case,mean_case,sd_case,n_control,mean_control,sd_control, data=metaInput, sm="SMD",title=geneName)
  return(metaData)
}


##' @title Extract the result of data with meta class
##' @description nothing
##'
##' @details nothing
##'
##' @param metaData A dataframe of meta input data
##' @return A data.frame of meta results for one gene
##' @export
##' @importFrom meta metacont forest
##'
metaRes <- function(metaData){
  if(!"meta" %in% class(metaData)) stop("metaRes only receive data with class attribution of meta")
  res <- data.frame(title=metaData$title,metaData)
  res[,c("fixed_TE","fixed_lower","fixed_upper","fixed_z","fixed_pvalue","random_TE","random_lower","random_upper","random_z","random_pvalue","I2","tao2","heter_pvalue")] <- NA
  res[1,c("fixed_TE","fixed_lower","fixed_upper","fixed_z","fixed_pvalue")] <- data.frame(summary(metaData)$fixed)[c("TE", "lower", "upper","statistic","p")]
  res[1,c("random_TE","random_lower","random_upper","random_z","random_pvalue")] <- data.frame(summary(metaData)$random)[c("TE", "lower", "upper","statistic","p")]
  res[1,c("I2","tao2","heter_pvalue")] <- c(metaData$I2*100,metaData$tau^2,metaData$pval.Q)
  return(res)
}


