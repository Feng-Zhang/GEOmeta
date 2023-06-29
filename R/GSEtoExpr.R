##' @title Download GSE expression matrix and convert probe to gene symbol
##' @description Given a GSE number, download related GSE data and convert probe to gene symbol, and return the gene matrix with gene symbol.
##'
##' @details nothing
##'
##' @param GSE A character, the number of GSE.
##' @param destdir A character, the path to download GSE related files.
##' @param annotSymbol A Boolean. The default is FALSE, do not annotate probe to gene Symbol.
##' @param getGPL if annoSymbol is TRUE, getGPL should be TRUE as well. A boolean defaulting to TRUE as to whether or not to download and include GPL information when getting a GSEMatrix file. You may want to set this to FALSE if you know that you are going to annotate your featureData using Bioconductor tools rather than relying on information provided through NCBI GEO. Download times can also be greatly reduced by specifying FALSE.
##' @return Save the files of expression matrix and phenotype information. The pattern of filename are GSE-GPL-phe.txt for phenotype and GSE-GPL-matrix.txt for genotype splilted by tab.
##' @export
##' @importFrom GEOquery getGEO Table
##' @importFrom Biobase pData exprs
##' @importFrom utils write.table

save_GSE <- function(GSE,destdir="temp-data",annotSymbol=FALSE,getGPL=FALSE){
  annotation <- NULL
  #GSE="GSE18508" GSE="GSE128562" GSE="GSE114517"
  if(!dir.exists(destdir)) dir.create(destdir)
  # download expresstion matrix
  gse <- getGEO(GSE, destdir = destdir, getGPL = getGPL)
  for(i in 1:length(gse)){
    eSet <- gse[[i]]
    GPL_id <- eSet@annotation
    phe_filename <- paste0(destdir,"/",GSE,"-",GPL_id,"-phe.txt")
    expr_filename <- paste0(destdir,"/",GSE,"-",GPL_id,"-matrix.txt")
    pdata <- pData(eSet)
    expr_mat <- exprs(eSet)
    if(file.exists(phe_filename) & file.exists(expr_filename)) return("Conversion Done!")
    GPL_data <- eSet@featureData@data
    if(annotSymbol){
      if(nrow(GPL_data)==0) stop("There is no GPL information for this GSE chip.")
      if(ncol(GPL_data)==0) GPL_data = Table(getGEO(GPL_id,getGPL=TRUE))
      probe_symbol <- get_probe_anno(GPL_id=GPL_id,GPL_data=GPL_data) #对探针进行注释
      expr_mat <- convert_probe_symbol_expr(expr_mat,probe_symbol) #把多个探针换成基因
    }
    write.table(pdata,file = phe_filename,sep="\t",quote = TRUE,row.names = FALSE)
    write.table(expr_mat,file = expr_filename,sep="\t",quote = TRUE,row.names = TRUE)
  }
  #invisible(exitStatus)
  return("Conversion Done!")
}

##' @title Convert ExpressionSet with probe to SummarizedExperiment with gene symbol
##' @description After download GSE data using getGEO, we use this function to convert an ExpressionSet with probe to SummarizedExperiment with gene symbol
##' @param eSet A ExpressionSet class
##' @return A SummarizedExperiment class. The pattern of filename are GSE-GPL-phe.txt for phenotype and GSE-GPL-matrix.txt for genotype splilted by tab.
##' @export
##' @importFrom GEOquery getGEO Table
##' @importFrom Biobase pData exprs
##' @importFrom SummarizedExperiment SummarizedExperiment
##' @importFrom S4Vectors DataFrame
convert_probe_se <- function(eSet){
  stopifnot("ExpressionSet" %in% class(eSet))
  if(nrow(exprs(eSet))==0) stop("There is no expression information for this GSE chip.")
  if(nrow(eSet@featureData@data)==0) stop("There is no GPL information for this GSE chip.")
  if(ncol(eSet@featureData@data)==0) eSet@featureData@data = Table(getGEO(eSet@annotation,getGPL=TRUE))
  probe_anno <- get_probe_anno(GPL_id=eSet@annotation,GPL_data=eSet@featureData@data) #对探针进行注释
  expr_mat <- convert_probe_symbol_expr(exprs(eSet),probe_anno) #把多个探针换成基因
  gene_anno <- probe_anno[match(rownames(expr_mat),probe_anno$probeID),]
  row.names(gene_anno) <- gene_anno$symbolID
  row.names(expr_mat) <- gene_anno$symbolID
  se <- SummarizedExperiment(assays = list(counts=expr_mat),
                             rowData = DataFrame(gene_anno),
                             colData = DataFrame(pData(eSet)),
                             metadata = list(experimentData=eSet@experimentData,GPL_id=eSet@annotation))
  return(se)
}


##' @title Get annotation information of probes
##' @description Given a GPL number or GPL data, return the data with first column of probes id and second column of gene symbol.
##'
##' @details nothing
##' @param GPL_id A character, the id number of GPL.
##' @param GPL_data A data.frame, the GPL annotation data from NCBI.
##' @return A data.frame.
##' @export
##' @importFrom GEOquery Table
##' @importFrom stringr str_trim str_split_fixed
get_probe_anno <- function(GPL_id="GPL9061",GPL_data=NULL){
  #有些芯片没有基因信息，因此很难获得。如GSE18508的平台为GPL9061，但是GEO就没有任何GPL信息。GSE10的平台为GPL4，但是没有基因信息。
  # obtain the header of gene symbol in GPL
  if(GPL_id %in% GPL_list$GPL){
    ids <- as.character(GPL_list[GPL_id,2:3]) #选择GPL_id的id and symbol
  } else {
    stop("The gene symbol header of the ",GPL_id," is not in our GPL_list! \n Please add GPL, ID,Gene symbol to the GPL_list variable.")
  }

  #download the data of GPL_id
  if(is.null(GPL_data)) {
    GPL_data = Table(getGEO(GPL_id,getGPL=TRUE))
    if(nrow(GPL_data)==0) stop("There is no GPL_id information for this GSE chip")
  }

  gene_anno <- GPL_data[,c(ids,setdiff(colnames(GPL_data),ids))]
  colnames(gene_anno)[1:2] <- c("probeID","symbolID")

  if(ids[2]=="gene_assignment"){
    symbols <- str_split_fixed(gene_anno$symbolID,"//",3)
    gene_anno$symbolID <- str_trim(symbols[,2])
  }

  row.names(gene_anno) <- gene_anno[,"probeID"]
  return(gene_anno)
}

##' @title Replace the probes id of GPL by gene symbol
##' @description Given a expression matrix and the GPL annotation with probeID and symbolID,
##' find the probe with maximum of total expression from same gene due to one gene have multiple porbes,
##' return the expression data with rowname of gene symbol. The row number of is commonly less than the row number of raw expression matrix.
##'
##' @details nothing
##' @param expr_mat matrix, the raw expression matrix with colname of individual ids and rowname of probe ids.
##' @param probe_anno data.frame with probeID and symbolID, and other information related probes.
##' @return matrix, the clean expression matrix with colname of individual ids and rowname of gene symbol.
convert_probe_symbol_expr <- function(expr_mat,probe_anno){
  print(paste0("The dim of raw expression matrix: number of row is ",nrow(expr_mat),", number of column is ",ncol(expr_mat)))
  tmp <- by(expr_mat,probe_anno[row.names(expr_mat),"symbolID"],function(x) rownames(x)[which.max(rowMeans(x))])
  probes <- as.character(tmp)
  expr_mat <- expr_mat[rownames(expr_mat) %in% probes ,]
  print(paste0("The dim of expression matrix with gene symbol: number of row is ",nrow(expr_mat),", number of column is ",ncol(expr_mat)))
  #rownames(expr_mat) <- probe_anno[match(rownames(expr_mat),probe_anno$probeID),2]
  return(expr_mat)
}


get_GPL_id = function(GSE, destdir = "temp-data"){
  annotation <- NULL
  if (!dir.exists(destdir)) dir.create(destdir)
  gse <- getGEO(GSE, destdir = destdir, getGPL = FALSE)
  gse_gpl <- NULL
  for (i in 1:length(gse)) {
    eSet <- gse[[i]]
    GPL <- eSet@annotation
    gse_gpl <- rbind(gse_gpl,c(GSE,GPL))
  }
  return(gse_gpl)
}

##' @title Get GPL id for multiple GSE
##' @param GSE A character, the number of GSE.
##' @param destdir A character, the path to download GSE related files.
##' @return A data frame with GSE and GPL.
##' @export
GSE_GPL <- function(GSE, destdir = "temp-data"){
  gse_gpl <- lapply(GSE, get_GPL_id, destdir)
  gse_gpl <- do.call(rbind,gse_gpl)
  colnames(gse_gpl) <- c("GSE","GPL")
  return(gse_gpl)
}


