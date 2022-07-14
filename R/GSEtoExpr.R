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

save_GSE <- function(GSE,destdir="temp-dir",annotSymbol=FALSE,getGPL=FALSE){
  annotation <- NULL
  #GSE="GSE18508" GSE="GSE128562" GSE="GSE114517"
  if(!dir.exists(destdir)) dir.create(destdir)
  # download expresstion matrix
  gse <- getGEO(GSE, destdir = destdir, getGPL = getGPL)
  for(i in 1:length(gse)){
    eSet <- gse[[i]]
    GPL <- eSet@annotation
    phe_filename <- paste0(destdir,"/",GSE,"-",GPL,"-phe.txt")
    expr_filename <- paste0(destdir,"/",GSE,"-",GPL,"-matrix.txt")
    pdata <- pData(eSet)
    expr_mat <- exprs(eSet)
    if(file.exists(phe_filename) & file.exists(expr_filename)) return("Conversion Done!")
    GPL_data <- eSet@featureData@data
    if(annotSymbol){
      if(nrow(GPL_data)==0) stop("There is no GPL information for this GSE chip.")
      if(ncol(GPL_data)==0) GPL_data = Table(getGEO(GPL,getGPL=TRUE))
      probe_symbol <- anno_probe(GPL=GPL,GPL_data=GPL_data) #对探针进行注释
      expr_mat <- convert_probe(expr_mat,probe_symbol) #把多个探针换成基因
    }
    expr_mat <- cbind(Symbol=row.names(expr_mat),expr_mat)
    write.table(pdata,file = phe_filename,sep="\t",quote = TRUE,row.names = FALSE)
    write.table(expr_mat,file = expr_filename,sep="\t",quote = TRUE,row.names = FALSE)
  }
  #invisible(exitStatus)
  return("Conversion Done!")
}

##' @title Get annotation information of probes
##' @description Given a GPL number or GPL data, return the data with first column of probes id and second column of gene symbol.
##'
##' @details nothing
##' @param GPL A character, the id number of GPL.
##' @param GPL_data A data.frame, the GPL annotation data from NCBI.
##' @return A data.frame.
##' @export
##' @importFrom GEOquery Table
##' @importFrom stringr str_trim str_split_fixed
anno_probe <- function(GPL="GPL9061",GPL_data=NA){
  #有些芯片没有基因信息，因此很难获得。如GSE18508的平台为GPL9061，但是GEO就没有任何GPL信息。GSE10的平台为GPL4，但是没有基因信息。
  # obtain the header of gene symbol in GPL
  if(any(GPL %in% GPLlist$GPL)){
    ids <- as.character(GPLlist[GPL,2:3]) #选择GPL的id and symbol
  } else {
    stop("The gene symbol header of the ",GPL," is not in our GPLlist! \n Please add GPL, ID,Gene symbol to GPLlist.")
  }

  # #download the data of GPL
  # if(is.null(nrow(GPL_data))) {
  #   GPL_data = Table(getGEO(GPL,getGPL=TRUE))
  #   if(nrow(GPL_data)==0) stop("There is no GPL information for this GSE chip")
  # }

  anno <- GPL_data[,ids]
  colnames(anno) <- c("probeID","symbolID")

  if(ids[2]=="gene_assignment"){
    symbols <- str_split_fixed(anno$symbolID,"//",3)
    anno$symbolID <- str_trim(symbols[,2])
  }
  row.names(anno) <- anno[,"probeID"]
  return(anno)
}

##' @title Replace the probes id of GPL by gene symbol
##' @description Given a expression matrix and the GPL annotation with probeID and symbolID,
##' find the probe with maximum of total expression from same gene due to one gene have multiple porbes,
##' return the expression data with rowname of gene symbol. The row number of is commonly less than the row number of raw expression matrix.
##'
##' @details nothing
##' @param expr_mat matrix, the raw expression matrix with colname of individual ids and rowname of probe ids.
##' @param probe_symbol data.frame with probeID and symbolID.
##' @return matrix, the clean expression matrix with colname of individual ids and rowname of gene symbol.
convert_probe <- function(expr_mat,probe_symbol){
  print(paste0("The dim of raw expression matrix: number of row is ",nrow(expr_mat),", number of column is ",ncol(expr_mat)))
  tmp <- by(expr_mat,probe_symbol[row.names(expr_mat),"symbolID"],function(x) rownames(x)[which.max(rowMeans(x))])
  probes <- as.character(tmp)
  expr_mat <- expr_mat[rownames(expr_mat) %in% probes ,]
  print(paste0("The dim of expression matrix with gene symbol: number of row is ",nrow(expr_mat),", number of column is ",ncol(expr_mat)))
  rownames(expr_mat) <- probe_symbol[match(rownames(expr_mat),probe_symbol$probeID),2]

  return(expr_mat)
}


get_GPL_id = function(GSE, destdir = "temp-dir"){
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
GSE_GPL <- function(GSE, destdir = "temp-dir"){
  gse_gpl <- lapply(GSE, get_GPL_id, destdir)
  gse_gpl <- do.call(rbind,gse_gpl)
  colnames(gse_gpl) <- c("GSE","GPL")
  return(gse_gpl)
}


