##' @title Download GSE expression matrix and convert probe to gene symbol
##' @description Given a GSE number, download related GSE data and convert probe to gene symbol, and return the gene matrix with gene symbol.
##'
##' @details nothing
##'
##' @param GSE a character vactor comparised of GSE and GPL number split by underline. Note if GSE was sequenced on only one platform, please just give GSE number.
##' @return download GSE expression matrix and GPL annotation,
##'         and return a matrix with colanmes of "GSE" and "GPL" and row names of GSE_GPL
##' @export
##' @importFrom GEOquery getGEO

GSEtoExpr = function(GSE,destdir="tmp"){
  #GSE="GSE18508" GSE="GSE128562"
  if(!file.exists(destdir)) dir.create(destdir)
  # download expresstion matrix
  gse = getGEO(GSE, destdir = destdir, getGPL = TRUE)
  for(i in 1:length(gse)){
    eSet = gse[[i]]
    GPL = eSet@annotation
    pheFileName = paste0(destdir,"/",GSE,"-",GPL,"-phe.txt")
    exprFileName = paste0(destdir,"/",GSE,"-",GPL,"-matrix.txt")
    pdata = pData(eSet)
    exprSet = exprs(eSet)
    if(file.exists(pheFileName) & file.exists(exprFileName)) invisible(0)
    GPLdata = eSet@featureData@data
    if(nrow(GPLdata)>0){
      probes_symbol = annoProbe(GPL=GPL,GPLdata=GPLdata) #对探针进行注释
      exprSet = probesToGene(exprSet,probes_symbol) #把多个探针换成基因
    } else stop("There is no GPL information for this GSE chip. The package can not ")

    write.table(pdata,file = pheFileName,sep="\t",quote = FALSE)
    write.table(exprSet,file = exprFileName,sep="\t",quote = FALSE)
  }
  invisible(0)
}

##' @title Get annotation information of probes
##' @description Given a GPL number or GPL data, return the data with first column of probes id and second column of gene symbol.
##'
##' @details nothing
annoProbe = function(GPL="GPL9061",GPLdata=NA){
  #有些芯片没有基因信息，因此很难获得。如GSE18508的平台为GPL9061，但是GEO就没有任何GPL信息。GSE10的平台为GPL4，但是没有基因信息。
  # obtain the header of gene symbol in GPL
  if(any(GPL %in% GPLlist$GPL)){
    ids = as.character(GPLlist[GPL,2:3]) #选择GPL的id and symbol
  } else {
    stop("The gene symbol header of the GPL is not in our GPLlist! \n Please add GPL, ID,Gene symbol to GPLlist.")
  }

  #download the data of GPL
  if(is.null(nrow(GPLdata))) {
    GPLdata = Table(getGEO(GPL,getGPL=TRUE))
    if(nrow(GPLdata)==0) stop("There is no GPL information for this GSE chip")
  }

  anno=GPLdata[,ids]
  colnames(anno)=c("probeID","symbolID")
  return(anno)
}

##' @title Replace the probes id of GPL by gene symbol
##' @description Given a expression matrix and the GPL annotation with probeID and symbolID,
##' find the probe with maximum of total expression from same gene due to one gene have multiple porbes,
##' return the expression data with rowname of gene symbol. The row number of is commonly less than the row number of raw expression matrix.
##'
##' @details nothing
probesToGene = function(exprSet,probes_symbol){
  print(paste0("The dim of raw expression matrix: number of row is ",nrow(exprSet),", number of column is ",ncol(exprSet)))
  tmp = by(exprSet,probes_symbol$symbolID,function(x) rownames(x)[which.max(rowMeans(x))])
  probes = as.character(tmp)
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  print(paste0("The dim of expression matrix wiht gene symbol: number of row is ",nrow(exprSet),", number of column is ",ncol(exprSet)))
  rownames(exprSet)=probes_symbol[match(rownames(exprSet),probes_symbol$probeID),2]

  return(exprSet)
}
