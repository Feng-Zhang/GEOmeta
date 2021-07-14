##' @title Merge two matrix by row names
##' @description Given two matrix with columns and row names, merge them by row names.
##' The function derives from merge function in base package, so it can receive arguments from merge.
##'
##' @details nothing
##'
##' @param x,y matrix, or objects to be coerced to one.
##' @param all logical;all = L is shorthand for all.x = L and all.y = L, where L is either TRUE or FALSE.
##' @param ... 	arguments to be passed to or from methods.
##' @return The merged matrix.
##' @export
##' @examples
##' x=matrix(1:9,nrow=3,dimnames = list(c("row1","row2","row3"),c("col1","col2","col3")))
##' y=matrix(1:8,nrow=4,dimnames = list(c("row2","row3","row4","row5"),c("col1","col5")))
##' mergeByRowname(x,y)
##' mergeByRowname(x,y,all.x=FALSE)
##' @importFrom utils modifyList

mergeByRowname = function(x,y,all=TRUE,...){
  dotargs = list(...)
  x = cbind(id=row.names(x),x)
  y = cbind(id=row.names(y),y)
  z = do.call("merge",modifyList(list(x=x,y=y,by="id",all=all),dotargs))
  row.names(z)=z$id
  z$id = NULL
  return(z)
}

##' @title Download multiple GSE expression matrix and convert simultaneously probe to gene symbol
##' @description Given a vector of GSE number, download related GSE data and convert probe to gene symbol, and save the gene matrix with gene symbol.
##'
##' @details combineExprs would serach expression matrix file with pattern of ^GSE.*-GPL.*-matrix.txt$ based on the GSE ids in phe argument.
##'
##' @param destdir A character, the path to save expression data, where file name of expression data has pattern ^GSE.*-GPL.*-matrix.txt$ .
##' @param GSEs A vector of GSE character. default to combine all expression data in destdir filefolder.
##' @return Save the files of expression matrix and phenotype information.
##' @export
##' @importFrom  stringr str_split_fixed
##' @importFrom  utils read.table
combineExprs = function(destdir="tmp",GSEs=NULL){
  allMatrix <- Symbol <- NULL
  exprFiles = dir(path=destdir,pattern = "^GSE.*-GPL.*-matrix.txt$")
  if(is.null(GSEs)){GSEs = unique(str_split_fixed(exprFiles,"-",2)[,1])}
  for(i in GSEs){
    GSEnames = exprFiles[str_detect(exprFiles,pattern = paste0("^",i))]
    expr = NULL
    for(j in GSEnames){
      temp = fread(paste0(pathTail(destdir),j))
      colnames(temp)[1]="Symbol"
      if(j==GSEnames[1])  expr = temp else {expr= merge(expr,temp,by="Symbol",all=TRUE) }
    }

    if(i==GSEs[1]) allMatrix=expr else{
      allMatrix = merge(allMatrix,expr,by="Symbol",all=TRUE)
    }
  }
  allMatrix = allMatrix[allMatrix$Symbol != c("")]
  return(allMatrix)
}

##' @title Convert delimiter to "/" for different system
##' @description Convert delimiter to "/" for different system, and end of "/"
##' @param  path A character for the path
##' @return A character for the path
##' @export
pathTail = function(path){
  path = normalizePath(path,winslash = "/")
  path = paste0(path,"/")
  return(path)
}
