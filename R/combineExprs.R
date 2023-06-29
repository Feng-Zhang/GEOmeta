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
##' merge_by_rowname(x,y)
##' merge_by_rowname(x,y,all.x=FALSE)
##' @importFrom utils modifyList

merge_by_rowname <- function(x,y,all=TRUE,...){
  dotargs <- list(...)
  x <- cbind(id=row.names(x),x)
  y <- cbind(id=row.names(y),y)
  z <- do.call("merge",modifyList(list(x=x,y=y,by="id",all=all),dotargs))
  row.names(z) <- z$id
  z$id <- NULL
  return(z)
}

##' @title Download multiple GSE expression matrix and convert simultaneously probe to gene symbol
##' @description Given a vector of GSE number, download related GSE data and convert probe to gene symbol, and save the gene matrix with gene symbol.
##'
##' @details combine_exprs would serach expression matrix file with pattern of ^GSE.*-GPL.*-matrix.txt$ based on the GSE ids in phe argument.
##'
##' @param destdir A character, the path to save expression data, where file name of expression data has pattern ^GSE.*-GPL.*-matrix.txt$ .
##' @param GSEs A vector of GSE character. default to combine all expression data in destdir filefolder.
##' @param is_duplicated A logical. default to duplicated the same features
##' @return Save the files of expression matrix and phenotype information.
##' @export
##' @importFrom stringr str_split_fixed
##' @importFrom utils read.table
##' @importFrom data.table fread
combine_exprs <- function(destdir="tmp",GSEs=NULL,is_duplicated=TRUE){
  res_matrix <- Symbol <- NULL
  expr_filenames <- dir(path=destdir,pattern = "^GSE.*-GPL.*-matrix.txt$")
  if(is.null(GSEs)){GSEs <- unique(str_split_fixed(expr_filenames,"-",2)[,1])}
  for(i in GSEs){
    GSE_names <- expr_filenames[str_detect(expr_filenames,pattern = paste0("^",i))]
    expr <- NULL
    for(j in GSE_names){
      temp <- fread(paste0(path_tail(destdir),j))
      colnames(temp)[1] <- "Symbol"
      if(is_duplicated) temp <- unique_feature(temp,key="Symbol")
      if(j==GSE_names[1])  expr <- temp else {expr= merge(expr,temp,by="Symbol",all=TRUE) }
    }

    if(i==GSEs[1]) res_matrix <- expr else{
      res_matrix  <-  merge(res_matrix,expr,by="Symbol",all=TRUE)
    }
  }
  res_matrix <- res_matrix[res_matrix$Symbol != c("")]
  return(res_matrix)
}

##' @title Convert delimiter to "/" for different system
##' @description Convert delimiter to "/" for different system, and end of "/"
##' @param  path A character for the path
##' @return A character for the path
##' @export
path_tail <- function(path){
  path <- normalizePath(path,winslash = "/")
  path <- paste0(path,"/")
  return(path)
}

##' @title Convert delimiter to "/" for different system
##' @description Convert delimiter to "/" for different system, and end of "/"
##' @param expr A data.table with Symbol and expression level
##' @param key A character for Symbol feature
##' @return A unique data.table
##' @export
unique_feature <- function(expr,key="Symbol"){
  tmp <- by(expr,expr[,key,with = FALSE],function(x) rownames(x)[which.max(rowMeans(x[,-key,with = FALSE]))])
  idx <- as.numeric(tmp)
  return(expr[idx,])
}
