##' @title Download multiple GSE expression matrix and convert simultaneously probe to gene symbol
##' @description Given a vector of GSE number, download related GSE data and convert probe to gene symbol, and save the gene matrix with gene symbol.
##'
##' @details nothing
##'
##' @param GSEs A character vector, the number of multiple GSEs.
##' @param destdir A character, the path to download GSE related files.
##' @return Save the files of expression matrix and phenotype information.
##' @export
##' @importFrom GEOquery getGEO
##' @importFrom Biobase pData exprs

GSEstoExprs = function(GSEs,destdir="tmp"){
  return(0)

}
