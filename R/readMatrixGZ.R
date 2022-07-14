# ##' @title Read expression matrix in compressed format
# ##' @description nothing
# ##'
# ##' @details nothing
# ##'
# ##' @param GSE a character of GSE number
# ##' @return a GSE expression dataframe
# ##' @export
# ##' @importFrom stringr str_detect
# ##' @importFrom data.table fread
#
#
# readMatrixGZ = function(GSE){
#
#
#   #GSE="GSE151839"
#   # get the number of redundant line
#   exprSet = readLines(paste(GSE,"_series_matrix.txt.gz",sep=""))
#   skipNum = which(str_detect(exprSet,"ID_REF"))-1
#
#   # read matrix
#   exprSet = try(fread(paste(GSE,"_series_matrix.txt.gz",sep=""),header = TRUE,skip=skipNum),silent=T)
#   if("try-error" %in% class(exprSet)){
#     exprSet= read.table(paste(GSE,"_series_matrix.txt.gz",sep=""),header=TRUE,skip=skipNum,sep="\t",comment.char = "!",row.names = "ID_REF")
#   } else {
#     probes = exprSet$ID_REF
#     exprSet = as.data.frame(exprSet[,-1])
#     row.names(exprSet) = probes
#   }
#
#   # replace null using 0
#   if(any(as.vector(exprSet)=="null",na.rm = T)){
#     temp=row.names(exprSet)
#     exprSet = apply(exprSet,2,function(x) {x[x=="null"]=0 ; as.numeric(x)})
#     exprSet = as.data.frame(exprSet)
#     row.names(exprSet) = temp
#   }
#   return(exprSet)
# }
#
