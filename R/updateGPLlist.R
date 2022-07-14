##' @title Add new annotation information of GPL to GPL_list
##' @description Given GPL, probeID, and symbolID from GPL data, update the GPL_list data.
##'
##' @details nothing
##' @param GPL A character, the number of GSE.
##' @param probeID A character, the header of probe from specific GPL platform.
##' @param symbolID A character, the header of gene symbol from specific GPL platform.
##' @param overwrite By default, use_data() will not overwrite existing files. If you really want to do so, set this to TRUE.
##' @return update the GPL_list.rda in data folder
##' @export
##' @examples
##' \dontrun{
##' GPL="GPL13648"
##' gplDat = Table(getGEO(GPL,getGPL=TRUE,destdir ="tmp"));colnames(gplDat)
##' updateAnno(GPL,probeID="ID",symbolID="Gene Symbol")
##' load("data/GPL_list.rda")
##' }
##' @importFrom usethis use_data

updateAnno <- function(GPL,probeID,symbolID,overwrite=FALSE){
  rawRowN <- nrow(GPL_list)
  if(!GPL %in%  GPL_list$GPL){
    temp <- c(GPL,probeID,symbolID)
    GPL_list <- rbind(GPL_list,temp)
    row.names(GPL_list)[GPL_list$GPL==GPL] <- GPL
    GPL_list <- GPL_list[order(GPL_list[,"GPL"]),]
    if(nrow(GPL_list)>rawRowN & overwrite){
      use_data(GPL_list,overwrite=overwrite)
    }
  }
  return(GPL_list)
}
