##' @title Add new annotation information of GPL to GPLlist
##' @description Given GPL, probeID, and symbolID from GPL data, update the GPLlist data.
##'
##' @details nothing
##' @param GPL A character, the number of GSE.
##' @param probeID A character, the header of probe from specific GPL platform.
##' @param symbolID A character, the header of gene symbol from specific GPL platform.
##' @return update the GPLlist.rda in data folder
##' @export
##' @examples
##' \dontrun{
##' GPL="GPL13648"
##' gplDat = Table(getGEO(GPL,getGPL=TRUE,destdir ="tmp"));colnames(gplDat)
##' updateAnno(GPL,probeID="ID",symbolID="Gene Symbol")
##' load("data/GPLlist.rda")
##' }
##' @importFrom usethis use_data

updateAnno = function(GPL,probeID,symbolID){
  rawRowN = nrow(GPLlist)
  if(!GPL %in%  GPLlist$GPL){
    temp=c(GPL,probeID,symbolID)
    GPLlist = rbind(GPLlist,temp)
    row.names(GPLlist)[GPLlist$GPL==GPL]=GPL
    GPLlist = GPLlist[order(GPLlist[,"GPL"]),]
    if(nrow(GPLlist)>rawRowN){
      use_data(GPLlist,overwrite=TRUE)
    }
  }

}
