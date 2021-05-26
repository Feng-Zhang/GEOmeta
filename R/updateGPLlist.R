##' @title Add new annotation information of GPL to GPLlist
##' @description Given GPL, probeID, and symbolID from GPL data, update the GPLlist data.
##'
##' @details nothing
##' @param GPL A character, the number of GSE.
##' @param probeID A character, the header of probe from specific GPL platform.
##' @param symbolID A character, the header of gene symbol from specific GPL platform.
##' @return
##' @export
##' @examples
##' \dontrun{
##' GPL="GPL10558"
##' gplDat = Table(getGEO(GPL,getGPL=TRUE));head(gplDat)
##' addAnno(GPL,probeID="ID",symbolID="ILMN_Gene")
##' }

addAnno = function(GPL,probeID,symbolID){
  rawRowN = nrow(GPLlist)
  if(!GPL %in%  GPLlist$GPL){
    temp=c(GPL,probeID,symbolID)
    GPLlist = rbind(GPLlist,temp)
    row.names(GPLlist)[GPLlist$GPL==GPL]=GPL
    if(nrow(GPLlist)>rawRowN){
      usethis::use_data(GPLlist,overwrite=TRUE)
    }
  }

}
