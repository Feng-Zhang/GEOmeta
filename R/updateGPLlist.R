##' @title Add new annotation information of GPL to GPLlist
##' @description Given GPL, probeID, and symbolID from GPL data, update the GPLlist data.
##'
##' @details nothing

addAnno = function(GPL,probeID,symbolID){
  # GPL="GPL10558"
  # gplDat = Table(getGEO(GPL,getGPL=TRUE));head(gplDat)
  # addAnno(GPL,probeID="ID",symbolID="ILMN_Gene")
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
