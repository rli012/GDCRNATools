##' @title Export network for Cytoscape
##' @description Export nodes and edges of ce network for \strong{Cytoscape} visualization
##' @param ceNetwork a dataframe generated from \code{\link{gdcCEAnalysis}}
##' @param net one of \code{'nodes'} and \code{'edges'}
##' @return A dataframe of nodes or edges
##' @export
##' @author Ruidong Li and Han Qu
gdcExportNetwork <- function(ceNetwork, net) {
  
  mirs <- unlist(strsplit(ceNetwork$miRNAs, ',', fixed=T))
  
  fromNode <- rep(c(ceNetwork$lncRNAs,ceNetwork$Genes), times=rep(as.numeric(ceNetwork$Counts),2))
  toNode <- rep(mirs, 2)
  altNode1Name <- ensembl2symbolFun(fromNode)
  edges <- data.frame(fromNode, toNode, altNode1Name, stringsAsFactors = F)
  
  #filter <- duplicated(edges)
  #edges <- edges[-filter,]
  edges <- unique(edges)
  
  nodeTable1 <- table(edges$fromNode)
  nodeTable2 <- table(edges$toNode)
  
  symbol <- c(ensembl2symbolFun(names(nodeTable1)), names(nodeTable2))
  
  type <- c(ifelse(names(nodeTable1) %in% ceNetwork$lncRNAs, 'lnc', 'pc'), rep('mir', length(nodeTable2)))
  
  
  nodes <- data.frame(gene=c(names(nodeTable1), names(nodeTable2)), symbol, type, 
                      numInteractions=as.numeric(c(nodeTable1, nodeTable2)))
  
  if (net=='nodes') {
    return (nodes)
  } else if (net=='edges') {
    return (edges)
  }
  
}
