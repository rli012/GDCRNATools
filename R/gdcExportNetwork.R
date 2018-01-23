##' @title Export network for Cytoscape
##' @description Export nodes and edges of ce network for 
##'     \strong{Cytoscape} visualization
##' @param ceNetwork a dataframe generated from \code{\link{gdcCEAnalysis}}
##' @param net one of \code{'nodes'} and \code{'edges'}
##' @return A dataframe of nodes or edges
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### ceRNA network analysis #######
##' ceOutput <- data.frame(lncRNAs=c('ENSG00000242125','ENSG00000242125',
##'                                 'ENSG00000245532'), 
##'                     Genes=c('ENSG00000043355','ENSG00000109586',
##'                                 'ENSG00000144355'), 
##'                     miRNAs=c('hsa-miR-340-5p','hsa-miR-340-5p',
##'                             'hsa-miR-320b,hsa-miR-320d,
##'                             hsa-miR-320c,hsa-miR-320a'),
##'                     Counts=c(1,1,4), stringsAsFactors=FALSE)
##' ####### Export edges #######
##' edges <- gdcExportNetwork(ceNetwork=ceOutput, net='edges')
##' 
##' ####### Export nodes #######
##' \dontrun{nodes <- gdcExportNetwork(ceNetwork=ceOutput, net='nodes')}
gdcExportNetwork <- function(ceNetwork, net) {

    mirs <- unlist(strsplit(ceNetwork$miRNAs, ',', fixed=TRUE))
    
    fromNode <- rep(c(ceNetwork$lncRNAs,ceNetwork$Genes), 
        times=rep(as.numeric(ceNetwork$Counts),2))
    toNode <- rep(mirs, 2)
    altNode1Name <- ensembl2symbolFun(fromNode)
    edges <- data.frame(fromNode, toNode, altNode1Name, 
        stringsAsFactors = FALSE)
    
    #filter <- duplicated(edges)
    #edges <- edges[-filter,]
    edges <- unique(edges)
    
    nodeTable1 <- table(edges$fromNode)
    nodeTable2 <- table(edges$toNode)
    
    symbol <- c(ensembl2symbolFun(names(nodeTable1)), names(nodeTable2))
    
    type <- c(ifelse(names(nodeTable1) %in% ceNetwork$lncRNAs, 'lnc', 'pc'), 
        rep('mir', length(nodeTable2)))
    
    
    nodes <- data.frame(gene=c(names(nodeTable1), names(nodeTable2)), 
        symbol, type, numInteractions=as.numeric(c(nodeTable1, nodeTable2)))
    
    if (net=='nodes') {
        return (nodes)
    } else if (net=='edges') {
        return (edges)
    }
}
