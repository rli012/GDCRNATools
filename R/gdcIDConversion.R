
#####################################
######     ID conversion       ######

ensembl2symbolFun <- function(ensemblID, info='symbol') {
    geneInfo <- biotype[match(ensemblID, biotype$ensemblID),]
    geneSymbol <- geneInfo$geneSymbol
    
    if (info=='symbol') {
        return (geneSymbol)
    } else if (info=='all') {
        return (geneInfo)
    }
}


symbol2ensemblFun <- function(symbol, info='ensemblID') {
    geneInfo <- biotype[match(symbol, biotype$geneSymbol),]
    ensemblID <- geneInfo$ensemblID
    
    if (info=='ensemblID') {
        return (ensemblID)
    } else if (info=='all') {
        return (geneInfo)
    }
}
