##' @title Functional enrichment analysis
##' @description Performs Gene Ontology (GO), Kyoto Encyclopedia of Genes 
##'     and Genomes (KEGG) pathway and Disease Ontology (DO) enrichment 
##'     analyses by \pkg{clusterProfiler} and \pkg{DOSE} packages
##' @param gene a vector of Ensembl gene id
##' @param simplify logical, specifying whether to remove redundant GO terms. 
##'     Default \code{simplify=TRUE}
##' @param level a numeric value, restrict the GO enrichment result at a 
##'     specific GO level. Default is \code{0}, which means all terms 
##'     should be returned
##' @importFrom clusterProfiler enrichGO
##' @importFrom clusterProfiler gofilter
##' @importFrom clusterProfiler simplify
##' @importFrom clusterProfiler enrichKEGG
##' @import org.Hs.eg.db
##' @importFrom DOSE enrichDO
##' @importFrom biomaRt useMart
##' @importFrom biomaRt getBM
##' @return A dataframe of enrichment analysis result containing 
##'     enriched terms, number of overlpped genes, p value of 
##'     hypergeometric test, fdr, fold of enrichment, Ensembl gene ids, 
##'     gene symbols, and functional categories, etc.
##' @references
##'     Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for 
##'     comparing biological themes among gene clusters. 
##'     Omics: a journal of integrative biology. 2012 May 1;16(5):284-7. \cr
##'     Yu G, Wang LG, Yan GR, He QY. DOSE: an R/Bioconductor package for 
##'     disease ontology semantic and enrichment analysis. Bioinformatics. 
##'     2014 Oct 17;31(4):608-9.
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### GO, KEGG, DO enrichment analysis #######
##' deg <- c('ENSG00000000938','ENSG00000000971','ENSG00000001036',
##'         'ENSG00000001084','ENSG00000001167','ENSG00000001460')
##' \dontrun{enrichOutput <- gdcEnrichAnalysis(gene=deg, simplify=TRUE)}
gdcEnrichAnalysis <- function(gene, simplify=TRUE, level=0) {

    message ('### This step may take a few minutes ###\n')
    
    goBP <- enrichGO(gene = gene,
        universe = biotype$ensemblID,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        keyType = 'ENSEMBL',
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.01,
        readable = FALSE)
    
    if (level != 0) {
        goBP <- gofilter(goBP, level=level)
    }
    
    if (simplify==TRUE) {
        goBP <- simplify(goBP, cutoff=0.7, by="p.adjust", select_fun=min)
    }
    
    message ('Step 1/5: BP analysis done!')
    
    goCC <- enrichGO(gene = gene,
        universe = biotype$ensemblID,
        OrgDb = org.Hs.eg.db,
        ont = "CC",
        keyType = 'ENSEMBL',
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.01,
        readable = FALSE)
    
    if (level != 0) {
        goCC <- gofilter(goCC, level=level)
    }
    
    if (simplify==TRUE) {
        goCC <- simplify(goCC, cutoff=0.7, by="p.adjust", select_fun=min)
    }
    
    message ('Step 2/5: CC analysis done!')
    
    goMF <- enrichGO(gene = gene,
        universe = biotype$ensemblID,
        OrgDb = org.Hs.eg.db,
        ont = "MF",
        keyType = 'ENSEMBL',
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.01,
        readable = FALSE)
    
    if (level != 0) {
        goMF <- gofilter(goMF, level=level)
    }
    
    if (simplify==TRUE) {
        goMF <- simplify(goMF, cutoff=0.7, by="p.adjust", select_fun=min)
    }
    
    message ('Step 3/5: MF analysis done!')
    
    genes <- biotype[match(gene, biotype$ensemblID),]
    genes <- genes[! is.na(genes$entrezgene),]
    universe <- biotype[!is.na(biotype$entrezgene),]
    
    #ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    #genes <- getBM(attributes = c('ensembl_gene_id','entrezgene'), 
    #    values=gene, filters='ensembl_gene_id', mart=ensembl)
    #genes <- genes[! is.na(genes$entrezgene),]
    #universe <- getBM(attributes = c('ensembl_gene_id','entrezgene'), 
    #    values=biotype$ensemblID, filters='ensembl_gene_id', mart=ensembl)
    
    kegg <- enrichKEGG(gene = as.character(genes$entrezgene),
        organism = 'hsa',
        universe = as.character(unique(
            universe$entrezgene[!is.na(universe$entrezgene)])),
        minGSSize = 10,
        maxGSSize = 500,
        pAdjustMethod = 'fdr',
        pvalueCutoff = 0.01)
    
    kegg <- data.frame(kegg@result)
    
    kegg$geneID <- unlist(lapply(kegg$geneID, function(v) 
        paste(genes$ensembl_gene_id[match(strsplit(v, '/', fixed=TRUE)[[1]],
            genes$entrezgene)], collapse = '/')))
    
    message ('Step 4/5: KEGG analysis done!')
    
    do <- enrichDO(gene = as.character(genes$entrezgene),
        universe = as.character(unique(
            universe$entrezgene[!is.na(universe$entrezgene)])),
        ont = "DO",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.01,
        readable = FALSE)
    
    do <- data.frame(do@result)
    
    do$geneID <- unlist(lapply(do$geneID, function(v) 
        paste(genes$ensembl_gene_id[match(strsplit(v, '/', fixed=TRUE)[[1]], 
            genes$entrezgene)], collapse = '/')))
    
    message ('Step 5/5: DO analysis done!')
    
    goBP <- organizeEnrichFun(data.frame(goBP@result))
    goCC <- organizeEnrichFun(data.frame(goCC@result))
    goMF <- organizeEnrichFun(data.frame(goMF@result))
    kegg <- organizeEnrichFun(kegg)
    do   <- organizeEnrichFun(do)
    
    enrichOutput <- data.frame(rbind(goBP, goCC, goMF, kegg, do))
    
    enrichOutput$Category <- rep(c('GO_BP','GO_CC','GO_MF','KEGG', 'DO'), 
        c(nrow(goBP),nrow(goCC),nrow(goMF),nrow(kegg),nrow(do)))
    
    return (enrichOutput)
}


###
organizeEnrichFun <- function(go) {

    Terms <- paste(go$ID, go$Description, sep='~')
    Counts <- go$Count
    
    GeneRatio <- go$GeneRatio
    BgRatio <- go$BgRatio
    
    pValue <- go$pvalue
    FDR <- go$p.adjust
    
    listTotal <- vapply(go$GeneRatio, function(v) 
        convertRatioFun(v, type='bg'), numeric(1))
    popHits <- vapply(go$BgRatio, function(v) 
        convertRatioFun(v, type='hit'), numeric(1))
    popTotal <- vapply(go$BgRatio, function(v) 
        convertRatioFun(v, type='bg'), numeric(1))
    
    foldEnrichment <- as.vector(Counts/listTotal*popTotal/popHits)
    
    geneID <- go$geneID
    geneSymbol <- unlist(lapply(strsplit(geneID, '/', fixed=TRUE), 
        function(v) paste(ensembl2symbolFun(v), collapse = '/')))
    
    goOutput <- data.frame(Terms, Counts, GeneRatio, BgRatio, pValue, FDR, 
        foldEnrichment, geneID, geneSymbol)
    
    return (goOutput)
}


###
convertRatioFun <- function(v, type='bg') {
    ratio <- strsplit(v, '/', fixed=TRUE)
    
    if (type=='bg') {
        num <- as.numeric(as.character(ratio[[1]][2]))
    } else if (type=='hit') {
        num <- as.numeric(as.character(ratio[[1]][1]))
    }
    
    return (num)
}
