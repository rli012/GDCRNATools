##' @title Correlation plot of two genes/miRNAs
##' @description Scatter plot showing the expression correlation 
##'     between two genes/miRNAs
##' @param gene1 an Ensembl gene id or miRBase v21 mature miRNA id
##' @param gene2 an Ensembl gene id or miRBase v21 mature miRNA id
##' @param rna.expr \code{\link[limma]{voom}} transformed expression data
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @return A scatter plot with line of best fit
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' genes <- c('ENSG00000000938','ENSG00000000971','ENSG00000001036',
##'         'ENSG00000001084','ENSG00000001167','ENSG00000001460')
##' 
##' samples <- c('TCGA-2F-A9KO-01', 'TCGA-2F-A9KP-01', 
##'             'TCGA-2F-A9KQ-01', 'TCGA-2F-A9KR-11', 
##'             'TCGA-2F-A9KT-11', 'TCGA-2F-A9KW-11')
##'              
##' metaMatrix <- data.frame(sample_type=rep(c('PrimaryTumor',
##'                         'SolidTissueNormal'),each=3),
##'                         sample=samples,
##'                         days_to_death=seq(100,600,100),
##'                         days_to_last_follow_up=rep(NA,6))
##' 
##' rnaExpr <- matrix(c(2.7,7.0,4.9,6.9,4.6,2.5,
##'                     0.5,2.5,5.7,6.5,4.9,3.8,
##'                     2.1,2.9,5.9,5.7,4.5,3.5,
##'                     2.7,5.9,4.5,5.8,5.2,3.0,
##'                     2.5,2.2,5.3,4.4,4.4,2.9,
##'                     2.4,3.8,6.2,3.8,3.8,4.2),6,6)
##' rownames(rnaExpr) <- genes
##' colnames(rnaExpr) <- samples
##' gdcCorPlot(gene1 = 'ENSG00000000938', 
##'         gene2    = 'ENSG00000001084',
##'         rna.expr = rnaExpr,
##'         metadata = metaMatrix)
gdcCorPlot <- function(gene1, gene2, rna.expr, metadata) {

    samples = intersect(colnames(rna.expr), metadata$sample)
    
    lncDa=rna.expr[gene1,samples]
    pcDa=rna.expr[gene2,samples]
    sampleType=as.factor(metadata$sample_type)
    
    x <- cor.test(x=lncDa, y=pcDa, alternative = 'greater')
    
    c <- format(x$estimate, digits=3)
    p <- format(x$p.value, digits=3)
    
    corDa <- data.frame(lncDa=rna.expr[gene1,], pcDa=rna.expr[gene2,], 
        sampleType=as.factor(metadata$sample_type))
    
    
    xpos <- (min(lncDa)+max(lncDa))/2
    ypos <- as.numeric(summary(pcDa)[6])+0.5
    
    
    ggplot(corDa, aes(x=corDa$lncDa, y=corDa$pcDa)) + 
        geom_point(aes(shape=sampleType, color=sampleType)) + 
        xlab(paste(gene1,' (',ensembl2symbolFun(gene1),')',sep='')) +
        ylab(paste(gene2,' (',ensembl2symbolFun(gene2),')',sep='')) + 
        geom_smooth(method="lm",se=FALSE, col='darkgreen', size=0.5) + 
        scale_colour_manual(breaks = sampleType, 
            values = c('chocolate1', 'blue')) +
        ggplot2::annotate("text", x = xpos, y = ypos, 
            label = paste('cor=', c, ', p=', p, sep=''), size = 5) +
        theme_bw()+theme(legend.title = element_blank(),
            legend.text = element_text(size=14),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour='white'),
            panel.background = element_blank(),
            axis.text = element_text(size=14),
            axis.title = element_text(size=14))
}
