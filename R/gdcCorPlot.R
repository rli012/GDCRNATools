##' @title Correlation plot of two genes/miRNAs
##' @description Scatter plot showing the expression correlation between two genes/miRNAs
##' @param gene1 an Ensembl gene id or miRBase v21 mature miRNA id
##' @param gene2 an Ensembl gene id or miRBase v21 mature miRNA id
##' @param rna.expr \code{\link[limma]{voom}} transformed expression data
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @return A scatter plot with line of best fit
##' @export
##' @author Ruidong Li and Han Qu
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
    scale_colour_manual(breaks = sampleType, values = c('chocolate1', 'blue')) +
    ggplot2::annotate("text", 
                      x = xpos, y = ypos, # x and y coordinates of the text
                      label = paste('cor=', c, ', p=', p, sep=''), size = 5) +
    theme_bw()+theme(legend.title = element_blank(),
                     legend.text = element_text(size=12),
                     axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_rect(colour='white'),
                     panel.background = element_blank(),
                     axis.text = element_text(size=10),
                     axis.title = element_text(size=14))
  
}
