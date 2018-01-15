##' @title Kaplan Meier plot
##' @description Plot Kaplan Meier survival curve
##' @param gene an Ensembl gene id
##' @param rna.expr \code{\link[limma]{voom}} transformed expression data
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @param sep a character string specifying which point should be used to 
##'   separate low-expression and high-expression groups. Possible values are \code{'1stQu'}, 
##'   \code{'mean'}, \code{'median'}, and \code{'3rdQu'}. Default is \code{'median'} 
##' @importFrom survival survfit
##' @importFrom survival survdiff
##' @importFrom survminer ggsurvplot
##' @return A plot of Kaplan Meier survival curve
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### KM plots #######
##' genes <- c('ENSG00000000938','ENSG00000000971','ENSG00000001036',
##'            'ENSG00000001084','ENSG00000001167','ENSG00000001460')
##' 
##' samples <- c('TCGA-2F-A9KO-01',	'TCGA-2F-A9KP-01',
##'              'TCGA-2F-A9KQ-01',	'TCGA-2F-A9KR-01',	
##'              'TCGA-2F-A9KT-01',	'TCGA-2F-A9KW-01')
##'              
##' metaMatrix <- data.frame(sample_type=rep('PrimaryTumor',6),
##'                          sample=samples,
##'                          days_to_death=seq(100,600,100),
##'                          days_to_last_follow_up=rep(NA,6))
##' rnaExpr <- matrix(c(2.7,7.0,4.9,6.9,4.6,2.5,
##'                     0.5,2.5,5.7,6.5,4.9,3.8,
##'                     2.1,2.9,5.9,5.7,4.5,3.5,
##'                     2.7,5.9,4.5,5.8,5.2,3.0,
##'                     2.5,2.2,5.3,4.4,4.4,2.9,
##'                     2.4,3.8,6.2,3.8,3.8,4.2),6,6)
##' rownames(rnaExpr) <- genes
##' colnames(rnaExpr) <- samples
##' gdcKMPlot(gene='ENSG00000000938', rna.expr=rnaExpr, metadata=metaMatrix, sep='median')
gdcKMPlot <- function(gene, rna.expr, metadata, sep='median') {
    metadata <- metadata[metadata$sample_type=='PrimaryTumor',]
    
    samples = intersect(colnames(rna.expr), metadata$sample)
    
    exprDa=rna.expr[gene,samples]
    
    if (sep=='1stQu') {
        thresh <- as.numeric(summary(exprDa)[2])
    } else if (sep=='median') {
        thresh <- as.numeric(summary(exprDa)[3])
    } else if (sep=='mean') {
        thresh <- as.numeric(summary(exprDa)[4])
    } else if (sep=='3rdQu') {
        thresh <- as.numeric(summary(exprDa)[5])
    }
    
    exprGroup <- exprDa > thresh
    
    clinicalDa=metadata[match(samples,metadata$sample),]
    
    daysToDeath <- as.numeric(clinicalDa$days_to_death)
    nonComplt <- is.na(daysToDeath)
    
    vitalStatus <- as.numeric(ifelse(nonComplt, 0, 1))
    daysToDeath[nonComplt] <- as.numeric(clinicalDa$days_to_last_follow_up[nonComplt])
    
    survDa <- data.frame(daysToDeath,vitalStatus, exprGroup)
    
    
    sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ exprGroup)
    pValue <- format(pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE),digits=3)
    
    fit <- survfit(Surv(daysToDeath, vitalStatus) ~ exprGroup, data=survDa)
    xpos = max(daysToDeath, na.rm=TRUE)/2
    ypos = 1.05
    
    ggsurvplot(fit, data=survDa, pval = FALSE, pval.coord = c(2200, 1),
        font.main = c(14, 'bold', 'blue'), conf.int = FALSE, legend = c(0.15, 0.2), 
        legend.labs = c('Low expression', 'High expression'),  legend.title='',
        xlab = 'Overall survival (days)', ylab = 'Survival probability',
        font.x = c(16), font.y = c(16), ylim=c(0,1.05), 
        ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour='black'),
            panel.background = element_blank(),
            legend.text = element_text(size=14))) +
        ggplot2::annotate("text", x = xpos, y = ypos,
            label = paste(ensembl2symbolFun(gene), ' (p=', pValue, ')', sep=''), size = 5.2)
  
}
