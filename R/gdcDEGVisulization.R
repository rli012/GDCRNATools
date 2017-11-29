##' @title Bar plot of differentially expressed genes/miRNAs
##' @description A bar plot showing the number of down-regulated and up-regulated DE genes/miRNAs of different biotypes
##' @param deg a dataframe generated from \code{\link{gdcDEReport}} containing DE genes/miRNAs ids, logFC, etc.
##' @param angle a numeric value specifying the angle of text on x-axis. Default is \code{0}
##' @param data.type one of \code{'RNAseq'} and \code{'miRNAs'}
##' @return A bar plot
##' @export
##' @author Ruidong Li and Han Qu
gdcBarPlot <- function(deg, angle=0, data.type) {
  
  if (data.type=='miRNAs') {
    down <- sum(deg$logFC < 0)
    up <- sum(deg$logFC > 0)
    
    d <- data.frame(geneClass = c('Up','Down'),
                    geneNums = c(up, down), 
                    Regulation = factor(c('Up-regulated','Down-regulated'),levels=c('Up-regulated','Down-regulated')))
    
    if (angle==0) {
      ggplot(data=d, aes(x=d$geneClass, y=d$geneNums, fill=Regulation)) + geom_bar(stat = 'identity') + 
        scale_x_discrete(limits=d$geneClass) + scale_fill_discrete(name = "") +
        ylab('No. of Differentially Expressed miRNAs') +xlab('') +
        #theme(axis.text.x = element_text(angle = angle)) + 
        theme_bw()+theme(axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour='white'),
                         panel.background = element_blank(),
                         axis.text.x = element_text(angle = angle,size=10),
                         axis.text.y=element_text(size=10))
    } else {
      ggplot(data=d, aes(x=d$geneClass, y=d$geneNums, fill=Regulation)) + geom_bar(stat = 'identity') + 
        scale_x_discrete(limits=d$geneClass) + scale_fill_discrete(name = "") +
        ylab('No. of Differentially Expressed miRNAs') +xlab('') +
        #theme(axis.text.x = element_text(angle = angle, vjust = 1, hjust=1)) +
        theme_bw()+theme(axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour='white'),
                         panel.background = element_blank(),
                         axis.text.x = element_text(angle = angle, vjust = 1, hjust=1,size=10),
                         axis.text.y=element_text(size=10))
    }
    
  } else if (data.type=='RNAseq') {
    gr <- list()
    for (tp in unique(deALL$group)) {
      gr[[tp]][['all']] <- sum(deALL$group==tp)
      gr[[tp]]['up'] <- sum(deALL$group==tp & deALL$logFC > 0) 
      gr[[tp]]['down'] <- sum(deALL$group==tp & deALL$logFC < 0)
      
    }
    
    names(gr)[names(gr) == 'protein_coding'] <- 'Protein coding'
    names(gr)[names(gr) == 'long_non_coding'] <- 'Long non-coding'
    names(gr)[names(gr) == 'pseudogene'] <- 'Pseudogene'
    names(gr)[names(gr) == 'ncRNA'] <- 'Other ncRNA'
    
    
    
    #names(gr)[names(gr) == 'IG'] <- 'Immunoglobulin gene'
    #names(gr)[names(gr) == 'TR'] <- 'T-cell receptor gene'
    #names(gr)[names(gr) == 'TEC'] <- 'To be experimentally confirmed'
    
    o <- order(unlist(lapply(gr, function(v) v[['all']])), decreasing=T)
    
    d <- data.frame(geneClass = rep(names(gr)[o], 2),
                    geneNums = c(unlist(lapply(gr, function(v) v[['up']]))[o],unlist(lapply(gr, function(v) v[['down']]))[o]), 
                    Regulation = factor(rep(c('Up-regulated','Down-regulated'), each=length(gr)), 
                                        levels=c('Up-regulated','Down-regulated')))
    
    if (angle==0) {
      ggplot(data=d, aes(x=d$geneClass, y=d$geneNums, fill=Regulation)) + geom_bar(stat = 'identity') + 
        scale_x_discrete(limits=d$geneClass[1:nrow(d)/2]) + scale_fill_discrete(name = "") +
        ylab('No. of Differentially Expressed Genes') +xlab('Type of genes') +
        #theme(axis.text.x = element_text(angle = angle)) +
        theme_bw()+theme(axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour='white'),
                         panel.background = element_blank(),
                         axis.text.x = element_text(angle = angle,size=10),
                         axis.text.y=element_text(size=10))
    } else {
      ggplot(data=d, aes(x=d$geneClass, y=d$geneNums, fill=Regulation)) + geom_bar(stat = 'identity') + 
        scale_x_discrete(limits=d$geneClass[1:nrow(d)/2]) + scale_fill_discrete(name = "") +
        ylab('No. of Differentially Expressed Genes') +xlab('Type of genes') +
        #theme(axis.text.x = element_text(angle = angle, vjust = 1, hjust=1)) +
        theme_bw()+theme(axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour='white'),
                         panel.background = element_blank(),
                         axis.text.x = element_text(angle = angle, vjust = 1, hjust=1,size=10),
                         axis.text.y=element_text(size=10))
    }
  }
}



##' @title Volcano plot of differentially expressed genes/miRNAs
##' @description A volcano plot showing differentially expressed genes/miRNAs
##' @param deg.all a dataframe generated from \code{\link{gdcDEAnalysis}} containing all genes of analysis no matter they are differentially expressed or not
##' @param fc a numeric value specifying the threshold of fold change
##' @param pval a nuemric value specifying the threshold of p value
##' @return A volcano plot
##' @export
##' @author Ruidong Li and Han Qu
gdcVolcanoPlot <- function(deg.all, fc=2, pval=0.01) {
  geneList <- deg.all
  geneList$threshold <- c()
  geneList$threshold[geneList$logFC>log(fc,2) & geneList$FDR<pval] <- 1
  geneList$threshold[geneList$logFC>=-log(fc,2) & geneList$logFC<=log(fc,2) | geneList$FDR>=pval] <- 2
  geneList$threshold[geneList$logFC < -log(fc,2) & geneList$FDR<pval] <- 3
  
  geneList$threshold <- as.factor(geneList$threshold)
  
  lim <- max(max(geneList$logFC), abs(min(geneList$logFC)))+0.5
  
  volcano <- ggplot(data=geneList, aes(x=geneList$logFC, y = -log10(geneList$FDR)))
  #volcano <- ggplot(data=geneList, aes(x=-log10(geneList$FDR), y = geneList$logFC))
  volcano+geom_point(aes(color=geneList$threshold), alpha=1, size=0.8) + xlab("log2(Fold Change)") + ylab("-log10(FDR)") +
    scale_colour_manual(breaks = geneList$threshold, values = c('red','black','green3')) + xlim(c(-lim,lim)) +
    geom_vline(xintercept = c(-log(fc,2),log(fc,2)), color='darkgreen', linetype=3) + geom_hline(yintercept = -log(pval,10), color='darkgreen',linetype=3)+
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour='black'),
          panel.background = element_blank()) +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
  
}



##' @title Heatmap of differentially expressed genes/miRNAs
##' @description A heatmap showing unsupervised hierarchical clustering of 
##'   DE genes/miRNAs by \code{\link[gplots]{heatmap.2}} in the \pkg{gplots} package
##' @param deg.id a vector of Ensembl gene ids or miRBase v21 mature miRNA ids
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @param rna.expr \code{\link[limma]{voom}} transformed expression data
##' @importFrom gplots heatmap.2
##' @importFrom gplots colorpanel
##' @importFrom gplots bluered
##' @return A heatmap with rows are DE genes/miRNAs and columns are samples.
##'   \emph{Solid Tissue Normal} samples are labeled with blue and \emph{Primary Tumor} samples are labeled with red
##' @export
##' @author Ruidong Li and Han Qu
gdcHeatmap <- function(deg.id, metadata, rna.expr) {
  
  degDa <- rna.expr[deg.id,]
  
  sampleCol <- ifelse(metadata$sample_type=='SolidTissueNormal', 'blue', 'red')
  
  lmat = rbind(c(4,3),c(2,1))
  lwid = c(2,4)
  lhei = c(1,5)
  
  heatmap.2(as.matrix(degDa), col=bluered(75), trace='none', cexCol=0.32, cexRow=0.1,
            dendrogram='both', srtCol=90, adjCol=c(0.8,0.15), density.info="none", labRow=NA,
            key.title=NA,na.color=NA,lwid=lwid, lhei=lhei,  margins =c(3,3), labCol=NA,
            key.xlab='Normalized intensity', scale='row', ColSideColors = sampleCol)
  
}
