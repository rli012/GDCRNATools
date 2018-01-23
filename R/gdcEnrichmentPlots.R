##' @title Plots for enrichment analysis
##' @description Bar plot and bubble plot for GO, KEGG, 
##'     and DO functional enrichment analysis
##' @param enrichment a dataframe generated from 
##'     \code{\link{gdcEnrichAnalysis}}
##' @param type type of the plot, should be one of \code{'bar'} 
##'     and \code{'bubble'}
##' @param category which category should be plotted. Possible values are 
##'     \code{'KEGG'}, \code{'GO'}, \code{'GO_BP'}, \code{'GO_CC'}, 
##'     \code{'GO_MF'}, and \code{'DO'}.
##'     Default is \code{'KEGG'}
##' @param num.terms number of terms to be plotted. Default is \code{10}
##' @param bar.color color of the bar plot. Default is \code{'black'}
##' @return A bar plot or bubble plot of functional enrichment analysis
##' @import ggplot2
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### Enrichment plots #######
##' enrichOutput<-data.frame(Terms=c('hsa05414~Dilated cardiomyopathy (DCM)',
##'                                 'hsa04510~Focal adhesion',
##'                                 'hsa05205~Proteoglycans in cancer'),
##'                             Category=rep('KEGG',3), 
##'                             FDR=c(0.001,0.002,0.003))
##' gdcEnrichPlot(enrichment=enrichOutput, type='bar', category='KEGG')
gdcEnrichPlot <- function(enrichment, type='bar', 
    category='KEGG', num.terms=10, bar.color='black') {

    goBP <- enrichment[enrichment$Category=='GO_BP',]
    goCC <- enrichment[enrichment$Category=='GO_CC',]
    goMF <- enrichment[enrichment$Category=='GO_MF',]
    kegg <- enrichment[enrichment$Category=='KEGG',]
    DO   <- enrichment[enrichment$Category=='DO',]
    
    if (nrow(goBP) > num.terms) {
        goBP <- goBP[seq_len(num.terms),]
    }
    
    if (nrow(goCC) > num.terms) {
        goCC <- goCC[seq_len(num.terms),]
    }
    
    if (nrow(goMF) > num.terms) {
        goMF <- goMF[seq_len(num.terms),]
    }
    
    if (nrow(kegg) > num.terms) {
        kegg <- kegg[seq_len(num.terms),]
    }
    
    if (nrow(DO) > num.terms) {
        DO <- DO[seq_len(num.terms),]
    }
    
    
    if (category == 'GO') {
        go <- rbind(goBP, goCC, goMF)
        
        if (type=='bubble') {
            go$Terms <- paste(go$Terms, '[', go$Category, ']', sep='')
            
            enrichBubblePlotFun(go)
        } else if (type=='bar') {
            enrichBarPlotFun(go, type='multiple')
        }
        
    } else if (category %in% c('GO_BP','GO_CC','GO_MF','KEGG', 'DO')) {
        goList <- list(goBP, goCC, goMF, kegg, DO)
        names(goList) <- c('GO_BP','GO_CC','GO_MF', 'KEGG', 'DO')
        go <- goList[[category]]
        
        if (type=='bubble') {
            enrichBubblePlotFun(go)
        } else if (type=='bar') {
            enrichBarPlotFun(go, bar.color=bar.color)
        }
        
    }

}


enrichBubblePlotFun <- function(kegg, pval=0.01) {
    keggBasic = ggplot(data=kegg, mapping=aes(x=Terms, 
        y=foldEnrichment,color=FDR,size=Counts))
    keggBasic+geom_point()+ coord_flip() + 
        scale_x_discrete(limits=rev(kegg$Terms)) + 
        scale_colour_gradientn(limits=c(0,pval),
            colors= c("red","yellow","green")) + 
        xlab('')+ylab('Fold enrichment') +
        guides(shape = guide_legend(order=1),
            colour = guide_colourbar(order=2)) + 
        theme_bw()+theme(axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour='black'),
            panel.background = element_blank()) +
        ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size=20)) +
        theme(axis.text.y=element_text(size=14), 
            axis.title=element_text(size=15),
            axis.text.x=element_text(size=14)) + 
        theme(legend.text = element_text(size = 14),
            legend.title = element_text(size = 14)) +
        theme(strip.text = element_text(size = 14), 
            legend.key.size = unit(0.8,'cm'))
}


enrichBarPlotFun <- function(kegg, type='single', bar.color='black') {
    
    if (type=='single') {
        keggBasic = ggplot(data=kegg, mapping=aes(x=Terms, y=-log(FDR,10)))
        keggBasic + geom_bar(stat='identity',fill=bar.color) + 
            scale_x_discrete(limits=rev(kegg$Terms)) + 
            ylim(0, max(-log(kegg$FDR,10))) + 
            theme(legend.title=element_blank())+ylab('-log10(FDR)')+
            xlab('') + coord_flip() + 
            theme_bw()+theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_rect(colour='white'),
                panel.background = element_blank()) +
            theme(axis.text.y=element_text(size=14), 
                axis.title=element_text(size=15),
                axis.text.x=element_text(size=14)) + 
            theme(legend.position = 'none')
        
    } else if (type=='multiple') {
        keggBasic = ggplot(data=kegg, mapping=aes(x=Terms, y=-log(FDR,10), 
            fill=Category))
        keggBasic + geom_bar(stat='identity') + 
            scale_x_discrete(limits=rev(kegg$Terms)) + 
            ylim(0, max(-log(kegg$FDR,10))) + 
            theme(legend.title=element_blank())+
            ylab('-log10(FDR)')+xlab('') + coord_flip() + 
            scale_fill_hue(name='',breaks=kegg$Category, 
                labels=kegg$Category) +
            theme_bw()+theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_rect(colour='white'),
                panel.background = element_blank()) +
            theme(axis.text=element_text(size=14), 
                axis.title=element_text(size=15)) + 
            theme(legend.text = element_text(size=16))
    }
}
