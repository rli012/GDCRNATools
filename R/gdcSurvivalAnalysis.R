##' @title Univariate survival analysis of multiple genes
##' @description 
##'     Univariate Cox Proportional-Hazards and Kaplan Meier 
##'     survival analysis of a vector of genes
##' @param gene a vector of Ensembl gene ids
##' @param rna.expr \code{\link[limma]{voom}} transformed expression data
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @param method method for survival analysis. Possible values are 
##'     \code{'coxph'} and \code{'KM'}. Default is \code{'coxph'}
##' @param sep which point should be used to separate low-expression 
##'     and high-expression groups for \code{method='KM'}.
##'     Possible values are \code{'1stQu'}, \code{'mean'}, 
##'     \code{'median'}, and \code{'3rdQu'}. Default is \code{'median'} 
##' @importFrom survival coxph
##' @importFrom survival Surv
##' @return A dataframe or numeric matrix of hazard ratio, 95\% 
##'     confidence interval, p value, and FDR
##' @references
##'     Therneau TM, Lumley T. Package ‘survival’. \cr
##'     Andersen PK, Gill RD. Cox's regression model for 
##'     counting processes: a large sample study. 
##'     The annals of statistics. 1982 Dec 1:1100-20. \cr
##'     Therneau TM, Grambsch PM. Extending the Cox model. 
##'     Edited by P. Bickel, P. Diggle, S. Fienberg, K. Krickeberg. 
##'     2000:51. \cr
##'     Harrington DP, Fleming TR. A class of rank test procedures 
##'     for censored survival data.Biometrika. 1982 Dec 1;69(3):553-66.
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' genes <- c('ENSG00000000938','ENSG00000000971','ENSG00000001036',
##'         'ENSG00000001084','ENSG00000001167','ENSG00000001460')
##' 
##' samples <- c('TCGA-2F-A9KO-01', 'TCGA-2F-A9KP-01',
##'             'TCGA-2F-A9KQ-01', 'TCGA-2F-A9KR-01',
##'             'TCGA-2F-A9KT-01', 'TCGA-2F-A9KW-01')
##'              
##' metaMatrix <- data.frame(sample_type=rep('PrimaryTumor',6),
##'                         sample=samples,
##'                         days_to_death=seq(100,600,100),
##'                         days_to_last_follow_up=rep(NA,6))
##' rnaExpr <- matrix(c(2.7,7.0,4.9,6.9,4.6,2.5,
##'                     0.5,2.5,5.7,6.5,4.9,3.8,
##'                     2.1,2.9,5.9,5.7,4.5,3.5,
##'                     2.7,5.9,4.5,5.8,5.2,3.0,
##'                     2.5,2.2,5.3,4.4,4.4,2.9,
##'                     2.4,3.8,6.2,3.8,3.8,4.2),6,6)
##' rownames(rnaExpr) <- genes
##' colnames(rnaExpr) <- samples
##' survOutput <- gdcSurvivalAnalysis(gene=genes,
##'     rna.expr=rnaExpr, metadata=metaMatrix)
gdcSurvivalAnalysis <- function(gene, rna.expr, metadata, 
    method='coxph', sep='median') {

    if (method=='coxph') {
        survOutput <- coxphTestFun(gene, rna.expr, metadata)
    } else if (method=='KM') {
        survOutput <- kmTestFun(gene, rna.expr, metadata, sep=sep)
    }
    
    return (survOutput)
}


###
coxphTestFun <- function(genes, rna.expr, metaMatrix) {

    metaMatrix <- metaMatrix[metaMatrix$sample_type=='PrimaryTumor',]
    
    samples = intersect(colnames(rna.expr), metaMatrix$sample)
    exprDa=rna.expr[genes,samples]
    
    clinicalDa=metaMatrix[match(samples,metaMatrix$sample),]
    daysToDeath <- as.numeric(clinicalDa$days_to_death)
    nonComplt <- is.na(daysToDeath)
    
    vitalStatus <- as.numeric(ifelse(nonComplt, 0, 1))
    daysToDeath[nonComplt] <- 
        as.numeric(clinicalDa$days_to_last_follow_up[nonComplt])
    
    coxphDEGs <- c()
    for (i in seq_len(nrow(exprDa))) {
        DEG <- unlist(exprDa[i,])
        coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ DEG)
        
        summcph <- summary(coxtest)
        coeffs <- c(summcph$coefficients[,1:2], summcph$conf.int[,3:4], 
            summcph$coefficients[,5])
        coxphDEGs <- rbind(coxphDEGs, coeffs)
        
    }
    rownames(coxphDEGs) <- rownames(exprDa)
    
    colnames(coxphDEGs) <- c('coef','HR','lower95','upper95','pValue')
    coxphDEGs <- data.frame(symbol=ensembl2symbolFun(rownames(exprDa)),
        coxphDEGs)
    #coxphDEGs$FDR <- p.adjust(coxphDEGs$pValue, method='fdr')
    
    #o <- order(coxphDEGs$pValue)
    #coxphDEGs <- coxphDEGs[o,]
    
    return (coxphDEGs)
}


###
kmTestFun <- function(genes, rna.expr, metaMatrix, sep='median') {
    metaMatrix <- metaMatrix[metaMatrix$sample_type=='PrimaryTumor',]
    
    samples = intersect(colnames(rna.expr), metaMatrix$sample)
    exprDa=rna.expr[genes,samples]
    
    clinicalDa=metaMatrix[match(samples,metaMatrix$sample),]
    daysToDeath <- as.numeric(clinicalDa$days_to_death)
    nonComplt <- is.na(daysToDeath)
    
    vitalStatus <- as.numeric(ifelse(nonComplt, 0, 1))
    daysToDeath[nonComplt] <- 
        as.numeric(clinicalDa$days_to_last_follow_up[nonComplt])
    
    kmDEGs <- c()
    for (i in seq_len(nrow(exprDa))) {
        DEG <- unlist(exprDa[i,])
        
        if (sep=='1stQu') {
            thresh <- as.numeric(summary(DEG)[2])
        } else if (sep=='median') {
            thresh <- as.numeric(summary(DEG)[3])
        } else if (sep=='mean') {
            thresh <- as.numeric(summary(DEG)[4])
        } else if (sep=='3rdQu') {
            thresh <- as.numeric(summary(DEG)[5])
        }
        
        exprGroup <- DEG > thresh
        
        sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ exprGroup)
        pValue <- format(pchisq(sdf$chisq, length(sdf$n)-1, 
            lower.tail = FALSE),digits=3)
        #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
        
        HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
        upper95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
        lower95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
        
        kmDEGs <- rbind(kmDEGs, c(HR, lower95, upper95, pValue))
        
    }
    
    rownames(kmDEGs) <- rownames(exprDa)
    colnames(kmDEGs) <- c('HR','lower95','upper95','pValue')
    kmDEGs <- data.frame(symbol=ensembl2symbolFun(rownames(exprDa)), kmDEGs)
    #kmDEGs$FDR <- p.adjust(kmDEGs$pValue, method='fdr')
    
    #o <- order(coxphDEGs$pValue)
    #coxphDEGs <- coxphDEGs[o,]
    
    return (kmDEGs)
}
