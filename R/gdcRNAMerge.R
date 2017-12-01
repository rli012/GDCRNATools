##' @title Merge RNA/miRNAs raw counts data
##' @description Merge raw counts data that is downloaded from GDC to a single expression matrix
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @param path path to downloaded files for merging
##' @param data.type one of \code{'RNAseq'} and \code{'miRNAs'}
##' @return A dataframe or numeric matrix of raw counts data with rows are genes or miRNAs and 
##'   columns are samples
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### Merge RNA expression data #######
##' metaMatrix <- gdcParseMetadata(project.id='TARGET-RT', data.type='RNAseq')
##' \dontrun{rnaExpr <- gdcRNAMerge(metadata=metaMatrix, path='RNAseq/', data.type='RNAseq')}
gdcRNAMerge <- function(metadata, path, data.type) {
    
    #if (endsWith(path, '/')) {
    #  path = substr(path, 1, nchar(path)-1)
    #}
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, fsep = .Platform$file.sep)
    
    if (data.type=='RNAseq') {
        cat ('############### Merging RNAseq data ################\n')
        cat ('### This step may take a few minutes ###\n')
        
        rnaMatrix <- do.call("cbind", lapply(filenames, function(fl) read.table(gzfile(fl))$V2))
        rownames(rnaMatrix) <- read.table(gzfile(filenames[1]))$V1
        rownames(rnaMatrix) <- unlist(lapply(strsplit(rownames(rnaMatrix), '.', fixed=TRUE), 
                                             function(gene) gene[1]))
        colnames(rnaMatrix) <- metadata$sample
        
        rnaMatrix <- rnaMatrix[biotype$ensemblID,]
        
        nSamples = ncol(rnaMatrix)
        nGenes = nrow(rnaMatrix)
        
        cat (paste('Number of samples: ', nSamples, '\n', sep=''))
        cat (paste('Number of genes: ', nGenes, '\n', sep=''))
        
        return (rnaMatrix)
    } else if (data.type=='miRNAs') {
        cat ('############### Merging miRNAs data ###############\n')
        
        mirMatrix <- lapply(filenames, function(fl) cleanMirFun(fl))
        #mirs <- sort(unique(names(unlist(mirMatrix))))
        mirs <- rownames(mirbase)
        mirMatrix <- do.call('cbind', lapply(mirMatrix, function(expr) expr[mirs]))
        
        rownames(mirMatrix) <- mirbase$v21[match(mirs,rownames(mirbase))]
        colnames(mirMatrix) <- metadata$sample
        
        mirMatrix[is.na(mirMatrix)] <- 0
        
        nSamples = ncol(mirMatrix)
        nGenes = nrow(mirMatrix)
        
        cat (paste('Number of samples: ', nSamples, '\n', sep=''))
        cat (paste('Number of miRNAs: ', nGenes, '\n', sep=''))
        
        return (mirMatrix)
    } else {
        return ('error !!!')
    }
}



cleanMirFun <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    
    mirs <- unlist(lapply(strsplit(expr$Group.1, ',', fixed=TRUE), function(mir) mir[2]))
    
    expr <- expr[,-1]
    names(expr) <- mirs
    #rownames(expr) <- mirs
    return(expr)
}
