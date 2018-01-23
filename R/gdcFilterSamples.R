##' @title Filter out duplicated samples
##' @description Filter out samples that are sequenced for two or more times
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @return A filtered dataframe of metadata without duplicated samples
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### Parse metadata by project id and data type #######
##' metaMatrix <- gdcParseMetadata(project.id='TARGET-RT', data.type='RNAseq')
##' metaMatrix <- gdcFilterDuplicate(metadata=metaMatrix)
gdcFilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
        metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
}


##' @title Filter out other type of samples
##' @description Filter out samples that are neither 
##'     \emph{Solid Tissue Normal} nor \emph{Primary Tumor}
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @return A filtered dataframe of metadata with \emph{Solid Tissue Normal} 
##'     and \emph{Primary Tumor} samples only
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### Parse metadata by project id and data type #######
##' metaMatrix <- gdcParseMetadata(project.id='TARGET-RT', data.type='RNAseq')
##' metaMatrix <- gdcFilterSampleType(metadata=metaMatrix)
gdcFilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
        c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
        metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
}
