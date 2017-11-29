##' @title Match samples in metadata and expression matrix
##' @description Check if samples in the metadata and expression data match
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @param rna.expr \code{\link[limma]{voom}} transformed expression data
##' @return A logical value. If \code{TRUE}, all the samples matched
##' @export
##' @author Ruidong Li and Han Qu
gdcMatchSamples <- function(metadata, rna.expr) {
  correctNum <- sum(metadata$sample==colnames(rna.expr))
  
  if (correctNum == nrow(metadata) & correctNum == ncol(rna.expr)) {
    cat ('Samples matched')
  }
}
