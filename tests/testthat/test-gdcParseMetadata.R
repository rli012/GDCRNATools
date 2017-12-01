library(GDCRNATools)

context("gdcParseMetadata")

test_that("Parse metadata automatically", {
    metaMatrix <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)
    expect_equal(unique(metaMatrix$project),"TCGA-CHOL")
})

