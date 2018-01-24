##' @title Shiny correlation plot
##' @description A simple \pkg{shiny} app to show scatter plot of correlations 
##'     between two genes/miRNAs on local web browser
##' @param gene1 a vector of Ensembl gene ids or miRBase v21 mature miRNA ids
##' @param gene2 a vector of Ensembl gene ids or miRBase v21 mature miRNA ids
##' @param rna.expr \code{\link[limma]{voom}} transformed expression data
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @return a local webpage for visualization of correlation plots
##' @import shiny
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
##' \dontrun{shinyCorPlot(gene1=genes[1:3], gene2=genes[4:5], rna.expr=rnaExpr, 
##'     metadata=metaMatrix)}
shinyCorPlot <- function (gene1, gene2, rna.expr,metadata) {
    ui <- fluidPage(
        headerPanel('Expression correlation of ceRNA pairs'),
        sidebarPanel(
            selectInput('xcol', 'Long non-coding gene', paste(gene1, 
                ' (', ensembl2symbolFun(gene1), ')', sep='')),
            selectInput('ycol', 'Protein coding gene', paste(gene2, 
                ' (', ensembl2symbolFun(gene2), ')', sep=''))
        ),
        mainPanel(
            plotOutput('plot1')
        )
    )
    
    server <- function(input, output, session) {
        output$plot1 <- renderPlot({
            gdcCorPlot(gene1 = strsplit(input$xcol,' (', fixed=TRUE)[[1]][1], 
                gene2 = strsplit(input$ycol,' (', fixed=TRUE)[[1]][1],
                rna.expr = rna.expr,
                metadata = metadata)
        })
    }
    shinyApp(ui = ui, server = server)
}




##' @title Shiny Kaplan Meier (KM) plot
##' @description A simple \pkg{shiny} app to show KM survival curves 
##'     on local web browser
##' @param gene a vector of Ensembl gene ids
##' @param rna.expr \code{\link[limma]{voom}} transformed expression data
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @return a local webpage for visualization of KM plots
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
##' \dontrun{shinyKMPlot(gene=genes, rna.expr=rnaExpr, 
##'     metadata=metaMatrix)}
shinyKMPlot <- function (gene, rna.expr, metadata) {
    sep=c('1stQu','median','mean','3rdQu')
    ui <- fluidPage(
        headerPanel('Kaplan Meier plot'),
        
        sidebarPanel(
            selectInput('xcol', 'Gene', paste(gene, ' (', 
                ensembl2symbolFun(gene), ')', sep='')),
            selectInput('ycol', 'Separator', sep, selected=sep[2])
        ),
        
        mainPanel(
            plotOutput('plot1')
        )
    )
    
    server <- function(input, output, session) {
        output$plot1 <- renderPlot({
            gdcKMPlot(gene = strsplit(input$xcol,' (', fixed=TRUE)[[1]][1], 
                sep = input$ycol,
                rna.expr = rna.expr,
                metadata = metadata)
        })
    }
    shinyApp(ui = ui, server = server)
}




##' @title Shiny pathview
##' @description A simple \pkg{shiny} app to show pathways genetrated by 
##'     \pkg{pathview} package on local web browser
##' @param gene a vector of numeric values (eg. fold change on log2 scale) 
##'     with names are Ensembl gene ids
##' @param pathways a vector of KEGG pathway ids
##' @param directory the folder to save pathway figures. 
##'     Default is the working directory
##' @return a local webpage for visualization of KEGG maps
##' @importFrom pathview pathview
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' genes <- c('ENSG00000000938','ENSG00000000971','ENSG00000001036',
##'         'ENSG00000001084','ENSG00000001167','ENSG00000001460')
##' pathways <- c("hsa05414~Dilated cardiomyopathy (DCM)",
##'             "hsa05410~Hypertrophic cardiomyopathy (HCM)",
##'             "hsa05412~Arrhythmogenic right ventricular cardiomyopathy",
##'             "hsa04512~ECM-receptor interaction",
##'             "hsa04510~Focal adhesion",
##'             "hsa04360~Axon guidance",
##'             "hsa04270~Vascular smooth muscle contraction",
##'             "hsa05205~Proteoglycans in cancer",
##'             "hsa04022~cGMP-PKG signaling pathway",
##'             "hsa00480~Glutathione metabolism")
##' \dontrun{shinyPathview(gene=genes, pathways=pathways)}
shinyPathview <- function(gene, pathways, directory='.') {

    if (! dir.exists(directory)) {
        dir.create(directory)
    }
    
    if ( ! endsWith(directory, '/')) {
        directory = paste(directory, '/', sep='')
    }
    
    ui <- fluidPage(
        headerPanel('Pathview'),
        sidebarPanel(
            selectInput('xcol', 'Pathway', pathways)
        ),
        mainPanel(
            plotOutput('plot1')
        )
    )
    
    server <- function(input, output, session) {
        output$plot1 <- renderImage({
            pathwayID <- strsplit(input$xcol, '~', fixed=TRUE)[[1]][1]
            
            outfile <- paste(directory, pathwayID, '.pathview.png', sep='')
            
            if (! file.exists(outfile)) {
                pathview(gene.data  = gene,
                    pathway.id = pathwayID,
                    species    = "hsa",
                    gene.idtype= 'ENSEMBL',
                    #limit      = list(gene=c(min(geneList),max(geneList)), 
                    #cpd=1),
                    limit      = list(gene=max(abs(gene)), cpd=1),
                    kegg.dir = directory)
                
                file.move(paste(pathwayID, '.pathview.png', sep=''), 
                    directory)
            }
            
            list(src = outfile)
        }, deleteFile = FALSE)
    }
    
    shinyApp(ui, server)
}
