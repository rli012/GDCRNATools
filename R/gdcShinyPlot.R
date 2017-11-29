##' @title Shiny correlation plot
##' @description A simple \pkg{shiny} app to show scatter plot of correlations between two genes/miRNAs on local web browser
##' @param gene1 a vector of Ensembl gene ids or miRBase v21 mature miRNA ids
##' @param gene2 a vector of Ensembl gene ids or miRBase v21 mature miRNA ids
##' @param rna.expr \code{\link[limma]{voom}} transformed expression data
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @import shiny
##' @export
##' @author Ruidong Li and Han Qu
shinyCorPlot <- function (gene1, gene2, rna.expr,metadata) {
  ui <- fluidPage(
    
    headerPanel('Expression correlation of ceRNA pairs'),
    
    sidebarPanel(
      selectInput('xcol', 'Long non-coding gene', paste(gene1, ' (', ensembl2symbolFun(gene1), ')', sep='')),
      selectInput('ycol', 'Protein coding gene', paste(gene2, ' (', ensembl2symbolFun(gene2), ')', sep=''))
    ),
    
    mainPanel(
      plotOutput('plot1')
    )
  )
  
  server <- function(input, output, session) {
    output$plot1 <- renderPlot({
      gdcCorPlot(gene1 = strsplit(input$xcol,' (', fixed=T)[[1]][1], 
                 gene2 = strsplit(input$ycol,' (', fixed=T)[[1]][1],
                 rna.expr = rna.expr,
                 metadata = metadata)
    })
  }
  shinyApp(ui = ui, server = server)
}




##' @title Shiny Kaplan Meier (KM) plot
##' @description A simple \pkg{shiny} app to show KM survival curves on local web browser
##' @param gene a vector of Ensembl gene ids
##' @param rna.expr \code{\link[limma]{voom}} transformed expression data
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @export
##' @author Ruidong Li and Han Qu
shinyKMPlot <- function (gene, rna.expr, metadata) {
  sep=c('1stQu','median','mean','3rdQu')
  
  ui <- fluidPage(
    
    headerPanel('Kaplan Meier plot'),
    
    sidebarPanel(
      selectInput('xcol', 'Gene', paste(gene, ' (', ensembl2symbolFun(gene), ')', sep='')),
      selectInput('ycol', 'Separator', sep, selected=sep[2])
    ),
    
    mainPanel(
      plotOutput('plot1')
    )
  )
  
  server <- function(input, output, session) {
    output$plot1 <- renderPlot({
      gdcKMPlot(gene = strsplit(input$xcol,' (', fixed=T)[[1]][1], 
                sep = input$ycol,
                rna.expr = rna.expr,
                metadata = metadata)
    })
  }
  shinyApp(ui = ui, server = server)
}




##' @title Shiny pathview
##' @description A simple \pkg{shiny} app to show pathways genetrated by \pkg{pathview} package on local web browser
##' @param gene a vector of numeric values (eg. fold change on log2 scale) with names are Ensembl gene ids
##' @param pathways a vector of KEGG pathway ids
##' @param directory the folder to save pathway figures. Default is the working directory
##' @export
##' @author Ruidong Li and Han Qu
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
      pathwayID <- strsplit(input$xcol, '~', fixed=T)[[1]][1]
      
      outfile <- paste(directory, pathwayID, '.pathview.png', sep='')
      
      if (! file.exists(outfile)) {
        pathview(gene.data  = gene,
                 pathway.id = pathwayID,
                 species    = "hsa",
                 gene.idtype= 'ENSEMBL',
                 #limit      = list(gene=c(min(geneList),max(geneList)), cpd=1),
                 limit      = list(gene=max(abs(gene)), cpd=1),
                 kegg.dir = directory)
        
        file.move(paste(pathwayID, '.pathview.png', sep=''), directory)
      }

      list(src = outfile)
    }, deleteFile = FALSE)
  }
  
  shinyApp(ui, server)
}



