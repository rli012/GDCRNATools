##' @title Download RNA data in GDC
##' @description Download gene expression quantification 
##'     and isoform expression quantification data from GDC 
##'     either by providing the manifest file or by sepcifying 
##'     the project id and data type
##' @param manifest menifest file that is downloaded from the GDC cart. 
##'     If provided, files whose UUIDs are in the manifest file 
##'     will be downloaded via gdc-client, otherwise, \code{project} and 
##'     \code{data.type} arguments should be provided to download 
##'     data automatically. Default is \code{NULL}
##' @param project.id project id in GDC
##' @param data.type one of \code{'RNAseq'} and \code{'miRNAs'}
##' @param directory the folder to save downloaded files. 
##'     Default is \code{'Data'}
##' @param write.manifest logical, whether to write out the manifest file
##' @param method method that is used to download data. Either 
##'     \code{'GenomicDataCommons'} which is a well established method 
##'     developed in the \pkg{GenomicDataCommons'} package, or alternatively 
##'     \code{'gdc-client'} which uses the \code{gdc-client} tool developed 
##'     by GDC. Default is \code{'gdc-client'}.
##' @importFrom GenomicDataCommons gdcdata
##' @importFrom DT datatable
##' @return Downloaded files in the specified directory
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### Download RNA data by menifest file #######
##' manifest <- 'RNAseq.manifest.txt'
##' \dontrun{gdcRNADownload(manifest=manifest)}
##'
##' ####### Download RNA data by project id and data type #######
##' project <- 'TCGA-PRAD'
##' \dontrun{gdcRNADownload(project.id=project, data.type='RNAseq')}
gdcRNADownload <- function(manifest=NULL, project.id, data.type, 
    directory='Data', write.manifest=FALSE, method='gdc-client') {

    if (! is.null(manifest)) {
        manifestDownloadFun(manifest=manifest,directory=directory)
        
    } else {
        
        url <- gdcGetURL(project.id=project.id, data.type=data.type)
        manifest <- read.table(paste(url, '&return_type=manifest', sep=''), 
            header=TRUE, stringsAsFactors=FALSE)
        
        systime <- gsub(' ', 'T', Sys.time())
        systime <- gsub(':', '-', systime)
        
        manifile <- paste(project.id, data.type, 'gdc_manifest', 
            systime, 'txt', sep='.')
        write.table(manifest, file=manifile, row.names=FALSE, 
            sep='\t', quote=FALSE)
        
        if (method=='GenomicDataCommons') {
            ex <- manifest$filename %in% dir(directory)
            nonex <- ! ex
            numFiles <- sum(ex)
            
            if(numFiles > 0) {
                message (paste('Already exists',numFiles,'files !',sep=' '))
                
                if (sum(nonex) > 0 ) {
                    message (paste('Download the other', 
                        sum(nonex), 'files !', sep=' '))
                    manifest <- manifest[nonex,]
                    fnames = lapply(manifest$id,gdcdata,
                        destination_dir=directory,overwrite=TRUE,
                        progress=TRUE)
                } else {
                    return(invisible())
                }
                
            } else {
                fnames = lapply(manifest$id,gdcdata,
                    destination_dir=directory,overwrite=TRUE,
                    progress=TRUE)
            }
            
        } else if (method=='gdc-client') {
            manifestDownloadFun(manifest=manifile,directory=directory)
        }
        
        if (write.manifest == FALSE) {
            invisible(file.remove(manifile))
        }
    }
}



###
downloadClientFun <- function(os) {
    if (os == 'Linux') {
        adress <- paste('https://gdc.cancer.gov/system/files/',
            'authenticated%20user/0/gdc-client_v1.3.0_Ubuntu14.04_x64.zip',
            sep='')
        download.file(adress, 
            destfile = './gdc-client_v1.3.0_Ubuntu14.04_x64.zip')
        unzip('./gdc-client_v1.3.0_Ubuntu14.04_x64.zip')
        
    } else if ( os== 'Windows') {
        adress <- paste('https://gdc.cancer.gov/system/files/',
            'authenticated%20user/0/gdc-client_v1.3.0_Windows_x64.zip',
            sep='')
        download.file(adress, 
            destfile = './gdc-client_v1.3.0_Windows_x64.zip')
        unzip('./gdc-client_v1.3.0_Windows_x64.zip')
        
    } else if (os == 'Darwin') {
        adress <- paste('https://gdc.cancer.gov/system/files/',
            'authenticated%20user/0/gdc-client_v1.3.0_OSX_x64.zip',
            sep='')
        download.file(adress, destfile = './gdc-client_v1.3.0_OSX_x64.zip')
        unzip('./gdc-client_v1.3.0_OSX_x64.zip')
    }
}


###
manifestDownloadFun <- function(manifest=manifest,directory) {

    ### download gdc-client
    if (! file.exists('gdc-client') & !file.exists('gdc-client.exe')) {
        downloadClientFun(Sys.info()[1])
    }
    
    Sys.chmod('gdc-client')
    
    manifestDa <- read.table(manifest, sep='\t', header=TRUE, 
        stringsAsFactors = FALSE)
    ex <- manifestDa$filename %in% dir(paste(directory, 
        dir(directory), sep='/'))
    nonex <- ! ex
    numFiles <- sum(ex)
    
    if(numFiles > 0) {
        message (paste('Already exists', numFiles, 'files !', sep=' '))
        
        if (sum(nonex) > 0 ) {
            message (paste('Download the other', 
                sum(nonex), 'files !', sep=' '))
            
            manifestDa <- manifestDa[nonex,]
            manifest <- paste(manifestDa$id, collapse =' ')
            system(paste('./gdc-client download ', manifest, sep=''))
        } else {
            return(invisible())
        }
        
        
    } else {
        system(paste('./gdc-client download -m ', manifest, sep=''))
    }
    
    
    #### move to the directory
    files <- manifestDa$id
    if (directory == 'Data') {
        if (! dir.exists('Data')) {
            dir.create('Data')
        }
    } else {
        if (! dir.exists(directory)) {
            dir.create(directory, recursive = TRUE)
        }
    }
    
    file.move(files, directory)
}


###
file.move <- function(files, directory) {
    file.copy(from=files, to=directory, recursive = TRUE)
    unlink(files, recursive=TRUE)
}


##############
gdcGetURL <- function(project.id, data.type) {
    urlAPI <- 'https://api.gdc.cancer.gov/files?'
    
    if (data.type=='RNAseq') {
        data.category <- 'Transcriptome Profiling'
        data.type <- 'Gene Expression Quantification'
        workflow.type <- 'HTSeq - Counts'
    } else if (data.type=='miRNAs') {
        data.category <- 'Transcriptome Profiling'
        data.type <- 'Isoform Expression Quantification'
        workflow.type <- 'BCGSC miRNA Profiling'
    } else if (data.type=='Clinical') {
        data.category <- 'Clinical'
        data.type <- 'Clinical Supplement'
        workflow.type <- NA
    } else if (data.type=='pre-miRNAs') {
        data.category <- 'Transcriptome Profiling'
        data.type <- 'miRNA Expression Quantification'
        workflow.type <- 'BCGSC miRNA Profiling'
    }
    
    project <- paste('{"op":"in","content":{"field":"cases.',
        'project.project_id","value":["', 
        project.id, '"]}}', sep='')
    dataCategory <- paste('{"op":"in","content":{"field":"files.', 
        'data_category","value":"', data.category, '"}}', sep='')
    dataType <- paste('{"op":"in","content":{"field":"files.data_type",',
        '"value":"', data.type, '"}}', sep='')
    workflowType <- paste('{"op":"in","content":{"field":"files.',
        'analysis.workflow_type","value":"', workflow.type, '"}}', sep='')
    
    
    if (is.na(workflow.type)) {
        dataFormat <- paste('{"op":"in","content":{"field":"files.',
            'data_format","value":"', 'BCR XML', '"}}', sep='')
        content <- paste(project, dataCategory, dataType, dataFormat, sep=',')
    } else {
        content <- paste(project, dataCategory, dataType, 
            workflowType, sep=',')
    }
    
    filters <- paste('filters=',URLencode(paste('{"op":"and","content":[', 
        content, ']}', sep='')),sep='')
    
    expand <- paste('analysis', 'analysis.input_files', 'associated_entities',
        'cases', 'cases.diagnoses','cases.diagnoses.treatments', 
        'cases.demographic', 'cases.project', 'cases.samples', 
        'cases.samples.portions', 'cases.samples.portions.analytes', 
        'cases.samples.portions.analytes.aliquots',
        'cases.samples.portions.slides', sep=',')
    
    expand <- paste('expand=', expand, sep='')
    
    payload <- paste(filters, 'pretty=true', 'format=JSON', 
        'size=10000', expand, sep='&')
    url <- paste(urlAPI, payload, sep='')
    
    return (url)
}