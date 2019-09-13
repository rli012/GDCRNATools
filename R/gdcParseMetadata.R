##' @title Parse metadata
##' @description Parse metadata either by providing the \emph{.json} 
##'     file that is downloaded from GDC cart or by parse metadata 
##'     automatically by providing the projct id and data type
##' @param metafile metadata file in \code{.json} format download 
##'     from GDC cart. If provided, the metadata will be parsed from 
##'     this file, otherwise, \code{project} and \code{data.type} arguments 
##'     should be provided to retrieve metadata automatically. 
##'     Default is \code{NULL}
##' @param project.id project id in GDC
##' @param data.type one of \code{'RNAseq'} and \code{'miRNAs'}
##' @param write.meta logical, whether to write the metadata to a 
##'     \code{.json} file
##' @importFrom rjson fromJSON
##' @importFrom jsonlite toJSON
##' @return A dataframe of metadata containing file_name, 
##'     sample_id, etc. as well as some basic clinical data
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### Merge RNA expression data #######
##' metaMatrix <- gdcParseMetadata(project.id='TARGET-RT', data.type='RNAseq')
gdcParseMetadata <- function(metafile=NULL, project.id, 
        data.type, write.meta=FALSE) {
    
    if (! is.null(metafile)) {
        metadata <- rjson::fromJSON(file=metafile)
    } else {
        url <- gdcGetURL(project.id = project.id, data.type = data.type)
        
        metadata <- rjson::fromJSON(file=url)
        metadata <- metadata$data$hits
        #keep <- unlist(lapply(metadata$data$hits, 
        #function(sam) sam$analysis$workflow_type %in% 
        #c('HTSeq - Counts', 'BCGSC miRNA Profiling')))
        #metadata <- metadata$data$hits[keep]
        
        if (write.meta==TRUE) {
            metafile <- jsonlite::toJSON(metadata, pretty=TRUE)
            
            systime <- gsub(' ', 'T', Sys.time())
            systime <- gsub(':', '-', systime)
            
            write(metafile, file=paste(project.id, data.type, 
                'metadata', systime, 'json', sep='.'))
        }
        
    }
    
    nSam <- length(metadata)
    
    file_name <- vapply(seq_len(nSam), 
        function(i) metadata[[i]]$file_name, 
        character(1))
    file_id <- vapply(seq_len(nSam), 
        function(i) metadata[[i]]$file_id, character(1))
    submitter_id <- vapply(seq_len(nSam), function(i) 
        metadata[[i]]$cases[[1]]$samples[[1]]$submitter_id,
        character(1))
    sample <- vapply(submitter_id, 
        function(v) substr(v, 1, nchar(v)-1), 
        character(1))
    entity_submitter_id <- vapply(seq_len(nSam), function(i) 
        metadata[[i]]$associated_entities[[1]]$entity_submitter_id,
        character(1))
    sample_type <- vapply(seq_len(nSam), function(i) 
        metadata[[i]]$cases[[1]]$samples[[1]]$sample_type,
        character(1))
    patient <- vapply(seq_len(nSam), function(i) 
        metadata[[i]]$cases[[1]]$submitter_id,
        character(1))
    
    
    gender <- null2naFun(vapply(seq_len(nSam), function(i) 
        chr2naFun(metadata[[i]]$cases[[1]]$demographic$gender), character(1)))
    #race <- sapply(1:length(metadata), function(i) 
    #metadata[[i]]$cases[[1]]$demographic$race)
    #ethnicity <- sapply(1:length(metadata), function(i) 
    #metadata[[i]]$cases[[1]]$demographic$ethnicity)
    
    project_id <- vapply(seq_len(nSam), function(i) 
        chr2naFun(metadata[[i]]$cases[[1]]$project$project_id), 
        character(1))
    
    tumor_stage <- null2naFun(vapply(seq_len(nSam), function(i) 
        chr2naFun(metadata[[i]]$cases[[1]]$diagnoses[[1]]$tumor_stage), 
        character(1)))
    tumor_grade <- null2naFun(vapply(seq_len(nSam), function(i) 
        chr2naFun(metadata[[i]]$cases[[1]]$diagnoses[[1]]$tumor_grade), 
        character(1)))
    
    
    age_at_diagnosis <- suppressWarnings(vapply(seq_len(nSam), function(i) 
        num2naFun(metadata[[i]]$cases[[1]]$diagnoses[[1]]$age_at_diagnosis), 
        numeric(1)))
    days_to_death <- suppressWarnings(vapply(seq_len(nSam), function(i) 
        num2naFun(metadata[[i]]$cases[[1]]$demographic$days_to_death), 
        numeric(1)))
    days_to_last_follow_up <- suppressWarnings(vapply(seq_len(nSam), 
        function(i) num2naFun(
            metadata[[i]]$cases[[1]]$diagnoses[[1]]$days_to_last_follow_up), 
        numeric(1)))
    
    vital_status <- null2naFun(vapply(seq_len(nSam), function(i) 
        chr2naFun(metadata[[i]]$cases[[1]]$demographic$vital_status), 
        character(1)))
    
    metaMatrix <- data.frame(file_name,file_id,patient,sample,submitter_id,
        entity_submitter_id, sample_type, gender,age_at_diagnosis,tumor_stage,
        tumor_grade,days_to_death, days_to_last_follow_up,vital_status, 
        project_id, stringsAsFactors = FALSE)
    
    
    metaMatrix <- metaMatrix[order(metaMatrix$submitter_id),]
    metaMatrix[metaMatrix=='not reported'] <- NA
    metaMatrix$sample_type <- gsub(' ', '', 
        metaMatrix$sample_type, fixed=TRUE)
    metaMatrix$tumor_stage <- gsub(' ', '', 
        metaMatrix$tumor_stage, fixed=TRUE)
    
    return (metaMatrix)
}



null2naFun <- function(v) {
    v[v=='NA'] <- NA
    return (v)
}



chr2naFun <- function(v) {
    if (is.null(v)) {
        return ('NA')
    } else {
        return (as.character(v))
    }
}

num2naFun <- function(v) {
    if (is.null(v)) {
        return (NA)
    } else if (is.numeric(v)) {
        return (as.numeric(as.character(v)))
    } else {
        return (as.character(v))
    }
}
