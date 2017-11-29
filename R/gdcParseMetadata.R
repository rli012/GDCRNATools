##' @title Parse metadata
##' @description Parse metadata saved in the \emph{.json} file that is downloaded from GDC cart 
##'   either by providing the metadata file or by the API method based on the \code{getGDCquery} function in \pkg{TCGAbiolinks} package
##' @param metafile metadata file in \code{.json} format download from GDC cart. If provided, the metadata will be parsed from this file,
##'   otherwise, \code{project} and \code{data.type} arguments should be provided to retrieve metadata via API based on the \code{getGDCquery}
##'   function in \pkg{TCGAbiolinks} package. Default is \code{NULL}
##' @param project.id project id in GDC
##' @param data.type one of \code{'RNAseq'} and \code{'miRNAs'}
##' @param write.meta logical, whether to write the metadata to a \code{.json} file
##' @importFrom rjson fromJSON
##' @importFrom jsonlite toJSON
##' @return A dataframe of metadata containing file_name, sample_id, etc. as well as some basic clinical data
##' @export
##' @author Ruidong Li and Han Qu
gdcParseMetadata <- function(metafile=NULL, project.id, data.type, write.meta=FALSE) {
  
  if (! is.null(metafile)) {
    metadata <- rjson::fromJSON(file=metafile)
  } else {
    url <- gdcGetURL(project.id = project.id, data.type = data.type)
    
    metadata <- rjson::fromJSON(file=url)
    metadata <- metadata$data$hits
    #keep <- unlist(lapply(metadata$data$hits, function(sam) sam$analysis$workflow_type %in% c('HTSeq - Counts', 'BCGSC miRNA Profiling')))
    #metadata <- metadata$data$hits[keep]

    if (write.meta==TRUE) {
      metafile <- jsonlite::toJSON(metadata, pretty=TRUE)
      
      systime <- gsub(' ', 'T', Sys.time())
      systime <- gsub(':', '-', systime)
      
      write(metafile, file=paste(project.id, data.type, 'metadata', systime, 'json', sep='.'))
    }
        
  }
  
  
  file_name <- sapply(1:length(metadata), function(i) metadata[[i]]$file_name)
  file_id <- sapply(1:length(metadata), function(i) metadata[[i]]$file_id)
  submitter_id <- sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$samples[[1]]$submitter_id)
  sample <- sapply(submitter_id, function(v) substr(v, 1, nchar(v)-1))
  entity_submitter_id <- sapply(1:length(metadata), function(i) metadata[[i]]$associated_entities[[1]]$entity_submitter_id)
  sample_type <- sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$samples[[1]]$sample_type)
  patient <- sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$submitter_id)
  
  
  gender <- null2naFun(sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$demographic$gender))
  #race <- sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$demographic$race)
  #ethnicity <- sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$demographic$ethnicity)
  
  project_id <- null2naFun(sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$project$project_id))
  
  tumor_stage <- null2naFun(sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$diagnoses[[1]]$tumor_stage))
  tumor_grade <- null2naFun(sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$diagnoses[[1]]$tumor_grade))
  
  
  age_at_diagnosis <- suppressWarnings(as.numeric(as.character(sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$diagnoses[[1]]$age_at_diagnosis))))
  days_to_death <- suppressWarnings(as.numeric(as.character(sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$diagnoses[[1]]$days_to_death))))
  days_to_last_follow_up <- suppressWarnings(as.numeric(as.character(sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$diagnoses[[1]]$days_to_last_follow_up))))
  vital_status <- null2naFun(sapply(1:length(metadata), function(i) metadata[[i]]$cases[[1]]$diagnoses[[1]]$vital_status))
  
  metaMatrix <- data.frame(file_name,file_id,patient,sample,submitter_id,entity_submitter_id,sample_type,
                           gender,age_at_diagnosis,tumor_stage,tumor_grade,days_to_death,days_to_last_follow_up,vital_status,
                           project_id, stringsAsFactors = F)
  
  
  metaMatrix <- metaMatrix[order(metaMatrix$submitter_id),]
  metaMatrix[metaMatrix=='not reported'] <- NA
  metaMatrix$sample_type <- gsub(' ', '', metaMatrix$sample_type, fixed=T)
  metaMatrix$tumor_stage <- gsub(' ', '', metaMatrix$tumor_stage, fixed=T)
  
  return (metaMatrix)
}



null2naFun <- function(v) {
  v[sapply(v, is.null)] <- NA
  return (unlist(v))
}


