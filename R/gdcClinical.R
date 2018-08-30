##' @title Download clinical data in GDC
##' @description Download clinical data in GDC 
##'     either by providing the manifest file or 
##'     specifying the project id and data type
##' @param manifest menifest file that is downloaded from the GDC cart. 
##'     If provided, files whose UUIDs are in the manifest file will be 
##'     downloaded via gdc-client, otherwise, \code{project id} argument 
##'     should be provided to download data automatically. 
##'     Default is \code{NULL}
##' @param project.id project id in GDC
##' @param directory the folder to save downloaded files. 
##'     Default is \code{'Clinical'}
##' @param write.manifest logical, whether to write out the manifest file
##' @param method method that is used to download data. Either 
##'     \code{'GenomicDataCommons'} which is a well established method 
##'     developed in the \pkg{GenomicDataCommons'} package, or alternatively 
##'     \code{'gdc-client'} which uses the \code{gdc-client} tool developed 
##'     by GDC. Default is \code{'gdc-client'}.
##' @return downloaded files in the specified directory
##' @export
##' @author Ruidong Li and Han Qu
##' @examples
##' ####### Download Clinical data by manifest file #######
##' manifest <- 'Clinical.manifest.txt'
##' \dontrun{gdcClinicalDownload(manifest  = manifest,
##'                     directory = 'Clinical')}
##'                    
##' ####### Download Clinical data by project id #######
##' project <- 'TCGA-PRAD'
##' \dontrun{gdcClinicalDownload(project.id     = project, 
##'                     write.manifest = TRUE,
##'                     directory      = 'Clinical')}
gdcClinicalDownload <- function(manifest=NULL, project.id, 
    directory='Clinical', write.manifest=FALSE, method='gdc-client') {

    data.type = 'Clinical'
    
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



##' @title Merge clinical data
##' @description Merge clinical data in \code{.xml} files 
##'     that are downloaded from GDC to a dataframe
##' @param path path to downloaded files for merging
##' @param key.info logical, whether to return the key clinical 
##'     information only. If \code{TRUE}, only clinical information 
##'     such as age, stage, grade, overall survial, etc. will be returned
##' @param organized logical, whether the clinical data have already
##'     been organized into a single folder (eg., data downloaded by the
##'     'GenomicDataCommons' method are already organized). 
##'     Default is \code{FALSE}.
##' @importFrom XML xmlParse
##' @importFrom XML xmlApply
##' @importFrom XML getNodeSet
##' @importFrom XML xmlValue
##' @importFrom XML xmlName
##' @return A dataframe of clinical data with rows are patients 
##'     and columns are clinical traits
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### Merge clinical data #######
##' path <- 'Clinical/'
##' \dontrun{clinicalDa <- gdcClinicalMerge(path=path, key.info=TRUE)}
gdcClinicalMerge <- function(path, key.info=TRUE, organized=FALSE) {

    options(stringsAsFactors = FALSE)
    
    #if (endsWith(path, '/')) {
    #  path = substr(path, 1, nchar(path)-1)
    #}
    
    message ('############### Merging Clinical data ###############\n')
    
    if (organized==TRUE) {
        filenames <- file.path(path, getXMLFile(path), 
            fsep = .Platform$file.sep)
    } else {
        folders <- file.path(path, dir(path), fsep = .Platform$file.sep)
        folders <- folders[dir.exists(folders)]
        filenames <- file.path(folders, unlist(lapply(folders, function(v) 
            getXMLFile(v))), fsep = .Platform$file.sep)
    }
    
    df <- lapply(filenames, function(fl) parseXMLFun(fl))
    traits <- unique(names(unlist(df)))
    xmlMatrix <- do.call('cbind', lapply(df, function(v) v[traits]))
    
    rownames(xmlMatrix) <- traits
    colnames(xmlMatrix) <- xmlMatrix['bcr_patient_barcode',]
    xmlMatrix[xmlMatrix==""]<- "NA"
    xmlMatrix<- data.frame(xmlMatrix)
    if (key.info== TRUE){
        line1<- xmlMatrix[c("age_at_initial_pathologic_diagnosis",
            "ethnicity", "gender", "race","clinical_stage", 
            "clinical_T","clinical_N", "clinical_M", "gleason_grading",
            "gleason_score", "primary_pattern","secondary_pattern",
            "tertiary_pattern","psa","psa_value","days_to_psa"),]
        
        
        ### days_to_death
        
        num3<- grep("^days_to_death", rownames(xmlMatrix))
        t3<- xmlMatrix[num3,]
        t3[t3=='NA']<-NA
        t3[is.na(t3)]<-"0"
        
        line2<- data.frame(t(apply(t3, 2, function(v) max(as.numeric(v)))))
        line2[line2<=0]<- NA
        rownames(line2)<- "days_to_death"
        
        
        ### days_to_last_followup
        
        num4<- grep("^days_to_last_followup", rownames(xmlMatrix))
        t4<- xmlMatrix[num4,]
        
        t4[t4=='NA']<-NA
        t4[is.na(t4)]<-"0"
        
        line3<- data.frame(t(apply(t4, 2, function(v) max(as.numeric(v)))))
        line3[line3<=0]<- NA
        rownames(line3)<- "days_to_last_followup"
        
        
        ### vital_status
        num11<- grep("^vital_status", rownames(xmlMatrix))
        t11<- xmlMatrix[num11,]
        
        t11[t11=='NA']<-NA
        t11[is.na(t11)]<-"0"
        
        line5<- data.frame(t(apply(t11, 2, max)))
        line5[line5=='0']<- NA
        rownames(line5)<- "vital_status"
        
        ### age_at_initial_pathologic_diagnosis
        
        line6<- xmlMatrix[c("initial_pathologic_diagnosis_method", 
            "lymphnodes_examined","number_of_lymphnodes_examined",
            "number_of_lymphnodes_positive_by_he","pathologic_categories", 
            "pathologic_stage","pathologic_T", "pathologic_M","pathologic_N",
            "new_tumor_event"),]
        
        ### days_to_new_tumor_event_after_initial_treatment
        
        
        if (xmlMatrix['disease_code',1]=='LAML') {
            #line7 <- data.frame(t(rep(NA, ncol(xmlMatrix))))
            #rownames(line7)<- "days_to_new_tumor_event_after_initial_treatment"
            
            cleantable<- rbind(line1, line2, line3, line5, line6)
            
            
        } else {
            num5<- grep("^days_to_new_tumor_event_after_initial_treatment", 
                        rownames(xmlMatrix))
            t5<- xmlMatrix[num5,]
            
            t5[t5=='NA']<-NA
            t5[is.na(t5)]<-"999999"
            
            line7<- data.frame(t(apply(t5, 2, function(v) min(as.numeric(v)))))
            line7[line7==999999 | line7<=0]<- NA
            
            rownames(line7)<- "days_to_new_tumor_event_after_initial_treatment"
            
            
            ### new_neoplasm_event_type
            
            num6<- grep("^new_neoplasm_event_type", rownames(xmlMatrix))
            t6<- xmlMatrix[num6,]
            # t6[is.na(t6)]<-"0"
            line8<- NULL
            
            for (i in seq_len(ncol(t6))) {
                t6.1<- t6[which(t6[,i] != "NA"),i]
                t6.1<- paste(t6.1,collapse=",")
                line8<- append(line8,t6.1)
            }
            line8<-data.frame(t(line8))
            line8[line8==""]<- "NA"
            rownames(line8)<-"new_neoplasm_event_type"
            colnames(line8) <- colnames(t6)
            
            
            ### new_tumor_event_after_initial_treatment
            
            num7<- grep("^new_tumor_event_after_initial_treatment", 
                        rownames(xmlMatrix))
            t7<- xmlMatrix[num7,]
            t7[t7=='NA']<-NA
            
            t7[is.na(t7)]<-"0"
            
            t7.1<- data.frame(t((colSums(t7=="YES"))))
            rownames(t7.1)<- "new_tumor_event_after_initial_treatment_yes"
            
            t7.2<- data.frame(t((colSums(t7=="NO"))))
            rownames(t7.2)<- "new_tumor_event_after_initial_treatment_no"
            
            line9<- rbind(t7.1,t7.2)
            
            
            if (xmlMatrix['disease_code',1] %in% c('DLBC','PCPG','TGCT')) {
                cleantable<- rbind(line1, line2, line3, line5, line6, line7, line8, line9)
            } else {
                #### additional_pharmaceutical_therapy
                num<- grep("^additional_pharmaceutical_therapy", rownames(xmlMatrix))
                t1<- xmlMatrix[num,]
                
                t1[t1=='NA']<-NA
                t1[is.na(t1)]<-"0"
                
                t1.1<- data.frame(t((colSums(t1=="YES"))))
                rownames(t1.1)<- "additional_pharmaceutical_therapy_yes"
                
                t1.2<- data.frame(t((colSums(t1=="NO"))))
                rownames(t1.2)<- "additional_pharmaceutical_therapy_no"
                
                line10<- rbind(t1.1,t1.2)
                
                
                ### additional_radiation_therapy
                
                num2<- grep("^additional_radiation_therapy", rownames(xmlMatrix))
                t2<- xmlMatrix[num2,]
                t2[t2=='NA']<-NA
                
                t2[is.na(t2)]<-"0"
                
                t2.1<- data.frame(t((colSums(t2=="YES"))))
                rownames(t2.1)<- "additional_radiation_therapy_yes"
                
                t2.2<- data.frame(t((colSums(t2=="NO"))))
                rownames(t2.2)<- "additional_radiation_therapy_no"
                
                line11<- rbind(t2.1,t2.2)
                
                ########## rbind #########
                
                cleantable<- rbind(line1, line2, line3, line5, line6, line7, 
                                   line8, line9, line10, line11)
                #names(cleantable)<- names(line3)
                #colnames(line7)<- names(line3)
                #cleantable<- rbind(cleantable, line3, line7)
                
            }
            
        }
        
        cleantable <- data.frame(t(cleantable), stringsAsFactors = FALSE)
        
        rownames(cleantable) <- gsub('.', '-', rownames(cleantable), 
            fixed=TRUE)
        
        filter <- grep('^NA',colnames(cleantable))
        
        if (length(filter) > 0) {
            cleantable <- cleantable[,-filter]
        }
        
        return(cleantable)
        
    } else{
        return(xmlMatrix)
    }
}


####
parseXMLFun <- function(fl) {
    doc<-xmlParse(file = fl)
    test1<- xmlApply(getNodeSet(doc, "//*"), xmlValue)
    test2<- xmlApply(getNodeSet(doc, "//*"), xmlName)
    names(test1) <- make.names(unlist(test2), unique = TRUE)
    test1 <- unlist(test1)[-c(1:2)]
    return (test1)
}



####
getXMLFile <- function(folder) {
    files <- dir(folder)
    xmlFile <- files[endsWith(files, '.xml')]
    return (xmlFile)
}


