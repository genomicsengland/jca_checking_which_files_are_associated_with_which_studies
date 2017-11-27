# Adapted from Marta Bleda and Dimitris Polychronopoulos
library(opencgaR); packageDescription ("opencgaR", fields = "Version") # "1.3.0"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.1.1.9000"
library(data.table); packageDescription ("data.table", fields = "Version") # "1.10.4"
library(janitor); packageDescription ("janitor", fields = "Version") # "0.3.0"
library(stringr); packageDescription ("stringr", fields = "Version") # "1.2.0"
library(getPass); packageDescription ("getPass", fields = "Version")

con <- initOpencgaR(host="https://opencgainternal.gel.zone/opencga/", version="v1")
con <- opencgaLogin(opencga = con, userid = "jambrose", passwd = getPass())

options(stringsAsFactors = FALSE)

#function to try and grab a file from openCGA using a pattern comprised of sample name and .bam or .bam.bai etc., or return "NA" if the file doesn't exist (i.e. doesn't have a name)
file_search_bam <- function(study, pattern){
  files <- NULL
  #use a uri pattern to ensure you pull out correct vcf, without this I was pulling out vcfs from other directories for some samples
  try(files <- fileClient(OpencgaR = con, action = "search", params = list(study = study, name = pattern, uri="~file:///genomes/by_date/", include = "name,uri,attributes.alignmentHeader.studyName")), silent = TRUE)
  if(!is.null(files$name)){
    files.df <- data.frame(name = files$name, uri = files$uri, study_name = files$attributes$alignmentHeader$studyName)
  }else{
    files.df <- NULL
  }
  return(files.df)
}

#function to try and grab a file from openCGA using a pattern comprised of sample name and .bam or .bam.bai etc., or return "NA" if the file doesn't exist (i.e. doesn't have a name)
file_search <- function(study, pattern){
  files <- NULL
  #use a uri pattern to ensure you pull out correct vcf, without this I was pulling out vcfs from other directories for some samples
  try(files <- fileClient(OpencgaR = con, action = "search", params = list(study = study, name = pattern, uri="~file:///genomes/by_date/", include = "name,uri")), silent = TRUE)
  #it is possible to check study name for GRCh38, for the bam only though
  #,attributes.alignmentHeader.studyName")), silent = TRUE)
  if(!is.null(files$name)){
    files.df <- data.frame(name = files$name, uri = files$uri)
  }else{
    files.df <- NULL
  }
  return(files.df)
}

#function to find file within list of files of one type (bam, bam.bai, etc.)
grep_bam_file <- function(df, pattern){
  if(length(grep(pattern, df$name)) > 0){
    name <- df[grep(pattern, df$name),"name"]
    #possibility to return uri, but currently not implemented
    uri <- df[grep(pattern, df$name),"uri"]
    study_name <- df[grep(pattern, df$name),"study_name"]
    grep_result <- data.frame(name = name, uri = uri, study_name = study_name)
  }else{
    grep_result <- data.frame(name = "NA", uri = "NA", study_name = "NA")
  }
  #this works to return name and uri if needs be
  return(grep_result)
  #return(grep_result$name)
}

#function to find file within list of files of one type (bam, bam.bai, etc.)
grep_file <- function(df, pattern){
  if(length(grep(pattern, df$name)) > 0){
    name <- df[grep(pattern, df$name),"name"]
    #possibility to return uri, but currently not implemented
    uri <- df[grep(pattern, df$name),"uri"]
    grep_result <- data.frame(name = name, uri = uri)
  }else{
    grep_result <- data.frame(name = "NA", uri = "NA")
  }
  #this works to return name and uri if needs be
  #return(grep_result)
  return(grep_result$name)
}

#get all cohorts from study 1000000038 to loop through
study <- 1000000038 #Cancer Program GRCh38 Somatic Genomes
cohorts <- cohortClient(OpencgaR = con, action = "search", params = list(study=study, include="name"))
tumour_bams <- file_search_bam(study, "~.bam$")
tumour_bam_bais <- file_search(study, "~.bam.bai$")
tumour_vcfs <- file_search(study, "~.vcf.gz$")
tumour_vcf_tbis <- file_search(study, "~.vcf.gz.tbi$")
study <- 1000000034
germline_bams <- file_search_bam(study, "~.bam$")
germline_bam_bais <- file_search(study, "~.bam.bai$")
germline_vcfs <- file_search(study, "~.vcf.gz$")
germline_vcf_tbis <- file_search(study, "~.vcf.gz.tbi$")

#create empty dataframe
cohort_all <- data.frame(cohort_id = character(0),
                         pair_name = character(0),
                         row.names = NULL)

#loop through all cohorts and check if they are ready for analysis and then get the files for the relevant samples
for(cohortId in cohorts$name){
  study <- "1000000038"
  #this grep ensures a correct cohort_id, which was needed to filter out "ALL"
  if(grepl("[0-9]{9,}_[0-9]{5,}", cohortId)){
    sampleInfo <- data.frame(sampleId = character(0), sampleType = character(0), LDPCode = character(0), disease = character(0), row.names = NULL)
    cohortAnnotationSet <- cohortAnnotationsetClient(con, cohort = cohortId, action = "search", params = list(study = study))
    if(unlist(cohortAnnotationSet$annotations[[1]]["value"][(cohortAnnotationSet$annotations[[1]]["name"]=="readyForAnalysis")]) == "TRUE"){
      matchedSamples <- cohortAnnotationSet$annotations[[1]]["value"][(cohortAnnotationSet$annotations[[1]]["name"]=="matchedSamples")]
      #go through each row of the matched samples to capture FF+GL and FFPE+GL, etc.
      for(i in 1:nrow(matchedSamples[[1]]["tumourSampleId"])){
        tumour_sample <- matchedSamples[[1]]["tumourSampleId"][i,]
        germline_sample <- matchedSamples[[1]]["germlineSampleId"][i,]
        pair_name <- paste0(tumour_sample, "_", germline_sample)
        
        #check there is something in the cohort_all dataframe before checking for multiple cohorts
        if(length(cohort_all$pair_name) > 0){
          #this is to check for two cohorts, something like "212000049_10000" and "212000049_10001"
          if(length(which(cohort_all$pair_name %in% pair_name)) > 0){
            #gets the bigger number, which should be the latest cohort
            if(cohortId > cohort_all$cohort_id[which(cohort_all$pair_name %in% pair_name)]){
              #remove the row from the overall dataframe before adding the new one
              cohort_all <- cohort_all[-which(cohort_all$pair_name %in% pair_name),]
              cohort_row <- data.frame(cohort_id = cohortId,
                                       pair_name = pair_name,
                                       row.names = NULL)
            }
          }else{
            cohort_row <- data.frame(cohort_id = cohortId,
                                     pair_name = pair_name,
                                     row.names = NULL)
          }
        }else{
          cohort_row <- data.frame(cohort_id = cohortId,
                                   pair_name = pair_name,
                                   row.names = NULL)
        }
        cohort_all <- rbind(cohort_all, cohort_row)
      }
    }else if(unlist(cohortAnnotationSet$annotations[[1]]["value"][(cohortAnnotationSet$annotations[[1]]["name"]=="readyForAnalysis")]) == "FALSE"){
        #not sure if we would ever get here - so stop and investigate!?!
        stop(cohortId)
    }
  }else{
      print(paste0("malformed cohortId :", cohortId))
  }
}

#get samples provided by Alona (see https://my.huddle.net/workspace/38447490/files/#/folder/44519107/list )
ready_to_upload <- read.csv(file = "/Users/johnambrose/Downloads/ready_to_upload.11-17-17.csv", header = TRUE)
ready_to_upload <- dplyr::select(ready_to_upload, PATIENT_ID, SAMPLE_WELL_ID, SAMPLE_TYPE)

#get those samples in ready_to_plot - need ready_to_upload for the patient ids
ready_to_plot <- read.csv(file = "/Users/johnambrose/Downloads/ready_to_plot.11-17-17.csv", header = TRUE)
ready_to_plot <- dplyr::select(ready_to_plot, WELL_ID)

ready_to_plot_with_patient_id <- dplyr::left_join(ready_to_plot, ready_to_upload, by = c("WELL_ID" = "SAMPLE_WELL_ID"))
ready_to_plot_with_patient_id <- ready_to_plot_with_patient_id[c("PATIENT_ID", "WELL_ID", "SAMPLE_TYPE")]

#separate out into FF/FFPE and GL to merge on the relevant column  in the cohort info
ready_to_plot_with_patient_id_tumour <- dplyr::filter(ready_to_plot_with_patient_id, SAMPLE_TYPE == "FF" | SAMPLE_TYPE == "FFPE")
ready_to_plot_with_patient_id_germline <- dplyr::filter(ready_to_plot_with_patient_id, SAMPLE_TYPE == "GL")

ready_to_plot_with_patient_id_tumour_and_germline <- dplyr::left_join(ready_to_plot_with_patient_id_tumour, ready_to_plot_with_patient_id_germline, by = "PATIENT_ID")
colnames(ready_to_plot_with_patient_id_tumour_and_germline) <- c("patient_id", "tumour_sample", "tumour_sample_type", "germline_sample", "germline_sample_type")

patient_all <- data.frame(pair_name = character(0),
                          tumour_sample = character(0),
                          germline_sample = character(0),
                          tumour_bam = character(0),
                          tumour_bam_bai = character(0),
                          tumour_som_vcf = character(0),
                          tumour_som_vcf_tbi = character(0),
                          tumour_sv_som_vcf = character(0),
                          tumour_sv_som_vcf_tbi = character(0),
                          germline_bam = character(0),
                          germline_bam_bai = character(0),
                          germline_vcf = character(0),
                          germline_vcf_tbi = character(0),
                          row.names = NULL)

for(i in 1:nrow(ready_to_plot_with_patient_id_tumour_and_germline)) {
  #print(paste0("Row is: ", i, " out of ", nrow(ready_to_plot_with_patient_id_tumour_and_germline), " rows"))
  tumour_sample <- ready_to_plot_with_patient_id_tumour_and_germline[i, c("tumour_sample")]
  germline_sample <- ready_to_plot_with_patient_id_tumour_and_germline[i, c("germline_sample")]
  pair_name <- paste0(tumour_sample, "_", germline_sample)
  
  tumour_bam <- grep_bam_file(tumour_bams, paste0("^", tumour_sample, "\\.bam$"))
  tumour_bam_bai <- grep_file(tumour_bam_bais, paste0("^", tumour_sample, "\\.bam.bai$"))
  tumour_som_vcf <- grep_file(tumour_vcfs, paste0("^", germline_sample, "_", tumour_sample, "\\.somatic.vcf.gz$"))
  tumour_som_vcf_tbi <- grep_file(tumour_vcf_tbis, paste0("^", germline_sample, "_", tumour_sample, "\\.somatic.vcf.gz.tbi$"))
  tumour_sv_som_vcf <- grep_file(tumour_vcfs, paste0("^", germline_sample, "_", tumour_sample, "\\.somatic.SV.vcf.gz$"))
  tumour_sv_som_vcf_tbi <- grep_file(tumour_vcf_tbis, paste0("^", germline_sample, "_", tumour_sample, "\\.somatic.SV.vcf.gz.tbi$"))
  
  germline_bam <- grep_bam_file(germline_bams, paste0("^", germline_sample, "\\.bam$"))
  germline_bam_bai <- grep_file(germline_bam_bais, paste0("^", germline_sample, "\\.bam.bai$"))
  germline_vcf <- grep_file(germline_vcfs, paste0("^", germline_sample, "\\.vcf.gz$"))
  germline_vcf_tbi <- grep_file(germline_vcf_tbis, paste0("^", germline_sample, "\\.vcf.gz.tbi$"))
  
  patient_row <- data.frame(pair_name = pair_name,
                            tumour_sample = tumour_sample,
                            germline_sample = germline_sample,
                            tumour_bam = tumour_bam,
                            tumour_bam_bai = tumour_bam_bai,
                            tumour_som_vcf = tumour_som_vcf,
                            tumour_som_vcf_tbi = tumour_som_vcf_tbi,
                            tumour_sv_som_vcf = tumour_sv_som_vcf,
                            tumour_sv_som_vcf_tbi = tumour_sv_som_vcf_tbi,
                            germline_bam = germline_bam,
                            germline_bam_bai = germline_bam_bai,
                            germline_vcf = germline_vcf,
                            germline_vcf_tbi = germline_vcf_tbi,
                            row.names = NULL)
  
  patient_all <- rbind(patient_all, patient_row)
}

patient_all_plus_cohort_check <- dplyr::left_join(patient_all, cohort_all, by = "pair_name")

#column order I want
cols <- c("pair_name",
          "cohort_id",
          "tumour_sample",
          "germline_sample",
          "tumour_bam.name",
          "tumour_bam.study_name",
          "tumour_bam_bai",
          "tumour_som_vcf",
          "tumour_som_vcf_tbi",
          "tumour_sv_som_vcf",
          "tumour_sv_som_vcf_tbi",
          "germline_bam.name",
          "germline_bam.study_name",
          "germline_bam_bai",
          "germline_vcf",
          "germline_vcf_tbi")
patient_all_plus_cohort_check <- patient_all_plus_cohort_check[cols]

#get rid of .x in tumour_sample.x and germline_sample.x
#cols_rename <- gsub(pattern = ".x", replacement = "", x = cols)
#colnames(patient_all_plus_cohort_check) <- cols_rename

write.csv(patient_all_plus_cohort_check, file = "patient_all_plus_cohort_check.csv", row.names = FALSE, quote = FALSE)
