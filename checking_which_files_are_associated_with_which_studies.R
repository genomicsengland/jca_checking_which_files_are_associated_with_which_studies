library(opencgaR); packageDescription ("opencgaR", fields = "Version") # "1.3.0"
library(getPass); packageDescription ("getPass", fields = "Version")

# Connect to openCGA
con <- initOpencgaR(host="https://opencgainternal.gel.zone/opencga/", version="v1")

# or the traditional way
con <- opencgaLogin(opencga = con, userid = "jambrose", passwd = getPass())

options(stringsAsFactors = F)

studies <- c("1000000038", "1000000034", "1000000036")
#studies <- c("1000000034")

samples_to_be_queried <- read.csv(file = "/Users/johnambrose/Downloads/ready_to_upload.11-12-17.csv", header = TRUE)

file_search <- function(study, pattern){
  files <- NULL
  try(files <- fileClient(OpencgaR = con, action = "search", params = list(study = study, name = pattern, include = "name")), silent = TRUE)
  if(!is.null(files$name)){return(files$name)}
  else{return("NA")}
}

filenames_df <- data.frame(names = character(0), study = character(0))
filenames_all_df <- data.frame(names = character(0), study = character(0))

for(study in studies){
  files_bam <- file_search(study, "~.bam$")
  files_bam_index <- file_search(study, "~.bam.bai$")
  files_vcf <- file_search(study, "~.vcf.gz$")
  files_vcf_index <- file_search(study, "~.vcf.gz.tbi$")
  filenames <- c(files_bam, files_bam_index, files_vcf, files_vcf_index)
  filenames_df <- data.frame(names = filenames, study = rep(study, length(filenames)))
  filenames_all_df <- rbind(filenames_all_df, filenames_df)
}

filenames_all_df <- unique(filenames_all_df)

filenames_subset_df <- data.frame(sampleId = character(0), file = character(0), study = character(0))

for(i in 1:nrow(samples_to_be_queried)){
  if(samples_to_be_queried[i,5] == "GL"){
    sample <- samples_to_be_queried[i,3]
    bam_indices <- grep(paste0("^", sample, "\\.bam"), filenames_all_df$names)
    bam_sample_name <- rep(sample, length(bam_indices))
    vcf_indices <- grep(paste0("^", sample, "\\.vcf.gz"), filenames_all_df$names)
    vcf_sample_name <- rep(sample, length(vcf_indices))
    filenames_subset_df.bam_row <- data.frame(sampleId = bam_sample_name, file = filenames_all_df$names[bam_indices], study = filenames_all_df$study[bam_indices])
    filenames_subset_df.vcf_row <- data.frame(sampleId = vcf_sample_name, file = filenames_all_df$names[vcf_indices], study = filenames_all_df$study[vcf_indices])
    filenames_subset_df <- rbind(filenames_subset_df, filenames_subset_df.bam_row, filenames_subset_df.vcf_row)
  }else if(samples_to_be_queried[i,5] == "FF" | samples_to_be_queried[i,5] == "FFPE"){
    sample <- samples_to_be_queried[i,3]
    bam_indices <- grep(paste0("^", sample, ".bam"), filenames_all_df$names)
    bam_sample_name <- rep(sample, length(bam_indices))
    vcf_indices <- grep(paste0("^LP[0-9]{7,}-DNA_[A-H]{1,}[0-9]{2,}_", sample, "\\.somatic\\.vcf\\.gz"), filenames_all_df$names)
    vcf_sample_name <- rep(sample, length(vcf_indices))
    vcf_sv_indices <- grep(paste0("^LP[0-9]{7,}-DNA_[A-H]{1,}[0-9]{2,}_", sample, "\\.somatic\\.SV\\.vcf\\.gz"), filenames_all_df$names)
    vcf_sv_sample_name <- rep(sample, length(vcf_sv_indices))
    filenames_subset_df.bam_row <- data.frame(sampleId = bam_sample_name, file = filenames_all_df$names[bam_indices], study = filenames_all_df$study[bam_indices])
    filenames_subset_df.vcf_row <- data.frame(sampleId = vcf_sample_name, file = filenames_all_df$names[vcf_indices], study = filenames_all_df$study[vcf_indices])
    filenames_subset_df.vcf_sv_row <- data.frame(sampleId = vcf_sv_sample_name, file = filenames_all_df$names[vcf_sv_indices], study = filenames_all_df$study[vcf_sv_indices])
    filenames_subset_df <- rbind(filenames_subset_df, filenames_subset_df.bam_row, filenames_subset_df.vcf_row, filenames_subset_df.vcf_sv_row)
  }
}

filenames_subset_df <- unique(filenames_subset_df)

#possibly inefficient, but now get cohort info from DDF data in 1000000038, in order to merge with the file info above
study <- "1000000038"
cohorts <- cohortClient(OpencgaR = con, action = "search", params = list(study=study, include="name"))

#set up empty data frame
sampleInfoOverall <- data.frame(sampleId=character(0), sampleType=character(0), LDPCode=character(0), disease=character(0), row.names = NULL)

for(cohortId in cohorts$name){
  #check that the cohortId is valid, hopefully catches cohortId "ALL"
  if(grepl("[0-9]{9,}_[0-9]{5,}", cohortId)){
    sampleInfo <- data.frame(sampleId=character(0), sampleType=character(0), LDPCode=character(0), disease=character(0), row.names = NULL)
    cohortAnnotationSet <- cohortAnnotationsetClient(con, cohort = cohortId, action = "search", params = list(study=study))
    if(unlist(cohortAnnotationSet$annotations[[1]]["value"][(cohortAnnotationSet$annotations[[1]]["name"]=="readyForAnalysis")]) == "TRUE"){
      matchedSamples <- unlist(cohortAnnotationSet$annotations[[1]]["value"][(cohortAnnotationSet$annotations[[1]]["name"]=="matchedSamples")])
      for(sample in matchedSamples){
        sampleAnnotationSet <- sampleAnnotationsetClient(con, sample = sample, action = "search", params = list(study=study))
        sampleType <- sampleAnnotationSet$name
        if(sampleType == "TumourSample"){
          disease <- unlist(sampleAnnotationSet$annotations[[1]]["value"][(sampleAnnotationSet$annotations[[1]]["name"]=="diseaseType")])
        }else if(sampleType == "GermlineSample"){
          disease <- "NA"
        }
        ldpcode <- unlist(sampleAnnotationSet$annotations[[1]]["value"][(sampleAnnotationSet$annotations[[1]]["name"]=="LDPCode")])
        sample.row <- data.frame(sampleId=sample, sampleType=sampleType, LDPCode=ldpcode, disease=disease, row.names = NULL)
        sampleInfo <- rbind(sampleInfo, sample.row)
      }
    }
    sampleInfoOverall <- rbind(sampleInfoOverall, sampleInfo)
  }else{
    print(paste0("malformed cohortId :", cohortId))
  }
}

files_and_ddf_merged <- dplyr::full_join(filenames_subset_df, sampleInfoOverall, by = "sampleId")

dim(files_and_ddf_merged[!is.na(files_and_ddf_merged$LDPCode),])
write.csv(files_and_ddf_merged, file = "/Users/johnambrose/Documents/checking_which_files_are_associated_with_which_studies/files_and_ddf_merged.output.csv", row.names = FALSE, quote = FALSE)

#write.csv(unique(filenames_subset_df), file = "/Users/johnambrose/Documents/checking_which_files_are_associated_with_which_studies/checking_which_study_files_are_associated_with.output.csv", row.names = FALSE, quote = FALSE)
