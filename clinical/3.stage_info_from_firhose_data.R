## get stage data from firehose downloaed data

data_path <- "/home/huff/project/TCGA_survival/data"

cancerlist <- readr::read_tsv(file.path(data_path,"cancer.list.no_combined.txt"),col_names = F)


# load and combine data ---------------------------------------------------
stage_data <- tibble::tibble()
for (cancer in cancerlist$X1) {
  tmp <- readr::read_tsv(file.path(data_path,cancer,paste("gdac.broadinstitute.org_",cancer,".Clinical_Pick_Tier1.Level_4.2016012800.0.0",sep=""),"extract_stage_info.txt")) %>%
    tidyr::gather(-bcr_patient_barcode, key="barcode",value="stage") %>%
    tidyr::spread(key="bcr_patient_barcode",value="stage") %>%
    dplyr::mutate(cancer_types=cancer)%>%
    tidyr::nest(-cancer_types,.key="stage")
  stage_data <- rbind(stage_data,tmp)
}

stage_data %>%
  readr::write_rds(file.path(data_path,"Pancan.Merge.clinical-STAGE.rds.gz"))