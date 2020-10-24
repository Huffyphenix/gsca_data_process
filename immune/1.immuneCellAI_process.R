library(magrittr)
# combine TIL from all cancers --------------------------------------------

immune_path <- file.path("/home/huff/data/GSCA/TIL")
file_list <- dir(immune_path)

for (file in file_list) {
  load(file.path(immune_path,file))
  cancer_name <- strsplit(file,"_")[[1]][1]
  if(cancer_name=="ACC"){
    predict_result %>%
      as.data.frame() %>%
      dplyr::as.tbl() %>%
      dplyr::mutate(cancer_types = strsplit(file,"_")[[1]][1],
                    aliquot = rownames(predict_result)) -> combine
  }else{
    predict_result %>%
      as.data.frame() %>%
      dplyr::as.tbl() %>%
      dplyr::mutate(cancer_types = strsplit(file,"_")[[1]][1],
                    aliquot = rownames(predict_result)) -> .tmp
    combine %>%
      rbind(.tmp) -> combine
  }
}

combine %>%
  dplyr::mutate(aliquot= gsub("\\.","-",aliquot)) %>%
  dplyr::mutate(barcode= substr(aliquot,1,16),
                sample_name= substr(aliquot,1,12)) %>%
  dplyr::select(cancer_types,aliquot,barcode,sample_name,everything()) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::rename("ImmuneCellAI"="data") -> pan33_ImmuneCellAI

pan33_ImmuneCellAI %>%
  readr::write_rds(file.path(immune_path,"pan33_ImmuneCellAI.rds.gz"))
