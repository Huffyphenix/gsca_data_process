############### tcga clinical info (survival stage subtype metastasis progression) process ###################

# Library -----------------------------------------------------------------

library(magrittr)

# load data ---------------------------------------------------------------

data_path <- file.path("/home/huff/data/GSCALite/TCGA")
gsca_path <- file.path("/home/huff/data/GSCA")

clinical <- readr::read_rds(file.path(data_path,"clinical","pancan34_clinical.rds.gz"))

# data process ------------------------------------------------------------


# stage -------------------------------------------------------------------

fn_stage_process <- function(cancer_types,.clinical){
  print(cancer_types)
  stage_index <- grep("stage",colnames(.clinical),value = T)
  if(length(stage_index)<1){
    return(tibble::tibble())
  } else{
    .clinical %>%
      dplyr::select(barcode,stage_index) %>%
      dplyr::rename(stage=stage_index) %>%
      dplyr::mutate(stage_type = stage_index) %>%
      dplyr::filter(!is.na(stage))
  }
}

clinical %>%
  dplyr::mutate(stage = purrr::map2(cancer_types,clinical,fn_stage_process)) %>%
  dplyr::select(-clinical) %>%
  dplyr::mutate(n=purrr::map(stage,.f=function(.x){nrow(.x)})) %>%
  tidyr::unnest(n) -> stage_info


# survival ----------------------------------------------------------------
# clinical %>%
#   dplyr::mutate(stage = purrr::map2(cancer_types,clinical,.f=function(.x,.y){
#     .os <- grep("^os",colnames(.y),value=TRUE)
#     .pfs <- grep("^pfs",colnames(.y),value=TRUE)
#     print(paste(.x,.os,.pfs,collapse = ","))
#   })) 
fn_survival_process <- function(cancer_types,.clinical){
  print(cancer_types)
  os_index <- grep("^os",colnames(.clinical),value = T)
  pfs_index <- grep("^pfs",colnames(.clinical),value = T)
  survival_index <- c(os_index,pfs_index)
  if(length(survival_index)<2){
    return(tibble::tibble())
  } else{
    .clinical %>%
      dplyr::select(barcode,survival_index) 
  }
}

clinical %>%
  dplyr::mutate(survival = purrr::map2(cancer_types,clinical,fn_survival_process)) %>%
  dplyr::select(-clinical) %>%
  dplyr::mutate(n=purrr::map(survival,.f=function(.x){nrow(.x)})) %>%
  tidyr::unnest(n) -> survival_info

# subtype ----------------------------------------------------------------
# clinical %>%
#   dplyr::mutate(stage = purrr::map2(cancer_types,clinical,.f=function(.x,.y){
#     .os <- grep("subtype",colnames(.y),ignore.case = TRUE,value=TRUE)
#     # .pfs <- grep("^pfs",colnames(.y),value=TRUE)
#     print(paste(.x,.os,collapse = ","))
#   }))
# fn_subtype_process <- function(cancer_types,.clinical){
#   print(cancer_types)
#   subtype_index <- grep("subtype",colnames(.clinical),value = T,ignore.case = TRUE)
#   if(length(subtype_index)<1){
#     return(tibble::tibble())
#   } else{
#     .clinical %>%
#       dplyr::select(barcode,subtype_index) 
#    }
# }
# 
# clinical %>%
#   dplyr::mutate(survival = purrr::map2(cancer_types,clinical,fn_subtype_process)) %>%
#   dplyr::select(-clinical) %>%
#   dplyr::mutate(n=purrr::map(survival,.f=function(.x){nrow(.x)})) %>%
#   tidyr::unnest(n) -> subtype_info

subtype_info <- readr::read_rds(file.path(data_path,"clinical","pancan34_clinical_subtype.rds.gz"))

subtype_info %>%
  dplyr::mutate(subtype=purrr::map(subtype,.f=function(.x){
    .x %>%
      dplyr::select(barcode,subtype)
  })) -> subtype_info

# combine -----------------------------------------------------------------

subtype_info %>%
  dplyr::full_join(survival_info,by="cancer_types") %>%
  dplyr::full_join(stage_info,by="cancer_types") -> clinical_stage_survival_subtype

clinical_stage_survival_subtype %>%
  readr::write_rds(file.path(gsca_path,"clinical","pancan34_clinical_stage_survival_subtype.rds.gz"),compress = "gz")

save.image(file.path("/home/huff/github/gsca_data_process/clinical/rda/1-clinical-info_stage_survival_subtype.rda"))
