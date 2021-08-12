# survival analysis of gene expression------------------------------------------

library(survival)
library(magrittr)

# path --------------------------------------------------------------------

data_path <- "/home/huff/data/GSCA"

survival_os <- readr::read_rds("/home/huff/project/TCGA_survival/data/Pancan.Merge.clinical-OS-Age-stage.rds.gz") %>%
  dplyr::rename("os"="clinical_data")
survival_pfs <- readr::read_rds("/home/huff/project/data/TCGA-survival-time/cell.2018.survival/TCGA_pancan_cancer_cell_survival_time.rds.gz") %>%
  dplyr::rename("cancer_types"="type","pfs"="data")

survival_os %>%
  dplyr::inner_join(survival_pfs,by="cancer_types") %>%
  dplyr::mutate(combine=purrr::map2(os,pfs,.f=function(.x,.y){
    .y %>%
      dplyr::rename("barcode"="bcr_patient_barcode","pfs_status"="PFS","pfs_days"="PFS.time") %>%
      dplyr::full_join(.x,by="barcode") %>%
      dplyr::mutate(os_days=as.numeric(OS), os_status=as.numeric(Status))
  })) %>%
  dplyr::select(cancer_types,combine) -> survival

survival %>% readr::write_rds(file.path(data_path,"clinical","pancan33_survival_age_stage_NEW.rds.gz"))
