
######################## snv survival ############################

library(survival)
library(magrittr)

# path --------------------------------------------------------------------

data_path <- "/home/huff/data/GSCA"
gsca_v2_path <- file.path("/home/huff/data/GSCA")
git_path <- "/home/huff/github/gsca_data_process/snv"

# load data ---------------------------------------------------------------

sample_info <- readr::read_rds(file.path("/home/huff/data/GSCA/sample_info.rds.gz")) %>%
  dplyr::mutate(barcode16=substr(barcode,1,16)) %>%
  dplyr::select(barcode16,cancer_types) %>%
  unique()

sample_info$cancer_types %>% unique() -> cancer_types

survival <- readr::read_rds(file.path(data_path,"clinical","pancan33_survival_age_stage_NEW.rds.gz")) %>%
  dplyr::mutate(survival=purrr::map(combine,.f=function(.x){
    .x %>%
      dplyr::select(barcode,pfs_status, pfs_days,os_days, os_status)
  })) %>%
  dplyr::select(-combine)

sample_with_snv <- readr::read_rds(file.path(data_path,"pancan33_sample_with_snv.rds.gz")) 

# functions ---------------------------------------------------------------
source("/home/huff/github/gsca_data_process/snv/4.snv_survival.funciton.R")

# calculation -------------------------------------------------------------
cluster <- multidplyr::new_cluster(5)
multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_assign(cluster, fn_survival_res=fn_survival_res)
multidplyr::cluster_assign(cluster, survival_group=survival_group)
multidplyr::cluster_assign(cluster, fn_cox_logp=fn_cox_logp)
multidplyr::cluster_assign(cluster, fn_survival=fn_survival)
multidplyr::cluster_assign(cluster, res_path=res_path)
multidplyr::cluster_assign(cluster, survival=survival)
multidplyr::cluster_assign(cluster, sample_with_snv=sample_with_snv)
multidplyr::cluster_assign(cluster, gsca_v2_path=gsca_v2_path)


tibble::tibble(cancer_types = cancer_types) %>%
  # dplyr::filter(!cancer_types %in% c("ACC","BLCA")) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(survival_res = purrr::map(cancer_types,.f=fn_survival_res,.survival=survival)) %>%
  dplyr::collect()-> pan33_snv_survival


pan33_snv_survival %>%
  readr::write_rds(file.path(gsca_v2_path,"snv","pan33_snv_survival_NEW_210813.rds.gz"),compress = "gz")

save.image(file.path(git_path,"3.snv_survival-GSCAv2.rda"))
parallel::stopCluster(cluster)