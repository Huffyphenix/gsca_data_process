
######################## cnv survival ############################

library(survival)
library(magrittr)

# path --------------------------------------------------------------------

data_path <- "/home/huff/data/GSCA"
gsca_v2_path <- file.path("/home/huff/data/GSCA")
git_path <- "/home/huff/github/gsca_data_process/cnv"

# load data ---------------------------------------------------------------

sample_info <- readr::read_rds(file.path("/home/huff/data/GSCA/sample_info.rds.gz")) %>%
  dplyr::mutate(barcode16=substr(barcode,1,16)) %>%
  dplyr::select(barcode16,cancer_types) %>%
  unique()

sample_info$cancer_types %>% unique() -> cancer_types

survival <- readr::read_rds(file.path(data_path,"clinical","pancan33_survival_age_stage_NEW.rds.gz")) %>%
  dplyr::mutate(survival=purrr::map(combine,.f=function(.x){
    .x %>%
      dplyr::select(sample_name=barcode,pfs_status, pfs_days,os_days, os_status)
  })) %>%
  dplyr::select(-combine)

cnv <- readr::read_rds(file.path(data_path,"cnv","pancan34_cnv_threshold.IdTrans.rds.gz"))

survival %>%
  dplyr::inner_join(cnv,by="cancer_types") -> combine_cnv_survival
# functions ---------------------------------------------------------------
source("/home/huff/github/gsca_data_process/cnv_new/4.cnv_survival.function.R")

# calculation -------------------------------------------------------------
filelist <- dir("/home/huff/data/GSCA/cnv/cancer_cnv_survival_210813") 
cancer_done <- c()
for (file in filelist) {
  cancer_done <- c(strsplit(file,"_")[[1]][1],cancer_done)
}
cluster <- multidplyr::new_cluster(4)
multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_assign(cluster, fn_survival_res=fn_survival_res)
multidplyr::cluster_assign(cluster, survival_group=survival_group)
multidplyr::cluster_assign(cluster, fn_cox_logp=fn_cox_logp)
multidplyr::cluster_assign(cluster, fn_survival=fn_survival)
multidplyr::cluster_assign(cluster, res_path=res_path)
multidplyr::cluster_assign(cluster, survival=survival)
multidplyr::cluster_assign(cluster, gsca_v2_path=gsca_v2_path)
multidplyr::cluster_assign(cluster, cnv_group=cnv_group)
multidplyr::cluster_assign(cluster, cancer_done=cancer_done)

combine_cnv_survival %>%
  dplyr::filter(!cancer_types %in% cancer_done) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(survival_res = purrr::pmap(list(cancer_types,survival,cnv),.f=fn_survival_res)) %>%  
  dplyr::collect() -> pan33_cnv_survival
#

tibble::tibble(done = list.files(file.path(res_path,"cancer_cnv_survival_210813"))) %>%
  dplyr::group_by(done) %>%
  dplyr::mutate(cancer_types = strsplit(done,"_")[[1]][1]) %>%
  dplyr::ungroup() -> done_cancers

cnv_survival<-tibble::tibble()
for (file in done_cancers$done) {
  done_cancers %>%
    dplyr::filter(done==file) -> tmp
  tmp_res <- readr::read_rds(file.path(res_path,"cancer_cnv_survival_210813",file)) %>%
    dplyr::mutate(cancer_types=tmp$cancer_types) %>%
    tidyr::nest(-cancer_types)
  cnv_survival<- rbind(cnv_survival,tmp_res)
}

cnv_survival %>%
  readr::write_rds(file.path(res_path,"pan33_cnv_survival_NEW210813.rds.gz"))

save.image(file.path(git_path,"3.cnv_survival-GSCAv2.rda"))
parallel::stopCluster(cluster)