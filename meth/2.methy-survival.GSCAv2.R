
# survival analysis of gene methylation------------------------------------------

library(survival)
library(magrittr)

# path --------------------------------------------------------------------

data_path <- "/home/huff/data/GSCA"
gsca_path <- file.path("/home/huff/data/GSCA")
git_path <- "/home/huff/github/gsca_data_process/meth"

survival <- readr::read_rds(file.path(data_path,"clinical","pancan33_survival_age_stage_NEW.rds.gz")) %>%
  dplyr::mutate(survival=purrr::map(combine,.f=function(.x){
    .x %>%
      dplyr::select(sample_name=barcode,pfs_status, pfs_days,os_days, os_status)
  })) %>%
  dplyr::select(-combine)


methy <- readr::read_rds(file.path(data_path,"methy","pancan33_meth.IdTrans.rds.gz"))

# cancer type done ------------------------------------------------------------
res_path <- file.path("/home/huff/data/GSCA","methy","survival_new20210812")

# tibble::tibble(done = list.files(res_path)) %>%
#   dplyr::group_by(done) %>%
#   dplyr::mutate(cancer_types = strsplit(done,"_")[[1]][1]) %>%
#   dplyr::ungroup() -> done_cancers


# functions ---------------------------------------------------------------
source("/home/huff/github/gsca_data_process/meth/2.function-methy-survival.GSCAv2.R")


# calculation -------------------------------------------------------------
cluster <- multidplyr::new_cluster(5)
multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_assign(cluster, fn_survival_res=fn_survival_res)
multidplyr::cluster_assign(cluster, fn_transform_samples=fn_transform_samples)
multidplyr::cluster_assign(cluster, survival_group=survival_group)
multidplyr::cluster_assign(cluster, fn_meth_group=fn_meth_group)
multidplyr::cluster_assign(cluster, fn_cox_logp=fn_cox_logp)
multidplyr::cluster_assign(cluster, fn_survival=fn_survival)
multidplyr::cluster_assign(cluster, res_path=res_path)
multidplyr::cluster_assign(cluster, survival=survival)

methy %>%
  # dplyr::filter(!cancer_types %in% done_cancers$cancer_types ) %>%
  dplyr::group_by(cancer_types) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(survival_res = purrr::pmap(list(cancer_types,methy),.f=fn_survival_res,.survival=survival)) %>%
  dplyr::collect() %>%
  dplyr::select(-data)-> pan33_methy_survival

tibble::tibble(done = list.files(res_path)) %>%
  dplyr::group_by(done) %>%
  dplyr::mutate(cancer_types = strsplit(done,"_")[[1]][1]) %>%
  dplyr::ungroup() -> done_cancers

methy_survival<-tibble::tibble()
for (file in done_cancers$done) {
  done_cancers %>%
    dplyr::filter(done==file) -> tmp
  tmp_res <- readr::read_rds(file.path(res_path,file)) %>%
    dplyr::mutate(cancer_types=tmp$cancer_types) %>%
    tidyr::nest(-cancer_types)
  methy_survival<- rbind(methy_survival,tmp_res)
}

methy_survival %>%
  readr::write_rds(file.path(gsca_path,"methy","pan33_methy_survival_NEW210812.rds.gz"))

save.image(file.path(git_path,"02.methy_survival.rda"))
parallel::stopCluster(cluster)
