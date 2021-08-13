
# survival analysis of gene expression------------------------------------------

library(survival)
library(magrittr)

# path --------------------------------------------------------------------

data_path <- "/home/huff/data/GSCA"
gsca_path <- file.path("/home/huff/data/GSCA")
git_path <- "/home/huff/github/gsca_data_process/expr"

survival_os <- readr::read_rds("/home/huff/project/TCGA_survival/data/Pancan.Merge.clinical-OS-Age-stage.rds.gz") %>%
  dplyr::rename("os"="clinical_data")
survival_pfs <- readr::read_rds("/home/huff/project/data/TCGA-survival-time/cell.2018.survival/TCGA_pancan_cancer_cell_survival_time.rds.gz") %>%
  dplyr::rename("cancer_types"="type","pfs"="data")

survival_os %>%
  dplyr::inner_join(survival_pfs,by="cancer_types") %>%
  dplyr::mutate(combine=purrr::map2(os,pfs,.f=function(.x,.y){
    .y %>%
      dplyr::rename("barcode"="bcr_patient_barcode","pfs_status"="PFS","pfs_days"="PFS.time") %>%
      dplyr::select(barcode,pfs_status,pfs_days) %>%
      dplyr::full_join(.x,by="barcode") %>%
      dplyr::mutate(os_days=as.numeric(OS), os_status=as.numeric(Status)) %>%
      dplyr::select(sample_name = barcode,os_days, os_status,pfs_status,pfs_days)
  })) %>%
  dplyr::select(cancer_types,combine) -> survival

expr <- readr::read_rds(file.path(data_path,"expr","pancan33_expr.IdTrans.rds.gz"))


gene_symbol <- readr::read_rds('/home/huff/data/GSCA/id/NCBI_id_in_TCGA-final.rds.gz')

# cancer type done ------------------------------------------------------------
res_path <- file.path("/home/huff/data/GSCA","expr","survival_new20210812")

tibble::tibble(done = list.files(res_path)) %>%
  dplyr::group_by(done) %>%
  dplyr::mutate(cancer_types = strsplit(done,"_")[[1]][1]) %>%
  dplyr::ungroup() -> done_cancers


# functions ---------------------------------------------------------------
source("/home/huff/github/gsca_data_process/expr/02.functions_survival_NEW_20210812.R")


# calculation -------------------------------------------------------------
cluster <- multidplyr::new_cluster(5)
multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_assign(cluster, fn_survival_res=fn_survival_res)
multidplyr::cluster_assign(cluster, fn_transform_samples=fn_transform_samples)
multidplyr::cluster_assign(cluster, survival_group=survival_group)
multidplyr::cluster_assign(cluster, fn_expr_group=fn_expr_group)
multidplyr::cluster_assign(cluster, fn_cox_logp=fn_cox_logp)
multidplyr::cluster_assign(cluster, fn_survival=fn_survival)
multidplyr::cluster_assign(cluster, res_path=res_path)
multidplyr::cluster_assign(cluster, gene_symbol=gene_symbol)
multidplyr::cluster_assign(cluster, survival=survival)

expr %>%
  dplyr::filter(!cancer_types %in% done_cancers$cancer_types ) %>%
  dplyr::group_by(cancer_types) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(survival_res = purrr::pmap(list(cancer_types,expr),.f=fn_survival_res,.survival=survival)) %>%
  dplyr::collect() %>%
  dplyr::select(-data)-> pan33_expr_survival

tibble::tibble(done = list.files(res_path)) %>%
  dplyr::group_by(done) %>%
  dplyr::mutate(cancer_types = strsplit(done,"_")[[1]][1]) %>%
  dplyr::ungroup() -> done_cancers

expr_survival<-tibble::tibble()
for (file in done_cancers$done) {
  done_cancers %>%
    dplyr::filter(done==file) -> tmp
  tmp_res <- readr::read_rds(file.path(res_path,file)) %>%
    dplyr::mutate(cancer_types=tmp$cancer_types) %>%
    tidyr::nest(-cancer_types)
  expr_survival<- rbind(expr_survival,tmp_res)
}

expr_survival %>%
  readr::write_rds(file.path(gsca_path,"expr","pan33_expr_survival_NEW210813.rds.gz"))

save.image(file.path(git_path,"02.gene_expr_survival.rda"))
parallel::stopCluster(cluster)
