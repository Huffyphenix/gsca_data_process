# Immune cell association with gene snv ----------------------------
# Method: wilcoxon test
# Author: Huffy
# Date: 2020-10-25

library(magrittr)
# path --------------------------------------------------------------------

gsca_v2_path <- "/home/huff/data/GSCA"
git_path <- "/home/huff/github/gsca_data_process/immune"

# load data ---------------------------------------------------------------

immune_cell_data <- readr::read_rds(file.path(gsca_v2_path,"TIL","pan33_ImmuneCellAI.rds.gz"))
# snv_data <- readr::read_rds(file.path(data_path,"snv","all_maf_data.IdTrans.maf.tsv.rds.gz"))

source(file.path(git_path,"4.snv-immune_functions.R"))

# snv data process ------------------------------------------------------------

# calculation -------------------------------------------------------------
immune_cnv_file_list <- grep("*snv_immune_wilcox.rds.gz",dir(file.path(gsca_v2_path,"TIL/snv_immune")),value = TRUE)

cancer_types_done <- purrr::map(immune_cnv_file_list,.f=function(.x){strsplit(.x,split = "\\.")[[1]][1]}) %>% unlist()

cluster <- multidplyr::new_cluster(7)
multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_assign(cluster, fn_immune_snv=fn_immune_snv)
multidplyr::cluster_assign(cluster, fn_res=fn_res)
multidplyr::cluster_assign(cluster, fn_wilcoxn=fn_wilcoxn)
multidplyr::cluster_assign(cluster, fn_immune_snv=fn_immune_snv)
multidplyr::cluster_assign(cluster, gsca_v2_path=gsca_v2_path)


immune_cell_data %>%
  dplyr::filter(!cancer_types %in% cancer_types_done) %>%
  dplyr::group_by(cancer_types) %>%
  # dplyr::filter(!cancer_types %in% c("ACC","BLCA")) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(wilcox_res = purrr::map2(cancer_types,ImmuneCellAI,.f=fn_immune_snv)) %>%
  dplyr::collect() %>%
  dplyr::select(-ImmuneCellAI)-> pan33_snv_immune


pan33_snv_immune %>%
  readr::write_rds(file.path(gsca_v2_path,"TIL","pan33_ImmuneCellAI_cor_geneSNV.rds.gz"),compress = "gz")

save.image(file.path(git_path,"4.snv-immune_association.rda"))
parallel::stopCluster(cluster)
# save image --------------------------------------------------------------

