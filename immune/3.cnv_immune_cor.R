
# Immune cell correlation with gene cnv ----------------------------
# Author: Huffy
# Date: 2020-10-22

library(magrittr)
# path --------------------------------------------------------------------

data_path <- "/home/huff/data/GSCA"
git_path <- "/home/huff/github/gsca_data_process/immune"
# load data ---------------------------------------------------------------

immune_cell_data <- readr::read_rds(file.path(data_path,"TIL","pan33_ImmuneCellAI.rds.gz"))
cnv_data <- readr::read_rds(file.path(data_path,"cnv","pancan34_cnv.IdTrans.rds.gz"))


# data combine ------------------------------------------------------------

immune_cell_data %>%
  dplyr::inner_join(cnv_data,by="cancer_types")-> combine_data


# function to get correlation ---------------------------------------------

fn_correlation <- function(.immune,.cnv,.cancers){
  message(glue::glue('Handling {.cancers} immune cnv correlation'))
  .cnv %>%
    dplyr::filter(!is.na(symbol)) %>%
    tidyr::gather(-entrez,-symbol,key="aliquot",value="cnv") %>%
    dplyr::mutate(barcode = substr(aliquot, start = 1, stop = 16)) %>%
    dplyr::select(-aliquot) %>%
    dplyr::distinct()-> .cnv
  .immune %>%
    tidyr::gather(-aliquot,-barcode,-sample_name,key="cell_type",value="TIL") %>%
    dplyr::select(-aliquot,-sample_name)-> .immune
  .combine <- .immune %>%
    dplyr::inner_join(.cnv,by="barcode")
  
  .combine %>%
    tidyr::nest(-symbol,-cell_type,-entrez) %>%
    dplyr::mutate(cor = purrr::map(data,.f=fn_spm)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cor)
}

fn_spm <- function(.data){
  broom::tidy(cor.test(.data$TIL,.data$cnv,method="kendall"))
}


# Get the results ---------------------------------------------------------
cluster <- multidplyr::new_cluster(10)

multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_assign(cluster, fn_correlation=fn_correlation)
multidplyr::cluster_assign(cluster, fn_spm=fn_spm)


combine_data %>%
  dplyr::group_by(cancer_types) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(cor=purrr::pmap(list(ImmuneCellAI,cnv,cancer_types),.f=fn_correlation)) %>%
  dplyr::collect() %>%
  dplyr::select(-ImmuneCellAI,-cnv) -> immune_cnv_correlation


# result output -----------------------------------------------------------

immune_cnv_correlation %>%
  readr::write_rds(file.path(data_path,"TIL","pan33_ImmuneCellAI_cor_geneExp.rds.gz"))

save.image(file = file.path(git_path,"rda","3.cnv_immune_cor.rda"))

parallel::stopCluster(cluster)