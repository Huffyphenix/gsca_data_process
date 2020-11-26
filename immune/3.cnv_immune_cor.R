
# Immune cell correlation with gene cnv ----------------------------
# Author: Huffy
# Date: 2020-10-25

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
    # head() %>%
    dplyr::mutate(cor = purrr::map(data,.f=fn_spm)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cor) -> .tmp
  .tmp %>%
    readr::write_rds(file.path(data_path,"TIL",paste(.cancers,"cnv_immune_cor.rds.gz",sep=".")),compress = "gz")
  .tmp
}

fn_spm <- function(.data){
  
  kendall = broom::tidy(
    tryCatch(
      cor.test(.data$TIL,.data$cnv,method="kendall"),
    error = function(e) {
      1
    })
    )
  if(ncol(kendall)>1){
    kendall %>%
    dplyr::mutate(cor=estimate) %>%   
    dplyr::mutate(fdr= p.adjust(p.value, method = "fdr")) %>%
    dplyr::mutate(logfdr=-log10(fdr)) %>%
    dplyr::mutate(logfdr=ifelse(logfdr>50,50,logfdr)) -> kendall
  } else {
    tibble::tibble(estimate=NA, statistic=NA, p.value=NA, method=NA,alternative=NA) -> kendall
  }
  
  pearson = broom::tidy(
    tryCatch(
      cor.test(.data$TIL,.data$cnv,method="pearson"),
      error = function(e) {
        1
      })
  )
  if(ncol(pearson)>1){
    pearson %>%
      dplyr::mutate(cor=estimate) %>%   
      dplyr::mutate(fdr= p.adjust(p.value, method = "fdr")) %>%
      dplyr::mutate(logfdr=-log10(fdr)) %>%
      dplyr::mutate(logfdr=ifelse(logfdr>50,50,logfdr)) %>%
      dplyr::select(-parameter,-conf.low,-conf.high)-> pearson
  } else {
    tibble::tibble(estimate=NA, statistic=NA, p.value=NA, method=NA,alternative=NA) -> pearson
  }
  
  spearman = broom::tidy(
    tryCatch(
      cor.test(.data$TIL,.data$cnv,method="spearman"),
      error = function(e) {
        1
      })
  )
  if(ncol(spearman)>1){
    spearman %>%
      dplyr::mutate(cor=estimate) %>%   
      dplyr::mutate(fdr= p.adjust(p.value, method = "fdr")) %>%
      dplyr::mutate(logfdr=-log10(fdr)) %>%
      dplyr::mutate(logfdr=ifelse(logfdr>50,50,logfdr)) -> spearman
  } else {
    tibble::tibble(estimate=NA, statistic=NA, p.value=NA, method=NA,alternative=NA) -> spearman
  }
  rbind(pearson,kendall) %>%
    rbind(spearman)
}


# Get the results ---------------------------------------------------------
# cluster <- multidplyr::new_cluster(33)
# 
# multidplyr::cluster_library(cluster,"magrittr")
# multidplyr::cluster_assign(cluster, fn_correlation=fn_correlation)
# multidplyr::cluster_assign(cluster, fn_spm=fn_spm)
# multidplyr::cluster_assign(cluster, data_path=data_path)
# 
# 
# combine_data %>%
#   dplyr::group_by(cancer_types) %>%
#   multidplyr::partition(cluster = cluster) %>%
#   dplyr::mutate(cor=purrr::pmap(list(ImmuneCellAI,cnv,cancer_types),.f=fn_correlation)) %>%
#   dplyr::collect() %>%
#   dplyr::select(-ImmuneCellAI,-cnv) -> immune_cnv_correlation
# 
# 
# # result output -----------------------------------------------------------
# 
# immune_cnv_correlation %>%
#   readr::write_rds(file.path(data_path,"TIL","pan33_ImmuneCellAI_cor_geneCNV.rds.gz"))
# 
# save.image(file = file.path(git_path,"rda","3.cnv_immune_cor.rda"))
# 
# parallel::stopCluster(cluster)

# read all cnv immune files -----------------------------------------------

immune_cnv_file_list <- grep("*cnv_immune_cor.rds.gz",dir(file.path(data_path,"TIL")),value = TRUE)

cancer_types_done <- purrr::map(immune_cnv_file_list,.f=function(.x){strsplit(.x,split = "\\.")[[1]][1]}) %>% unlist()

cluster <- multidplyr::new_cluster(1)

multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_assign(cluster, fn_correlation=fn_correlation)
multidplyr::cluster_assign(cluster, fn_spm=fn_spm)
multidplyr::cluster_assign(cluster, data_path=data_path)


combine_data %>%
  dplyr::filter(!cancer_types %in% cancer_types_done) %>%
  dplyr::group_by(cancer_types) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(cor=purrr::pmap(list(ImmuneCellAI,cnv,cancer_types),.f=fn_correlation)) %>%
  dplyr::collect() %>%
  dplyr::select(-ImmuneCellAI,-cnv) -> immune_cnv_correlation_restcancers


# result output -----------------------------------------------------------

# immune_cnv_correlation %>%
#   readr::write_rds(file.path(data_path,"TIL","pan33_ImmuneCellAI_cor_geneCNV.rds.gz"))
# 
# save.image(file = file.path(git_path,"rda","3.cnv_immune_cor.rda"))

parallel::stopCluster(cluster)