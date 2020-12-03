
# Immune cell correlation with gene expression ----------------------------
# Author: Huffy
# Date: 2020-10-22

library(magrittr)
# path --------------------------------------------------------------------

data_path <- "/home/huff/data/GSCA"
git_path <- "/home/huff/github/gsca_data_process/immune"
# load data ---------------------------------------------------------------

immune_cell_data <- readr::read_rds(file.path(data_path,"TIL","pan33_ImmuneCellAI.rds.gz"))
expr_data <- readr::read_rds(file.path(data_path,"expr","pancan33_expr.IdTrans.rds.gz"))


# data combine ------------------------------------------------------------

immune_cell_data %>%
  dplyr::inner_join(expr_data,by="cancer_types")-> combine_data


# function to get correlation ---------------------------------------------

fn_correlation <- function(.immune,.expr,.cancers){
  message(glue::glue('Handling {.cancers} immune expr correlation'))
  message(class(.expr))
  .expr %>%
    dplyr::filter(!is.na(symbol)) %>%
    tidyr::gather(-entrez_id,-symbol,key="aliquot",value="exp") -> .expr
  .immune %>%
    tidyr::gather(-aliquot,-barcode,-sample_name,key="cell_type",value="TIL") -> .immune
  .combine <- .immune %>%
    dplyr::inner_join(.expr,by="aliquot")
  
  .combine %>%
    tidyr::nest(-symbol,-cell_type,-entrez_id) %>%
    dplyr::mutate(cor = purrr::map(data,.f=fn_spm)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cor) %>%
    dplyr::mutate(cor=estimate)-> tmp
  tmp%>%
    readr::write_rds(file.path(data_path,"TIL/expr_immune",paste(.cancers,"exp_immune_seaprmancor.rds.gz",sep=".")),compress = "gz")
  return(tmp)
}

fn_spm <- function(.data){
  broom::tidy(cor.test(.data$TIL,.data$exp,method="spearman"))
}


# Get the results ---------------------------------------------------------
cancer_done <- list(dir(file.path(data_path,"TIL/expr_immune")) )%>%
  purrr::pmap(.f=function(.x){strsplit(x = .x, split = '\\.')[[1]][1]}) %>% unlist()

cluster <- multidplyr::new_cluster(3)
multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_assign(cluster, fn_correlation=fn_correlation)
multidplyr::cluster_assign(cluster, fn_spm=fn_spm)
multidplyr::cluster_assign(cluster, git_path=git_path)
multidplyr::cluster_assign(cluster, data_path=data_path)

combine_data %>%
  dplyr::filter(!cancer_types %in% cancer_done) %>%
  dplyr::group_by(cancer_types) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(cor=purrr::pmap(list(ImmuneCellAI,expr,cancer_types),.f=fn_correlation)) %>%
  dplyr::collect() %>%
  dplyr::select(-ImmuneCellAI,-expr) -> immune_expr_correlation


# result output -----------------------------------------------------------

immune_expr_correlation %>%
  readr::write_rds(file.path(data_path,"TIL","pan33_ImmuneCellAI_spearmancor_geneExp.rds.gz"))

save.image(file = file.path(git_path,"rda","2.immune_exor_spearmancor.rda"))
parallel::stopCluster(cluster)