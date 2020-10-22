
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
    tidyr::unnest(cor)
}

fn_spm <- function(.data){
  broom::tidy(cor.test(.data$TIL,.data$exp,method="kendall"))
}


# Get the results ---------------------------------------------------------

combine_data %>%
  dplyr::mutate(cor=purrr::pmap(list(ImmuneCellAI,expr,cancer_types),.f=fn_correlation)) %>%
  dplyr::select(-ImmuneCellAI,-expr) -> immune_expr_correlation


# result output -----------------------------------------------------------

immune_expr_correlation %>%
  readr::write_rds(file.path(data_path,"TIL","pan33_ImmuneCellAI_cor_geneExp.rds.gz"))

save.image(file = file.path(git_path,"rda","2.immune_exor_cor.rda"))
