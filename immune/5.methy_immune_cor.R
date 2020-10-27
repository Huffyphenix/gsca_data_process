
# Immune cell correlation with gene methy ----------------------------
# Author: Huffy
# Date: 2020-10-25

library(magrittr)
# path --------------------------------------------------------------------

data_path <- "/home/huff/data/GSCA"
git_path <- "/home/huff/github/gsca_data_process/immune"
# load data ---------------------------------------------------------------

immune_cell_data <- readr::read_rds(file.path(data_path,"TIL","pan33_ImmuneCellAI.rds.gz"))
methy_data <- readr::read_rds(file.path(data_path,"methy","pancan33_meth.IdTrans.rds.gz"))


# data combine ------------------------------------------------------------

immune_cell_data %>%
  dplyr::inner_join(methy_data,by="cancer_types")-> combine_data


# function to get correlation ---------------------------------------------
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
}

fn_correlation <- function(.immune,.methy,.cancers){
  message(glue::glue('Handling {.cancers} immune methy correlation'))
  .methy %>%
    dplyr::filter(!is.na(symbol)) %>%
    tidyr::gather(-entrez,-symbol,-gene,key="barcode",value="methy") %>%
    dplyr::filter(substr(barcode,14,15)!= "11") %>%
    dplyr::mutate(sample=substr(barcode,1,12)) %>%
    dplyr::select(-barcode) %>%
    dplyr::distinct() -> .methy
  .immune %>%
    tidyr::gather(-aliquot,-barcode,-sample_name,key="cell_type",value="TIL") %>%
    dplyr::select(-aliquot) %>%
    dplyr::filter(substr(barcode,14,15)!= "11") %>%
    dplyr::rename(sample=sample_name) %>%
    dplyr::select(-barcode) -> .immune
  
  .combine <- .immune %>%
    dplyr::inner_join(.methy,by="sample")
  
  .combine %>% 
    tidyr::nest(-symbol,-cell_type,-entrez,-gene) %>%
    dplyr::mutate(cor = purrr::map(data,.f=fn_spm)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cor) %>%
    dplyr::mutate(cor=estimate) %>%
    dplyr::mutate(fdr= p.adjust(p.value, method = "fdr")) %>%
    dplyr::mutate(logfdr=-log10(fdr)) %>%
    dplyr::mutate(logfdr=ifelse(logfdr>50,50,logfdr))

  }

fn_spm <- function(.data){
  broom::tidy(cor.test(.data$TIL,.data$methy,method="kendall"))
}


# Get the results ---------------------------------------------------------
cluster <- multidplyr::new_cluster(10)

multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_assign(cluster, fn_correlation=fn_correlation)
multidplyr::cluster_assign(cluster, fn_spm=fn_spm)


combine_data %>%
  dplyr::group_by(cancer_types) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(cor=purrr::pmap(list(ImmuneCellAI,methy,cancer_types),.f=fn_correlation)) %>%
  dplyr::collect() %>%
  dplyr::select(-ImmuneCellAI,-methy) -> immune_methy_correlation


# result output -----------------------------------------------------------

immune_methy_correlation %>%
  readr::write_rds(file.path(data_path,"TIL","pan33_ImmuneCellAI_cor_genemethy.rds.gz"))

save.image(file = file.path(git_path,"rda","3.methy_immune_cor.rda"))

parallel::stopCluster(cluster)
