#!/usr/bin/Rscript

library(magrittr)
# load snv data  ----------------------------------------------------------
config <-list(database="/data/GSCALite")
snv <- readr::read_rds(file.path(config$database, "TCGA","snv", "pancan33_snv.rds.gz"))

# get all gene list -------------------------------------------------------

# all_gene_list <- snv$snv[[1]]["symbol"][1,]


# funtions ----------------------------------------------------------------


fn_get_percent <- function(cancer_types, filter_snv){
  print(cancer_types)
  n <- length(filter_snv) - 1
  filter_snv %>%
    dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.na(.), 0, .))) -> .d
  .d %>% 
    tidyr::gather(key = barcode, value = count, -symbol) %>% 
    dplyr::mutate(samples = ifelse(count > 0, 1, 0)) %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::summarise(sm_count = sum(count)) %>% 
    dplyr::mutate(per = sm_count / n) -> .d_count
  
  tibble::tibble(cancer_types = cancer_types, n = n, mut_count = list(.d_count))
}


# run ---------------------------------------------------------------------

cl <-34 
cluster <- multidplyr::create_cluster(core=cl)
snv %>% 
#  head() %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_get_percent", fn_get_percent)  %>%
  dplyr::mutate(res = purrr::map2(cancer_types, snv, fn_get_percent)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::select(-cancer_types, -snv) %>% 
  tidyr::unnest(res) -> gene_list_snv_count
parallel::stopCluster(cluster)

# data output  ------------------------------------------------------------

gene_list_snv_count %>% 
  readr::write_rds(path = file.path(config$database, "TCGA","snv",".rds_snv_all_gene_snv_count.rds.gz"), compress = "gz")



