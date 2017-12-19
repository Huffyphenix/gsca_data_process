#!/usr/bin/Rscript

library(magrittr)

# data input --------------------------------------------------------------
cnv <- readr::read_rds(file.path("/data/GSCALite", "TCGA","cnv","pancan34_cnv_threshold.rds.gz"))




all_gene_list<- cnv$cnv[[1]]["symbol"]

# function get exclusive cnv -------------------------------------------------------
fn_cnv_exclusive <- function(V1, V2, .data, cancer_types){
  # .data <- filter_cnv
  # V1 <- 'TP53'
  # V2 <- 'EZH2'
  .data %>% 
    dplyr::filter(symbol %in% c(V1, V2)) %>% 
    tidyr::gather(key = barcode, value = gistic, -symbol) %>% 
    tidyr::spread(key = symbol, value = gistic) %>% 
    dplyr::select(-barcode) -> .d
  .g_name <- colnames(.d)
  # colnames(.d) <- c("A", "B")
  name <- paste(c(cancer_types, .g_name), collapse = "_")
  .d %>% 
    dplyr::filter_all(.vars_predicate = dplyr::all_vars(. == 0)) %>% 
    nrow() -> nn
  .d %>%
    dplyr::filter_all(.vars_predicate = dplyr::all_vars(. != 0)) %>% 
    nrow()-> aa
  .d %>% 
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. != 0)) %>% 
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == 0)) -> .d_an
  
  sum(.d_an %>% dplyr::pull(1) != 0) -> an
  sum(.d_an %>% dplyr::pull(2) != 0) -> na
  c(nn = nn, an = an, na = na, aa = aa) %>% 
    cometExactTest::comet_exact_test(mutmatplot = F) -> p_val
  
  tibble::tibble(te = name, nn = nn, an = an ,na = na, aa = aa, p_val = p_val)
}
fn_cnv_mutal_exclusive <- function(cancer_types, cnv,cluster){
  # cancer_types <- te$cancer_types
  # filter_cnv <- te$filter_cnv[[1]]
  cnv %>% 
    dplyr::pull(symbol) %>% 
    combn(m = 2) %>% 
    t() %>% 
    dplyr::as_data_frame() -> .gene_pairs
  
  .gene_pairs %>%
    multidplyr::partition(cluster = cluster) %>%
    multidplyr::cluster_library("magrittr") %>%
    multidplyr::cluster_library("tidyverse") %>%
    multidplyr::cluster_assign_value("fn_cnv_exclusive", fn_cnv_exclusive) %>%
    multidplyr::cluster_assign_value("cnv", cnv) %>%
    multidplyr::cluster_assign_value("cancer_types", cancer_types) %>%
    dplyr::mutate(rs = purrr::map2(V1, V2, .f = fn_cnv_exclusive, .data = cnv, cancer_types = cancer_types)) %>% 
    dplyr::collect() %>%
    dplyr::as_tibble() %>%
    dplyr::ungroup() %>%
    dplyr::select(-PARTITION_ID) %>%
    dplyr::select(rs) %>% 
    tidyr::unnest() %>% 
    tidyr::separate(col = te, into = c('cancer_types', 'g1', 'g2')) -> .gene_pairs_pval
  
  .gene_pairs_pval %>% 
    dplyr::mutate(fdr = p.adjust(p_val, method = 'fdr'))
}

# get exclusive cnv  -----------------------------------------------------
library(parallel)

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(core=cl)

cnv %>%
  purrr::pmap(.f = fn_cnv_mutal_exclusive,cluster=cluster) %>%
  dplyr::bind_rows() -> mutual_exclusive

parallel::stopCluster(cluster)
readr::write_rds(mutual_exclusive,"/data/GSCALite/TCGA/cnv/pancan34_exclusive_cnv.rds.gz",compress = 'gz')
