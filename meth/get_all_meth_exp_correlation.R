#!/usr/bin/Rscript
library(magrittr)
# load meth data ----------------------------------------------------------

methy <- readr::read_rds(file.path("/data/TCGA/TCGA_data/pancan33_meth.rds.gz"))
expr <- readr::read_rds(file.path("/data/TCGA/TCGA_data/pancan33_expr.rds.gz"))
colnames(methy)[2] <-"filter_methy"
colnames(expr)[2] <-"filter_expr"


# expr %>%
#   dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
#   dplyr::select(-expr) -> gene_list_expr
# cnv %>%
#   dplyr::mutate(filter_cnv = purrr::map(cnv, filter_gene_list, gene_list = gene_list)) %>%
#   dplyr::select(-cnv) -> gene_list_cnv
# methy %>%
#   dplyr::mutate(filter_methy = purrr::map(filter_methy, filter_gene_list, gene_list = gene_list)) -> gene_list_methy
# functions ---------------------------------------------------------------

filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}
fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
}
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
}

fn_transform_meth <- function(.d){
  # .d <- te$filter_methy[[1]]
  .d %>% 
    # dplyr::select(-2) %>%
    tidyr::gather(key = barcode, value = value, -symbol,-gene) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type != "11") %>% 
    dplyr::select(-type) %>% 
    dplyr::mutate(sample = fun_barcode(barcode)) %>%
    dplyr::distinct(symbol, sample, .keep_all = T) %>% 
    dplyr::select(-barcode,-gene)
}

# methy %>%
#   head() %>%
#   dplyr::mutate(filter_methy = purrr::map(filter_methy, filter_gene_list, gene_list = gene_list)) -> gene_list_methy

# gene_list_methy %>%
#   dplyr::mutate(filter_methy = purrr::map(.x = filter_methy, .f = fn_transform_meth)) %>%
#   tidyr::unnest() %>%
#   dplyr::filter(! is.na(value)) %>%
#   dplyr::rename(b_val = value) -> methy_df

cl<-15
cluster <- multidplyr::create_cluster(cl)
methy %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type)  %>%
  multidplyr::cluster_assign_value("fn_transform_meth", fn_transform_meth)  %>%
  dplyr::mutate(filter_methy = purrr::map(.x =filter_methy, .f = fn_transform_meth)) %>% 
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>% 
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::filter(! is.na(value)) %>% 
  dplyr::rename(b_val = value) -> methy_df


fn_transform_exp <- function(.d){
  # .d <- te$filter_methy[[1]]
  .d %>% 
    # dplyr::select(-2) %>%
    tidyr::gather(key = barcode, value = value, -symbol,-entrez_id) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type != "11") %>% 
    dplyr::select(-type) %>% 
    dplyr::mutate(sample = fun_barcode(barcode)) %>%
    dplyr::distinct(symbol, sample, .keep_all = T) %>% 
    dplyr::select(-barcode,-entrez_id)
}

# expr %>%
#   head() %>%
#   dplyr::mutate(filter_expr = purrr::map(filter_expr, filter_gene_list, gene_list = gene_list)) -> gene_list_expr
# gene_list_expr %>%
#   dplyr::mutate(filter_expr = purrr::map(.x = filter_expr, .f = fn_transform_exp)) %>%
#   # dplyr::collect() %>%
#   # dplyr::ungroup() %>%
#   tidyr::unnest() %>%
#   # dplyr::select(-PARTITION_ID) %>%
#   dplyr::filter(! is.na(value)) %>%
#   dplyr::rename(expr = value) %>%
#   dplyr::mutate(expr = log2(expr + 1)) ->expr_df

expr %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type)  %>%
  multidplyr::cluster_assign_value("fn_transform_exp", fn_transform_exp)  %>%
  dplyr::mutate(filter_expr = purrr::map(.x = filter_expr, .f = fn_transform_exp)) %>% 
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::select(-PARTITION_ID) %>%
  dplyr::filter(! is.na(value)) %>%
  dplyr::rename(expr = value) %>%
  dplyr::mutate(expr = log2(expr + 1)) ->expr_df

  


# gene_list_cnv %>%
#   dplyr::mutate(filter_cnv = purrr::map(.x = filter_cnv, .f = fn_transform)) %>% 
#   tidyr::unnest() %>% 
#   dplyr::filter(! is.na(value)) %>% 
#   dplyr::rename(cnv = value) -> cnv_df



# expression and methylation ----------------------------------------------
fn_get_cor<-function(data){
  data %>%
    dplyr::group_by(symbol) %>%
    dplyr::do(
      tryCatch(
        broom::tidy(
          cor.test(.$expr,.$b_val,method = c("pearson")))
        )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(spm=estimate) %>%
    dplyr::mutate(fdr= p.adjust(p.value, method = "fdr")) %>%
    dplyr::filter(fdr<=0.05) %>%
    dplyr::mutate(logfdr=-log10(fdr)) %>%
    dplyr::mutate(logfdr=ifelse(fdr>50,50,logfdr)) %>%
    dplyr::select(symbol,spm,logfdr) ->.out
  return(.out)
}

df_expr_methy <-
  methy_df %>% 
  dplyr::inner_join(expr_df, by = c("cancer_types", "symbol", "sample")) %>%
  tidyr::nest(symbol,b_val,sample,expr,.key="data")
df_expr_methy %>%
  dplyr::mutate(spm=purrr::map(data,fn_get_cor)) %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-data)->df_expr_methy_cor

df_expr_methy %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_get_cor", fn_get_cor)  %>%
  dplyr::mutate(spm=purrr::map(data,fn_get_cor)) %>%
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-data) %>%
  dplyr::select(-PARTITION_ID) ->df_expr_methy_cor
parallel::stopCluster(cluster)

# 
# df_expr_methy %>% 
#   dplyr::group_by(cancer_types, symbol) %>% 
#   dplyr::do(broom::glance(lm(expr ~ b_val, data = .))) %>% 
#   dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::select(cancer_types, symbol, ars = adj.r.squared, fdr) %>%
#   dplyr::filter(fdr<=0.05) %>%
#   tidyr::nest(symbol,ars,fdr,.key="cnv_exp") -> df_expr_methy_cor

df_expr_methy_cor %>%
  readr::write_rds("/data/GSCALite/TCGA/meth/pancan34_all_gene_exp-cor-meth.rds.gz",compress="gz")
# expression and cnv ------------------------------------------------------




# expression and cnv and methylation --------------------------------------

# df_expr_cnv_methy <- 
#   df_expr_methy %>% 
#   dplyr::inner_join(cnv_df, by = c("cancer_types", "symbol", "sample"))
# df_expr_cnv_methy %>% 
#   # dplyr::filter(cancer_types == "BRCA", symbol == "ATG5") %>%
#   dplyr::group_by(cancer_types, symbol) %>% 
#   dplyr::do(broom::glance(lm(expr ~ cnv + b_val, data = .))) %>% 
#   dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::select(cancer_types, symbol, ars = adj.r.squared, p.value, fdr) -> df_expr_cnv_methy_cor
# 
# df_expr_cnv_methy_cor %>% 
#   ggplot(aes(x = ars, y = -log10(fdr))) +
#   geom_point()
