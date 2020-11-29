#!/usr/bin/Rscript
library(magrittr)

# load meth data ----------------------------------------------------------
data_path <- "/home/huff/data/GSCA"
methy <- readr::read_rds(file.path(data_path,"methy","pancan33_meth.IdTrans.rds.gz"))

expr <- readr::read_rds(file.path(data_path,"expr","pancan33_expr.IdTrans.rds.gz"))

colnames(methy)[2] <-"filter_methy"
colnames(expr)[2] <-"filter_expr"


# functions ---------------------------------------------------------------

filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}
fun_sample <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
}
fun_barcode16 <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 15
  )
}
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 1)
}

fn_transform_meth <- function(.d){
  # .d <- te$filter_methy[[1]]
  .d %>% 
    # dplyr::filter(symbol %in% c("TP53","CBX2","EZH2","A2M"))%>% 
    dplyr::mutate(entrez=as.numeric(entrez)) %>%
    tidyr::gather(key = barcode, value = value, -symbol,-gene,-entrez) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type != "1") %>% 
    dplyr::select(-type) %>% 
    dplyr::mutate(sample = fun_sample(barcode)) %>%
    dplyr::mutate(barcode16 = fun_barcode16(barcode)) %>%
    dplyr::distinct(symbol, sample, .keep_all = T) %>%
    dplyr::select(-barcode)
}

fn_transform_exp <- function(.d){
  .d %>% 
    # dplyr::filter(symbol %in% c("TP53","CBX2","EZH2","A2M"))%>% 
    dplyr::filter(!is.na(symbol)) %>% 
    dplyr::rename(entrez=entrez_id) %>%
    dplyr::mutate(entrez=as.numeric(entrez)) %>%
    tidyr::gather(key = barcode, value = value, -symbol,-entrez) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type != "11") %>% 
    dplyr::select(-type) %>% 
    dplyr::mutate(sample = fun_sample(barcode)) %>%
    dplyr::mutate(barcode16 = fun_barcode16(barcode)) %>%
    dplyr::distinct(symbol, sample, .keep_all = T) %>%
    dplyr::select(-barcode)
}

# correlation calculation --------
fn_get_cor<-function(data){
  data %>%
    dplyr::do(
      tryCatch(
        broom::tidy(
          cor.test(.$expr,.$b_val,method = c("spearman")))
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(spm=estimate) %>%
    dplyr::select(spm,p.value) ->.out
  fdr <- p.adjust(.out$p.value, method = "fdr")
  .out %>%
    dplyr::mutate(fdr = fdr) %>%
    dplyr::mutate(logfdr=-log10(fdr)) %>%
    dplyr::mutate(logfdr=ifelse(fdr>50,50,logfdr))%>%
    dplyr::mutate(logfdr=ifelse(logfdr=="Inf",50,logfdr))  -> .out
  return(.out)
}

cluster <- multidplyr::new_cluster(5)
multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_assign(cluster, fun_barcode16=fun_barcode16)
multidplyr::cluster_assign(cluster, fun_sample=fun_sample)
multidplyr::cluster_assign(cluster, fun_tn_type=fun_tn_type)
multidplyr::cluster_assign(cluster, fn_transform_meth=fn_transform_meth)

multidplyr::cluster_assign(cluster, fn_transform_exp=fn_transform_exp)
multidplyr::cluster_assign(cluster, fn_get_cor=fn_get_cor)
methy %>% 
  dplyr::group_by(cancer_types) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(filter_methy = purrr::map(.x =filter_methy, .f = fn_transform_meth)) %>% 
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  tidyr::unnest(cols = c(filter_methy)) %>% 
  dplyr::filter(! is.na(value)) %>% 
  dplyr::rename(b_val = value) -> methy_df

expr %>% 
  dplyr::group_by(cancer_types) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(filter_expr = purrr::map(.x = filter_expr, .f = fn_transform_exp)) %>% 
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  tidyr::unnest(cols = c(filter_expr)) %>%
  dplyr::filter(! is.na(value)) %>%
  dplyr::rename(expr = value) %>%
  dplyr::mutate(expr = log2(expr + 1)) ->expr_df

# expression and methylation ----------------------------------------------

df_expr_methy <-
  methy_df %>% 
  dplyr::inner_join(expr_df, by = c("cancer_types", "symbol","entrez","sample","barcode16")) %>%
  tidyr::nest(data = c(barcode16, sample, b_val, expr))

df_expr_methy %>%
  dplyr::group_by(cancer_types,symbol) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(spm=purrr::map(data,fn_get_cor)) %>%
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-data) %>%
  tidyr::unnest(cols = c(spm)) ->df_expr_methy_cor

df_expr_methy_cor %>%
  readr::write_rds(file.path(data_path,"methy","pancan34_exp-cor-meth_GSCAv2.rds.gz"),compress="gz")

parallel::stopCluster(cluster)





