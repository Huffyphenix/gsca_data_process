

library(magrittr)
# load meth data ----------------------------------------------------------
cnv <- readr::read_rds(file.path("/data/TCGA/TCGA_data/pancan34_cnv.rds.gz")) 
expr <- readr::read_rds(file.path("/data/TCGA/TCGA_data/pancan33_expr.rds.gz"))
colnames(cnv)[2] <-"filter_cnv"
colnames(expr)[2] <-"filter_expr"

# expr %>%
#   dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
#   dplyr::select(-expr) -> gene_list_expr
# cnv %>%
#   dplyr::mutate(filter_cnv = purrr::map(cnv, filter_gene_list, gene_list = gene_list)) %>%
#   dplyr::select(-cnv) -> gene_list_cnv



# fuctions ----------------------------------------------------------------

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

fn_transform <- function(.d){
  # .d <- te$filter_methy[[1]]
  .d %>% 
    # dplyr::select(-2) %>%
    tidyr::gather(key = barcode, value = value, -symbol) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type != "11") %>% 
    dplyr::select(-type) %>% 
    dplyr::mutate(sample = fun_barcode(barcode)) %>%
    dplyr::distinct(symbol, sample, .keep_all = T) %>% 
    dplyr::select(-barcode)
}

# data processing ---------------------------------------------------------

cl<-33
cluster <- multidplyr::create_cluster(cl)
cnv %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type)  %>%
  multidplyr::cluster_assign_value("fn_transform", fn_transform)  %>%
  dplyr::mutate(filter_cnv = purrr::map(.x = filter_cnv, .f = fn_transform)) %>% 
  dplyr::collect() %>%
  tidyr::unnest() %>% 
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::filter(! is.na(value)) %>% 
  dplyr::rename(cnv = value) -> cnv_df
# parallel::stopCluster(cluster)

# cl<-33
# cluster <- multidplyr::create_cluster(cl)
expr %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type)  %>%
  multidplyr::cluster_assign_value("fn_transform_exp", fn_transform_exp)  %>%
  dplyr::mutate(filter_expr = purrr::map(.x = filter_expr, .f = fn_transform_exp)) %>% 
  dplyr::collect() %>%
  tidyr::unnest() %>% 
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::filter(! is.na(value)) %>% 
  dplyr::rename(expr = value) %>% 
  dplyr::mutate(expr = log2(expr + 1)) -> expr_df

# calculate ---------------------------------------------------------------
fn_get_cor<-function(data){
  data %>%
    dplyr::group_by(symbol) %>%
    dplyr::do(
      tryCatch(
        broom::tidy(
          cor.test(.$expr,.$cnv,method = c("pearson")))
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(spm=estimate) %>%
    dplyr::mutate(fdr= p.adjust(p.value, method = "fdr")) %>%
    dplyr::filter(fdr<=0.05) %>%
    dplyr::mutate(logfdr=-log10(fdr)) %>%
    dplyr::mutate(logfdr=ifelse(fdr>50,50,fdr)) %>%
    dplyr::select(symbol,spm,logfdr) ->.out
  return(.out)
}

df_expr_cnv <- 
  cnv_df %>% 
  dplyr::inner_join(expr_df, by = c("cancer_types", "symbol", "sample")) %>%
  tidyr::nest(symbol,cnv,sample,expr,.key="data")

df_expr_cnv %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_get_cor", fn_get_cor)  %>%
  dplyr::mutate(spm=purrr::map(data,fn_get_cor)) %>%
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-data) %>%
  dplyr::select(-PARTITION_ID) ->df_expr_cnv_cor
parallel::stopCluster(cluster)


df_expr_cnv_cor %>%
  readr::write_rds("/data/GSCALite/TCGA/cnv/pancan34_all_gene_exp-cor-cnv.rds.gz",compress="gz")

