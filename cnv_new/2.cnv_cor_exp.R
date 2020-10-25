
# correlation between cnv and gene expression -----------------------------

# add some gene in present results and delete some genes from previous results #


# library -----------------------------------------------------------------

library(magrittr)

library(parallel)

# config path -------------------------------------------------------------

git_path <- "/home/huff/github/gsca_data_process"
gsca_v2_path <- "/home/huff/data/GSCA"


# load data ---------------------------------------------------------------

expr <- readr::read_rds(file.path(gsca_v2_path,"expr","pancan33_expr.IdTrans.rds.gz"))
cnv <- readr::read_rds(file.path(gsca_v2_path,"cnv","pancan34_cnv.IdTrans.rds.gz")) 

# add genes to do correlation ---------------------------------------------

# load(file = file.path(git_path,"cnv_new/rda/1.cnv_symbol_process.rda"))

# functions ----------------------------------------------------------------

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

fn_transform_exp <- function(.d,cancer){
  
  message(glue::glue('exp tranform: {cancer}'))
  # .d <- te$filter_methy[[1]]
  .d %>% 
    # dplyr::select(-2) %>%
    dplyr::rename("entrez"="entrez_id") %>%
    dplyr::filter(!is.na(symbol)) %>%
    dplyr::mutate(entrez=as.character(entrez)) %>%
    tidyr::gather(key = barcode, value = value, -symbol,-entrez) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type != "11") %>% 
    dplyr::select(-type) %>% 
    dplyr::mutate(sample = fun_barcode(barcode)) %>%
    dplyr::distinct(symbol, sample, .keep_all = T) %>% 
    dplyr::select(-barcode)
}

fn_transform_cnv <- function(.d,cancer){
  
  message(glue::glue('cnv tranform: {cancer}'))
  .d %>% 
    # dplyr::mutate(entrez=as.character(entrez)) %>%
    dplyr::filter(!is.na(symbol)) %>%
    dplyr::mutate(entrez=as.character(entrez)) %>%
    tidyr::gather(key = barcode, value = value, -symbol,-entrez) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type != "11") %>% 
    dplyr::select(-type) %>% 
    dplyr::mutate(sample = fun_barcode(barcode)) %>%
    dplyr::distinct(symbol, sample, .keep_all = T) %>% 
    dplyr::select(-barcode)
}

# data processing ---------------------------------------------------------

# cluster <- multidplyr::new_cluster(10)
# multidplyr::cluster_library(cluster,"magrittr")
# multidplyr::cluster_assign(cluster, fun_barcode=fun_barcode)
# multidplyr::cluster_assign(cluster,fun_tn_type=fun_tn_type)
# multidplyr::cluster_assign(cluster,fn_transform_cnv=fn_transform_cnv)  

cnv %>% 
  dplyr::group_by(cancer_types) %>%
  # multidplyr::partition(cluster) %>%
  dplyr::mutate(cnv = purrr::map2(cnv,cancer_types, .f = fn_transform_cnv)) %>% 
  # dplyr::collect() %>%
  tidyr::unnest() %>% 
  dplyr::ungroup() %>%
  # dplyr::select(-PARTITION_ID) %>%
  dplyr::filter(! is.na(value)) %>% 
  dplyr::rename(cnv = value) -> cnv_df
# parallel::stopCluster(cluster)

# cl<-33
# cluster <- multidplyr::create_cluster(cl)

# multidplyr::cluster_library(cluster,"magrittr")
# multidplyr::cluster_assign_value(cluster,fun_barcode=fun_barcode)
# multidplyr::cluster_assign_value(cluster,fun_tn_type=fun_tn_type)
# multidplyr::cluster_assign_value(cluster,fn_transform_exp=fn_transform_exp)
  
expr %>% 
  dplyr::group_by(cancer_types) %>%
  # multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(expr = purrr::map2(expr,cancer_types, .f = fn_transform_exp)) %>% 
  # dplyr::collect() %>%
  tidyr::unnest() %>% 
  dplyr::ungroup() %>%
  dplyr::filter(! is.na(value)) %>% 
  dplyr::rename(expr = value) -> expr_df

# calculate ---------------------------------------------------------------
fn_get_cor<-function(data,cancer){
  message(glue::glue('correlation: {cancer}'))
  data %>%
    dplyr::group_by(symbol) %>%
    dplyr::do(
      tryCatch(
        broom::tidy(
          cor.test(.$expr,.$cnv,method = c("spearman")))
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(spm=estimate) %>%
    dplyr::mutate(fdr= p.adjust(p.value, method = "fdr")) %>%
    dplyr::mutate(logfdr=-log10(fdr)) %>%
    dplyr::mutate(logfdr=ifelse(logfdr>50,50,logfdr)) %>%
    dplyr::select(symbol,spm,fdr,logfdr) ->.out
  return(.out)
}

df_expr_cnv <- 
  cnv_df %>% 
  dplyr::inner_join(expr_df, by = c("cancer_types", "entrez", "symbol","sample")) %>%
  tidyr::nest(entrez,symbol,cnv,sample,expr,.key="data")


# multidplyr::cluster_library(cluster,"magrittr")
# multidplyr::cluster_assign_value(cluster,fn_get_cor=fn_get_cor)

df_expr_cnv %>%
  dplyr::group_by(cancer_types) %>%
  # multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(spm=purrr::map2(data,cancer_types,fn_get_cor)) %>%
  # dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-data) ->df_expr_cnv_cor

df_expr_cnv_cor %>%
  readr::write_rds(file.path(gsca_v2_path,"cnv","pancan34_all_exp-cor-cnv.NEW.IdTrans.rds.gz"),compress="gz")

save.image(file.path(git_path,"cnv_new","rda","2.cnv_cor_exp.rda"))

# parallel::stopCluster(cluster)