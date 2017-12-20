library(magrittr)
# get data ----------------------------------------------------------------

rppa_score <- readr::read_rds(file.path("/data/TCGA/TCGA_data/pancan32_rppa_score.rds.gz"))
expr <- readr::read_rds(file.path("/data/TCGA/TCGA_data/pancan33_expr.rds.gz"))

colnames(expr)[2] <-"filter_expr"

# functions ---------------------------------------------------------------
# get expression median cluster
fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
} #get short barcode from long barcode
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
}
fun_median_cluster <- function(cancer_types, filter_expr){
  cancer <- cancer_types
  filter_expr %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type == "01") %>%
    dplyr::mutate(barcode = fun_barcode(barcode)) %>% 
    dplyr::select(-type) %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::mutate(cluster = ifelse(expr >= median(expr), 2, 1)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(cancer_types = cancer) %>% 
    tidyr::nest(-cancer_types, .key = median_cluster)
}

# merge rppa score and expression median cluster
fn_cluster_rppa <- function(median_cluster, rppa){
  median_cluster %>% 
    dplyr::filter(!is.na(expr)) %>% 
    dplyr::mutate(cluster = as.factor(cluster)) %>% 
    dplyr::mutate(cluster = plyr::revalue(cluster, replace = c("1" = "down", "2" = "up"))) %>% 
    dplyr::inner_join(rppa, by = "barcode") %>% 
    dplyr::distinct() -> merged_clean
}

# t test used to get the difference between high and low expression 
fn_test <- function(merged_clean, cancer_types){
  print(cancer_types)
  merged_clean %>%
    dplyr::group_by(symbol, pathway) %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          t.test(score ~ cluster, data = .),
          error = function(e){1},
          warning = function(e){1})
      )
    ) %>% 
    dplyr::select(symbol, pathway, p.value) %>%
    dplyr::group_by(symbol) %>% 
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
    # dplyr::mutate(bfi = p.adjust(p.value, method = "bonferroni")) %>% 
    dplyr::ungroup() -> clean_pval
  #9761, 11901
  merged_clean %>% 
    tidyr::spread(key = cluster, value = score) %>% 
    dplyr::group_by(symbol, pathway) %>% 
    dplyr::summarise(diff = mean(up, na.rm = T) - mean(down, na.rm = T)) %>% 
    dplyr::ungroup() -> clean_diff

  clean_pval %>% 
    dplyr::inner_join(clean_diff, by = c("symbol", "pathway")) #%>%
    # dplyr::filter(is.na(p.value)) %>%
    # # dplyr::mutate(pathway = plyr::revalue(pathway, pathway_replace)) %>% 
    # dplyr::mutate(class = ifelse(fdr <= 0.05 & diff > 0, "Activation", "None")) %>% 
    # dplyr::mutate(class = ifelse(fdr <= 0.05 & diff < 0, "Inhibition", class))
}

fn_ai_n <- function(symbol, data){
  print(symbol)
  data %>% 
    dplyr::group_by(pathway, class) %>% 
    dplyr::count() %>% 
    dplyr::ungroup() %>% 
    tidyr::spread(key = class, value = n) %>% 
    dplyr::mutate(Activation = ifelse(rep(tibble::has_name(., "Activation"), 10), Activation, 0)) %>% 
    dplyr::mutate(Inhibition = ifelse(rep(tibble::has_name(., "Inhibition"), 10), Inhibition, 0)) %>% 
    tidyr::replace_na(replace = list(Activation = 0, Inhibition = 0, None = 0)) %>% 
    dplyr::select(pathway, Activation, Inhibition, None)
}
# calculation -------------------------------------------------------------
pathway_replace <- c(
  "PI3KAKT"="PI3K/AKT",
  "RASMAPK"="RAS/MAPK",
  "TSCmTOR"="TSC/mTOR",
  "CellCycle"="Cell Cycle",
  "DNADamage"="DNA Damage Response"
)
expr %>% 
  purrr::pmap(fun_median_cluster) %>% 
  dplyr::bind_rows() -> expr_median_cluster

expr_median_cluster %>% 
  dplyr::inner_join(rppa_score, by = "cancer_types") -> gene_rppa 

cl <- 15
cluster <- multidplyr::create_cluster(cl)
gene_rppa %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_cluster_rppa", fn_cluster_rppa)  %>%
  multidplyr::cluster_assign_value("fn_test", fn_test) %>%
  dplyr::mutate(merged_clean = purrr::map2(median_cluster, rppa, fn_cluster_rppa)) %>% 
  dplyr::select(cancer_types, merged_clean) %>% 
  dplyr::mutate(diff_pval = purrr::map2(merged_clean, cancer_types, fn_test)) %>%
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>%
  dplyr::select(-merged_clean) %>%
  tidyr::unnest() -> gene_rppa_pval 
parallel::stopCluster(cluster)

gene_rppa_pval %>%
  dplyr::filter(!is.na(p.value)) %>% 
  dplyr::mutate(pathway = plyr::revalue(pathway, pathway_replace)) %>% 
  dplyr::mutate(class = ifelse(fdr < 0.05 & diff > 0, "Activation", "None")) %>% 
  dplyr::mutate(class = ifelse(fdr < 0.05 & diff < 0, "Inhibition", class)) -> gene_rppa_sig_pval_class

gene_rppa_sig_pval_class %>%
  dplyr::select("symbol", "pathway", "cancer_types", "class", "p.value", "fdr") %>%
  tidyr::nest(-symbol) %>%
  dplyr::mutate(st = purrr::map2(symbol, data, .f = fn_ai_n)) -> gene_ai_n


# output ------------------------------------------------------------------


gene_ai_n %>%
  readr::write_rds("/data/GSCALite/TCGA/rppa/pan32_gene_activate.inhibit_pathway_data.rds.gz",compress = "gz")

gene_ai_n %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(a=Activation/32,i=Inhibition/32,n=None/32) %>%
  dplyr::select(symbol, pathway, a, i, n) %>%
  tidyr::nest(-symbol) ->gene_percent

gene_percent %>%
  readr::write_rds("/data/GSCALite/TCGA/rppa/pan32_gene_A-I-N_percent.rds.gz",compress = "gz")