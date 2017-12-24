library(magrittr)
config <- list()
config$database <- c("/data/GSCALite/TCGA/rppa")

gene_rppa_pval <- readr::read_rds(file.path(config$database,"pan32_gene_A-I-N_percent.rds.gz"))

# functions ---------------------------------------------------------------

fn_ai_n <- function(symbol, data){
  print(symbol)
  data %>% 
    dplyr::group_by(pathway, class) %>% 
    dplyr::count() %>% 
    dplyr::ungroup() %>% 
    tidyr::spread(key = class, value = n) ->.tmp
  nrow(.tmp) ->n
  .tmp %>%
    dplyr::mutate(Activation = ifelse(rep(tibble::has_name(., "Activation"), n), Activation, 0)) %>% 
    dplyr::mutate(Inhibition = ifelse(rep(tibble::has_name(., "Inhibition"), n), Inhibition, 0)) %>% 
    tidyr::replace_na(replace = list(Activation = 0, Inhibition = 0)) %>% #, None = 0
    dplyr::select(pathway, Activation, Inhibition)#, None
}
# calculation -------------------------------------------------------------
pathway_replace <- c(
  "PI3KAKT"="PI3K/AKT",
  "RASMAPK"="RAS/MAPK",
  "TSCmTOR"="TSC/mTOR",
  "CellCycle"="Cell Cycle",
  "DNADamage"="DNA Damage Response"
)


gene_rppa_pval %>%
  tidyr::unnest() %>%
  tidyr::drop_na() %>%
   # dplyr::filter(symbol %in% c("PRAMEF10","ZZZ3")) %>%
  # dplyr::filter(fdr <= 0.05) %>%
  dplyr::mutate(pathway = plyr::revalue(pathway, pathway_replace)) %>%
  dplyr::mutate(class = ifelse(fdr <= 0.05 & diff > 0, "Activation", "None")) %>%
  dplyr::mutate(class = ifelse(fdr <= 0.05 & diff < 0, "Inhibition", class)) %>%
  dplyr::filter(class!="None") -> gene_rppa_sig_pval_class


cl <- 10
cluster <- multidplyr::create_cluster(cl)
gene_rppa_sig_pval_class %>%
  # dplyr::filter(symbol=="PRAMEF10") %>%
  # head(10000) %>%
  dplyr::select("symbol", "pathway", "cancer_types", "class", "p.value", "fdr") %>%
  tidyr::nest(-symbol) %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_ai_n", fn_ai_n)  %>%
  dplyr::mutate(st = purrr::map2(symbol, data, .f = fn_ai_n)) %>%
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  dplyr::select(-data,-PARTITION_ID) -> gene_ai_n
parallel::stopCluster(cluster)

gene_ai_n %>%
  tidyr::unnest() %>%
  dplyr::mutate(a=Activation/32,i=Inhibition/32) %>%
  dplyr::mutate(n=1-a-i) %>%
  dplyr::select(symbol, pathway, a, i, n) %>%
  tidyr::nest(-symbol) ->gene_percent

gene_percent %>%
  readr::write_rds("/data/GSCALite/TCGA/rppa/pan32_gene_activate.inhibit_pathway_percent.rds.gz",compress = "gz")