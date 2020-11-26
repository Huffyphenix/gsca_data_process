
# survival analysis of gene expression------------------------------------------

library(survival)
library(magrittr)

# path --------------------------------------------------------------------

data_path <- "/home/huff/data/GSCA"
gsca_path <- file.path("/home/huff/data/GSCA")
git_path <- "/home/huff/github/gsca_data_process/expr"

survival <- readr::read_rds(file.path(data_path,"clinical","pancan34_clinical_stage_survival_subtype.rds.gz")) %>%
  dplyr::select(cancer_types,survival)

expr <- readr::read_rds(file.path(data_path,"expr","pancan33_expr.IdTrans.rds.gz"))


gene_symbol <- readr::read_rds('/home/huff/data/GSCA/id/NCBI_id_in_TCGA-final.rds.gz')
# combine data ------------------------------------------------------------
expr %>%
  dplyr::mutate(gene_expr = purrr::map(expr,.f=function(.x){
    .x %>%
      dplyr::filter(symbol %in% gene_symbol$symbol) %>% 
      dplyr::rename(entrez = entrez_id) %>% 
      dplyr::mutate(entrez = as.numeric(entrez)) %>% 
      dplyr::group_by(entrez, symbol) %>% 
      tidyr::nest()
  })) %>%
  dplyr::select(-expr) %>%
  tidyr::unnest() -> expr_by_gene

expr_sur_file_list <- grep("*_survival.exp.rds.gz",dir(file.path(data_path,"expr/cancer_gene_survival_separate")),value = TRUE)
data_done <- purrr::pmap(list(expr_sur_file_list),.f=function(.x){strsplit(.x,split = "_survival")[[1]][1]}) %>% unlist() 
expr_by_gene %>%
  dplyr::mutate(cancer_symbol = paste(cancer_types,symbol,sep="_")) %>%
  dplyr::filter(!cancer_symbol %in% data_done) -> expr_by_gene_not_done
# functions ---------------------------------------------------------------
source("/home/huff/github/gsca_data_process/expr/02.functions_survival.R")


# calculation -------------------------------------------------------------
cluster <- multidplyr::new_cluster(5)
multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_assign(cluster, fn_survival_res=fn_survival_res)
multidplyr::cluster_assign(cluster, fn_transform_samples=fn_transform_samples)
multidplyr::cluster_assign(cluster, survival_group=survival_group)
multidplyr::cluster_assign(cluster, fn_expr_group=fn_expr_group)
multidplyr::cluster_assign(cluster, fn_cox_logp=fn_cox_logp)
multidplyr::cluster_assign(cluster, fn_survival=fn_survival)
multidplyr::cluster_assign(cluster, res_path=res_path)
multidplyr::cluster_assign(cluster, gene_symbol=gene_symbol)
multidplyr::cluster_assign(cluster, survival=survival)



expr_by_gene_not_done %>%
  dplyr::group_by(cancer_types,symbol) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(survival_res = purrr::pmap(list(cancer_types,data,symbol),.f=fn_survival_res,.survival=survival)) %>%
  dplyr::collect() %>%
  dplyr::select(-data)-> pan33_expr_survival

# pan33_expr_survival %>%
#   tidyr::nest(-cancer_types) %>%
#   readr::write_rds(file.path(gsca_path,"expr","pan33_expr_survival.rds.gz"),compress = "gz")

save.image(file.path(git_path,"02.gene_expr_survival.rda"))
parallel::stopCluster(cluster)
