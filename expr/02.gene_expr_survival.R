
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

survival %>%
  dplyr::inner_join(expr,by="cancer_types") -> combine_data

# functions ---------------------------------------------------------------
source("/home/huff/github/gsca_data_process/expr/02.functions_survival.R")


# calculation -------------------------------------------------------------

combine_data %>%
  dplyr::mutate(survival_res = purrr::pmap(list(cancer_types,survival,expr),.f=fn_survival_res)) -> pan33_expr_survival

pan33_expr_survival %>%
  readr::write_rds(file.path(gsca_path,"expr","pan33_expr_survival.rds.gz"),compress = "gz")

save.image(file.path(git_path,"02.gene_expr_survival.rda"))
