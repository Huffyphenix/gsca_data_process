
# methy files trans id ------------------------------------------------------


# symbol process of methy data ----------------------------------------------

# Library -----------------------------------------------------------------

library(magrittr)

# Load data ---------------------------------------------------------------
rda_path <- "/home/huff/github/gsca_data_process/trans_id/rda/"
gscalite_path <- "/home/huff/data/GSCALite/TCGA"
gsca_v2_path <- "/home/huff/data/GSCA"

# read methy gene info ------------------------------------------------

methy_symbol_search_symbol <- readr::read_rds(file.path(rda_path,'methy_symbol_search_symbol_final.rds.gz'))
methy_symbol_entrez <- methy_symbol_search_symbol %>%
  dplyr::select(methy_symbol,entrez,symbol)

####################################################################
# methy files ---------------------------------------------------------------


# read methy  ---------------------------------------------------------------

methy_data <- readr::read_rds(file.path(gscalite_path,"meth", "pancan33_meth.rds.gz"))

# translate methy  ----------------------------------------------------------
methy_data %>%
  dplyr::mutate(methy_trans = purrr::map(methy,.f=function(.x){
    methy_symbol_entrez %>%
      dplyr::inner_join(.x %>% dplyr::rename("methy_symbol"="symbol"),by="methy_symbol") %>%
      dplyr::select(-methy_symbol)
  })) %>%
  dplyr::select(-methy) %>%
  dplyr::rename("methy"="methy_trans") -> pancan33_meth.rds.gz.IdTrans

pancan33_meth.rds.gz.IdTrans %>%
  readr::write_rds(file.path(gsca_v2_path,"methy","pancan33_meth.IdTrans.rds.gz"),compress="gz")

# read methy diff ---------------------------------------------------------------

pan33_allgene_methy_diff.rds.gz <- readr::read_rds(file.path(gscalite_path,"meth", "pan33_allgene_methy_diff.rds.gz"))

# translate methy  ----------------------------------------------------------
pan33_allgene_methy_diff.rds.gz %>%
  dplyr::mutate(methy_trans = purrr::map(methy_comparison,.f=function(.x){
    methy_symbol_entrez %>%
      dplyr::inner_join(.x %>% dplyr::rename("methy_symbol"="symbol"),by="methy_symbol") %>%
      dplyr::select(-methy_symbol)
  })) %>%
  dplyr::select(-methy_comparison) %>%
  dplyr::rename("methy"="methy_trans") -> pan33_allgene_methy_diff.rds.gz.IdTrans

pan33_allgene_methy_diff.rds.gz.IdTrans %>%
  readr::write_rds(file.path(gsca_v2_path,"methy","pan14_allgene_methy_diff.IdTrans.rds.gz"),compress="gz")

# read methy survival ---------------------------------------------------------------

pancan32_meth_survival_genelist.rds.gz <- readr::read_rds(file.path(gscalite_path,"meth", "pancan32_meth_survival_genelist.rds.gz"))

# translate methy  ----------------------------------------------------------
pancan32_meth_survival_genelist.rds.gz %>%
  dplyr::mutate(methy_trans = purrr::map(diff_pval,.f=function(.x){
    methy_symbol_entrez %>%
      dplyr::inner_join(.x %>% dplyr::rename("methy_symbol"="symbol"),by="methy_symbol") %>%
      dplyr::select(-methy_symbol)
  })) %>%
  dplyr::select(-diff_pval) %>%
  dplyr::rename("methy"="methy_trans") -> pancan32_meth_survival_genelist.rds.gz.IdTrans

pancan32_meth_survival_genelist.rds.gz.IdTrans %>%
  readr::write_rds(file.path(gsca_v2_path,"methy","pancan32_meth_survival.IdTrans.rds.gz"),compress="gz")

# read methy exp cor ---------------------------------------------------------------

pancan34_all_gene_exp_cor_meth.rds.gz <- readr::read_rds(file.path(gscalite_path,"meth", "pancan34_all_gene_exp-cor-meth.rds.gz"))

# translate methy  ----------------------------------------------------------
pancan34_all_gene_exp_cor_meth.rds.gz %>%
  dplyr::mutate(methy_trans = purrr::map(spm,.f=function(.x){
    methy_symbol_entrez %>%
      dplyr::inner_join(.x %>% dplyr::rename("methy_symbol"="symbol"),by="methy_symbol") %>%
      dplyr::select(-methy_symbol)
  })) %>%
  dplyr::select(-spm) %>%
  dplyr::rename("methy"="methy_trans") -> pancan34_all_gene_exp_cor_meth.rds.gz.IdTrans

pancan34_all_gene_exp_cor_meth.rds.gz.IdTrans %>%
  readr::write_rds(file.path(gsca_v2_path,"methy","pancan34_all_gene_exp_cor_meth.IdTrans.rds.gz"),compress="gz")


# save image --------------------------------------------------------------

save.image(file.path(rda_path,"7.methy_trans_id.rda"))
