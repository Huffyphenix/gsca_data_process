##################################################################
#++++++++Tanslate symbol in TCGA mutation into NCBI symbol +++++++++++
library(magrittr)
# GSCALite path -----------------------------------------------------------

gscalite_path <- "/home/huff/data/GSCALite/TCGA"
gsca_v2_path <- "/home/huff/data/GSCA"


# read cnv gene info ------------------------------------------------

cnv_symbol_search_symbol <- readr::read_rds(file.path(gsca_v2_path,"id","cnv_symbol_search_symbol_final.rds.gz"))

####################################################################
# cnv files ---------------------------------------------------------------


# read cnv  ---------------------------------------------------------------

pancan34_cnv.rds.gz <- readr::read_rds(file.path(gscalite_path,"cnv","pancan34_cnv.rds.gz")) 


# translate cnv  ----------------------------------------------------------

pancan34_cnv.rds.gz %>%
  dplyr::mutate(cnv_trans = purrr::map(cnv,.f=function(.x){
    .x %>%
      dplyr::rename("cnvsymbol"="symbol") %>%
      dplyr::inner_join(cnv_symbol_search_symbol,by="cnvsymbol") %>%
      dplyr::select(-cnvsymbol) %>%
      dplyr::rename("symbol"="NCBI_sym")
  })) %>%
  dplyr::select(-cnv) %>%
  dplyr::rename("cnv"="cnv_trans") -> pancan34_cnv.rds.gz.IdTrans

pancan33_expr.rds.gz.IdTrans %>% readr::write_rds(file.path(gsca_v2_path,"cnv","pancan33_expr.IdTrans.rds.gz"))

# read cnv threshold ------------------------------------------------------
pancan34_cnv_threshold.rds.gz <- readr::read_rds(file.path(gscalite_path,"cnv","pancan34_cnv_threshold.rds.gz"))


# read cnv percentage -----------------------------------------------------

# read cnv gene correlation -----------------------------------------------


# translate SNV -----------------------------------------------------------


# read SNV ---------------------------------------------------------------------

pancan33_snv.rds.gz <- readr::read_rds(file.path(gscalite_path,"snv","pancan33_snv.rds.gz"))



# read maf ----------------------------------------------------------------

syn_mutation_syn7824274_mc3_public.pass.simplification.maf.rds.gz <- readr::read_rds(file.path(gscalite_path,"snv","syn_mutation_syn7824274_mc3_public.pass.simplification.maf.rds.gz"))

syn_mutation_syn7824274_mc3_public.pass.maf.rds.gz <- readr::read_rds(file.path(gscalite_path,"snv","syn_mutation_syn7824274_mc3_public.pass.maf.rds.gz"))


