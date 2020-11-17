##################################################################
#++++++++Tanslate symbol in TCGA mutation into NCBI symbol +++++++++++
library(magrittr)

# GSCALite path -----------------------------------------------------------

gscalite_path <- "/home/huff/data/GSCALite/TCGA"
gsca_v2_path <- "/home/huff/data/GSCA"


# read cnv gene info ------------------------------------------------

cnv_symbol_search_symbol <- readr::read_rds(file.path(gsca_v2_path,"id","cnv_symbol_search_symbol_final.rds.gz"))
cnv_symbol_entrez <- cnv_symbol_search_symbol %>%
  dplyr::select(cnvsymbol,entrez,symbol)
####################################################################
# cnv files ---------------------------------------------------------------


# read cnv  ---------------------------------------------------------------

pancan34_cnv.rds.gz <- readr::read_rds(file.path(gscalite_path,"cnv","pancan34_cnv.rds.gz")) 


# translate cnv  ----------------------------------------------------------

pancan34_cnv.rds.gz %>%
  dplyr::mutate(cnv_trans = purrr::map(cnv,.f=function(.x){
    cnv_symbol_entrez %>%
      dplyr::inner_join(.x %>% dplyr::rename("cnvsymbol"="symbol"),by="cnvsymbol") %>%
      dplyr::select(-cnvsymbol)
  })) %>%
  dplyr::select(-cnv) %>%
  dplyr::rename("cnv"="cnv_trans") -> pancan34_cnv.rds.gz.IdTrans

pancan34_cnv.rds.gz.IdTrans %>% readr::write_rds(file.path(gsca_v2_path,"cnv","pancan34_cnv.IdTrans.rds.gz"))


# read cnv threshold ------------------------------------------------------
pancan34_cnv_threshold.rds.gz <- readr::read_rds(file.path(gscalite_path,"cnv","pancan34_cnv_threshold.rds.gz"))


# translate cnv threshold -------------------------------------------------

pancan34_cnv_threshold.rds.gz %>%
  dplyr::mutate(cnv_trans = purrr::map(cnv,.f=function(.x){
    cnv_symbol_entrez %>%
      dplyr::inner_join(.x %>% dplyr::rename("cnvsymbol"="symbol"),by="cnvsymbol") %>%
      dplyr::select(-cnvsymbol)
  })) %>%
  dplyr::select(-cnv) %>%
  dplyr::rename("cnv"="cnv_trans") -> pancan34_cnv_threshold.rds.gz.IdTrans
s
pancan34_cnv_threshold.rds.gz.IdTrans  %>% readr::write_rds(file.path(gsca_v2_path,"cnv","pancan34_cnv_threshold.IdTrans.rds.gz"))

# read cnv percentage -----------------------------------------------------

pancan34_cnv_percent.rds.gz <- readr::read_rds(file.path(gscalite_path,"cnv","pancan34_cnv_percent.rds.gz"))

# translate cnv percentage -----------------------------------------------------

pancan34_cnv_percent.rds.gz %>%
  dplyr::mutate(cnv_trans = purrr::map(cnv,.f=function(.x){
    cnv_symbol_entrez %>%
      dplyr::inner_join(.x %>% dplyr::rename("cnvsymbol"="symbol"),by="cnvsymbol") %>%
      dplyr::select(-cnvsymbol)
  })) %>%
  dplyr::select(-cnv) %>%
  dplyr::rename("cnv"="cnv_trans") -> pancan34_cnv_percent.rds.gz.IdTrans

pancan34_cnv_percent.rds.gz.IdTrans  %>% readr::write_rds(file.path(gsca_v2_path,"cnv","pancan34_cnv_percent.IdTrans.rds.gz"))

# read cnv gene correlation -----------------------------------------------
# need a new correlation results

