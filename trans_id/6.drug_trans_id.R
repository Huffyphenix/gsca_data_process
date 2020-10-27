
# durg trans id -----------------------------------------------------------


library(magrittr)

# GSCALite path -----------------------------------------------------------

gscalite_path <- "/home/huff/data/GSCALite"
gsca_v2_path <- "/home/huff/data/GSCA"


# read gene info ------------------------------------------------

search_symbol <- readr::read_tsv(file.path(gsca_v2_path,"id","NCBI_id_in.with_TCGAsym.tsv")) %>%
  dplyr::select(entrez_id, NCBI_sym, TCGA_sym)


# ctrp drug data ---------------------------------------------------------------

ctrp_exp_spearman.rds.gz <- readr::read_rds(file.path(gscalite_path,"Drug","ctrp_exp_spearman.rds.gz"))


# translate ctrp gene id --------------------------------------------------
search_symbol %>%
  dplyr::inner_join(ctrp_exp_spearman.rds.gz %>%
                      dplyr::rename(TCGA_sym=symbol),by="TCGA_sym") %>%
  dplyr::rename(symbol=NCBI_sym) %>%
  dplyr::select(-TCGA_sym) -> ctrp_exp_spearman.rds.gz.IdTrans

ctrp_exp_spearman.rds.gz.IdTrans %>%
  readr::write_rds(file.path(gsca_v2_path,"drug","ctrp_exp_spearman.IdTrans.rds.gz"),compress = "gz")

# ctrp drug data ---------------------------------------------------------------

gdsc_exp_spearman.rds.gz <- readr::read_rds(file.path(gscalite_path,"Drug","gdsc_exp_spearman.rds.gz"))


# translate gdsc gene id --------------------------------------------------
search_symbol %>%
  dplyr::inner_join(gdsc_exp_spearman.rds.gz %>%
                      dplyr::rename(TCGA_sym=symbol),by="TCGA_sym")%>%
  dplyr::rename(symbol=NCBI_sym) %>%
  dplyr::select(-TCGA_sym)  -> gdsc_exp_spearman.rds.gz.IdTrans

gdsc_exp_spearman.rds.gz.IdTrans %>%
  readr::write_rds(file.path(gsca_v2_path,"drug","gdsc_exp_spearman.IdTrans.rds.gz"),compress = "gz")

