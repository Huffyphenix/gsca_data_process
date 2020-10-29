##################################################################
#++++++++Tanslate symbol in rppa data into NCBI symbol +++++++++++

library(magrittr)
# GSCALite path -----------------------------------------------------------

gscalite_path <- "/home/huff/data/GSCALite"
gsca_v2_path <- "/home/huff/data/GSCA"

# read search names  -----------------------------------------------------------
NCBI_id_in.with_TCGAsym <- readr::read_tsv(file.path(gsca_v2_path,"id","NCBI_id_in.with_TCGAsym.tsv")) %>%
  dplyr::filter(!is.na(NCBI_sym)) %>%
  dplyr::select(entrez_id,TCGA_sym,NCBI_sym)


# load rppa percent data ----------------------------------------------------------

pan32_gene_activate.inhibit_pathway_percent.rds.gz <- readr::read_rds(file.path(gscalite_path,"TCGA","rppa","pan32_gene_activate.inhibit_pathway_percent.rds.gz"))

# translate rppa percent id -----------------------------------------------
pan32_gene_activate.inhibit_pathway_percent.rds.gz %>%
  dplyr::rename(TCGA_sym=symbol) %>%
  dplyr::inner_join(NCBI_id_in.with_TCGAsym,by="TCGA_sym") %>%
  dplyr::rename(entrez=entrez_id,symbol=NCBI_sym) %>%
  dplyr::select(-TCGA_sym) %>%
  dplyr::select(entrez,symbol,data) -> pan32_gene_activate.inhibit_pathway_percent.rds.gz.IdTrans

pan32_gene_activate.inhibit_pathway_percent.rds.gz.IdTrans %>%
  readr::write_rds(file.path(gsca_v2_path,"rppa","pan32_gene_activate.inhibit_pathway_percent.IdTrans.rds.gz"))

# load rppa percent data ----------------------------------------------------------

pan32_gene_AIN_sig_pval_class.siplification.rds.gz <- readr::read_rds(file.path(gscalite_path,"TCGA","rppa","pan32_gene_A-I-N_sig_pval_class.siplification.rds.gz"))

# translate rppa percent id -----------------------------------------------
pan32_gene_AIN_sig_pval_class.siplification.rds.gz %>%
  dplyr::mutate(data=purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::rename(TCGA_sym=symbol) %>%
      dplyr::inner_join(NCBI_id_in.with_TCGAsym,by="TCGA_sym") %>%
      dplyr::rename(entrez=entrez_id,symbol=NCBI_sym) %>%  
      dplyr::select(-TCGA_sym) %>%
      dplyr::select(entrez,symbol,everything())
  })) -> pan32_gene_AIN_sig_pval_class.siplification.rds.gz.IdTrans

pan32_gene_AIN_sig_pval_class.siplification.rds.gz %>%
  readr::write_rds(file.path(gsca_v2_path,"rppa","pan32_gene_AIN_sig_pval_class.siplification.IdTrans.rds.gz"))
