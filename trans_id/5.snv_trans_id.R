
# snv files trans id ------------------------------------------------------


# symbol process of snv data ----------------------------------------------

# Library -----------------------------------------------------------------

library(magrittr)

# Load data ---------------------------------------------------------------
rda_path <- "/home/huff/github/gsca_data_process/trans_id/rda/"
gscalite_path <- "/home/huff/data/GSCALite/TCGA"
gsca_v2_path <- "/home/huff/data/GSCA"

# read cnv gene info ------------------------------------------------

snv_symbol_search_symbol <- readr::read_rds('trans_id/rda/snv_symbol_search_symbol_final.rds.gz')
snv_symbol_entrez <- snv_symbol_search_symbol %>%
  dplyr::select(snv_symbol,entrez,symbol)

####################################################################
# snv files ---------------------------------------------------------------


# read snv  ---------------------------------------------------------------

mc3_pass <- readr::read_rds(file.path(gscalite_path,"snv", "01-snv_mutation_mc3_public.pass.filtered_maf.rds.gz"))
snv_data <- mc3_pass@data

# translate cnv  ----------------------------------------------------------

snv_data %>%
  dplyr::as.tbl() %>%
  dplyr::rename("snv_symbol"="Hugo_Symbol") %>%
  dplyr::inner_join(snv_symbol_entrez,by="snv_symbol") %>%
  dplyr::select(-snv_symbol,-Entrez_Gene_Id) %>%
  dplyr::mutate(barcode = substr(Tumor_Sample_Barcode,1,14)) %>%
  dplyr::rename(Hugo_Symbol=symbol) -> snv_data.IdTrans

snv_data.IdTrans %>% dplyr::select(Cancer_Types,Tumor_Sample_Barcode) %>% unique() -> snv_clinical.IdTrans


pancan34_cnv.rds.gz.IdTrans %>% readr::write_rds(file.path(gsca_v2_path,"cnv","pancan34_cnv.IdTrans.rds.gz"))

mc3_pass.IdTrans <- read.maf(maf=snv_data.IdTrans,
                             clinicalData = snv_clinical.IdTrans)
mc3_pass.IdTrans %>%
  readr::write_rds(file.path(gsca_v2_path,"snv","mc3_pass.IdTrans.maf.rds.gz"),compress="gz")

# snv count ---------------------------------------------------------------

gene_list_snv_count.new.new.rds.gz <- readr::read_rds(file.path(gscalite_path,"snv","gene_list_snv_count.new.new.rds.gz"))

# translate snv count  ----------------------------------------------------

gene_list_snv_count.new.new.rds.gz %>%
  dplyr::mutate(snv_trans = purrr::map(mut_count,.f=function(.x){
    snv_symbol_entrez %>%
      dplyr::inner_join(.x %>% dplyr::rename("snv_symbol"="symbol"),by="snv_symbol") %>%
      dplyr::select(-snv_symbol)
  })) %>%
  dplyr::select(-mut_count) %>%
  dplyr::rename("snv"="snv_trans") -> gene_list_snv_count.rds.gz.IdTrans

gene_list_snv_count.rds.gz.IdTrans %>% readr::write_rds(file.path(gsca_v2_path,"snv","gene_list_snv_count.IdTrans.rds.gz"))

# snv survival ---------------------------------------------------------------

pancan32_snv_survival_genelist_sig_pval.rds.gz <- readr::read_rds(file.path(gscalite_path,"snv","pancan32_snv_survival_genelist_sig_pval.rds.gz"))

# translate snv count  ----------------------------------------------------

pancan32_snv_survival_genelist_sig_pval.rds.gz %>%
  dplyr::mutate(snv_trans = purrr::map(diff_pval,.f=function(.x){
    snv_symbol_entrez %>%
      dplyr::inner_join(.x %>% dplyr::rename("snv_symbol"="symbol"),by="snv_symbol") %>%
      dplyr::select(-snv_symbol)
  })) %>%
  dplyr::select(-diff_pval) %>%
  dplyr::rename("snv_survival"="snv_trans") -> pancan32_snv_survival_genelist_sig_pval.rds.gz.IdTrans

pancan32_snv_survival_genelist_sig_pval.rds.gz.IdTrans %>% readr::write_rds(file.path(gsca_v2_path,"snv","pancan32_snv_survival_genelist_sig_pval.IdTrans.rds.gz"))

