# Immune cell association with gene snv ----------------------------
# Method: wilcoxon test
# Author: Huffy
# Date: 2020-10-25

library(magrittr)
# path --------------------------------------------------------------------

data_path <- "/home/huff/data/GSCA"
git_path <- "/home/huff/github/gsca_data_process/immune"
# load data ---------------------------------------------------------------

immune_cell_data <- readr::read_rds(file.path(data_path,"TIL","pan33_ImmuneCellAI.rds.gz"))
mc3_pass <- readr::read_rds(file.path(data_path,"snv","mc3_pass.IdTrans.maf.rds.gz"))
snv_data <- mc3_pass@data %>%
  dplyr::rename(symbol=Hugo_Symbol,cancer_types=Cancer_Types)


# snv data process ------------------------------------------------------------
snv_data %>%
  dplyr::as.tbl() %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(data = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::group_by(symbol,Tumor_Sample_Barcode,sample, entrez, barcode) %>%
      tidyr::nest()
  })) %>%
  tidyr::unnest(data) %>%
  dplyr::mutate(n = purrr::map(data,.f=function(.x){nrow(.x)})) %>%
  tidyr::unnest(n) %>%
  dplyr::mutate(barcode=substr(Tumor_Sample_Barcode,1,16))-> snv_data_count

snv_data_count %>%
  dplyr::select(cancer_types,Tumor_Sample_Barcode,barcode,symbol,entrez,n) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::rename(snv=data) -> snv_data_for_calculation

# data combine ------------------------------------------------------------

immune_cell_data %>%
  dplyr::inner_join(snv_data_for_calculation,by="cancer_types") -> combine_data



# function ----------------------------------------------------------------
fn_wilcoxon <- function(.cancer,.immune,.snv){
  message(glue::glue('wilcoxon: {.cancer}'))
  
  .immune %>%
    tidyr::gather(-aliquot,-barcode,-sample_name,value="TIL",key="cell_type")%>%
    dplyr::select(-aliquot,-sample_name) -> .immune
  .snv %>%
    dplyr::filter(!is.na(symbol)) %>%
    dplyr::distinct()-> .snv  
  
  .combine <- .immune %>%
    dplyr::inner_join(.snv,by="barcode")
  
  .combine %>%
    tidyr::nest(-symbol,-cell_type,-entrez) %>%
    dplyr::mutate(cor = purrr::map(data,.f=fn_spm)) 
}

# calculation -------------------------------------------------------------

combine_data %>%
  dplyr::mutate(wilcoxon_test = purrr::pmap(list(cancer_types,ImmuneCellAI,snv),.f=fn_wilcoxon))
# save image --------------------------------------------------------------

rm(mc3_pass)
