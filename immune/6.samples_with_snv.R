############### samples with snv data ################
gsca_v2_path <- "/home/huff/data/GSCA"

library(magrittr)

immune_cell_data <- readr::read_rds(file.path(gsca_v2_path,"TIL","pan33_ImmuneCellAI.rds.gz"))



fn_immune_snv <- function(.cancer_types,.immune){
  .snv_data <- readr::read_rds(file.path(gsca_v2_path,"snv","sub_cancer_maf_tsv",paste(.cancer_types,"maf_data.IdTrans.tsv.rds.gz",sep = "_")))%>%
    dplyr::mutate(group = ifelse(Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Splice_Site","Frame_Shift_Del","In_Frame_Del","In_Frame_Ins"), "2_mutated","1_nonmutated")) %>%
    dplyr::mutate(barcode=substr(Tumor_Sample_Barcode,1,16)) %>%
    dplyr::select(Hugo_Symbol,entrez,barcode,group) 
  
  .snv_data$barcode %>%
    unique() -> sample_with_snv
  
  sample_with_snv
}


immune_cell_data %>%
  dplyr::mutate(sample_with_snv = purrr::map2(cancer_types,ImmuneCellAI,.f=fn_immune_snv)) %>%
  dplyr::select(-ImmuneCellAI ) -> pancan33_sample_with_snv

pancan33_sample_with_snv %>%
  readr::write_rds(file.path(gsca_v2_path,"TIL","pancan33_sample_with_snv.rds.gz"),compress="gz")