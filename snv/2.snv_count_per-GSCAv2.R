
library(magrittr)
######################## snv count and percentage ############################################
# data path ---------------------------------------------------------------

gscalite_path <- "/home/huff/data/GSCALite/TCGA"
gsca_v2_path <- "/home/huff/data/GSCA"
rda_path <- "/home/huff/github/gsca_data_process/snv"
# load data ---------------------------------------------------------------

sample_info <- readr::read_rds(file.path("/home/huff/data/GSCA/sample_info.rds.gz")) %>%
  dplyr::mutate(barcode16=substr(barcode,1,16)) %>%
  dplyr::select(barcode16,cancer_types) %>%
  unique()

sample_info$cancer_types %>% unique() -> cancer_types

# get snv count -----------------------------------------------------------

fn_get_count <- function(mutation_type){
  mutation_type %>%
    dplyr::filter(mutation == "Mut") %>%
    dplyr::select(Tumor_Sample_Barcode) %>%
    unique() %>% nrow() -> mut_n
  data.frame(mut_n = mut_n)
}

snv_count_combine <- tibble::tibble()
for (cancer_type in cancer_types) {
  print(cancer_type)
  snv_data <- readr::read_rds(file.path(gsca_v2_path,"snv","sub_cancer_maf_tsv",paste(cancer_type,"maf_data.IdTrans.tsv.rds.gz",sep = "_")))
  snv_data %>%
    dplyr::select(Tumor_Sample_Barcode) %>%
    unique() %>%
    nrow() -> cancer_sample_n
  snv_data %>% 
    dplyr::select(Hugo_Symbol,entrez,Variant_Classification,VARIANT_CLASS,Tumor_Sample_Barcode,barcode16,cancer_types) %>%
    dplyr::mutate(mutation = "Mut") %>%
    tidyr::nest(-c(Hugo_Symbol,VARIANT_CLASS,entrez)) %>%
    dplyr::group_by(Hugo_Symbol,VARIANT_CLASS) %>%
    dplyr::mutate(count = purrr::map(data,fn_get_count)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::mutate(mut_sum = sum(mut_n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cancer_types=cancer_type,cancer_sample=cancer_sample_n) %>%
    dplyr::mutate(per = 100*mut_sum/cancer_sample) %>%
    dplyr::rename("symbol" = "Hugo_Symbol")-> cancer.mut_count 
  Variant_type_not_in <- setdiff(c("deletion", "insertion", "SNV", "substitution"),unique(cancer.mut_count$VARIANT_CLASS))
  if(length(Variant_type_not_in) > 0){
    for (type in Variant_type_not_in) {
      cancer.mut_count %>%
        head(1) %>%
        dplyr::mutate(VARIANT_CLASS=type,mut_n=NA, mut_sum=NA,per=NA) %>%
        rbind(cancer.mut_count) -> cancer.mut_count
    }
  } else {
    cancer.mut_count -> cancer.mut_count
  }
  cancer.mut_count %>% 
    tidyr::spread(key=VARIANT_CLASS,value=mut_n)%>%
    dplyr::filter(!is.na(mut_sum)) %>%
    dplyr::mutate(deletion = ifelse(is.na(deletion),0,deletion),
                  insertion = ifelse(is.na(insertion),0,insertion),
                  SNV = ifelse(is.na(SNV),0,SNV),
                  substitution = ifelse(is.na(substitution),0,substitution)) %>%
    tidyr::nest(-cancer_types,-cancer_sample)  -> cancer.mut_count
  
  cancer.mut_count %>%
    readr::write_rds(file.path(gsca_v2_path,"snv/cancer_snv_count",paste(cancer_type,"snv_count.rds.gz",sep="_")),compress = "gz")
  if (nrow(snv_count_combine)<1) {
    snv_count_combine <- cancer.mut_count
  } else {
    snv_count_combine %>%
      rbind(cancer.mut_count) -> snv_count_combine
  }
}

snv_count_combine %>%
  readr::write_rds(file.path(gsca_v2_path,"snv","all_maf.snv_count_per.rds.gz"))

save.image(file.path(rda_path,"2.snv_count_per-GSCAv2.rda"))


