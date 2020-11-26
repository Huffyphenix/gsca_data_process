
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
    dplyr::filter(mutation == "EffectiveMut") %>%
    dplyr::select(Tumor_Sample_Barcode) %>%
    unique() %>% nrow() -> mut_n
 
  mutation_type %>%
    dplyr::filter(mutation == "NonEffectiveMut") %>%
    dplyr::select(Tumor_Sample_Barcode) %>%
    unique() %>% nrow() -> nmut_n

  tibble::tibble(EffectiveMut_n = mut_n,
                 NonEffectiveMut_n=nmut_n) 
}

fn_mutate_type <- function(mutation_type){
  mutation_type %>%
    dplyr::filter(mutation == "EffectiveMut") %>%
    dplyr::select(Variant_Classification) %>%
    table() %>%
    as.data.frame() %>%
    dplyr::as.tbl() %>%
    dplyr::rename(Variant_Classification=".") %>%
    tidyr::nest(EffectiveMutType=everything()) %>%
    dplyr::mutate(x="x") -> EffectiveMutType
  mutation_type %>%
    dplyr::filter(mutation == "NonEffectiveMut") %>%
    dplyr::select(Variant_Classification) %>%
    table() %>%
    as.data.frame() %>%
    dplyr::as.tbl() %>%
    dplyr::rename(Variant_Classification=".") %>%
    tidyr::nest(NonEffectiveMutType=everything()) %>%
    dplyr::mutate(x="x")-> NonEffectiveMutType
  EffectiveMutType %>%
    dplyr::inner_join(NonEffectiveMutType,by="x") %>%
    dplyr::select(-x)
}
oncoplot_effective_mutation <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Splice_Site","Frame_Shift_Del","In_Frame_Del","In_Frame_Ins")
excluded_mut_type <- c("Silent", "Intron", "IGR", "3'UTR", "5'UTR", "3'Flank", "5'Flank") # https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats


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
    dplyr::mutate(mutation = ifelse(Variant_Classification %in% oncoplot_effective_mutation,"EffectiveMut","NonEffectiveMut")) -> snv_mut_diagnos
  # overall mutation rate: sample level
  snv_mut_diagnos %>%
    tidyr::nest(data = c(Variant_Classification, VARIANT_CLASS, Tumor_Sample_Barcode, barcode16, cancer_types, mutation)) %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::mutate(count = purrr::map(data,fn_get_count)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cols = c(count)) %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::rename(EffectiveMut_sample = EffectiveMut_n) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cancer_types=cancer_type,cancer_sample=cancer_sample_n) %>%
    dplyr::mutate(per = 100*EffectiveMut_sample/cancer_sample) %>%
    dplyr::rename("symbol" = "Hugo_Symbol")-> cancer.mut_count 
  snv_mut_diagnos %>%
    tidyr::nest(data = c(Variant_Classification, Tumor_Sample_Barcode, barcode16, cancer_types, mutation)) %>%
    dplyr::group_by(Hugo_Symbol,VARIANT_CLASS) %>%
    dplyr::mutate(count = purrr::map(data,fn_get_count)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cols = c(count)) %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::ungroup() %>%
    dplyr::rename("symbol" = "Hugo_Symbol")-> cancer.mut_count.type
  
  # Variant_type_not_in <- setdiff(c("deletion", "insertion", "SNV", "substitution"),unique(cancer.mut_count.type$VARIANT_CLASS))
  # if(length(Variant_type_not_in) > 0){
  #   for (type in Variant_type_not_in) {
  #     cancer.mut_count.type %>%
  #       head(1) %>%
  #       dplyr::mutate(VARIANT_CLASS=type,mut_site_n=NA) %>%
  #       rbind(cancer.mut_count.type) -> cancer.mut_count.type
  #   }
  # } else {
  #   cancer.mut_count.type -> cancer.mut_count.type
  # }
  cancer.mut_count.type %>% 
    tidyr::nest(VARIANT_SITE_CLASS = c(VARIANT_CLASS, EffectiveMut_n, NonEffectiveMut_n)) -> cancer.mut_count.type
  
  cancer.mut_count %>%
    dplyr::inner_join(cancer.mut_count.type,by=c("symbol","entrez")) %>%
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


