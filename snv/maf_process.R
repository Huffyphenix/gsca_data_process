
library(magrittr)
library(finalfit)
library(maftools)


# data path ---------------------------------------------------------------

gscalite_path <- "/home/huff/data/GSCALite/TCGA"
gsca_v2_path <- "/home/huff/data/GSCA"

# load data ---------------------------------------------------------------

syn7824274.mc3.public.maf <- readr::read_tsv(file.path("/home/huff/data/GSCALite/TCGA/snv/syn-mutation-syn7824274-mc3-public.maf"))



# maf data explore --------------------------------------------------------

syn7824274.mc3.public.maf %>%
  ff_glimpse(levels_cut = 20) -> syn7824274.mc3.public.maf_summary


syn7824274.mc3.public.maf$Variant_Classification %>% table()
syn7824274.mc3.public.maf$Variant_Type %>% table()
syn7824274.mc3.public.maf %>% dplyr::select(Variant_Classification ,Variant_Type) %>% table()

# add cancer type  --------------------------------------------------------


clinical_info <- readr::read_rds(file.path("/home/huff/data/GSCALite/TCGA/clinical/pancan34_clinical.rds.gz"))
clinical_info %>%
  dplyr::mutate(clinical = purrr::map(clinical,.f=function(.x){
    .x %>%
      dplyr::select(barcode)
  })) %>%
  tidyr::unnest(clinical) -> clinical_samples


sample_info <- readr::read_rds(file.path("/home/huff/data/GSCA/sample_info.rds.gz")) %>%
  dplyr::mutate(barcode16=substr(barcode,1,16)) %>%
  dplyr::select(barcode16,cancer_types) %>%
  unique()

syn7824274.mc3.public.maf %>%
  dplyr::mutate(barcode16=substr(Tumor_Sample_Barcode,1,16)) %>%
  dplyr::inner_join(sample_info,by="barcode16") -> syn7824274.mc3.public.maf.cancer_type

syn7824274.mc3.public.maf.cancer_type$Tumor_Sample_Barcode %>% unique() %>% length() # 10234

# id trans ----------------------------------------------------------------

ncbi_name<- readr::read_rds(file = '/home/huff/github/GSCA/data/rda/ncbi_name.rds.gz')
load(file = file.path("/home/huff/github/GSCA/data/","rda",'01-gene-symbols.rda'))

# snv symbol
snv_symbol <- syn7824274.mc3.public.maf.cancer_type %>% 
  dplyr::select(Hugo_Symbol) %>%
  unique() %>%
  dplyr::rename(snv_symbol="Hugo_Symbol")

search_symbol %>%
  dplyr::filter(symbol %in% snv_symbol$snv_symbol) -> snv_symbol_in_search_symbol # 17994

snv_symbol %>%
  dplyr::filter(!snv_symbol %in% search_symbol$symbol) -> snv_symbol_NOT_in_search_symbol # 3335

# find snv symbol not in search symbol from ncbi

snv_symbol_NOT_in_search_symbol %>%
  dplyr::mutate(new_data = purrr::map(.x = snv_symbol, .f = function(.x) {
    .x_lower <- stringr::str_to_lower(string = .x)
    ncbi_name %>%
      dplyr::filter(purrr::map_lgl(.x = Synonyms, .f = function(.x) {
        .x_lower %in% .x
      }))
  })) ->  snv_symbol_NOT_in_search_symbol_IN_ncbi # 3335

snv_symbol_NOT_in_search_symbol_IN_ncbi %>% 
  dplyr::filter(purrr::map_lgl(.x = new_data, .f = function(.x){nrow(.x) != 0})) -> 
  snv_symbol_NOT_in_search_symbol_IN_ncbi_filter # 2393


snv_symbol_NOT_in_search_symbol_IN_ncbi_filter %>% 
  dplyr::filter(purrr::map_lgl(.x = new_data, .f = function(.x){nrow(.x) > 1})) %>%  # 74
  tidyr::unnest(cols = new_data) %>% 
  dplyr::filter(snv_symbol == symbol) ->
  snv_symbol_NOT_in_search_symbol_IN_ncbi_filter_2unique # 30


snv_symbol_NOT_in_search_symbol_IN_ncbi_filter %>% 
  dplyr::filter(purrr::map_lgl(.x = new_data, .f = function(.x){nrow(.x) == 1})) %>% # 1649
  tidyr::unnest(cols = new_data) %>% 
  dplyr::bind_rows(snv_symbol_NOT_in_search_symbol_IN_ncbi_filter_2unique) %>% 
  dplyr::select(-Synonyms) ->
  snv_symbol_NOT_in_search_symbol_IN_ncbi_filter_retain # 2340

snv_symbol_in_search_symbol %>% 
  tibble::add_column('snv_symbol' = snv_symbol_in_search_symbol$symbol, .before = 1) %>% 
  dplyr::bind_rows(snv_symbol_NOT_in_search_symbol_IN_ncbi_filter_retain) %>% 
  dplyr::distinct(entrez, symbol, .keep_all = T) ->
  snv_symbol_search_symbol_final

snv_symbol_search_symbol_final %>%
  readr::write_rds(file = '/home/huff/github/GSCA/data/rda/snv_symbol_search_symbol_final-new.rds.gz', compress = 'gz')


# maf files trans id ------------------------------------------------------
snv_symbol_search_symbol_final %>%
  dplyr::select(snv_symbol,entrez,symbol) -> snv_symbol_entrez

syn7824274.mc3.public.maf.cancer_type %>%
  dplyr::rename("snv_symbol"="Hugo_Symbol") %>%
  dplyr::inner_join(snv_symbol_entrez,by="snv_symbol") %>%
  dplyr::select(-snv_symbol,-Entrez_Gene_Id) %>%
  dplyr::mutate(barcode14 = substr(Tumor_Sample_Barcode,1,14)) %>%
  dplyr::rename(Hugo_Symbol=symbol) -> snv_data.IdTrans

snv_data.IdTrans %>% dplyr::select(cancer_types,Tumor_Sample_Barcode) %>% unique() -> snv_clinical.IdTrans
allmaf.IdTrans <- read.maf(maf=snv_data.IdTrans,
                             clinicalData = snv_clinical.IdTrans)
allmaf.IdTrans %>%
  readr::write_rds(file.path(gsca_v2_path,"snv","all_maf_data.IdTrans.maf.rds.gz"),compress="gz")

# cancer specific ---------------------------------------------------------
cancers_for_subset <- unique(snv_data.IdTrans$cancer_types)
for (cancer in cancers_for_subset) {
  snv_clinical.IdTrans %>%
    dplyr::filter(cancer_types %in% cancer)-> subset_clinical
  snv_data.IdTrans %>%
    dplyr::filter(cancer_types %in% cancer) -> subset_snv
  read.maf(maf = subset_snv,
           clinicalData = subset_clinical) %>%
    readr::write_rds(file.path(gsca_v2_path,"snv","sub_cancer_maf",paste(cancer,"maf_data.IdTrans.maf.rds.gz",sep = "_")))
}