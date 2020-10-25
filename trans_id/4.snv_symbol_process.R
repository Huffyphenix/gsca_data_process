
# symbol process of snv data ----------------------------------------------

# Library -----------------------------------------------------------------

library(magrittr)

# Load data ---------------------------------------------------------------
rda_path <- "/home/huff/github/gsca_data_process/trans_id/rda/"
gscalite_path <- "/home/huff/data/GSCALite/TCGA"
gsca_v2_path <- "/home/huff/data/GSCA"

load(file = file.path("/home/huff/github/GSCA/data/","rda",'01-gene-symbols.rda'))

ncbi_name<- readr::read_rds(file = '/home/huff/github/GSCA/data/rda/ncbi_name.rds.gz')


# load snv data -----------------------------------------------------------

mc3_pass <- readr::read_rds(file.path(gscalite_path,"snv", "01-snv_mutation_mc3_public.pass.filtered_maf.rds.gz"))
# "gene_list_snv_count.new.new200215.rds.gz"
# "01-snv_mutation_mc3_public.pass.filtered_maf.rds.gz"
# "pancan32_snv_survival_genelist_sig_pval.rds.gz"

# snv symbol --------------------------------------------------------------
snv_symbol <- mc3_pass %>% 
  maftools::getGeneSummary() %>% 
  .$Hugo_Symbol %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::rename(snv_symbol=".")
  

search_symbol %>%
  dplyr::filter(symbol %in% snv_symbol$snv_symbol) -> snv_symbol_in_search_symbol # 17567

snv_symbol %>%
  dplyr::filter(!snv_symbol %in% search_symbol$symbol) -> snv_symbol_NOT_in_search_symbol # 2063

# find snv symbol not in search symbol from ncbi

snv_symbol_NOT_in_search_symbol %>%
  dplyr::mutate(new_data = purrr::map(.x = snv_symbol, .f = function(.x) {
    .x_lower <- stringr::str_to_lower(string = .x)
    ncbi_name %>% 
      dplyr::filter(purrr::map_lgl(.x = Synonyms, .f = function(.x) {
        .x_lower %in% .x
      }))
  })) ->  snv_symbol_NOT_in_search_symbol_IN_ncbi # 2063

snv_symbol_NOT_in_search_symbol_IN_ncbi %>% 
  dplyr::filter(purrr::map_lgl(.x = new_data, .f = function(.x){nrow(.x) != 0})) -> 
  snv_symbol_NOT_in_search_symbol_IN_ncbi_filter # 1723


snv_symbol_NOT_in_search_symbol_IN_ncbi_filter %>% 
  dplyr::filter(purrr::map_lgl(.x = new_data, .f = function(.x){nrow(.x) > 1})) %>%  # 74
  tidyr::unnest(cols = new_data) %>% 
  dplyr::filter(snv_symbol == symbol) ->
  snv_symbol_NOT_in_search_symbol_IN_ncbi_filter_2unique # 25


snv_symbol_NOT_in_search_symbol_IN_ncbi_filter %>% 
  dplyr::filter(purrr::map_lgl(.x = new_data, .f = function(.x){nrow(.x) == 1})) %>% # 1649
  tidyr::unnest(cols = new_data) %>% 
  dplyr::bind_rows(snv_symbol_NOT_in_search_symbol_IN_ncbi_filter_2unique) %>% 
  dplyr::select(-Synonyms) ->
  snv_symbol_NOT_in_search_symbol_IN_ncbi_filter_retain # 1674

snv_symbol_in_search_symbol %>% 
  tibble::add_column('snvsymbol' = snv_symbol_in_search_symbol$symbol, .before = 1) %>% 
  dplyr::bind_rows(snv_symbol_NOT_in_search_symbol_IN_ncbi_filter_retain) %>% 
  dplyr::distinct(entrez, symbol, .keep_all = T) ->
  snv_symbol_search_symbol_final

snv_symbol_search_symbol_final %>%
  readr::write_rds(file = 'trans_id/rda/snv_symbol_search_symbol_final.rds.gz', compress = 'gz')
snv_symbol_search_symbol_final %>%
  readr::write_rds(file = '/home/huff/github/GSCA/data/rda/snv_symbol_search_symbol_final.rds.gz', compress = 'gz')
# save image --------------------------------------------------------------


save.image(file.path(rda_path,"4.snv_symbol_process.rda"))
