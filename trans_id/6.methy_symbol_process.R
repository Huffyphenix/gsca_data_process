
# symbol process of methylation data ----------------------------------------------

# Library -----------------------------------------------------------------

library(magrittr)

# Load data ---------------------------------------------------------------
rda_path <- "/home/huff/github/gsca_data_process/trans_id/rda"
gscalite_path <- "/home/huff/data/GSCALite/TCGA"
gsca_v2_path <- "/home/huff/data/GSCA"

load(file = file.path("/home/huff/github/GSCA/data/","rda",'01-gene-symbols.rda'))

ncbi_name<- readr::read_rds(file = '/home/huff/github/GSCA/data/rda/ncbi_name.rds.gz')


# load methy data -----------------------------------------------------------

methy_data <- readr::read_rds(file.path(gscalite_path,"meth", "pancan33_meth.rds.gz"))

# methy symbol --------------------------------------------------------------
methy_symbol <- methy_data %>% 
  tidyr::unnest() %>%
  dplyr::select(methy_symbol=symbol) %>%
  unique() 


search_symbol %>%
  dplyr::filter(symbol %in% methy_symbol$methy_symbol) -> methy_symbol_in_search_symbol # 16,887 

methy_symbol %>%
  dplyr::filter(!methy_symbol %in% search_symbol$symbol) -> methy_symbol_NOT_in_search_symbol # 654

# find methy symbol not in search symbol from ncbi

methy_symbol_NOT_in_search_symbol %>%
  dplyr::mutate(new_data = purrr::map(.x = methy_symbol, .f = function(.x) {
    .x_lower <- stringr::str_to_lower(string = .x)
    ncbi_name %>% 
      dplyr::filter(purrr::map_lgl(.x = Synonyms, .f = function(.x) {
        .x_lower %in% .x
      }))
  })) ->  methy_symbol_NOT_in_search_symbol_IN_ncbi # 654

methy_symbol_NOT_in_search_symbol_IN_ncbi %>% 
  dplyr::filter(purrr::map_lgl(.x = new_data, .f = function(.x){nrow(.x) != 0})) -> 
  methy_symbol_NOT_in_search_symbol_IN_ncbi_filter # 592


methy_symbol_NOT_in_search_symbol_IN_ncbi_filter %>% 
  dplyr::filter(purrr::map_lgl(.x = new_data, .f = function(.x){nrow(.x) > 1})) %>%  # 74
  tidyr::unnest(cols = new_data) %>% 
  dplyr::filter(methy_symbol == symbol) ->
  methy_symbol_NOT_in_search_symbol_IN_ncbi_filter_2unique # 0


methy_symbol_NOT_in_search_symbol_IN_ncbi_filter %>% 
  dplyr::filter(purrr::map_lgl(.x = new_data, .f = function(.x){nrow(.x) == 1})) %>% # 568
  tidyr::unnest(cols = new_data) %>% 
  dplyr::bind_rows(methy_symbol_NOT_in_search_symbol_IN_ncbi_filter_2unique) %>% 
  dplyr::select(-Synonyms) ->
  methy_symbol_NOT_in_search_symbol_IN_ncbi_filter_retain # 568

methy_symbol_in_search_symbol %>% 
  tibble::add_column('methy_symbol' = methy_symbol_in_search_symbol$symbol, .before = 1) %>% 
  dplyr::bind_rows(methy_symbol_NOT_in_search_symbol_IN_ncbi_filter_retain) %>% 
  dplyr::distinct(entrez, symbol, .keep_all = T) ->
  methy_symbol_search_symbol_final # 17456

methy_symbol_search_symbol_final %>%
  readr::write_rds(file = file.path(rda_path,'methy_symbol_search_symbol_final.rds.gz'), compress = 'gz')
# save image --------------------------------------------------------------


save.image(file.path(rda_path,"6.methy_symbol_process.rda"))
# load(file.path(rda_path,"4.methy_symbol_process.rda"))
