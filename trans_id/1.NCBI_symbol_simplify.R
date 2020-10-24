
# ncbi gene info simplify -------------------------------------------------

readr::read_tsv(file.path("/home/huff/data/gene_info/homo_sapines/NCBI/Homo_sapiens.gene_info")) %>%
  dplyr::select(GeneID,Symbol,dbXrefs,Synonyms,description,type_of_gene,Full_name_from_nomenclature_authority,Modification_date) %>%
  tidyr::separate(dbXrefs,into = c("other","Ensembl"),sep = "Ensembl:") %>%
  dplyr::select(-other) %>%
  dplyr::mutate(Ensembl=purrr::map(Ensembl,.f=function(.x){
    stringr::str_trim(strsplit(.x,"\\|")[[1]][1],side = "both")
  })) %>%
  tidyr::unnest() -> ncbi_symbol_info_9606

ncbi_symbol_info_9606 %>% 
  readr::write_tsv(file.path("/home/huff/data/gene_info/homo_sapines/NCBI/Homo_sapiens.gene_info.simplify"))
