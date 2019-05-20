##########################
# id corresponding between TCGA symbol, GETx symbol by using entrez id

library(magrittr)
# data path ---------------------------------------------------------------

basic_path <- "/home/huff/project"

ncbi_symbol_info_9606 <- readr::read_tsv(file.path(basic_path,"data/id_corresponding/Homo_NCBIgeneid2symbol","Homo_sapiens.gene_info.9606.no_blank"))
TCGA_RNAseq_gene_info <- readr::read_tsv(file.path(basic_path,"TCGA_lungcancer_data/RNAseq/all_sample_exp","rna_exp.confirm")) %>%
  dplyr::select(gene_id) %>%
  dplyr::mutate(n=1:20531) %>%
  tidyr::separate(gene_id,c("Symbol","GeneID"),sep="\\|")

GTEx_gene_info <- readr::read_rds(file.path(basic_path,"data/GSCALite", "GTEx", "expression", "gtex_gene_mean_exp.rds.gz")) %>%
  head(1) %>%
  dplyr::select(Mean) %>%
  tidyr::unnest() %>%
  dplyr::select(symbol) %>%
  dplyr::rename("Symbol"="symbol")

# id process --------------------------------------------------------------

fn_get_synonyms <- function(.x){
  .x %>% 
    stringr::str_split(pattern = "\\|", simplify = TRUE) %>%
    .[1, ] %>%
    stringr::str_trim(side = "both") -> .ss
  tibble::tibble(alias = .ss)
}

ncbi_symbol_info_9606 %>%
  dplyr::select(GeneID,Symbol,Synonyms) %>%
  tidyr::nest(-GeneID,-Symbol) %>%
  dplyr::mutate(alias = purrr::map(data,fn_get_synonyms)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ncbi_symbol_info_9606.alias
  
ncbi_symbol_info_9606.alias %>%
  dplyr::mutate(alias = ifelse(alias == "-",Symbol,alias)) -> ncbi_symbol_info_9606.alias

# exist replecates ? YES------------------------------------------------------

ncbi_symbol_info_9606.alias %>%
  dplyr::group_by(alias) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::arrange(desc(n))


# merge ncbi gene info and TCGa gene info ---------------------------------

TCGA_RNAseq_gene_info %>%
  dplyr::mutate(GeneID=as.numeric(GeneID)) %>%
  dplyr::inner_join(ncbi_symbol_info_9606.alias,by="GeneID") %>%
  dplyr::rename("TCGA_sym"="Symbol.x","NCBI_sym"="Symbol.y") %>%
  dplyr::select(-n) %>%
  dplyr::mutate(TCGA_sym = ifelse(TCGA_sym == "?", NCBI_sym,TCGA_sym)) -> gene_info.TCGA.NCBI

# merge ncbi gene info and GTEX gene info ---------------------------------
GTEx_gene_info %>%
  dplyr::inner_join(ncbi_symbol_info_9606.alias,by = "Symbol") %>%
  dplyr::mutate(NCBI_sym = Symbol) %>%
  dplyr::rename("GTEX_sym" = "Symbol") -> gene_info.TCGA.GTEXsym

GTEx_gene_info %>%
  dplyr::filter(!Symbol %in% gene_info.TCGA.GTEXsym$GTEX_sym) -> GTEx_gene.alias
ncbi_symbol_info_9606.alias %>%
  dplyr::filter(alias %in% GTEx_gene.alias$Symbol) %>%
  dplyr::rename("NCBI_sym" = "Symbol") %>%
  dplyr::mutate(GTEX_sym = alias) %>%
  dplyr::select(GTEX_sym,GeneID, alias, NCBI_sym)-> gene_info.TCGA.GTEXalias
gene_info.TCGA.GTEXsym %>%
  rbind(gene_info.TCGA.GTEXalias) -> gene_info.GTEX.TCGA

gene_info.GTEX.TCGA %>%
  dplyr::full_join(gene_info.TCGA.NCBI,by=c("GeneID","NCBI_sym","alias")) -> gene_info.merge

gene_info.merge %>%
  readr::write_rds(file.path(basic_path,"data/id_corresponding","id_correspond_between_NCBI_TCGA.rds.gz"),compress="gz")
