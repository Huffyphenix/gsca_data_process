##################################################################
#++++++++Tanslate symbol in TCGA expr into NCBI symbol +++++++++++


# read ncbi id  -----------------------------------------------------------

Homo_sapiens.gene_info.simplify <- readr::read_tsv(file.path("/home/huff/data/gene_info/homo_sapines/NCBI/Homo_sapiens.gene_info.simplify")) %>%
  dplyr::mutate(GeneID=as.character(GeneID))

Homo_sapiens.gene_info.simplify %>%
  dplyr::select(GeneID,Symbol) %>%
  dplyr::rename("entrez_id"="GeneID") -> Homo_sapiens.gene_info.simplify.symbolId

# read tcga expression file -----------------------------------------------

readr::read_rds(file.path("/home/liucj/shiny-data/GSCALite/TCGA/expr/pancan33_expr.rds.gz")) -> pancan33_expr.rds.gz


# translate ---------------------------------------------------------------

pancan33_expr.rds.gz %>%
  dplyr::mutate(expr_trans = purrr::map(expr,.f=function(.x){
    Homo_sapiens.gene_info.simplify.symbolId %>%
      dplyr::right_join(.x,by="entrez_id") %>%
      dplyr::select(-symbol) %>%
      dplyr::rename("symbol"="Symbol")
  })) %>%
  dplyr::select(-expr) %>%
  dplyr::rename("expr"="expr_trans") -> pancan33_expr.rds.gz.IdTrans

pancan33_expr.rds.gz.IdTrans %>% readr::write_rds(file.path("/home/huff/data/GSCA/pancan33_expr.rds.IdTrans.gz"))
