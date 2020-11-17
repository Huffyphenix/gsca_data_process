##################################################################
#++++++++Tanslate symbol in TCGA expr into NCBI symbol +++++++++++


# GSCALite path -----------------------------------------------------------

gscalite_path <- "/home/huff/data/GSCALite"
gsca_v2_path <- "/home/huff/data/GSCA"

# read ncbi id  -----------------------------------------------------------

Homo_sapiens.gene_info.simplify <- readr::read_tsv(file.path("/home/huff/data/gene_info/homo_sapines/NCBI/Homo_sapiens.gene_info.simplify")) %>%
  dplyr::mutate(entrez_id=as.character(GeneID)) %>%
  dplyr::select(-GeneID)

Homo_sapiens.gene_info.simplify %>%
  dplyr::select(entrez_id,Symbol) -> Homo_sapiens.gene_info.simplify.symbolId

# read tcga expression file -----------------------------------------------

readr::read_rds(file.path(gscalite_path,"TCGA/expr/pancan33_expr.rds.gz")) -> pancan33_expr.rds.gz


# translate tcga expression file--------------------------------------------

pancan33_expr.rds.gz %>%
  dplyr::mutate(expr_trans = purrr::map(expr,.f=function(.x){
    Homo_sapiens.gene_info.simplify.symbolId %>%
      dplyr::right_join(.x,by="entrez_id") %>%
      dplyr::select(-symbol) %>%
      dplyr::rename("symbol"="Symbol")
  })) %>%
  dplyr::select(-expr) %>%
  dplyr::rename("expr"="expr_trans") -> pancan33_expr.rds.gz.IdTrans

pancan33_expr.rds.gz.IdTrans %>% readr::write_rds(file.path(gsca_v2_path,"pancan33_expr.IdTrans.rds.gz"))


# get corresponding id of TCGA in NCBI ------------------------------------

pancan33_expr.rds.gz$expr[[1]] %>%
  dplyr::select(entrez_id,symbol) %>%
  dplyr::left_join(Homo_sapiens.gene_info.simplify,by=c("entrez_id")) %>%
  dplyr::rename("TCGA_sym"="symbol","NCBI_sym"="Symbol") %>%
  dplyr::mutate(searchname_ncbi=tolower(NCBI_sym)) -> NCBI_id_in_TCGA

NCBI_id_in_TCGA %>%
  dplyr::filter(!is.na(NCBI_sym)) %>%
  dplyr::select(entrez=entrez_id,symbol=NCBI_sym,description,biotype=type_of_gene,ensembl=Ensembl,searchname=searchname_ncbi) %>%
  readr::write_rds(file.path(gsca_v2_path,"/NCBI_id_in_TCGA-final.rds.gz"),compress = "gz")

NCBI_id_in_TCGA %>%
  readr::write_tsv(file.path(gsca_v2_path,"NCBI_id_in.with_TCGAsym.rds.gz"))
# read TCGA DE file ----------------------------------------------------------------

pancan14_expr_fc_pval.rds.gz <- readr::read_rds(file.path(gscalite_path,"TCGA/expr/pancan14_expr_fc_pval.rds.gz"))


# translate TCGA DE file --------------------------------------------------
NCBI_id_in_TCGA %>% 
  dplyr::select(entrez_id,TCGA_sym,NCBI_sym) %>%
  dplyr::inner_join(pancan14_expr_fc_pval.rds.gz %>%
                      dplyr::filter(symbol!="?") %>%  # delete illegal value
                      dplyr::rename("TCGA_sym"="symbol"),by="TCGA_sym") %>%
  dplyr::rename("symbol"="NCBI_sym","entrez"="entrez_id") %>%
  dplyr::select(-TCGA_sym) -> pancan14_expr_fc_pval.IdTrans

pancan14_expr_fc_pval.IdTrans %>% readr::write_rds(file.path(gsca_v2_path,"pancan14_expr_fc_pval.IdTrans.rds.gz"),compress = "gz")



# read TCGA expr_subtype file ---------------------------------------------

expr_subtype.rds.gz <- readr::read_rds(file.path(gscalite_path,"TCGA/expr/expr_subtype.rds.gz"))


# translate TCGA expr subtype ---------------------------------------------

NCBI_id_in_TCGA %>% 
  dplyr::select(entrez_id,TCGA_sym,NCBI_sym) %>%
  dplyr::inner_join(expr_subtype.rds.gz %>%
                      dplyr::filter(symbol!="?") %>%  # delete illegal value
                      dplyr::rename("TCGA_sym"="symbol"),by="TCGA_sym") %>%
  dplyr::rename("symbol"="NCBI_sym","entrez"="entrez_id") %>%
  dplyr::select(-TCGA_sym) -> expr_subtype.IdTrans

expr_subtype.IdTrans %>% readr::write_rds(file.path(gsca_v2_path,"expr_subtype.IdTrans.rds.gz"),compress = "gz")


# read TCGA gene expr with survival ---------------------------------------

expr_survival.rds.gz <- readr::read_rds(file.path(gscalite_path,"TCGA/expr","expr_survival.rds.gz"))


# tanslate TCGA gene expr with survival -----------------------------------

NCBI_id_in_TCGA %>% 
  dplyr::select(entrez_id,TCGA_sym,NCBI_sym) %>%
  dplyr::inner_join(expr_survival.rds.gz %>%
                      dplyr::filter(symbol!="?") %>%  # delete illegal value
                      dplyr::rename("TCGA_sym"="symbol"),by="TCGA_sym") %>%
  dplyr::rename("symbol"="NCBI_sym","entrez"="entrez_id") %>%
  dplyr::select(-TCGA_sym) -> expr_survival.IdTrans

expr_survival.IdTrans %>% readr::write_rds(file.path(gsca_v2_path,"expr_survival.IdTrans.rds.gz"),compress = "gz")

# read GTEx mean expression file -----------------------------------------------

gtex_gene_mean_exp.rds.gz <- readr::read_rds(file.path(gscalite_path,"GTEx/expression","gtex_gene_mean_exp.rds.gz"))

gtex_gene_tmp_annotation_phenotype_v7.rds.gz <- readr::read_rds(file.path(gscalite_path,"GTEx/expression","gtex_gene_tmp_annotation_phenotype_v7.rds.gz"))

# translate GTEx mean expr ------------------------------------------------

gtex_gene_mean_exp.rds.gz %>%
  dplyr::mutate(Mean_trans = purrr::map(Mean,.f=function(.x){
    NCBI_id_in_TCGA %>%
      dplyr::select(symbol=NCBI_sym) %>%
      dplyr::filter(!is.na(symbol)) %>%
      dplyr::inner_join(.x,by="symbol")})
    ) %>%
  dplyr::select(-Mean) -> gtex_gene_mean_exp.IdTrans

gtex_gene_mean_exp.IdTrans %>% readr::write_rds(file.path(gsca_v2_path,"gtex_gene_mean_exp.IdTrans.rds.gz"),compress = "gz")


# read GTEx eqtl ----------------------------------------------------------

GTEx_egene.merged.tissue.rds.gz <- readr::read_rds(file.path(gscalite_path,"GTEx/eqtl","GTEx_egene.merged.tissue.rds.gz"))      

# translate GTEx mean expr ------------------------------------------------

NCBI_id_in_TCGA %>%
      dplyr::select(symbol=NCBI_sym) %>%
      dplyr::filter(!is.na(symbol)) %>%
      dplyr::inner_join(GTEx_egene.merged.tissue.rds.gz %>%
                          dplyr::rename("symbol"="gene_name"),by="symbol") -> GTEx_egene.merged.tissue.IdTrans

gtex_gene_mean_exp.IdTrans %>% readr::write_rds(file.path(gsca_v2_path,"GTEx_egene.merged.tissue.IdTrans.rds.gz"),compress = "gz")

