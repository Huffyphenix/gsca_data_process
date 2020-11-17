
# check overlap between old and new genes in cnv files --------------------


library(magrittr)
cnv_symbol_final <- readr::read_rds("/home/huff/github/GSCA/data/rda/cnv_symbol_search_symbol_final.rds.gz")

cnv_symbol_before <- readr::read_rds("/home/huff/data/GSCALite/TCGA/cnv/pancan34_cnv.rds.gz")$cnv[[1]]$symbol
 
expr_symbol_final <- readr::read_tsv("/home/huff/data/GSCA/id/NCBI_id_in.with_TCGAsym.tsv") %>%
  dplyr::filter(!is.na(NCBI_sym))

expr_symbol_before <- readr::read_rds("/home/huff/data/GSCALite/TCGA/expr/pancan33_expr.rds.gz")$expr[[1]]$symbol

cnv_symbol_before %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::filter(. %in% expr_symbol_before) %>%
  dplyr::rename("cnvsymbol"=".") -> cnv_expr_before #  17951

cnv_symbol_final %>%
  dplyr::filter(symbol %in% expr_symbol_final$NCBI_sym) -> cnv_expr_final # 19122

cnv_expr_final %>%
  dplyr::inner_join(cnv_expr_before,by="cnvsymbol")  -> final_before_overlap # 17843

cnv_expr_final %>%
  dplyr::filter(! cnvsymbol %in% cnv_expr_before$cnvsymbol) ->  final_not_in_before# 1279, genes should be added in present cor results

cnv_expr_before %>%
  dplyr::filter(! cnvsymbol %in% cnv_expr_final$cnvsymbol) -> before_not_in_final # 108, genes should be delete in previous results

load("/home/huff/github/GSCA/data/rda/07-cnv-symbol.rda")
before_not_in_final %>%
  dplyr::filter(cnvsymbol %in% cnv_symbol_no_search_symbol_ncbi_filter$cnvsymbol) -> before_not_in_final_symbolAliasConfused # 100, symbol alias confused genes

before_not_in_final %>%
  dplyr::filter(! cnvsymbol %in% cnv_symbol_no_search_symbol_ncbi_filter$cnvsymbol) -> before_not_in_final_Missing # 8, not human genes


before_not_in_final -> genes_should_delete_in_previous_res
final_not_in_before -> genes_should_add_in_present_res

save.image("cnv_new/rda/1.cnv_symbol_process.rda")
