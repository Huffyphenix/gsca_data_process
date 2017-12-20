
# methy survival p value filter -------------------------------------------

meth_survival <- readr::read_rds(file.path(config$database, "TCGA", "meth", "pancan32_meth_survival_genelist_sig_pval.rds.gz"))
meth_survival %>%
  tidyr::unnest() %>%
  dplyr::filter(logRankP<=0.05) %>%
  dplyr::mutate(log10logrankP=-log10(logRankP)) %>%
  dplyr::select(cancer_types,symbol,Hyper_worse,log10logrankP) %>%
  tidyr::nest(symbol,Hyper_worse,log10logrankP,.key="diff_pval") -> meth_survival_p0.05

meth_survival_p0.05 %>%
  readr::write_rds("/data/GSCALite/TCGA/meth/pancan32_meth_survival_genelist_sig_pval0.05.rds.gz",compress = "gz")

# methy diff cancer filter ------------------------------------------------

meth_diff <- readr::read_rds(file.path(config$database, "TCGA", "meth", "pan33_allgene_methy_diff.rds.gz"))

c("CHOL","OV","ACC","THYM","UCS","SKCM","GBM","LGG","LAML","UVM","MESO","DLBC","TGCT","SARC","READ","STAD","CESC","KICH","PCPG") ->cancer.no.diff

meth_diff %>%
dplyr::filter(! cancer_types %in% cancer.no.diff) %>%
  tidyr::unnest() %>%
  dplyr::select(-gene,-p_val) %>%
  tidyr::nest(symbol,diff,direction,fdr,.key=methy_comparison)->meth_diff.simplification

meth_diff.simplification %>%
  readr::write_rds("/data/GSCALite/TCGA/meth/pan33_allgene_methy_diff.simplification.rds.gz",compress = "gz")
