#!/usr/bin/Rscript
library(magrittr)
# load meth data ----------------------------------------------------------
out_path <- "/data/GSCALite/TCGA"


# load meth data ----------------------------------------------------------
methy <- readr::read_rds(file.path("/data/TCGA/TCGA_data/pancan33_meth.rds.gz"))

# load clinical data
clinical <- readr::read_rds(path = file.path("/data/TCGA/TCGA_data/", "pancan34_clinical.rds.gz"))

# merge clinical and snv
meth_clinical <-
  methy %>%
  # dplyr::filter(cancer_types %in% cancer_pairs$cancer_types) %>%
  dplyr::inner_join(clinical, by = "cancer_types")

# functions ---------------------------------------------------------------

fun_tn_type <- function(.b) {
  type <- .b %>%
    stringr::str_split(pattern = "-", simplify = T) %>%
    .[, 4] %>%
    stringr::str_sub(1, 2)
} # get tumor and normal info

fun_barcode <- function(.b) {
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
} # get short barcode from long barcode


fun_expr_survival_merge <- function(meth, clinical) {
  # merge clinical and expr data
  meth %>%
    # head() %>%
    dplyr::select(-gene) %>%         
    tidyr::gather(key = barcode, value = meth, -symbol) %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarise(mid=median(meth)) ->mid_mean
  meth %>%
    # head() %>%
    dplyr::select(-gene) %>%
    tidyr::gather(key = barcode, value = meth, -symbol) %>%
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type != "11") %>% 
    dplyr::select(-type) %>% 
    dplyr::mutate(barcode = fun_barcode(barcode)) %>%
    dplyr::select(symbol, barcode, meth) %>%
    dplyr::left_join(mid_mean,by="symbol") %>%
    dplyr::mutate(group=ifelse(meth>mid,"Up","Down")) %>%
    dplyr::select(symbol,barcode,group) %>%
    # dplyr::mutate(meth = as.character(meth)) %>%
    # dplyr::mutate(snv=ifelse(is.na(snv),"non mut",snv)) %>%
    tidyr::drop_na() %>%
    dplyr::inner_join(clinical, by = "barcode") %>%
    dplyr::select(symbol, barcode, group, time = os_days, status = os_status) %>%
    dplyr::filter(!is.na(time), time > 0, !is.na(status)) %>%
    dplyr::mutate(status = plyr::revalue(status, replace = c("Alive" = 0, "Dead" = 1))) %>%
    dplyr::mutate(status = as.numeric(status)) %>%
    # dplyr::mutate(expr = log2(expr + 1)) %>%
    # tidyr::drop_na(expr) %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(group = as.character(group)) %>%
    dplyr::ungroup() -> expr_clinical_ready
}

fun_clinical_test <- function(expr_clinical_ready, cancer_types) {
  if (nrow(expr_clinical_ready) < 1) {
    return(tibble::tibble())
  }
  # if(expr_clinical_ready %>% dplyr::filter(group=="mut") %>% nrow() <nrow(expr_clinical_ready)*0.02){return(tibble::tibble())}
  
  print(cancer_types)
  # print(1)
  expr_clinical_ready %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarise(up = sum(group == "Up"), down = sum(group == "Down")) %>%
    dplyr::filter_all(.vars_predicate = dplyr::all_vars(.>10)) %>%
    dplyr::pull(symbol) -> ready_gene
  
  if (ready_gene %>% length() == 0) {
    return(tibble::tibble())
  } else {
    expr_clinical_ready %>%
      dplyr::filter(symbol %in% ready_gene) %>%
      dplyr::group_by(symbol) %>%
      dplyr::do(
        cox = broom::tidy(
          tryCatch(
            survival::coxph(survival::Surv(time, status) ~ group, data = ., na.action = na.exclude),
            error = function(e) {
              1
            }
          )
        ),
        logrank = broom::tidy(
          tryCatch(
            1 - pchisq(survival::survdiff(survival::Surv(time, status) ~ group, data = ., na.action = na.exclude)$chisq, df = length(levels(as.factor(.$group))) - 1),
            error = function(e) {
              1
            }
          )
        )
      ) %>%
      dplyr::ungroup() %>%
      tidyr::unnest() %>%
      # dplyr::filter(p.value < 0.05) %>%
      dplyr::mutate(logRankP = x) %>%
      dplyr::mutate(coxP = p.value) %>%
      dplyr::select(symbol, estimate, coxP, logRankP) %>%
      dplyr::mutate(Hyper_worse = ifelse(estimate > 0,  "High","Low")) -> d
    
    # d %>%
    #   dplyr::select(symbol, coxP) %>%
    #   purrr::pwalk(fun_draw_survival, cancer_types = cancer_types, expr_clinical_ready = expr_clinical_ready)
    
    return(d)
  }
}

cl <- 15
cluster <- multidplyr::create_cluster(core = cl)
meth_clinical %>%
  # dplyr::filter(cancer_types=="LUAD") %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  # multidplyr::cluster_library("ggplot2") %>%
  multidplyr::cluster_assign_value("fun_clinical_test", fun_clinical_test) %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type) %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode) %>%
  multidplyr::cluster_assign_value("fun_expr_survival_merge", fun_expr_survival_merge) %>%
  # multidplyr::cluster_assign_value("ready_gene", ready_gene) %>%
  dplyr::mutate(merged_clean = purrr::map2(methy, clinical, fun_expr_survival_merge)) %>%
  dplyr::select(-methy, -clinical) %>%
  dplyr::mutate(diff_pval = purrr::map2(merged_clean, cancer_types, fun_clinical_test)) %>%
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-merged_clean) %>%
  dplyr::select(-PARTITION_ID) -> expr_clinical_sig_pval
parallel::stopCluster(cluster)

expr_clinical_sig_pval %>%
  readr::write_rds(file.path("/data/GSCALite/TCGA/meth", "pancan32_meth_survival_genelist_sig_pval.rds.gz"), compress = "gz")
