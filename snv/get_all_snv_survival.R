# snv with survival
library(magrittr)
library(ggplot2)

out_path <- "/data/GSCALite/TCGA"
snv_path <- file.path(out_path, "snv")
tcga_path <- "/project/huff/huff/immune_checkpoint/data/TCGA_data"



# load cnv and gene list
snv <- readr::read_rds(file.path(tcga_path, "pancan33_snv.rds.gz"))


# load clinical data
clinical <- readr::read_rds(path = file.path(tcga_path, "pancan34_clinical.rds.gz"))

# merge clinical and snv
snv_clinical <-
  snv %>%
  # dplyr::filter(cancer_types %in% cancer_pairs$cancer_types) %>%
  dplyr::inner_join(clinical, by = "cancer_types")

fun_barcode <- function(.b) {
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
} # get short barcode from long barcode
fun_tn_type <- function(.b) {
  type <- .b %>%
    stringr::str_split(pattern = "-", simplify = T) %>%
    .[, 4] %>%
    stringr::str_sub(1, 2)
} # get tumor and normal info

fun_expr_survival_merge <- function(filter_snv, clinical) {
  # merge clinical and expr data
  filter_snv %>%
    # dplyr::select(-entrez_id) %>%
    tidyr::gather(key = barcode, value = snv, -symbol) %>%
    dplyr::mutate(barcode = fun_barcode(barcode)) %>%
    dplyr::select(symbol, barcode, snv) %>%
    dplyr::mutate(snv = as.character(snv)) %>%
    # dplyr::mutate(snv=ifelse(is.na(snv),"non mut",snv)) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(snv = plyr::revalue(snv, replace = c("0" = "non mut"))) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(snv_class = ifelse(snv != "non mut", "mut", snv)) %>%
    dplyr::select(symbol, barcode, snv_class) %>%
    dplyr::inner_join(clinical, by = "barcode") %>%
    dplyr::select(symbol, barcode, snv_class,time = os_days, status = os_status) %>%
    dplyr::filter(!is.na(time), time > 0, !is.na(status)) %>%
    dplyr::mutate(status = plyr::revalue(status, replace = c("Alive" = 0, "Dead" = 1))) %>%
    dplyr::mutate(status = as.numeric(status)) %>%
    # dplyr::mutate(expr = log2(expr + 1)) %>%
    # tidyr::drop_na(expr) %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(snv_class = as.character(snv_class)) %>%
    dplyr::mutate(group = snv_class) %>%
    dplyr::ungroup() -> expr_clinical_ready
}

fun_draw_survival <- function(symbol, coxP, cancer_types, expr_clinical_ready) {
  # symbol="CD96"
  # p.value=0.05
  # cancer_types="SKCM"
  gene <- symbol
  p_val <- signif(-log10(coxP), digits = 3)
  fig_name <- paste(cancer_types, gene, "png", sep = ".")
  # print(fig_name)
  .d <-
    expr_clinical_ready %>%
    dplyr::filter(symbol == gene)

  .d_diff <- survival::survdiff(survival::Surv(time, status) ~ group, data = .d)

  kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)

  if (kmp > 0.05) {
    return(NA)
  } else {
    fit_x <- survival::survfit(survival::Surv(time, status) ~ group, data = .d, na.action = na.exclude)
    # pdf(file = file.path(snv_path,"snv_survival",fig_name),width = 480,height = 480)
    survminer::ggsurvplot(
      fit_x, data = .d, pval = T, pval.method = T,
      tables.height = 0.2,
      # tables.theme = theme_cleantable(),
      title = paste(paste(cancer_types, gene, sep = "-"), "Cox P =", signif(coxP, 2), ", Logrank P = ", signif(kmp, 2)),
      xlab = "Survival in days",
      ylab = "Probability of survival"
    )
    # dev.off()
    ggsave(filename = fig_name, device = "png", path = file.path("/data/GSCALite/TCGA", "snv", "snv_survival"), width = 6, height = 6)
  }
}
fun_clinical_test <- function(expr_clinical_ready, cancer_types) {
  if (nrow(expr_clinical_ready) < 1) {
    return(tibble::tibble())
  }
  # if(expr_clinical_ready %>% dplyr::filter(group=="mut") %>% nrow() <nrow(expr_clinical_ready)*0.02){return(tibble::tibble())}

  print(cancer_types)
  expr_clinical_ready %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarise(sm_mut = sum(snv_class == "mut"), sm_nm = sum(snv_class == "non mut")) %>%
    dplyr::mutate(mut_per = sm_mut / (sm_mut + sm_nm))  -> gene_mut_per
  if (gene_mut_per$sm_mut[1]+gene_mut_per$sm_nm[1] > 50) {
    gene_mut_per %>%
      dplyr::filter(mut_per >= 0.05) %>% # num of samle >50,only do survival analysis for genes who have mutated in more than 10% samples
      dplyr::pull(symbol) -> ready_gene
  } else {
    gene_mut_per %>%
      dplyr::filter(sm_mut >= 3) %>% # num of samle <50,only do survival analysis for genes who have mutated in more than 3 samples
      dplyr::pull(symbol) -> ready_gene
  }
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
      dplyr::mutate(worse = ifelse(estimate > 0, "Low", "High")) -> d
    
    # d %>%
    #   dplyr::select(symbol, coxP) %>%
    #   purrr::pwalk(fun_draw_survival, cancer_types = cancer_types, expr_clinical_ready = expr_clinical_ready)

    return(d)
  }
}

cl <- 15
cluster <- multidplyr::create_cluster(core = cl)
snv_clinical %>%
  # dplyr::filter(cancer_types=="LUAD") %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("ggplot2") %>%
  multidplyr::cluster_assign_value("fun_clinical_test", fun_clinical_test) %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode) %>%
  # multidplyr::cluster_assign_value("fun_draw_survival", fun_draw_survival) %>%
  multidplyr::cluster_assign_value("fun_expr_survival_merge", fun_expr_survival_merge) %>%
  dplyr::mutate(merged_clean = purrr::map2(snv, clinical, fun_expr_survival_merge)) %>%
  dplyr::select(-snv, -clinical) %>%
  dplyr::mutate(diff_pval = purrr::map2(merged_clean, cancer_types, fun_clinical_test)) %>%
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-merged_clean) %>%
  dplyr::select(-PARTITION_ID)  -> expr_clinical_sig_pval
parallel::stopCluster(cluster)

expr_clinical_sig_pval %>%
  readr::write_rds(file.path("/data/GSCALite/TCGA/snv", "pancan32_snv_survival_genelist_sig_pval.rds.gz"), compress = "gz")
