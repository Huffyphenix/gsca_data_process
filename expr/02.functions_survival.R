library(magrittr)
library(survival)
library(dplyr)

# transform samples barcode -----------------------------------------------

fn_transform_samples <- function(.x) {
  .x %>% 
    tidyr::gather(key = 'aliquot', value = 'expr') ->
    .y
  
  .y %>% 
    dplyr::mutate(tmp = substr(x = aliquot, start = 14, stop = 14)) %>% 
    dplyr::mutate(type = ifelse(tmp == '0', 'tumor', 'normal')) %>% 
    dplyr::mutate(barcode = substr(x = aliquot, start = 1, stop = 16)) %>% 
    dplyr::mutate(sample_name = substr(x = barcode, start = 1, stop = 12)) 
}


# expr group for survival -------------------------------------------------

survival_group <- tibble::tibble(type=c("os","pfs"),
                                 time=c("os_days","pfs_days"),
                                 status=c("os_status","pfs_status"))

fn_expr_group <- function(.data,.cutoff){
  .data %>%
    dplyr::mutate(group=ifelse(expr>quantile(expr,.cutoff)[[1]],"Higher expr.", "Lower expr."))
}

fn_cox_logp <- function(.d){
  .d %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(n=n()) %>%
    dplyr::select(group,n) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n>5) %>%
    .$group %>% unique() %>% length() -> len_group
  if(len_group==2){
    .d_diff <- survival::survdiff(survival::Surv(time, status) ~ group, data = .d, na.action = na.exclude)
    kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
    coxp_categorical <- broom::tidy(survival::coxph(survival::Surv(time, status) ~ group, data = .d, na.action = na.exclude))$p.value
  } else {
    kmp <- NA
    coxp_categorical <- NA
  }
  
  coxp_continus <- broom::tidy(survival::coxph(survival::Surv(time, status) ~ expr, data = .d, na.action = na.exclude))$p.value

  tibble::tibble(logrankp=kmp,coxp_categorical=coxp_categorical,coxp_continus=coxp_continus)
}

fn_survival <- function(.data,.cutoff,sur_type){
  .data %>%
    fn_expr_group(.cutoff = .cutoff) -> .data_grouped
  
  survival_group %>%
    dplyr::filter(type==sur_type) -> sur_type_do
  
  .data_grouped %>%
    dplyr::select(sample_name,group,time=sur_type_do$time,status=sur_type_do$status,expr) %>%
    fn_cox_logp()
}

fn_survival_res <- function(.cancer_types,.survival,.expr){
  print(.cancer_types)
  if (length(grep("pfs",colnames(.survival)))>0) {
    .survival %>%
      dplyr::rename("sample_name"="barcode") %>%
      dplyr::mutate(os_status=purrr::map(os_status,.f=function(.x){
        if(!is.na(.x)){
          ifelse(.x=="Dead",1,0)
        } else {
          NA
        }
      })) %>%
      dplyr::mutate(pfs_status =purrr::map(pfs_status,.f=function(.x){
        if(!is.na(.x)){
          ifelse(.x=="progression",1,0)
        } else {
          NA
        }
      })) %>%
      tidyr::unnest() -> .survival
  }else {
    .survival %>%
      dplyr::rename("sample_name"="barcode") %>%
      dplyr::mutate(os_status=purrr::map(os_status,.f=function(.x){
        if(!is.na(.x)){
          ifelse(.x=="Dead",1,0)
        } else {
          NA
        }
      })) %>%
      tidyr::unnest() -> .survival
  }
  
  .expr %>%
    dplyr::filter(symbol %in% gene_symbol$symbol) %>% 
    dplyr::rename(entrez = entrez_id) %>% 
    dplyr::mutate(entrez = as.numeric(entrez)) %>% 
    dplyr::group_by(entrez, symbol) %>% 
    tidyr::nest() -> .expr
  
  .expr %>%
    dplyr::mutate(data = purrr::map(data,.f=fn_transform_samples)) -> .expr 
  
  .expr %>%
    dplyr::mutate(data=purrr::map(data,.f=function(.x){
      .x %>%
        dplyr::filter(type=="tumor") %>%
        dplyr::filter(!is.na(expr)) %>%
        dplyr::inner_join(.survival,by="sample_name")
    })) -> .combine
  # os survival ----
  .combine %>%
    dplyr::mutate(cutoff=0.5,sur_type="os") %>%
    dplyr::mutate(survival_res = purrr::pmap(list(data,cutoff,sur_type),fn_survival)) %>%
    dplyr::select(-data) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(survival_res)-> os_median
  .combine %>%
    dplyr::mutate(cutoff=0.75,sur_type="os") %>%
    dplyr::mutate(survival_res = purrr::pmap(list(data,cutoff,sur_type),fn_survival)) %>%
    dplyr::select(-data) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(survival_res) -> os_upquantile
  .combine %>%
    dplyr::mutate(cutoff=0.25,sur_type="os") %>%
    dplyr::mutate(survival_res = purrr::pmap(list(data,cutoff,sur_type),fn_survival)) %>%
    dplyr::select(-data) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(survival_res) -> os_lowquantile
  # pfs survival -----
  if (length(grep("pfs",colnames(.survival)))>0) {
    .combine %>%
      dplyr::mutate(cutoff=0.5,sur_type="pfs") %>%
      dplyr::mutate(survival_res = purrr::pmap(list(data,cutoff,sur_type),fn_survival)) %>%
      dplyr::select(-data) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(survival_res) -> pfs_median
    .combine %>%
      dplyr::mutate(cutoff=0.75,sur_type="pfs") %>%
      dplyr::mutate(survival_res = purrr::pmap(list(data,cutoff,sur_type),fn_survival)) %>%
      dplyr::select(-data) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(survival_res) -> pfs_upquantile
    .combine %>%
      dplyr::mutate(cutoff=0.25,sur_type="pfs") %>%
      dplyr::mutate(survival_res = purrr::pmap(list(data,cutoff,sur_type),fn_survival)) %>%
      dplyr::select(-data) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(survival_res) -> pfs_lowquantile
  } else {
    .combine %>%
      dplyr::mutate(cutoff=0.5,sur_type="pfs") %>%
      dplyr::mutate(survival_res = purrr::map(data,.f=function(.x){
        tibble::tibble(logrankp=NA, coxp_categorical=NA, coxp_continus=NA)
      })) %>%
      dplyr::select(-data) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(survival_res)-> pfs_median
    .combine %>%
      dplyr::mutate(cutoff=0.75,sur_type="pfs") %>%
      dplyr::mutate(survival_res = purrr::map(data,.f=function(.x){
        tibble::tibble(logrankp=NA, coxp_categorical=NA, coxp_continus=NA)
      })) %>%
      dplyr::select(-data) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(survival_res)-> pfs_upquantile
    .combine %>%
      dplyr::mutate(cutoff=0.25,sur_type="pfs") %>%
      dplyr::mutate(survival_res = purrr::map(data,.f=function(.x){
        tibble::tibble(logrankp=NA, coxp_categorical=NA, coxp_continus=NA)
      })) %>%
      dplyr::select(-data) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(survival_res)-> pfs_lowquantile
  }
  
  os_median %>%
    rbind(os_lowquantile) %>%
    rbind(os_upquantile) %>%
    rbind(pfs_lowquantile) %>%
    rbind(pfs_median) %>%
    rbind(pfs_upquantile) %>%
    tidyr::nest(survival=c(cutoff, sur_type, logrankp, coxp_categorical, coxp_continus))
}
