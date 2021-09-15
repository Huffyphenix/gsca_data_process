library(magrittr)
library(survival)
library(dplyr)


res_path <- file.path("/home/huff/data/GSCA","methy","DSS-DFI_survival_210914")
# transform samples barcode -----------------------------------------------

fn_transform_samples <- function(.x) {
  .x %>% 
    tidyr::gather(-entrez,-symbol,-gene,key = 'aliquot', value = 'meth') ->
    .y
  
  .y %>% 
    dplyr::mutate(tmp = substr(x = aliquot, start = 14, stop = 14)) %>% 
    dplyr::mutate(type = ifelse(tmp == '0', 'tumor', 'normal')) %>% 
    dplyr::mutate(sample_name = substr(x = aliquot, start = 1, stop = 12)) %>%
    dplyr::group_by(entrez, symbol,  gene,sample_name) %>%
    dplyr::mutate(meth = mean(meth)) %>%
    dplyr::select(-aliquot ) %>%
    unique()
}


# meth group for survival -------------------------------------------------

survival_group <- tibble::tibble(type=c("dss","dfi"),
                                 time=c("dss_days","dfi_days"),
                                 status=c("dss_status","dfi_status"))

fn_meth_group <- function(.data,.cutoff){
  .data %>%
    dplyr::mutate(group=ifelse(meth>quantile(meth,.cutoff)[[1]],"Higher meth.", "Lower meth."))
}

fn_cox_logp <- function(.d){
  .d %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::filter(time!="#N/A") %>%
    dplyr::filter(status %in% c("1","0")) %>%
    dplyr::mutate(time=as.numeric(time), status=as.numeric(status)) %>%
    dplyr::filter(!is.na(meth)) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(n=dplyr::n()) %>%
    dplyr::select(group,n) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n>=2) %>%
    .$group %>% unique() %>% length() -> len_group
  .d %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::filter(time!="#N/A") %>%
    dplyr::filter(status %in% c("1","0")) %>%
    dplyr::mutate(time=as.numeric(time), status=as.numeric(status)) ->.d
  if(len_group==2){
    kmp <- tryCatch(
      1 - pchisq(survival::survdiff(survival::Surv(time, status) ~ group, data = .d, na.action = na.exclude)$chisq, df = len_group - 1),
      error = function(e) {1}
    )
    
    cox_categorical <- tryCatch(
      broom::tidy(survival::coxph(survival::Surv(time, status) ~ group, data = .d, na.action = na.exclude)),
      error = function(e) {1}
    )
    if (!is.numeric(cox_categorical)) {
      coxp_categorical <- cox_categorical$p.value
      hr_categorical <- exp(cox_categorical$estimate)
    } else {
      coxp_categorical <- NA
      hr_categorical <- NA
    }
    
    cox_continus <- tryCatch(
      broom::tidy(survival::coxph(survival::Surv(time, status) ~ meth, data = .d, na.action = na.exclude)),
      error = function(e) {1}
    )
    if (!is.numeric(cox_continus)) {
      coxp_continus <- cox_continus$p.value
      hr_continus <- exp(cox_continus$estimate)
    } else {
      coxp_continus <- NA
      hr_continus <- NA
    }
    if(!is.na(hr_categorical)){
      if(is.numeric(hr_categorical)){
        if(hr_categorical>1){
          higher_risk_of_death <- "Lower meth."
        }else if (hr_categorical<1){
          higher_risk_of_death <- "Higher meth."
        } else{
          higher_risk_of_death <- "Unable to judge"
        }
      } else{
        higher_risk_of_death <- "Not applicable"
      }
    }else{
      higher_risk_of_death <- "Not applicable"
    }
    
  } else {
    kmp<-NA
    coxp_categorical<-NA
    coxp_continus<-NA
    hr_continus <- NA
    hr_categorical<- NA
    higher_risk_of_death <- NA
  }
  tibble::tibble(logrankp=kmp,
                 coxp_categorical=coxp_categorical,
                 hr_categorical=hr_categorical,
                 coxp_continus=coxp_continus,
                 hr_continus=hr_continus,
                 higher_risk_of_death=higher_risk_of_death)
}

fn_survival <- function(.data,.cutoff,sur_type){
  .data %>%
    fn_meth_group(.cutoff = .cutoff) -> .data_grouped
  
  survival_group %>%
    dplyr::filter(type==sur_type) -> sur_type_do
  
  .data_grouped %>%
    dplyr::select(sample_name,group,time=sur_type_do$time,status=sur_type_do$status,meth) %>%
    fn_cox_logp()
}

fn_survival_res <- function(.cancer_types,.meth,.survival){
  print(.cancer_types)
  .survival %>%
    dplyr::filter(cancer_types == .cancer_types) %>%
    tidyr::unnest(cols = c(survival )) ->.survival_f
  
  fn_transform_samples(.meth) -> .meth_t
  
  .meth_t %>%
    dplyr::filter(type=="tumor") %>%
    dplyr::filter(!is.na(meth)) %>%
    dplyr::inner_join(.survival_f,by="sample_name") %>%
    dplyr::group_by(entrez, symbol, gene) %>%
    tidyr::nest()-> .combine
  
  # dss survival ----
  if (length(grep("dss",colnames(.survival_f)))>0) {
    .combine %>%
      dplyr::mutate(res = purrr::map(data,fn_survival,.cutoff=0.5,sur_type="dss")) %>%
      dplyr::select(-data) %>%
      tidyr::unnest() %>%
      dplyr::mutate(sur_type="dss")-> dss_median
  } else {
    print("no dss data")
    tibble::tibble() -> dss_median
  }
  
  
  # dfi survival -----
  if (length(grep("dfi",colnames(.survival_f)))>0) {
    #.survival_f %>% dplyr::filter(!is.na(.survival_f)) %>% nrow() -> .n_dfi
    #if(.n_dfi>0){
      .combine %>%
        dplyr::mutate(res = purrr::map(data,fn_survival,.cutoff=0.5,sur_type="dfi")) %>%
        dplyr::select(-data) %>%
        tidyr::unnest()  %>%
        dplyr::mutate(sur_type="dfi")-> dfi_median
    #}else{
    #  tibble::tibble(logrankp=NA, hr_categorical=NA,coxp_categorical=NA, coxp_continus=NA,hr_continus=NA,higher_risk_of_death=NA,cutoff=0.5,sur_type="dfi") -> dfi_median
    #}
  } else {
    tibble::tibble(logrankp=NA, hr_categorical=NA,coxp_categorical=NA, coxp_continus=NA,hr_continus=NA,higher_risk_of_death=NA,cutoff=0.5,sur_type="dfi") -> dfi_median
  }
  
  dss_median %>%
    rbind(dfi_median) -> tmp
  tmp %>%
    readr::write_rds(file.path(res_path,paste(.cancer_types,"DSS-DFI_survival.exp.rds.gz",sep="_")),compress = "gz")
  tmp
}
