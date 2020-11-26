library(magrittr)
library(survival)
library(dplyr)


res_path <- file.path("/home/huff/data/GSCA","snv")



# expr group for survival -------------------------------------------------

survival_group <- tibble::tibble(type=c("os","pfs"),
                                 time=c("os_days","pfs_days"),
                                 status=c("os_status","pfs_status"))

fn_cox_logp <- function(.d){
  .d %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::filter(!is.na(time)) %>%
    dplyr::filter(!is.na(status)) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(n=dplyr::n()) %>%
    dplyr::select(group,n) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n>5) %>%
    .$group %>% unique() %>% length() -> len_group
  if(!is.na(len_group)){
    if(len_group==2){
      kmp <- tryCatch(
        1 - pchisq(survival::survdiff(survival::Surv(time, status) ~ group, data = .d, na.action = na.exclude)$chisq, df = len_group - 1),
        error = function(e) {1}
      )
      
      coxp <- tryCatch(
        broom::tidy(survival::coxph(survival::Surv(time, status) ~ group, data = .d, na.action = na.exclude)),
        error = function(e) {1}
      )
      if (!is.numeric(coxp)) {
        cox_p <- coxp$p.value
        hr <- exp(coxp$estimate)
      } else {
        cox_p <- 1
        hr <- 1
      }
      
      if(!is.na(hr)){
        if(hr>1){
          higher_risk_of_death <- "Mutated"
        }else if (hr<1){
          higher_risk_of_death <- "Non-mutated"
        }else {
          higher_risk_of_death <- NA
        }
      } else {
        higher_risk_of_death <- NA
      }
      
    } else {
      kmp<-1
      cox_p<-1
      hr <- 1
      higher_risk_of_death <- NA
    }
    tibble::tibble(logrankp=kmp,cox_p=cox_p,hr=hr,higher_risk_of_death=higher_risk_of_death)
  } else {
    tibble::tibble(logrankp=1,cox_p=1,hr=1,higher_risk_of_death=NA)
  }
}

fn_survival <- function(.data,sur_type){
  
  survival_group %>%
    dplyr::filter(type==sur_type) -> sur_type_do
  
  .data %>%
    dplyr::select(sample_name,group,time=sur_type_do$time,status=sur_type_do$status) %>%
    fn_cox_logp()
}

fn_survival_res <- function(.cancer_types,.survival){
  .survival %>%
    dplyr::filter(cancer_types == .cancer_types) %>%
    tidyr::unnest(c(survival)) ->.survival
  
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
      tidyr::unnest(c(os_status, pfs_status)) -> .survival
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
      tidyr::unnest(c(os_status)) -> .survival
  }
  oncoplot_effective_mutation <- 
  .snv_data <- readr::read_rds(file.path(gsca_v2_path,"snv","sub_cancer_maf_tsv",paste(.cancer_types,"maf_data.IdTrans.tsv.rds.gz",sep = "_")))%>%
    dplyr::mutate(group = ifelse(Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Splice_Site","Frame_Shift_Del","In_Frame_Del","In_Frame_Ins"), "2_mutated","1_nonmutated")) %>%
    dplyr::mutate(sample_name=substr(Tumor_Sample_Barcode,1,12)) %>%
    dplyr::select(Hugo_Symbol,entrez,sample_name) 
  
  .snv_data$sample_name %>%
    unique() -> sample_with_snv
  .survival %>%
    dplyr::filter(sample_name %in% sample_with_snv) -> .survival
  
  .snv_data %>% 
    tidyr::nest(data=c(sample_name,group)) %>%
    dplyr::mutate(combine_data = purrr::map(data,.f=function(.x){
      .x %>%
        dplyr::right_join(.survival,by="sample_name") %>%
        dplyr::mutate(group = ifelse(is.na(group),"1_nonmutated",group))
    })) %>%
    dplyr::select(-data) -> .combine
  
  # os survival ----
  
  .combine %>%
    dplyr::mutate(survival_res = purrr::map(combine_data,fn_survival,sur_type="os")) %>%
    dplyr::select(-combine_data) %>%
    tidyr::unnest(cols = c(survival_res)) %>%
    dplyr::mutate(sur_type="os")-> os
  
  # pfs survival -----
  if (length(grep("pfs",colnames(.survival)))>0) {
    .combine %>% 
      dplyr::mutate(survival_res = purrr::map(combine_data,fn_survival,sur_type="pfs")) %>%
      dplyr::select(-combine_data) %>%
      tidyr::unnest(cols = c(survival_res)) %>%
      dplyr::mutate(sur_type="pfs") -> pfs
    
    os %>%
      rbind(pfs) -> tmp
  } else {
    print("no pfs data")
    os -> tmp
  }
  
  tmp %>%
    readr::write_rds(file.path(res_path,"cancer_snv_survival",paste(.cancer_types,"survival.snv.rds.gz",sep="_")),compress = "gz")
  tmp
}
