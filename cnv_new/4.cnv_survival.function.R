library(magrittr)
library(survival)
library(dplyr)


res_path <- file.path("/home/huff/data/GSCA","cnv")

cnv_group <- tibble::tibble(cnv=c(-2,-1,0,1,2),
                            group_detail=c("Homo. dele.","Hete. dele.","WT","Hete. amp.","Homo. amp."),
                            group=c("2Dele.","2Dele.","1WT","3Amp.","3Amp."))

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
    if(len_group>=2){
      kmp <- tryCatch(
        1 - pchisq(survival::survdiff(survival::Surv(time, status) ~ group, data = .d, na.action = na.exclude)$chisq, df = len_group - 1),
        error = function(e) {1}
      )
    } else {
      kmp<-NA
    }
    tibble::tibble(logrankp=kmp)
  } else {
    tibble::tibble(logrankp=NA)
  }
}

fn_survival <- function(.data,sur_type){
  
  survival_group %>%
    dplyr::filter(type==sur_type) -> sur_type_do
  
  .data %>%
    dplyr::select(sample_name,group,time=sur_type_do$time,status=sur_type_do$status) %>%
    fn_cox_logp()
}

fn_survival_res <- function(.cancer_types,.survival,.cnv){
  
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
  .cnv_data <- .cnv %>%
    tidyr::gather(-entrez,-symbol,key="barcode",value="cnv") %>%
    dplyr::inner_join(cnv_group,by="cnv") %>%
    dplyr::filter(substr(barcode,14,14)=="0") %>%
    dplyr::mutate(sample_name = substr(barcode,1,12))
  
  .cnv_data %>%
    dplyr::inner_join(.survival,by="sample_name") %>%
    tidyr::nest(-entrez,-symbol) -> .combine
  
  # os survival ----
  
  .combine %>% 
    dplyr::mutate(survival_res = purrr::map(data,fn_survival,sur_type="os")) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cols = c(survival_res)) %>%
    dplyr::mutate(sur_type="os")-> os
  
  # pfs survival -----
  if (length(grep("pfs",colnames(.survival)))>0) {
    .combine %>% 
      dplyr::mutate(survival_res = purrr::map(data,fn_survival,sur_type="pfs")) %>%
      dplyr::select(-data) %>%
      tidyr::unnest(cols = c(survival_res)) %>%
      dplyr::mutate(sur_type="pfs") -> pfs
    
    os %>%
      rbind(pfs) -> tmp
  } else {
    print("no pfs data")
    os -> tmp
  }
  
  tmp %>%
    readr::write_rds(file.path(res_path,"cancer_cnv_survival",paste(.cancer_types,"survival.cnv.rds.gz",sep="_")),compress = "gz")
  tmp
}
