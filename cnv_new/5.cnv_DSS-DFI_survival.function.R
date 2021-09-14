library(magrittr)
library(survival)
library(dplyr)


res_path <- file.path("/home/huff/data/GSCA","cnv")

cnv_group <- tibble::tibble(cnv=c(-2,-1,0,1,2),
                            group_detail=c("Homo. dele.","Hete. dele.","WT","Hete. amp.","Homo. amp."),
                            group=c("2Dele.","2Dele.","1WT","3Amp.","3Amp."))

# expr group for survival -------------------------------------------------

survival_group <- tibble::tibble(type=c("dss","dfi"),
                                 time=c("dss_days","dfi_days"),
                                 status=c("dss_status","dfi_status"))

fn_cox_logp <- function(.d){
  .d %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::filter(time!="#N/A") %>%
    dplyr::filter(status %in% c("1","0")) %>%
    dplyr::mutate(time=as.numeric(time), status=as.numeric(status)) %>%
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
  if(!is.na(len_group)){
    if(len_group>=2){
      kmp <- tryCatch(
        1 - pchisq(survival::survdiff(survival::Surv(time, status) ~ group, data = .d, na.action = na.exclude)$chisq, df = len_group - 1),
        error = function(e) {1}
      )
      
      coxp <- tryCatch(
        broom::tidy(survival::coxph(survival::Surv(time, status) ~ group, data = .d, na.action = na.exclude)),
        error = function(e) {1}
      )
      if (!is.numeric(coxp)) {
        cox_p <- coxp
      } else {
        cox_p <- tibble::tibble(term=NA,estimate=NA, std.error=NA,statistic=NA, p.value=NA)
      }
    } else {
      kmp<-NA
      cox_p <- tibble::tibble(term=NA,estimate=NA, std.error=NA,statistic=NA, p.value=NA)
    }
    tibble::tibble(logrankp=kmp,coxp=tidyr::nest(cox_p))
    
  } else {
    tibble::tibble(logrankp=NA,coxp=tidyr::nest("not avaliable"))
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
  print(.cancer_types)
  .cnv_data <- .cnv %>%
    tidyr::gather(-entrez,-symbol,key="barcode",value="cnv") %>%
    dplyr::inner_join(cnv_group,by="cnv") %>%
    dplyr::filter(substr(barcode,14,14)=="0") %>%
    dplyr::mutate(sample_name = substr(barcode,1,12))
  
  .cnv_data %>%
    dplyr::inner_join(.survival,by="sample_name") %>%
    tidyr::nest(-entrez,-symbol) -> .combine
  
  # dss survival ----
  if (length(grep("dss",colnames(.survival)))>0) {
    .combine %>% 
      dplyr::mutate(survival_res = purrr::map(data,fn_survival,sur_type="dss")) %>%
      dplyr::select(-data) %>%
      tidyr::unnest(cols = c(survival_res)) %>%
      dplyr::mutate(sur_type="dss")-> dss
  }else{
    print("no dss data")
    tibble::tibble() -> dss
  }
  
  
  # dfi survival -----
  if (length(grep("dfi",colnames(.survival)))>0) {
    .combine %>% 
      dplyr::mutate(survival_res = purrr::map(data,fn_survival,sur_type="dfi")) %>%
      dplyr::select(-data) %>%
      tidyr::unnest(cols = c(survival_res)) %>%
      dplyr::mutate(sur_type="dfi") -> dfi
    
    dss %>%
      rbind(dfi) -> tmp
  } else {
    print("no dfi data")
    dss -> tmp
  }
  
  tmp %>%
    readr::write_rds(file.path(res_path,"cancer_cnv_DSS-DFI_survival_210914",paste(.cancer_types,"DSS-DFI_survival.cnv.rds.gz",sep="_")),compress = "gz")
  tmp
}
