
# clinical stage with expr ------------------------------------------------

library(magrittr)
library(dplyr)

# load data ---------------------------------------------------------------

data_path <- file.path("/home/huff/data/GSCALite/TCGA")
gsca_path <- file.path("/home/huff/data/GSCA")
git_path <- "/home/huff/github/gsca_data_process/expr"

clinical <- readr::read_rds(file.path(data_path,"clinical","pancan34_clinical.rds.gz"))

expr <- readr::read_rds(file.path(gsca_path,"expr","pancan33_expr.IdTrans.rds.gz"))

# combine data ------------------------------------------------------------
fn_stage_process <- function(cancer_types,.clinical){
  print(cancer_types)
  stage_index <- grep("stage",colnames(.clinical),value = T)
  if(length(stage_index)<1){
    return(tibble::tibble())
  } else{
    .clinical %>%
      dplyr::select(barcode,stage_index) %>%
      dplyr::rename(stage=stage_index) %>%
      dplyr::mutate(stage_type = stage_index) %>%
      dplyr::filter(!is.na(stage))
  }
}

clinical %>%
  dplyr::mutate(stage = purrr::map2(cancer_types,clinical,fn_stage_process)) %>%
  dplyr::select(-clinical) %>%
  dplyr::mutate(n=purrr::map(stage,.f=function(.x){nrow(.x)})) %>%
  tidyr::unnest(n) %>%
  dplyr::filter(n>0) -> stage_info

expr %>%
  dplyr::inner_join(stage_info,by="cancer_types") -> combine


# functions ---------------------------------------------------------------

fn_expr_process <- function(.expr) {
  .expr %>%
    dplyr::filter(!is.na(symbol)) %>%
    tidyr::gather(-entrez_id,-symbol,key="aliquot",value="expr") %>%
    dplyr::mutate(barcode = substr(aliquot,1,12)) %>%
    dplyr::mutate(type = substr(aliquot,14,14)) %>%
    dplyr::filter(type == "0") %>%
    dplyr::select(entrez=entrez_id,symbol,barcode,expr) %>%
    dplyr::distinct()
}

fn_oneway <- function(.x) {
  .x$stage %>% unique() %>% length() -> .n_stage
  if(.n_stage > 2){
    .res <- broom::tidy(oneway.test(log2(expr + 1) ~ stage, data = .x))
  }else{
    .res <- broom::tidy(wilcox.test(log2(expr + 1) ~ stage, data = .x))
  }
  .x %>%
    tidyr::nest(data=c(barcode, expr, l)) %>%
    dplyr::mutate(mean=purrr::map(data,.f=function(.x){mean(.x$expr)})) %>%
    dplyr::mutate(n=purrr::map(data,.f=function(.x){nrow(.x)})) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cols=c(mean,n)) %>%
    dplyr::mutate(mean.n = paste(signif(mean,2),n,sep = "/"),
                  stage = paste(stage,"(mean/n)")) %>%
    dplyr::select(-mean,-n) %>%
    tidyr::spread(key="stage",value="mean.n") -> .mean
  .res %>%
    dplyr::select(method,pval=p.value) %>%
    dplyr::mutate(fdr=p.adjust(pval, method = "fdr")) %>%
    cbind(.mean) 
}
fn_test <- function(.cancer_types,.expr,.stage,n) {
  print(.cancer_types)
  .expr %>% 
    fn_expr_process() -> .expr_process
  .stage %>%
    dplyr::inner_join(.expr_process,by="barcode") %>%
    dplyr::distinct() %>%
    dplyr::group_by(symbol,entrez,stage) %>%
    dplyr::mutate(l = n()) %>%
    dplyr::ungroup()-> .combine
  
   #filter out stages with less than 5 samples in one of stage.
    .combine %>%
      dplyr::filter(l>=5) %>%
      dplyr::select(stage) %>%
      dplyr::distinct() -> .stage_more5
    
    if(nrow(.stage_more5)<2){ #filter out cancers with less than 2 stages
      return(tibble::tibble())
    } else {
    .combine %>% 
      dplyr::filter(stage %in% .stage_more5$stage) %>%
      tidyr::drop_na(expr) %>%
      tidyr::nest(data = c(barcode, stage, expr, l)) %>% 
      dplyr::mutate(test = purrr::map(data,.f=fn_oneway)) %>%
      dplyr::select(-data) %>%
      tidyr::unnest(cols=c(test))
  }
}

combine %>%
  dplyr::mutate(res = purrr::pmap(list(cancer_types,expr,stage,n),.f=fn_test)) %>%
  dplyr::select(cancer_types,res) -> stage_results

stage_results %>%
  readr::write_rds(file.path(gsca_path,"expr","pan25_expr_with_stage.rds.gz"),compress = "gz")

save.image(file.path(git_path,"rda","gene_exp_with_stage.rda"))