
# clinical stage with expr ------------------------------------------------

library(magrittr)
library(dplyr)
library(doParallel)

# load data ---------------------------------------------------------------

data_path <- file.path("/home/huff/data/GSCALite/TCGA")
gsca_path <- file.path("/home/huff/data/GSCA")
git_path <- "/home/huff/github/gsca_data_process/expr"

stage <- readr::read_rds(file.path(gsca_path,"clinical","Pancan.Merge.clinical-STAGE.rds.gz")) 

stage %>%
  dplyr::mutate(count= purrr::map(stage,.f=function(.x){
    .x %>%
      dplyr::filter(!is.na(pathologic_stage)) %>%
      nrow() ->.n_pathological_stage
    .x %>%
      dplyr::filter(!is.na(clinical_stage)) %>%
      nrow() ->.n_clinical_stage
    .x %>%
      dplyr::filter(!is.na(igcccg_stage)) %>%
      nrow() ->.n_igcccg_stage
    .x %>%
      dplyr::filter(!is.na(masaoka_stage)) %>%
      nrow() ->.n_masaoka_stage
    tibble::tibble(pathological_stage=.n_pathological_stage,clinical_stage=.n_clinical_stage,igcccg_stage=.n_igcccg_stage,masaoka_stage=.n_masaoka_stage)
  })) %>%
  dplyr::select(-stage) %>%
  tidyr::unnest() -> stage.count

stage.count %>%
  dplyr::filter(masaoka_stage>0) -> cancers.masaoka_stage

expr <- readr::read_rds(file.path(gsca_path,"expr","pancan33_expr.IdTrans.rds.gz"))

# combine data ------------------------------------------------------------
stage %>%
  dplyr::filter(cancer_types%in% cancers.masaoka_stage$cancer_types) %>%
  tidyr::unnest() %>%
  .$masaoka_stage %>%
  unique() %>% sort()
stages_included <- tibble::tibble(stage=c("Stage II","Stage II",
                                          "Stage III",
                                          "Stage IV","Stage IV"),
                                  stage_raw=c("iia","iib",
                                              "iii",
                                              "iva","ivb"))
fn_stage_process <- function(cancer_types,.clinical){
  print(cancer_types)
  stage_index <- grep("stage",colnames(.clinical),value = T)
  if(length(stage_index)<1){
    return(tibble::tibble())
  } else{
    .clinical %>%
      dplyr::rename("Stage"=stage_index) %>%
      dplyr::filter(!is.na(Stage)) %>%
      dplyr::rename(stage_raw=Stage) %>%
      dplyr::inner_join(stages_included,by="stage_raw") %>%
      dplyr::select(-stage_raw) %>%
      dplyr::mutate(stage_type = stage_index) %>%
      dplyr::filter(!is.na(stage))
  }
}

stage %>%
  dplyr::filter(cancer_types%in% cancers.masaoka_stage$cancer_types) %>%
  dplyr::mutate(stage = purrr::map(stage,.f=function(.x){
    .x %>%
      dplyr::mutate(barcode=toupper(barcode)) %>%
      dplyr::select(barcode,masaoka_stage)
  })) %>%
  dplyr::mutate(stage = purrr::map2(cancer_types,stage,fn_stage_process)) %>%
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
    # dplyr::mutate(fdr=p.adjust(pval, method = "fdr")) %>%
    cbind(.mean) 
}
fn_test <- function(.cancer_types) {
  print(.cancer_types)
  combine %>%
    dplyr::filter(cancer_types == .cancer_types) -> .tmp
  .tmp$expr[[1]]-> .expr
  .tmp$stage[[1]]-> .stage
  .tmp$n[[1]]-> n
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
      tidyr::unnest(cols=c(test)) %>%
      dplyr::mutate(cancer_types=.cancer_types)-> .tmp
    .tmp %>%
      dplyr::mutate(fdr=p.adjust(.tmp$pval, method = "fdr"))
  }
}

# fn_parallel_start <- function(n_cores = 5) {
#   n_detected_cores <- parallel::detectCores()
#   global_cluster <<- parallel::makeForkCluster(nnodes = n_cores)
#   doParallel::registerDoParallel(cl = global_cluster)
# }
# fn_parallel_stop <- function() {
#   parallel::stopCluster(cl = global_cluster)
#   foreach::registerDoSEQ()
# }

# calculation ---------------------------------------------------------------------
# fn_parallel_start(n_cores = 5)
stage_results <- foreach(i = combine$cancer_types, .packages = c('magrittr')) %dopar% {
  fn_test(i)
}
# fn_parallel_stop()

# res arrangement ---------------------------------------------------------

for (i in 1:length(stage_results)) {
  if(i==1){
    stage_results[[i]] %>%
      tidyr::nest(-cancer_types)-> stage_results.arrange
  }else{
    if(nrow(stage_results[[i]])>0){
      stage_results[[i]]  %>%
        tidyr::nest(-cancer_types) %>%
        rbind(stage_results.arrange)-> stage_results.arrange
    }else{
      stage_results.arrange -> stage_results.arrange
    }
    
  }
}
# combine %>%
# dplyr::mutate(res = purrr::pmap(list(cancer_types,expr,stage,n),.f=fn_test)) %>%
# dplyr::select(cancer_types,res) -> stage_results

stage_results.arrange %>%
  readr::write_rds(file.path(gsca_path,"expr","TGCT_expr_with_masaokaStage-210813.rds.gz"),compress = "gz")

save.image(file.path(git_path,"rda","gene_exp_with_masaokaStage--210813.rda"))