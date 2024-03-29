library(magrittr)
library(doParallel)

# data path---------------------------------------------------------------

gsca_path <- file.path("/home/huff/data/GSCA")

# load data ---------------------------------------------------------------
clinical <- readr::read_rds(file.path(gsca_path,"clinical","pancan34_clinical_stage_survival_subtype.rds.gz"))

all_expr <- readr::read_rds(file = '/home/huff/data/GSCA/expr/pancan33_expr.IdTrans.rds.gz')


# combine data ------------------------------------------------------------

clinical %>%
  dplyr::filter(n.x>0) %>%
  dplyr::select(cancer_types,subtype) %>%
  dplyr::inner_join(all_expr,by="cancer_types") -> combine_data

# function ----------------------------------------------------------------
fn_test <- function(.cancer){
  combine_data %>%
    dplyr::filter(cancer_types==.cancer) -> .filter
  .filter$expr[[1]] %>%
    tidyr::gather(-entrez_id,-symbol,key="sample_name",value="exp") %>%
    dplyr::mutate(barcode = substr(sample_name,1,12),type=substr(sample_name,14,14)) %>%
    dplyr::filter(type=="0") %>%
    tidyr::drop_na(exp) -> .exp_gather
  .exp_gather %>%
    dplyr::inner_join(.filter$subtype[[1]],by="barcode") %>%
    dplyr::select(entrez_id, symbol, exp, subtype) %>%
    tidyr::nest(-entrez_id,-symbol)-> .combine
  
  .combine %>%
    dplyr::mutate(mean=purrr::map(data,.f=function(.x){
      .x %>%
        dplyr::group_by(subtype) %>%
        dplyr::mutate(meanExp=mean(exp)) %>%
        dplyr::mutate(n=dplyr::n()) %>%
        dplyr::select(-exp) %>%
        unique()
    })) %>%
    dplyr::select(entrez_id,symbol,mean) %>%
    tidyr::unnest() %>%
    dplyr::mutate(cancer_types=.cancer)
}

fn_parallel_start <- function(n_cores = 11) {
  n_detected_cores <- parallel::detectCores()
  global_cluster <<- parallel::makeForkCluster(nnodes = n_cores)
  doParallel::registerDoParallel(cl = global_cluster)
}
fn_parallel_stop <- function() {
  parallel::stopCluster(cl = global_cluster)
  foreach::registerDoSEQ()
}


# calculation ---------------------------------------------------------------------
fn_parallel_start(n_cores = 5)
exp_subtype_test <- foreach(i = combine_data$cancer_types, .packages = c('magrittr')) %dopar% {
  fn_test(i)
}
fn_parallel_stop()


# res arrangement ---------------------------------------------------------

# for (i in 1:length(exp_subtype_test)) {
#   if(i==1){
#     exp_subtype_test[[i]] %>%
#       dplyr::select(entrez_id, symbol, p.value, method, cancer_types)-> exp_subtype_test.arrange
#   }else{
#     if(nrow(exp_subtype_test[[i]])>0){
#       exp_subtype_test[[i]] %>%
#         dplyr::select(entrez_id, symbol, p.value, method, cancer_types) %>%
#         rbind(exp_subtype_test.arrange)-> exp_subtype_test.arrange
#     }else{
#       exp_subtype_test.arrange -> exp_subtype_test.arrange
#     }
#     
#   }
# }
# 
# exp_subtype_test.arrange %>%
#   tidyr::nest(-cancer_types) -> exp_subtype_test.arrange.nest

# save --------------------------------------------------------------------

exp_subtype_test %>%
  readr::write_rds(file.path(gsca_path,"expr","meanexpr_in_subtypes.NEW.IdTrans.rds.gz"),compress = "gz")
