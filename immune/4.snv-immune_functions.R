fn_wilcoxn <- function(.data){
  wilcox <- tryCatch(
    broom::tidy(wilcox.test(immune ~ group, data = .data)),
    error = function(e) {1}
  )
  .data %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(immune_mean = mean(immune)) %>%
    dplyr::select(group,immune_mean) %>%
    unique() %>%
    tidyr::spread(key="group",value="immune_mean") -> meanimmune
  
  logfc <- log2(meanimmune$`2_mutated`/(meanimmune$`1_nonmutated`+0.00001))
  if(logfc>0){
    Higher_immune_group <- "Mutated"
  }else if(logfc<0){
    Higher_immune_group <- "WT"
  }else{
    Higher_immune_group <- "NA"
  }
  if(!is.na(wilcox)){
    if (!is.numeric(wilcox)) {
      wilcox %>% dplyr::mutate(logfc=logfc,Higher_immune_group=Higher_immune_group)
    } else {
      tibble::tibble()
    }
  }else {
    tibble::tibble()
  }
}
fn_res <- function(.combine_data){
  .combine_data %>%
    tidyr::gather(-barcode,-group,-aliquot,-sample_name,key="cell_type",value="immune") %>%
    dplyr::group_by(cell_type) %>%
    tidyr::nest() -> for_calculation
  for_calculation %>%
    dplyr::mutate(res = purrr::map(data,fn_wilcoxn)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(res)
}

fn_immune_snv <- function(.cancer_types,.immune){
  message(glue::glue('Handling {.cancer_types} snv-immune'))
  
    .snv_data <- readr::read_rds(file.path(gsca_v2_path,"snv","sub_cancer_maf_tsv",paste(.cancer_types,"maf_data.IdTrans.tsv.rds.gz",sep = "_")))%>%
    dplyr::mutate(group = ifelse(Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Splice_Site","Frame_Shift_Del","In_Frame_Del","In_Frame_Ins"), "2_mutated","1_nonmutated")) %>%
    dplyr::mutate(barcode=substr(Tumor_Sample_Barcode,1,16)) %>%
    dplyr::select(Hugo_Symbol,entrez,barcode,group) 
  
  .snv_data$barcode %>%
    unique() -> sample_with_snv
  .immune %>%
    dplyr::filter(barcode %in% sample_with_snv) -> .immune
  
  .snv_data %>% 
    tidyr::nest(data=c(barcode,group)) %>%
    dplyr::mutate(combine_data = purrr::map(data,.f=function(.x){
      .x %>%
        dplyr::right_join(.immune,by="barcode") %>%
        dplyr::mutate(group = ifelse(is.na(group),"1_nonmutated",group))
    })) %>%
    dplyr::select(-data) -> .combine
  
  # os survival ----
  
  .combine %>%
    dplyr::mutate(res = purrr::map(combine_data,fn_res)) %>%
    dplyr::select(-combine_data) %>%
    tidyr::unnest(cols = c(res)) -> res
  message(glue::glue('complete {.cancer_types} snv-immune'))
  res %>%
    readr::write_rds(file.path(gsca_v2_path,"TIL","snv_immune",paste(.cancer_types,"snv_immune_wilcox.rds.gz",sep=".")))
  return(res)
}
