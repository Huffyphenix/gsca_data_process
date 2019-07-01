library("magrittr")

# sample info prepare -----------------------------------------------------
readr::read_rds("/home/huff/project/data/TCGA/sample_info.rds.gz") ->sample_info #.1
readr::read_rds("/data/TCGA/TCGA_data/sample_info.rds.gz") ->sample_info #.3
# sample_info %>% dplyr::select(-barcode) ->sample_info.filter
# colnames(sample_info.filter)[1]="Tumor_Sample_Barcode"
sample_info %>%
  dplyr::mutate(sample = substr(barcode,1,15)) %>%
  dplyr::select(-barcode) %>%
  unique() %>%
  dplyr::rename(Cancer_Types="cancer_types") -> sample_info.unique

# maf prepare -------------------------------------------------------------
readr::read_rds("/home/huff/project/data/GSCALite/TCGA/snv/syn_mutation_syn7824274_mc3_public.pass.simplification.maf.rds.gz") ->pan_maf_filter
pan_maf_filter %>%
  dplyr::mutate(sample=substr(Tumor_Sample_Barcode,1,15)) ->pan_maf_filter.sample

# get snv of each samples -------------------------------------------------

fn_get_sample_count <- function(mutation_type){
  mutation_type %>% 
    dplyr::filter(mutation == "Mut") %>%
    nrow() -> mut_n
  mutation_type %>% 
    dplyr::filter(mutation == "Nmut") %>%
    nrow() -> nmut_n
  tibble::tibble(mut_n = mut_n) 
}

pan_maf_filter.sample %>%
  dplyr::select(Hugo_Symbol,Variant_Classification,Variant_Type,Tumor_Sample_Barcode ,sample) %>%
  dplyr::mutate(mutation =  "Mut") %>%
  tidyr::nest(-c(Hugo_Symbol,Tumor_Sample_Barcode)) -> pan_maf_filter.sample.for_PARTITION_ID

cl <-10 
cluster <- multidplyr::create_cluster(core=cl)
pan_maf_filter.sample.for_PARTITION_ID %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_get_sample_count", fn_get_sample_count)  %>%
  dplyr::group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>%
  dplyr::mutate(count = purrr::map(data,fn_get_sample_count)) %>%
  dplyr::collect() %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::mutate(sample = substr(Tumor_Sample_Barcode,1,15)) %>%
  dplyr::rename("symbol" = "Hugo_Symbol") -> pan_maf_filter.sample.mut_count
parallel::stopCluster(cluster)


# mutation burden ---------------------------------------------------------

pan_maf_filter.sample.mut_count%>%
  dplyr::left_join(sample_info.unique,by="sample") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(mutation_durden = sum(mut_n)) %>%
  dplyr::select(-symbol,-mut_n) %>%
  unique() %>%
  dplyr::ungroup() -> pan_maf_filter.sample.mut_burden
pan_maf_filter.sample.mut_burden %>%
  readr::write_rds(file.path("/home/huff/project/data/TCGA/snv","pancan33_maf_syn7824274.sample.mut_burden.rds.gz"),compress = "gz")

fn_spread <- function(.x){
  .x %>%
    tidyr::spread(key=sample,value = mut_n) %>%
    tidyr::gather(-symbol,key="sample",value="mut_n") %>%
    dplyr::mutate(mut_n = ifelse(is.na(mut_n),0,mut_n)) %>%
    tidyr::spread(key=sample,value = mut_n)
}

pan_maf_filter.sample.mut_count %>%
  dplyr::left_join(sample_info.unique,by="sample") %>%
  dplyr::select(-Tumor_Sample_Barcode) %>%
  tidyr::nest(-Cancer_Types,.key="snv") %>% 
  dplyr::group_by(Cancer_Types) %>%
  dplyr::mutate(snv = purrr::map(snv,fn_spread)) %>%
  readr::write_rds(file.path("/home/huff/project/data/TCGA/snv/","NEW-pancan33_snv_from_syn7824274-allMutationType.rds.gz"),compress = "gz")
