library("magrittr")

# sample info prepare -----------------------------------------------------
readr::read_rds("/data/TCGA/TCGA_data/sample_info.rds.gz") ->sample_info
# sample_info %>% dplyr::select(-barcode) ->sample_info.filter
# colnames(sample_info.filter)[1]="Tumor_Sample_Barcode"
sample_info %>%
  dplyr::select(-barcode) %>%
  unique() -> sample_info.unique


# maf prepare -------------------------------------------------------------
# readr::read_tsv("/data/TCGA/TCGA_data/syn_mutation_syn7824274_mc3_public.pass.maf") -> pan_maf
# pan_maf %>% dplyr::select(Hugo_Symbol,Entrez_Gene_Id,Center,NCBI_Build,Chromosome,
#                           Start_Position,End_Position,Strand,Variant_Classification,
#                           Variant_Type,Reference_Allele,Tumor_Seq_Allele1,
#                           Tumor_Seq_Allele2,Tumor_Sample_Barcode) ->pan_maf_filter
# readr::write_rds(pan_maf_filter,path="/data/GSCALite/TCGA/snv/syn_mutation_syn7824274_mc3_public.pass.simplification.maf.rds.gz",compress = "gz")
readr::read_rds("/data/GSCALite/TCGA/snv/syn_mutation_syn7824274_mc3_public.pass.simplification.maf.rds.gz") ->pan_maf_filter
######??????????????????????????? barcode don't match
pan_maf_filter %>%
  dplyr::mutate(sample=substr(Tumor_Sample_Barcode,1,12)) ->pan_maf_filter.sample

# get tbs info for cancer type file ---------------------------------------

pan_maf_filter.sample %>%
  dplyr::select(Tumor_Sample_Barcode,sample) %>%
  unique()->pan_maf_filter.tbs

pan_maf_filter.tbs %>%
  dplyr::inner_join(sample_info.unique,by="sample") %>%
  dplyr::select(-sample)->sample_info.tbs
readr::write_rds(sample_info.tbs,path="/data/GSCALite/TCGA/snv/syn_mutation_syn7824274_mc3_public.pass.clinical.tsv.rds.gz",compress = "gz")
# give cancer type infomation to maf file ---------------------------------

fun_filter_cancer_type <-function(id){
  # data<-x1
  # id<-"z_1"
  sample_info.unique %>%
    dplyr::filter(sample %in% id) %>%
    dplyr::select(cancer_types) %>%
    unique() %>% 
    dplyr::pull(cancer_types) %>%
    as.character() ->a
  return(a)
}



cl <- 32
cluster <- multidplyr::create_cluster(core = cl)
pan_maf_filter.sample %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fun_filter_cancer_type", fun_filter_cancer_type) %>%
  multidplyr::cluster_assign_value("sample_info.unique",sample_info.unique) %>%
  dplyr::mutate(cancer_types=purrr::map(sample,fun_filter_cancer_type)) %>%
  dplyr::mutate(cancer_types=unlist(cancer_types)) %>%
  dplyr::filter(cancer_types!="character(0)") %>%
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) -> pan_maf_filter.cancer
parallel::stopCluster(cluster)



# get maf -----------------------------------------------------------------
maftools::read.maf(maf=pan_maf_filter.cancer,clinicalData=sample_info.tbs) ->maf_filter

# output ------------------------------------------------------------------



readr::write_rds(maf_filter,"/data/GSCALite/TCGA/snv/snv_mutation_mc3_public.pass.filtered_maf.rds.gz",compress="gz")
