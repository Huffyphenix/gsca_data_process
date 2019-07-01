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
# readr::read_tsv("/data/TCGA/TCGA_data/syn_mutation_syn7824274_mc3_public.pass.maf") -> pan_maf
# pan_maf %>% dplyr::select(Hugo_Symbol,Entrez_Gene_Id,Center,NCBI_Build,Chromosome,
#                           Start_Position,End_Position,Strand,Variant_Classification,
#                           Variant_Type,Reference_Allele,Tumor_Seq_Allele1,
#                           Tumor_Seq_Allele2,Tumor_Sample_Barcode) ->pan_maf_filter
# readr::write_rds(pan_maf_filter,path="/data/GSCALite/TCGA/snv/syn_mutation_syn7824274_mc3_public.pass.simplification.maf.rds.gz",compress = "gz")
# readr::read_rds("/data/shiny-data/GSCALite/TCGA/snv/syn_mutation_syn7824274_mc3_public.pass.simplification.maf.rds.gz") ->pan_maf_filter
readr::read_rds("/home/huff/project/data/GSCALite/TCGA/snv/syn_mutation_syn7824274_mc3_public.pass.simplification.maf.rds.gz") ->pan_maf_filter
pan_maf_filter %>%
  dplyr::mutate(sample=substr(Tumor_Sample_Barcode,1,15)) ->pan_maf_filter.sample

#############################################
# give cancer type infomation to maf file ---------------------------------

# fun_filter_cancer_type <-function(id){
#   # data<-x1
#   # id<-"z_1"
#   sample_info.unique %>%
#     dplyr::filter(sample %in% id) %>%
#     dplyr::select(Cancer_Types) %>%
#     unique() %>% 
#     dplyr::pull(Cancer_Types) %>%
#     as.character() ->a
#   return(a)
# }
# cl <- 15
# cluster <- multidplyr::create_cluster(core = cl)
# pan_maf_filter.sample %>%
#   multidplyr::partition(cluster = cluster) %>%
#   multidplyr::cluster_library("magrittr") %>%
#   multidplyr::cluster_assign_value("fun_filter_cancer_type", fun_filter_cancer_type) %>%
#   multidplyr::cluster_assign_value("sample_info.unique",sample_info.unique) %>%
#   dplyr::mutate(Cancer_Types=purrr::map(sample,fun_filter_cancer_type)) %>%
#   dplyr::mutate(Cancer_Types=unlist(Cancer_Types)) %>%
#   dplyr::filter(Cancer_Types!="character(0)") %>%
#   dplyr::collect() %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-PARTITION_ID) -> pan_maf_filter.cancer
# parallel::stopCluster(cluster)
# 
# pan_maf_filter.cancer
#############################################


#############################################
# calculation of the SNV persentage ---------------------------------------
oncoplot_effective_mutation <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Splice_Site","Frame_Shift_Del","In_Frame_Del","In_Frame_Ins")

# fn_get_count <- function(mutation_type){
#   mutation_type %>% 
#     dplyr::filter(mutation == "Mut") %>%
#     dplyr::select(Tumor_Sample_Barcode) %>%
#     unique() %>% nrow() -> mut_n
#   mutation_type %>% 
#     dplyr::filter(mutation == "Nmut") %>%
#     dplyr::select(Tumor_Sample_Barcode) %>%
#     unique() %>% nrow() -> nmut_n
#   data.frame(mut_n = mut_n)
# }

# pan_maf_filter.cancer %>%
#   dplyr::select(Hugo_Symbol,Variant_Classification,Variant_Type,Tumor_Sample_Barcode ,sample,Cancer_Types) %>%
#   dplyr::mutate(mutation = ifelse(Variant_Classification %in% oncoplot_effective_mutation, "Mut", "Nmut")) %>%
#   tidyr::nest(-c(Hugo_Symbol,Cancer_Types)) %>%
#   dplyr::group_by(Hugo_Symbol,Cancer_Types) %>%
#   dplyr::mutate(count = purrr::map(data,fn_get_count)) %>%
#   dplyr::select(-data) %>%
#   tidyr::unnest() %>%
#   tidyr::nest(-Cancer_Types) -> pan_maf_filter.cancer.mut_count
#  
# readr::read_tsv("/data/shiny-data/GSCALite/TCGA/snv/sample_count_SNV.tsv") %>%
#   dplyr::rename("Cancer_Types" = "cancer_types") -> cancer_count
# 
# pan_maf_filter.cancer.mut_count %>%
#   dplyr::inner_join(cancer_count,by="Cancer_Types") %>%
#   tidyr::unnest() %>%
#   dplyr::mutate(per = mut_n/n) %>%
#   dplyr::rename("sm_count" = "mut_n", "symbol" = "Hugo_Symbol" , "cancer_types" = "Cancer_Types") %>%
#   tidyr::nest(-cancer_types,-n) %>%
#   dplyr::rename("mut_count" = "data") -> gene_list_snv_count.new
# 
# gene_list_snv_count.new %>% 
#   readr::write_rds(path = file.path("/data/shiny-data/GSCALite/TCGA/snv/",".rds_snv_all_gene_snv_count-new.rds.gz"), compress = "gz")

#############################################




#############################################
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
  dplyr::mutate(mutation = ifelse(Variant_Classification %in% oncoplot_effective_mutation, "Mut", "Nmut")) %>%
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
  dplyr::mutate(sample = substr(Tumor_Sample_Barcode,1,12)) %>%
  dplyr::rename("symbol" = "Hugo_Symbol") -> pan_maf_filter.sample.mut_count
parallel::stopCluster(cluster)

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
  readr::write_rds(file.path("/data/shiny-data/GSCALite/TCGA/snv/","pancan33_snv_from_syn7824274_spread.rds.gz"),compress = "gz")


# readr::read_rds("/data/shiny-data/GSCALite/TCGA/snv/pancan33_snv.rds.gz") -> pancan33_snv
# fn_gather <- function(.x,.y){
#   print(.x)
#   .y %>%
#     tidyr::gather(-symbol,key="sample",value=n_mut_old)
# }
# pancan33_snv %>%
#   dplyr::mutate(gather = purrr::map2(cancer_types,snv,fn_gather)) %>%
#   dplyr::select(-snv) %>%
#   tidyr::unnest() -> pan_maf_filter.sample.mut_count.old
# 
# pan_maf_filter.sample.mut_count.old %>%
#   dplyr::full_join(pan_maf_filter.sample.mut_count,by=c("symbol","sample")) %>%
#   dplyr::mutate(mut_n = ifelse(is.na(mut_n),0,mut_n)) %>%
#   dplyr::select(-n_mut_old,-Tumor_Sample_Barcode) -> pan_maf_filter.sample.mut_count.oldnew.combine
# 
# fn_spread <- function(.x){
#   .x %>%
#     tidyr::spread(key=sample,value = mut_n)
# }
# pan_maf_filter.sample.mut_count.oldnew.combine %>%
#   tidyr::nest(-cancer_types) %>%
#   dplyr::group_by(cancer_types) %>%
#   dplyr::mutate(data = purrr::map(data,.f=fn_spread)) -> pan_maf_filter.sample.mut_count.oldnew.combine.spread

# pan_maf_filter.sample.mut_count.oldnew.combine.spread %>%
#   readr::write_rds(file.path("/data/shiny-data/GSCALite/TCGA/snv/","pancan33_snv_from_syn7824274.rds.gz"), compress = "gz")
#   
# pan_maf_filter.sample.mut_count.oldnew.combine.spread %>%
#   readr::write_rds(file.path("/data/shiny-data/GSCALite/TCGA/snv/","pancan33_snv_from_syn7824274.rds.gz"), compress = "gz")
