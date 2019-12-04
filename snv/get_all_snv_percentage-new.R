library("magrittr")

# sample info prepare -----------------------------------------------------
readr::read_rds("/home/huff/project/data/TCGA/sample_info.rds.gz") ->sample_info #.1
readr::read_rds("/home/huff/project/data/GSCALite/TCGA/snv/syn_mutation_syn7824274_mc3_public.pass.simplification.maf.rds.gz") ->pan_maf_filter
pan_maf_filter %>%
  dplyr::mutate(sample=substr(Tumor_Sample_Barcode,1,15)) ->pan_maf_filter.sample

sample_info %>%
  dplyr::mutate(sample = substr(barcode,1,15)) %>%
  dplyr::select(-barcode) %>%
  unique() %>%
  dplyr::rename(Cancer_Types="cancer_types") -> sample_info.unique

pan_maf_filter.sample %>%
  dplyr::left_join(sample_info.unique, by="sample")  -> pan_maf_filter.cancer

# wrong cancer count, not run --------
# readr::read_tsv("/data/shiny-data/GSCALite/TCGA/snv/sample_count_SNV.tsv") %>%
#   dplyr::rename("Cancer_Types" = "cancer_types") -> cancer_count
pan_maf_filter.cancer %>%
  dplyr::select(Cancer_Types,Tumor_Sample_Barcode) %>%
  unique() %>%
  .$Cancer_Types %>%
  table() %>%
  t() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::select(-Var1) %>%
  dplyr::rename("cancer_types" = ".", "n" = "Freq") -> cancer_count


gene_list_snv_count.new <- readr::read_rds("/home/huff/project/data/TCGA/snv/.rds_snv_all_gene_snv_count-new.rds.gz")

pan_maf_filter.cancer.mut_count %>%
  dplyr::inner_join(cancer_count,by="Cancer_Types") %>%
  tidyr::unnest() %>%
  dplyr::mutate(per = mut_n/n) %>%
  dplyr::rename("sm_count" = "mut_n", "symbol" = "Hugo_Symbol" , "cancer_types" = "Cancer_Types") %>%
  tidyr::nest(-cancer_types,-n) %>%
  dplyr::rename("mut_count" = "data") -> gene_list_snv_count.new

gene_list_snv_count.new %>%
  dplyr::select(-n) %>%
  dplyr::inner_join(cancer_count,by="cancer_types") %>%
  tidyr::unnest() %>%
  dplyr::mutate(per = sm_count/n) %>%
  tidyr::nest(-cancer_types,-n) %>%
  dplyr::rename("mut_count" = "data") -> gene_list_snv_count.new.new

gene_list_snv_count.new.new %>%
  readr::write_rds(file.path("/home/huff/project/data/TCGA/snv/gene_list_snv_count.new.new.rds.gz"),compress = "gz")
