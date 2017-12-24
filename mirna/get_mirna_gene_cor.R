
library(magrittr)
# config ------------------------------------------------------------------

config <- list()
config$database <- c("/data/TCGA/TCGA_data")
config$out <- c("/data/GSCALite/")

# load data ---------------------------------------------------------------

mirna_exp<- readr::read_rds(file.path(config$database,"pancan33_mirna_expr.rds.gz"))
gene_exp <- readr::read_rds(file.path(config$database,"pancan33_expr.rds.gz"))
mirna2gene <- readr::read_rds(file.path("/data/GSCALite/TCGA/mirna/miRNA2gene.all.nest.rds.gz"))

# merge
fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 16
  )
}

fn_mirna_exp <- function(mirna_exp){
  mirna_exp %>%
    dplyr::select(-gene) %>%
    tidyr::gather(-name,key="sample",value="mirna_exp") %>%
    dplyr::mutate(sample=fun_barcode(sample)) %>%
    dplyr::rename(mirna=name)->mirna_exp.p
  return(mirna_exp.p)
}
fn_gene_exp <- function(gene_exp){
  gene_exp %>%
    dplyr::select(-entrez_id) %>%
    dplyr::filter(symbol!="?") %>%
    tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
    dplyr::mutate(sample=fun_barcode(sample)) -> gene_exp.p
  return(gene_exp.p)
}
gene_exp %>%
  # dplyr::filter(cancer_types=="CESC") %>%
  dplyr::mutate(gene_exp=purrr::map(expr,fn_gene_exp)) %>%
  dplyr::select(-expr) -> gene_exp.gather
mirna_exp %>%
  # dplyr::filter(cancer_types=="CESC") %>%
  dplyr::mutate(mirna_exp=purrr::map(mirna,fn_mirna_exp)) %>%
  dplyr::select(-mirna) -> mirna_exp.gather

fn_get_spm_a <-function(cancer,gene_e) {
  print(cancer)
  mirna_exp.gather %>%
    dplyr::filter(cancer_types %in% cancer) %>%
    tidyr::unnest() -> mirna_ready
  
  gene_e %>%
    tidyr::nest(-symbol) %>%
    dplyr::group_by(symbol) %>%
    # dplyr::filter(symbol %in% c("ABI2")) %>%
    dplyr::mutate(spm=purrr::map2(symbol,data,fn_get_spm_b)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
  dplyr::ungroup()  -> .out
return(.out)
}
fn_get_spm_b <- function(s,data){
  .out<-data.frame(mirna="test") 
  # print(s)
  mirna2gene %>%
    dplyr::filter(symbol %in% s) ->.symbol_mirna_data
  if(nrow(.symbol_mirna_data)==1){
    .symbol_mirna_data %>%
      tidyr::unnest() %>%
      dplyr::select(mirna) %>%
      dplyr::inner_join(mirna_ready,by="mirna") ->m_data_ready
  for(m in m_data_ready$mirna %>% unique()){
    # print(m)
    m_data_ready %>%
      dplyr::filter(mirna %in% m) %>%
      dplyr::inner_join(data,by="sample") ->tmp
      tryCatch(cor.test(tmp$gene_exp,tmp$mirna_exp,method = c("pearson")),
               warning =function(e) 2 ,
               error=function(e) 1)->tmp.person
    if(length(tmp.person)!=1){
      if(tmp.person$p.value<0.05 && tmp.person$estimate <= -0.3){
        data.frame(mirna=m) ->tmp.out
        rbind(.out,tmp.out) ->.out
      }
      # tmp.person$estimate ->cor
      # tmp.person$p.value -> p
      # data.frame(mirna=m,gene=s,cor=cor,p=p) ->tmp.out
      # rbind(.out,tmp.out) ->.out
    }
  }
  }
  .out$mirna <- .out$mirna %>% as.character()
  return(.out[-1,])
}

cl<- 10
cluster <- multidplyr::create_cluster(cl)
gene_exp.gather %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_get_spm_b", fn_get_spm_b)  %>%
  multidplyr::cluster_assign_value("fn_get_spm_a", fn_get_spm_a)  %>%
  multidplyr::cluster_assign_value("mirna2gene", mirna2gene)  %>%
  multidplyr::cluster_assign_value("mirna_exp.gather",mirna_exp.gather)  %>%
  multidplyr::cluster_assign_value("mirna_ready",mirna_ready)  %>%
  dplyr::mutate(mirna_result=purrr::map2(cancer_types,gene_exp,fn_get_spm_a)) %>%
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>%
  dplyr::select(-gene_exp) -> pan33_gene_cor_mirna
parallel::stopCluster(cluster)

pan33_gene_cor_mirna %>%
  readr::write_rds(file.path(config$out,"TCGA","mirna","pan33_gene_cor_with_mirna.rds.gz"),compress="gz")

 # save.image(file = "/project/huff/huff/github/gsca_data_process/mirna/.RData")
