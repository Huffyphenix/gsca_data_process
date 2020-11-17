library(maftools)
library(magrittr)
data_path <- "/home/huff/data/GSCA/mutation/snv"

pan_maf %>%
  getClinicalData() %>%
  .$Cancer_Types %>%
  unique() -> cancer_types

for (cancer in cancer_types) {
  query <- as.expression(paste0("Cancer_Types %in% c('", cancer_types, "')"))
  subset_maf <- subsetMaf(maf = pan_maf, query = query)
  filename <- paste(cancer,".maf.rds.gz",sep="")
  subset_maf %>%
    readr::write_rds(file.path(data_path,"subset_maf",filename),compress = "gz")
}
