
# Library -----------------------------------------------------------------

library(magrittr)

# data path ---------------------------------------------------------------
data_path <- "/home/huff/data/GSCA/TIL"

# Load expr survival----------------------------------------------------------------

file_list <- grep("*.cnv_immune_cor.rds.gz",dir(file.path(data_path)),value = TRUE)
cnv_immune <- tibble::tibble()
for (file in file_list) {
  .cnv_immune <- readr::read_rds(file.path(data_path,file)) %>%
    dplyr::mutate(cancer_types= strsplit(file,split = "\\.")[[1]][1])
  if(nrow(cnv_immune)<1){
    cnv_immune<-.cnv_immune
  } else {
    rbind(cnv_immune,.cnv_immune) ->cnv_immune
  }
}
cnv_immune %>%
  dplyr::group_by(cancer_types) %>%
  tidyr::nest() %>%
  readr::write_rds(file.path(data_path,"pan33_ImmuneCellAI_cor_geneCNV.rds.gz"),compress = "gz")
