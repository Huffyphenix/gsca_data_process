# packages prepare --------------------------------------------------------

library(magrittr)
library(ggplot2)

# load data ---------------------------------------------------------------

drug_sample <- readr::read_tsv("S:\\study\\GSCALite\\data\\summary_data\\drug_log2.txt") %>%
  tidyr::gather(-Database,key="Data",value="Log2(Counts)")

drug_sample_raw <- readr::read_tsv("S:\\study\\GSCALite\\data\\summary_data\\drug.txt") %>%
  tidyr::gather(-Database,key="Data",value="Counts")

# path  -------------------------------------------------------------------

out_path <- c("S:/study/GSCALite/data/summary_data")

# plot --------------------------------------------------------------------

rank <- drug_sample %>%
  dplyr::group_by(Data) %>%
  dplyr::mutate(sum=sum(`Log2(Counts)`)) %>%
  dplyr::arrange(sum) %>%
  dplyr::ungroup() %>%
  dplyr::select(Data,sum) %>%
  unique()
text <- drug_sample_raw %>%
  dplyr::mutate(ypos=ifelse(Database=="CTRP",15,5))

drug_sample %>%
  dplyr::mutate(pos=ifelse(Database=="CTRP",))

drug_sample %>%
  ggplot(aes(x=Data,y = `Log2(Counts)`,fill=Database)) +
  geom_col() +
  scale_x_discrete(limits = rank$Data) +
  # scale_fill_brewer(palette = "Accent") +
  scale_fill_manual(values=RColorBrewer::brewer.pal(3,"Pastel2")[1:2]) +
  geom_text(label=text$Counts,y=text$ypos) +
  coord_polar("y") +
  theme(
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) 

ggsave(file.path(out_path,"drug_sample_info-circle.pdf"),device = "pdf")
