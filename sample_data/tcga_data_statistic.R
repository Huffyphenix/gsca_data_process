
# packages prepare --------------------------------------------------------

library(magrittr)
library(ggplot2)

# load data ---------------------------------------------------------------

tcga_sample <- readr::read_tsv("S:\\study\\GSCALite\\data\\summary_data\\projects-table.tsv") %>%
  dplyr::select(Project,Cases,Seq,Exp,SNV,CNV,Meth,Clinical) %>%
  tidyr::gather(-Project,-Cases,key="Data_type",value="Counts")


# path  -------------------------------------------------------------------

out_path <- c("S:/study/GSCALite/data/summary_data")

# plot --------------------------------------------------------------------
tcga_sample %>%
  tidyr::spread(Data_type, Counts) %>%
  dplyr::mutate(sum=Clinical+CNV+Exp+Meth+Seq+SNV) %>%
  tidyr::gather(-Project,-Cases,-sum,key="Data_type",value="Counts") %>%
  dplyr::mutate(Percent=Counts/sum) %>%
  dplyr::mutate(`Cancer Type`=paste(Project," (n=", Cases, ")", sep = "")) %>%
  dplyr::mutate(ypos=ifelse(Data_type=="Clinical",0.9,0.1)) %>%
  dplyr::mutate(ypos=ifelse(Data_type=="CNV",0.75,ypos)) %>%
  dplyr::mutate(ypos=ifelse(Data_type=="Exp",0.6,ypos)) %>%
  dplyr::mutate(ypos=ifelse(Data_type=="Meth",0.42,ypos)) %>%
  dplyr::mutate(ypos=ifelse(Data_type=="Seq",0.25,ypos)) %>%
  dplyr::mutate(ypos=ifelse(Data_type=="SNV",0.05,ypos)) -> plot_ready

plot_ready %>%
  ggplot(mapping= aes(x=`Cancer Type`,y=Percent,fill=Data_type)) +
  geom_col(position = "stack")  +
  geom_text(aes(label = Counts,y=ypos),size = 3) +
  theme(
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  coord_flip()->p;p

ggsave(file.path(out_path,"tcga_sample_info.pdf"),device = "pdf")


# circle plot -------------------------------------------------------------

data.frame(Project=as.character(plot_ready$Project[1:33]),
           Cases=plot_ready$Cases[1:33],
           sum=plot_ready$sum[1:33],
           Data_type=as.character("SS"),
           Counts=plot_ready$Counts[1:33],
           Percent=0.2,
           `Cancer Type`=as.character(plot_ready$`Cancer Type`[1:33]),
           ypos=0) %>%
  dplyr::as_tibble() %>%
  dplyr::rename("Cancer Type"="Cancer.Type")-> white.plot
white.plot$Project <- as.character(white.plot$Project)
white.plot$Data_type <- as.character(white.plot$Data_type)
white.plot$`Cancer Type` <- as.character(white.plot$`Cancer Type`)

text <- rbind(plot_ready, white.plot) %>%
  dplyr::select(Project,Cases) %>%
  unique()
rbind(plot_ready, white.plot) -> circle_plot_ready
circle_plot_ready$Cases[34:231] <- NA


circle_plot_ready %>%
  ggplot(mapping= aes(x=Project,y=Percent,fill=Data_type)) +
  geom_col(position = "stack")  +
  geom_text(aes(label = Cases),y=1.2,size = 3) +
  scale_fill_brewer(palette = "Set3") +
  theme(
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  coord_polar("x")->p;p
ggsave(file.path(out_path,"tcga_sample_info-circle.pdf"),device = "pdf")
