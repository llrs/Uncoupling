library("omicade4")
library("ggplot2")


A <- readRDS("data_out/RGCCA_uncoupling_data.RDS")
meta <- A$meta

out <- mcia(A[1:2])
data_plot <- cbind(out$mcoa$SynVar, meta)
ggplot(data_plot) +
  geom_point(aes(SynVar1, SynVar2, col = state, shape = Location), size = 5) +
  scale_color_brewer(type = "qual") +
  theme_minimal() +
  labs(title = "MCIA on data from Uncoupling...") 
ggsave("Figures/mcia_uncoupling.png")
