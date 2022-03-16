library("omicade4")
library("ggplot2")


A <- readRDS("data_out/RGCCA_uncoupling_data.RDS")
meta <- A$meta
theme_update(strip.background = element_blank(), 
             panel.grid.minor = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank())
out <- mcia(A[1:2])
data_plot <- cbind(out$mcoa$SynVar, meta)
ggplot(data_plot) +
  geom_point(aes(SynVar1, SynVar2, col = state, shape = Location), size = 5) +
  scale_color_brewer(type = "qual") +
  theme_minimal() +
  labs(title = 'MCIA on data from "Uncoupling..."', col = "State")
ggsave("Figures/mcia_uncoupling.png")
ggplot(data_plot) +
  geom_point(aes(SynVar1, SynVar2, col = Type, shape = Location), size = 5) +
  scale_color_brewer(type = "qual") +
  theme_minimal() +
  labs(title = 'MCIA on data from "Uncoupling..."', col = "Diagnosis")
ggsave("Figures/mcia_uncoupling_type.png")

library("pROC")
roc_loc <- multiclass.roc(response = meta$Location,
               predictor = data_plot$SynVar1,
               levels  = unique(meta$Location))
roc_type <- multiclass.roc(response = meta$Type,
               predictor = data_plot$SynVar2,
               levels  = unique(meta$Type))
