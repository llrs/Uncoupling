library("ggplot2")
library("UpSetR")
library("integration")

A <- readRDS("data_out/RGCCA_uncoupling_data.RDS")
meta <- A$meta

model0 <- readRDS("data_out/model0.RDS")
model1.1 <- readRDS("data_out/model1.1.RDS")
model1.2 <- readRDS("data_out/model1.2.RDS")
model2.1 <- readRDS("data_out/model2.1.RDS")
model2.2 <- readRDS("data_out/model2.2.RDS")
# superblock
# superblock.2

models <- list("0" = model0, "1.1" = model1.1, "1.2" = model1.2, 
               "2.1" = model2.1, "2.2" = model2.2)

m0GE <- tidyer(model0$Y[[1]], "0", "GE")
m1.1GE <- tidyer(model1.1$Y[[1]], "1.1", "GE")
m1.2GE <- tidyer(model1.2$Y[[1]], "1.2", "GE")
m2.1GE <- tidyer(model2.1$Y[[1]], "2.1", "GE")
m2.2GE <- tidyer(model2.2$Y[[1]], "2.2", "GE")
# msGE <- tidyer(superblock$Y[[1]], "superblock", "GE")
# ms.2GE <- tidyer(superblock.2$Y[[1]], "superblock.2", "GE")

m0M <- tidyer(model0$Y[[2]], "0", "M")
m1.1M <- tidyer(model1.1$Y[[2]], "1.1", "M")
m1.2M <- tidyer(model1.2$Y[[2]], "1.2", "M")
m2.1M <- tidyer(model2.1$Y[[2]], "2.1", "M")
m2.2M <- tidyer(model2.2$Y[[2]], "2.2", "M")
# msM <- tidyer(superblock$Y[[2]], "superblock", "M")
# ms.2M <- tidyer(superblock.2$Y[[2]], "superblock.2", "M")



df <- rbind(
  merge(m0GE, m0M, all = TRUE),
  merge(m1.1GE, m1.1M, all = TRUE),
  merge(m1.2GE, m1.2M, all = TRUE),
  merge(m2.1GE, m2.1M, all = TRUE),
  merge(m2.2GE, m2.2M, all = TRUE)
  
)

df <- merge(df, meta, by.x = "Rownames", by.y = "ID", all = TRUE)
saveRDS(df, "data_out/models_summary.RDS")

theme_set(theme_bw())
theme_update(strip.background = element_blank(), 
             panel.grid.minor = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank())
df %>% 
  filter(!grepl(" i$", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = Type, shape = Type)) +
  facet_wrap(~Model, scales = "free", nrow = 2) + 
  labs(title = "Samples by model", x = "Transcriptome", y = "Microbiome")  +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.25)) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.25))
df %>% 
  filter(!grepl(" i$", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = Location, shape = Location)) +
  facet_wrap(~Model, scales = "free", nrow = 2) + 
  labs(title = "Samples by model", x = "Transcriptome", y = "Microbiome") +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.25)) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.25))
