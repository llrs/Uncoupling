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
superblock
superblock.2


m0GE <- tidyer(model0$Y[[1]], "0", "GE")