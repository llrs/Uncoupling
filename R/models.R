library("inteRmodel")
library("integration")
library("RGCCA")
library("ggplot2")
library("dplyr")

A <- readRDS("data_out/RGCCA_uncoupling_data.RDS")
meta <- A$meta

# Prepare the data for the models ####
At <- lapply(A[1:2], t)
At2 <- clean_unvariable(At)

# Tau estimate
# Calculated on the Juajo's server
shrinkage <- c(RNAseq = 0.486223918802408, micro = 0.938776330169865)
# shrinkage[[2]] <- tau.estimate(A[[2]]) 
# Model 0 ####
# Prepare the data:
Ab <- lapply(At, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
C <- matrix(c(0, 1, 1, 0), ncol = 2, dimnames = list(names(Ab), names(Ab)))
model0 <- sgcca(Ab, scale = FALSE, scheme = "centroid", C =  C, bias = TRUE,
                c1 = shrinkage)
model0 <- improve(model0, names(Ab))
saveRDS(model0, "data_out/model0.RDS")

#  Metadata to dummy variables ####
Time <- model_RGCCA(meta, "AgeAtDateOfSampling")
Demographic <- model_RGCCA(meta, c("Gender", "Type", "state"))
Location <- model_RGCCA(meta, c("Location"))
meta2 <- model_RGCCA(meta, c("Gender", "Type", "state", "Location", "AgeAtDateOfSampling"))

# Models 1 ####
A1 <- At2
A1$meta <- meta2
A1b <- lapply(A1, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
 
m1.1 <- diag(3)
diag(m1.1) <- 0
dimnames(m1.1) <- list(names(A1b), names(A1b))
m1.1 <- subSymm(m1.1, "RNAseq", "meta", 1)
m1.1 <- subSymm(m1.1, "micro", "meta", 1)

model1.1 <- sgcca(A1b, scale = FALSE, scheme = "centroid", C =  m1.1, bias = TRUE,
                c1 = c(shrinkage, 1))
model1.1 <- improve(model1.1, names(A1b))
saveRDS(model1.1, "data_out/model1.1.RDS")

# # Search all models
designs <- weight_design(weights = 11, size = 3)
keep <- vapply(designs, correct, logical(1L))
designsc <- designs[keep]
out_model <- search_model(A = A1b, c1 = c(shrinkage, 1), scheme = "centroid",
                          scale = FALSE, verbose = FALSE,
                          ncomp = rep(1, length(A1b)),
                          bias = TRUE,
                          nWeights = 11,
                          BPPARAM = BiocParallel::MulticoreParam(workers = 6)) # To end up with .1 intervals
saveRDS(out_model, "data_out/uncoupling_models1.RDS")
out_model <- readRDS("data_out/uncoupling_models1.RDS")

C <- diag(3)
diag(C) <- 0
columns <- grep("var", colnames(out_model))
m1.2 <- symm(C, out_model[which.max(out_model$AVE_inner), columns])
model1.2 <- sgcca(A1b, C = m1.2, scale = FALSE, c1 = c(shrinkage, 1))
model1.2 <- improve(model1.2, names(A1b))
saveRDS(model1.2, "data_out/model1.2.RDS")

model1.2_plot <- cbind(RNAseq = model1.2$Y[[1]], model1.2$Y[[2]],  
                       model1.2$Y[[3]], meta)
colnames(model1.2_plot)[1:2] <- c("RNAseq", "micro")
ggplot(model1.2_plot) +
  geom_point(aes(RNAseq, micro, col = Gender, shape = Gender))
ggplot(model1.2_plot) +
  geom_point(aes(RNAseq, micro, col = Type, shape = Type))
ggplot(model1.2_plot) +
  geom_point(aes(RNAseq, micro, col = Location, shape = Location))
ggplot(model1.2_plot) +
  geom_point(aes(micro, comp1, col = Type, shape = Location))

# Analyze the best model in deep:

# Models 2 ####
A2 <- A1[1:2]
A2$Location <- Location
A2$Demographic <- Demographic
A2$Time <- Time

A2b <- lapply(A2, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))

# Model 2.1 ####
C <- diag(5)
diag(C) <- 0
dimnames(C) <- list(names(A2), names(A2))
m2.1 <- subSymm(C, "RNAseq", "Location", 1)
m2.1 <- subSymm(m2.1, "micro", "Demographic", 1)
m2.1 <- subSymm(m2.1, "micro", "Location", 0.5)
m2.1 <- subSymm(m2.1, "Demographic", "Time", 1)
model2.1 <- sgcca(A = A2b, C = m2.1, scale = FALSE, c1 = c(shrinkage, 1, 1, 1))
model2.1 <- improve(model2.1, names(A2))
saveRDS(model2.1, "data_out/model2.1.RDS")

# * Prepare some random designs ####
designs <- weight_design(weights = 3, size = 5)
keep <- vapply(designs, integration::correct, logical(1L))
designs <- designs[keep]

# Random subsample of 10% of the trials
# Store all AVEs in the path so that it can be confirmed that it is the max value
set.seed(4672679)
s <- sample(designs, size = min(length(designs)*.1, 5000))
# Perform the sgcca on these samples
testing <- function(x, type, ...) {
  result.sgcca <- RGCCA::sgcca(C = x,
                               scheme = type,
                               verbose = FALSE,
                               scale = FALSE,
                               ...)
  analyze(result.sgcca)
}
# Estimated time of three days with designs and about 1 hour with the sample of 1000
out <- sapply(s, testing, type = "centroid", A = A2b, c1 = c(shrinkage, 1, 1, 1), USE.NAMES = FALSE)
out2 <- as.data.frame(t(out))
saveRDS(out2, "data_out/sample_model3_boot.RDS")
out2 <- readRDS("data_out/sample_model3_boot.RDS")


columns <- grep("var", colnames(out2))
model <- symm(C, out2[which.max(out2$AVE_inner), columns])
model <- subSymm(model, 4, 5, 1)
w <- which(lower.tri(model) & model != 0)
d <- weight_design(11, 5, w)

out3 <- sapply(d, testing, type = "centroid", A = A2b, c1 = c(shrinkage, 1, 1, 1), USE.NAMES = FALSE)
saveRDS(out3, "data_out/sample_model3_refined.RDS")

out3 <- readRDS("data_out/sample_model3_refined.RDS")
a <- vapply(out3, is, class = "numeric", logical(1L))
out4 <- simplify2array(out3[, a])

out <- rbind(out2, t(out3))



columns <- grep("var", colnames(out))
m2.2 <- symm(C, out[which.max(out$AVE_inner), columns])
dimnames(m2.2) <- list(names(A2b), names(A2b))
model2.2 <- sgcca(A = A2b, scale = FALSE, scheme =  "centroid", c1 = c(shrinkage, 1, 1, 1),  C = m2.2)
model2.2 <- improve(model2.2, names(A2b))
saveRDS(model2.2, "data_out/model2.2.RDS")

model2.2_Y <- cbind(RNAseq = model2.2$Y[[1]], model2.2$Y[[2]], 
                  model2.2$Y[[3]],
                  model2.2$Y[[4]],
                  model2.2$Y[[5]],
                  meta)
colnames(model2.2_Y)[1:5] <- colnames(model)
colnames(model2.2_Y)[3] <- "Loc"
ggplot(model2.2_Y) +
  geom_point(aes(RNAseq, micro, col = Gender, shape = Gender))
ggplot(model2.2_Y) +
  geom_point(aes(RNAseq, micro, col = Type, shape = Type))
ggplot(model2.2_Y) +
  geom_point(aes(RNAseq, micro, col = Location, shape = Location))
ggplot(model2.2_Y) +
  geom_point(aes(Demographic, Time, col = Type))

meta %>% 
  count(Location, Type, sort = TRUE)
# Not possible to run due to the number of combinations of the designs
# out_model <- search_model(A = A2b, c1 = c(shrinkage, 1, 1, 1), scheme = "centroid",
#                           scale = FALSE, verbose = FALSE,
#                           ncomp = rep(1, length(A2b)),
#                           bias = TRUE,
#                           nWeights = 11,
#                           BPPARAM = BiocParallel::MulticoreParam(workers = 6))

# Search all models
saveRDS(out_model, "data_out/uncoupling_models2.RDS")

# Analyze the best model in deep:


