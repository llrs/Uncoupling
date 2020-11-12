library("inteRmodel")
library("integration")
library("RGCCA")
library("fastDummies")

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

#  Metadata to dummy variables ####
convert2dummy <- function(meta, names) {
  meta2 <- dummy_cols(meta[, names, drop = FALSE], names, 
                      remove_first_dummy = TRUE)
  new_names <- colnames(meta2)[!colnames(meta2) %in% names]
  meta2[, new_names, drop = FALSE]
}
Time <- meta[, c("AgeAtDateOfSampling"), drop = FALSE]
Demographic <- convert2dummy(meta, c("Gender", "Type", "state"))
Location <- convert2dummy(meta, c("Location"))
meta2 <- convert2dummy(meta, c("Gender", "Type", "state", "Location"))

# Models 1 ####
designs <- weight_design(weights = 11, size = 3)
keep <- vapply(designs, RGCCA::correct, logical(1L))
designsc <- designs[keep]

A1 <- At2
A1$meta <- cbind(meta2, meta$AgeAtDateOfSampling)
A1b <- lapply(A1, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))

# Search all models
out_model <- search_model(A = A1b, c1 = c(shrinkage, 1), scheme = "centroid",
                          scale = FALSE, verbose = FALSE,
                          ncomp = rep(1, length(A1b)),
                          bias = TRUE, 
                          nWeights = 11) # To end up with .1 intervals
saveRDS(out_model, "data_out/uncoupling_models1.RDS")

# Analyze the best model in deep:

# Models 2 ####
A2 <- A1[1:2]
A2$Location <- Location
A2$Demographic <- Demographic
A2$Time <- Time

A2b <- lapply(A2, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))

out_model <- search_model(A = A2b, c1 = c(shrinkage, 1, 1, 1), scheme = "centroid",
                          scale = FALSE, verbose = FALSE,
                          ncomp = rep(1, length(A2b)),
                          bias = TRUE, 
                          nWeights = 11)

# Search all models
saveRDS(out_model, "data_out/uncoupling_models2.RDS")

# Analyze the best model in deep: