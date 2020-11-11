library("RGCCA")
library("ggplot2")
library("patchwork")
library("integration")
theme_set(theme_minimal())

# Prepare data ####
A <- readRDS("data_out/RGCCA_uncoupling_data.RDS")
meta <- A$meta


A <- A[1:2]
A <- lapply(A, t)
A <- clean_unvariable(A)

y <- model_RGCCA(meta, c("Gender", "Type", "state", "Location"))
A$y <- y

A0 <- A
A <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))

# Prepare designs
nweight <- 11
C_list2 <- weight_design(11, 3)
keep <- vapply(C_list2, correct, logical(1L))
design <- C_list2[keep]
# Design matrix:
# 0.0 var12 var13
# var12 0.0  var23
# var13 var23  0.0

shrinkage <- c(RNAseq = 0.486223918802408, micro = 0.938776330169865, y = 1)

testing <- function(x, type) {
  result.sgcca <- RGCCA::sgcca(A, x, c1 = shrinkage, ncomp = c(1, 1, 1), 
                               scheme = type, verbose = FALSE, scale = FALSE)
  
  analyze(result.sgcca)
}

# Schemes

# * Centroid  ####
out <- sapply(design, testing, type = "centroid", USE.NAMES = FALSE)
out2 <- t(out)
def <- as.data.frame(out2)
def$weights <- as.factor(def$weights)
offset <- is.na(def$AVE_inner)
centroid_weights <- droplevels(def[!offset, ])
saveRDS(centroid_weights, file = "centroid_weights.RDS")