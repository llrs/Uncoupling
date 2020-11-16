library("RGCCA")
library("ggplot2")
library("inteRmodel")
library("integration")



A <- readRDS("data_out/RGCCA_uncoupling_data.RDS")
meta <- A$meta

# Prepare the data for the models ####
At <- lapply(A[1:2], t)
At2 <- clean_unvariable(At)

At2[[3]] = At2[[3]]
B = At2
B[[3]] = Reduce("cbind", At2[1:(length(At2)-1)])
# superblock
B[[4]] = A[[3]]
names(B) <- c("GE", "CGH", "Superblock", "Loc")
# Prepare some preliminary design/weights
C <-  matrix(0,ncol = 4, nrow = 4)
dimnames(C) <- list(names(B), names(B))


nb_block <- 4
#length(B)
D <- matrix(0, nb_block, nb_block)
D[, nb_block - 1] <- 1
D[nb_block - 1, ] <- 1
diag(D) <- 0
shrinkage <- c(RNAseq = 0.486223918802408, micro = 0.938776330169865, y = 0.382540494844239, 1)
sgcca.glioma <- sgcca(B, D, c1 = shrinkage,
                     ncomp = c(rep(2, (length() - 1)),1),
                     scheme = "centroid", scale = TRUE,
                     verbose = FALSE)
saveRDS(sgcca.glioma, "superblock_sgcca.RDS")
glioma.superblock <- cbind.data.frame(sgcca.glioma$Y[[3]], Loc, Samples = rownames(sgcca.glioma$Y[[3]]))
glioma.GE <- cbind.data.frame(sgcca.glioma$Y[[1]], Loc, Samples = rownames(sgcca.glioma$Y[[3]]))
ggplot(glioma.superblock) +
  geom_text(aes(comp1, comp2, col = Loc, label = Samples)) +
  ggtitle("Original design", subtitle = "Superblock")
ggplot(glioma.GE) +
  geom_text(aes(comp1, comp2, col = Loc, label = Samples)) +
  ggtitle("Original design", subtitle = "GE")
sgcca.glioma$AVE$AVE_inner

nweights <- 4
weight <- seq(from = 0, to = 1, length.out = nweights)
C_list2 <- weight_design(nweights, 4)
C_list3 <- c(C_list2, 
             lapply(C_list2, function(x){x[1, 1] <- weight[2];x}),
             lapply(C_list2, function(x){x[1, 1] <- weight[3];x}),
             lapply(C_list2, function(x){x[1, 1] <- weight[4];x}))
keep <- check_design(C_list3)
C_list3 <- C_list3[keep]

keep <- vapply(C_list3, correct, logical(1L))
design <- C_list3[keep]
testing <- function(x, type) {
  result.sgcca <- RGCCA::sgcca(B, x, 
                               c1 = c(.071, .2, 1/sqrt(NCOL(B[[3]])) + .2, 1), 
                               ncomp = c(1, 1, 1), 
                               scheme = type, verbose = FALSE)
  
  out <- analyze(result.sgcca)
  c(out, "var11" = result.sgcca$C[1, 1])
}


out <- sapply(design, testing, type = "centroid", USE.NAMES = FALSE)
out2 <- t(out)
def <- as.data.frame(out2)
def$weights <- as.factor(def$weights)
offset <- is.na(def$AVE_inner)
def <- droplevels(def[!offset, ])
saveRDS(def, file = "data_out/superblock_interactions.RDS")