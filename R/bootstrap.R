library("inteRmodel")

# Prepare data ####
A <- readRDS("data_out/RGCCA_uncoupling_data.RDS")
meta <- A$meta

# Prepare index ####

boots <- 10
index <- vector("list", length = boots)
for (i in seq_len(boots)) {
  index[[i]] <- sample(nrow(meta), replace = TRUE)
}

# Prepare the data for the models ####
At <- lapply(A[1:2], t)
At2 <- clean_unvariable(At)
Ab <- lapply(At2, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))

# Tau estimate
# Calculated on the Juanjo's server
shrinkage <- c(RNAseq = 0.486223918802408, micro = 0.938776330169865)
# * Model 0 ####
m0 <- matrix(c(0, 1, 1, 0), ncol = 2, dimnames = list(names(Ab), names(Ab)))


# * Models 1 ####
A1 <- At2
meta2 <- model_RGCCA(meta, c("Gender", "Type", "state", "Location", "AgeAtDateOfSampling"))
A1$meta <- meta2
A1b <- lapply(A1, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
model1.2 <- readRDS("data_out/model1.2.RDS")
m1.2 <- model1.2$C

# * Models 2 ####
A2 <- A1[1:2]

# ** Metadata to dummy variables ####
Time <- model_RGCCA(meta, "AgeAtDateOfSampling")
Demographic <- model_RGCCA(meta, c("Gender", "Type", "state"))
Location <- model_RGCCA(meta, c("Location"))

A2$Location <- Location
A2$Demographic <- Demographic
A2$Time <- Time
A2b <- lapply(A2, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
model2.2 <- readRDS("data_out/model2.2.RDS")
m2.2 <- model2.2$C


index <- boot_index(nrow(meta), 1000)

boot0 <- boot_index_sgcca(index, A = Ab, C = m0, c1 = shrinkage, scale = FALSE)
saveRDS(boot0, "data_out/boot0.RDS")
boot1.2 <- boot_index_sgcca(index, A = A1b, C = m1.2, c1 = c(shrinkage, 1), scale = FALSE)
saveRDS(boot1.2, "data_out/boot1.2.RDS")
boot2.2 <- boot_index_sgcca(index, A = A2b, C = m2.2, c1 = c(shrinkage, 1, 1, 1), scale = FALSE)
saveRDS(boot2.2, "data_out/boot2.2.RDS")
