library("RGCCA")
library("integration")
library("inteRmodel")
library("ggplot2")
library("dplyr")
library("patchwork")

A <- readRDS("data_out/RGCCA_uncoupling_data.RDS")
meta <- A$meta


A <- A[1:2]
A <- lapply(A, t)
A <- clean_unvariable(A)

y <- model_RGCCA(meta, c("Gender", "Type", "state", "Location"))
A$y <- y

(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))

taus <- lapply(min_shrinkage, seq, to = 1, length.out = 10)
# ts <- lapply(A, tau.estimate) # On the computer
ts <- c(RNAseq = 0.486223918802408, micro = 0.938776330169865, y = 0.382540494844239)
taus[["RNAseq"]] <- c(taus[["RNAseq"]], ts["RNAseq"])
taus[["micro"]] <- c(taus[["micro"]], ts["micro"])
taus.combn <- expand.grid(taus)

C <-  matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), ncol = 3, nrow = 3,
             dimnames = list(names(A), names(A)))

sgcca_eval <- function(x, scheme) {
  result.sgcca <- sgcca(A, C, c1 = x, scheme = scheme, verbose = FALSE)
  
  analyze(result.sgcca)
}

# Centroid ####
out <- apply(taus.combn, 1, sgcca_eval, scheme = "centroid")
centroid <- cbind.data.frame(t(out), taus.combn)
saveRDS(centroid, "data_out/centroid_tau.RDS")


plot_tau2AVE <- function(data, title) {
  limits <- c(0.4, 0.7)
  p1 <- data %>%
    group_by(GE, y) %>%
    summarize(AVE_inner = median(AVE_inner)) %>%
    ungroup() %>%
    ggplot() +
    geom_point(aes(GE, y, col = AVE_inner)) +
    scale_color_viridis_c(limits = limits) +
    theme_minimal() +
    labs(col = "inner AVE") +
    scale_y_continuous(position = "right") +
    theme(axis.title.y.right = element_text(angle = 0, vjust = 0.5, hjust = 1))
  
  p2 <- data %>%
    group_by(CGH, y) %>%
    summarize(AVE_inner = median(AVE_inner)) %>%
    ggplot() +
    geom_point(aes(CGH, y, col = AVE_inner)) +
    scale_color_viridis_c(limits = limits) +
    theme_minimal() +
    labs(col = "inner AVE") +
    theme(axis.title.y = element_blank())
  p1 + p2 +  plot_layout(guides = 'collect')  +
    plot_annotation(tag_levels = "A", title = title)
  
}

plot_tau2AVE(centroid, "Inner AVE depending on tau")

#  factorial ####
out <- apply(taus.combn, 1, sgcca_eval, scheme = "factorial")
factorial <- cbind.data.frame(t(out), taus.combn)
saveRDS(factorial, "data_out/factorial_tau.RDS")

plot_tau2AVE(factorial, "Factorial scheme")

# horst ####
out <- apply(taus.combn, 1, sgcca_eval, scheme = "horst")
horst <- cbind.data.frame(t(out), taus.combn)
saveRDS(horst, "data_out/horst_tau.RDS")

plot_tau2AVE(horst, "Horst scheme")