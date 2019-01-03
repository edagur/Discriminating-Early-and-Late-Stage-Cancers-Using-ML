library(grid)
library(pheatmap)

threshold <- 0.01

draw_colnames_90 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_90",
                  ns=asNamespace("pheatmap"))

result_path <- "./results_1_vs_234"
replication_count <- 100
methods <- c("KNN", "LR","MLP")
weights_list <- rep(NA,replication_count)
data_path <- "./primary_data"

cohort <- "TCGA-TGCT"
load(sprintf("%s/%s.RData", data_path, cohort))
TCGA$mrna <- TCGA$mrna[,unique(colnames(TCGA$mrna))]
common_patients <- intersect(rownames(TCGA$clinical)[which(is.na(TCGA$clinical$pathologic_stage) == FALSE)], rownames(TCGA$mrna))

X <- log2(TCGA$mrna[common_patients,] + 1)
y <- rep(NA, length(common_patients))

y[TCGA$clinical[common_patients, "pathologic_stage"] %in% c("Stage I",  "Stage IA",  "Stage IB",  "Stage IC")] <- +1
y[TCGA$clinical[common_patients, "pathologic_stage"] %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC",
                                                            "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                                            "Stage IV",  "Stage IVA",  "Stage IVB",  "Stage IVC")] <- -1

valid_patients <- which(is.na(y) == FALSE)
valid_features <- as.numeric(which(apply(X[valid_patients,], 2, sd) != 0))
X <- X[valid_patients, valid_features]
y <- y[valid_patients]

negative_indices <- which(y == -1)
positive_indices <- which(y == +1)

alpha_set <- seq(from = 0, to = 1, by = 0.1)
fold_count <- 4
train_ratio <- 0.8

if (pathway != "none") {
  pathways <- read_pathways(pathway)
  gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
  X <- X[, which(colnames(X) %in% gene_names)]
}

set.seed(1606 * replication)
train_negative_indices <- sample(negative_indices, ceiling(train_ratio * length(negative_indices)))
train_positive_indices <- sample(positive_indices, ceiling(train_ratio * length(positive_indices)))

auroc_matrix <- matrix(NA, nrow = fold_count, ncol = length(alpha_set), dimnames = list(1:fold_count, sprintf("%g", alpha_set)))

negative_allocation <- sample(rep(1:fold_count, ceiling(length(train_negative_indices) / fold_count)), length(train_negative_indices))
positive_allocation <- sample(rep(1:fold_count, ceiling(length(train_positive_indices) / fold_count)), length(train_positive_indices))

train_indices <- c(train_negative_indices, train_positive_indices)
test_indices <- setdiff(1:length(y), train_indices)

X_train <- X[train_indices,]
X_test <- X[test_indices,]
X_train <- scale(X_train)
X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
X_train[is.na(X_train)] <- 0
X_test[is.na(X_test)] <- 0
invalid_indices <- which(apply(X_train, MARGIN = 2, sd) == 0)
counts <- rep(0,length(weights))
names(counts) <- colnames(X_train)
  for (replication in 1:replication_count) {
    load(file = sprintf("%s/%s/logistic_regression_%s_weights_%d_result.RData", result_path, cohort, "none", replication))
    
    counts[order(weights, decreasing = TRUE)[1:250]] <- counts[order(weights, decreasing = TRUE)[1:250]] + 1
    counts[order(weights, decreasing = FALSE)[1:250]] <- counts[order(weights, decreasing = FALSE)[1:250]] + 1
  }

#picking the best 5 weights for this cancer type
counts[order(counts, decreasing = TRUE)[1:5]]

  