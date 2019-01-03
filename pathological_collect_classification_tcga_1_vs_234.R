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

read_pathways <- function(name) {
  symbols_lines <- read.table(sprintf("../msigdb/%s.gmt", name), header = FALSE, sep = ",", stringsAsFactor = FALSE)
  pathways <- vector("list", nrow(symbols_lines))
  for (line in 1:nrow(symbols_lines)) {
    symbols_entries <- strsplit(symbols_lines[line, 1], "\t")
    pathways[[line]]$name <- symbols_entries[[1]][1]
    pathways[[line]]$link <- symbols_entries[[1]][2]
    pathways[[line]]$symbols <- sort(symbols_entries[[1]][-2:-1])
  }
  return(pathways)
}

cohorts <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH",
             "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
             "TCGA-PAAD", "TCGA-READ", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA")

result_path <- "./results_1_vs_234"
replication_count <- 100
methods <- c("KNN", "LR","MLP")

perfs <- array(0, c(length(cohorts), length(methods), replication_count))
rownames(perfs) <- cohorts
colnames(perfs) <- methods


for (cohort in cohorts) {
  
  
  aurocs_knn_hallmark <- rep(NA, replication_count)
  for (replication in 1:replication_count) {
    load(file = sprintf("%s/%s/knn_%s_measure_AUROC_replication_%d_result.RData", result_path, cohort, "none", replication))
    aurocs_knn_hallmark[replication] <- result$AUROC
  }

  aurocs_logistic_regression_hallmark <- rep(NA, replication_count)
  for (replication in 1:replication_count) {
    load(file = sprintf("%s/%s/logistic_%s_measure_AUROC_replication_%d_result.RData", result_path, cohort, "none", replication))
    aurocs_logistic_regression_hallmark[replication] <- result$AUROC
  }
  aurocs_mlp_hallmark <- rep(NA, replication_count)
  for (replication in 1:replication_count) {
    load(file = sprintf("%s/%s/multilayer_perceptron_%s_measure_AUROC_replication_%d_result.RData", result_path, cohort, "none", replication))
    aurocs_mlp_hallmark[replication] <- result$AUROC
  }

  perfs[cohort, "KNN",] <- aurocs_knn_hallmark
  perfs[cohort, "LR",] <- aurocs_logistic_regression_hallmark
  perfs[cohort, "MLP",] <- aurocs_mlp_hallmark
 
}

colMeans(apply(perfs, MARGIN = c(1, 2), FUN = mean))

methods <- c("KNN","LR","MLP")

method_colors <- c("#cb181d", "#238b45","#2171b5")
pdf(file = "pathological_performance_comparison_1_vs_234.pdf", width = 8.5, height = 5.1)
par(oma = c(4, 1, 0, 0.5), cex = 1, lwd = 1, cex.axis = 1.0, cex.main = 1.25)
layout(matrix(1:15, 3, 5, byrow = TRUE))
for (cohort in cohorts) {
  cohort_data <- t(perfs[cohort, methods,])
  limits <- c(min(cohort_data, 0.5), max(cohort_data))
  limits <- c(floor(10 * limits[1]), ceiling(10 * limits[2]))
  boxplot(cohort_data, notch = TRUE, outline = FALSE, horizontal = FALSE,
          xaxt = "n", yaxt = "n", main = cohort, las = 1,
          ylim = c(limits[1] / 10, limits[2] / 10 + (limits[2] - limits[1]) / 10 * 0.30), ylab = "", par = par(mar = c(0, 3, 2, 0)))
  for (method in 2:length(methods)) {
    test_result <- t.test(cohort_data[,1], cohort_data[,method], alternative = "two.sided", paired = TRUE)
    lines(c(1, method), c(limits[2] / 10 + (2 * method - 3) * (limits[2] - limits[1]) / 10 * 0.07, limits[2] / 10 + (2 * method - 3) * (limits[2] - limits[1]) / 10 * 0.07), type = "l")
    lines(c(1, 1), c(limits[2] / 10 + (2 * method - 3) * (limits[2] - limits[1]) / 10 * 0.07 - (limits[2] - limits[1]) / 10 * 0.02, limits[2] / 10 + (2 * method - 3) * (limits[2] - limits[1]) / 10 * 0.07), type = "l")
    lines(c(method, method), c(limits[2] / 10 + (2 * method - 3) * (limits[2] - limits[1]) / 10 * 0.07 - (limits[2] - limits[1]) / 10 * 0.02, limits[2] / 10 + (2 * method - 3) * (limits[2] - limits[1]) / 10 * 0.07), type = "l")
    text((1 + method) / 2, y = limits[2] / 10 + (method - 1) * (limits[2] - limits[1]) / 10 * 0.14, labels = ifelse(test_result$p.value < 1e-3, "p < 1e-3", sprintf("p = %0.3f", test_result$p.value)), col = ifelse(test_result$p.value < 0.05, ifelse(test_result$statistic > 0, method_colors[1], method_colors[method]), "black"), font = ifelse(test_result$p.value < 0.05, 2, 1))
  }
  abline(h = 0.5, lty = 2, col = "orange", lwd = 1.5)
  for (m in 1:length(methods)) {
    method <- methods[m]
    points(x = m + 0.10 * rnorm(replication_count), cex = 0.625,
           y = perfs[cohort, method,], pch = 19, col = method_colors[m]) 
  }
  if (cohort >= "TCGA-PAAD") {
    axis(side = 1, at = 1:length(methods), labels = methods, las = 2)
  }
  if (cohort == "TCGA-BRCA" | cohort == "TCGA-KIRC" | cohort == "TCGA-PAAD") {
    mtext(side = 2, at = mean(limits) / 10, text = "AUROC", line = 2.5, cex = 0.8)
  }
  axis(side = 2, at = (limits[1]:limits[2]) / 10, labels = sprintf("%.1f", (limits[1]:limits[2]) / 10), las = 1)
  dev.next()
}

combined_result <- data.frame()
for (c in 1:length(cohorts)) {
  combined_result <- rbind(combined_result, as.data.frame(cbind(cohorts[c], 1:replication_count, t(perfs[c,,]))))
}
colnames(combined_result)[1:2] <- c("Cohort", "Replication")
write.csv(combined_result, file = "pathological_combined_result_1_vs_234.csv", row.names = FALSE)
dev.off()
M <- apply(perfs, c(1, 2), mean)
D <- apply(perfs, c(1, 2), sd)
R <- cbind(cohorts, matrix(paste0(sprintf("%.3f", M), "\\pm", sprintf("%.3f", D)), length(cohorts), 2))
colnames(R) <- c("Cohort", "KNN", "LR","MLP")

paste0(sapply(1:nrow(R), function(r) paste0(R[r,], collapse = " & ")), collapse = " \\ ")
