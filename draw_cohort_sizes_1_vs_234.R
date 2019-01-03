cohort_sizes <- c("BRCA" = 1067,
                  "COAD" = 443,
                  "ESCA" = 142,
                  "HNSC" = 429,
                  "KICH" = 65,
                  "KIRC" = 527,
                  "KIRP" = 259,
                  "LIHC" = 347,
                  "LUAD" = 505,
                  "LUSC" = 497,
                  "PAAD" = 174,
                  "READ" = 156,
                  "STAD" = 352,
                  "TGCT" = 81,
                  "THCA" = 500)

cohort_line_colors <- list("ACC" = c(255, 165, 0),
                           "BLCA" = c(128, 0, 128),
                           "BRCA" = c(255, 192, 203),
                           "CESC" = c(0, 128, 128),
                           "CHOL" = c(76, 187, 23),
                           "COAD" = c(0, 0, 139),
                           "DLBC" = c(191, 255, 0),
                           "ESCA" = c(204, 204, 255),
                           "GBM" = c(128, 128, 128),
                           "HNSC" = c(140, 0, 26),
                           "KICH" = c(255, 165, 0),
                           "KIRC" = c(255, 165, 0),
                           "KIRP" = c(255, 165, 0),
                           "LAML" = c(255, 165, 0),
                           "LGG" = c(128, 128, 128),
                           "LIHC" = c(80, 200, 120),
                           "LUAD" = c(0, 0, 0),
                           "LUSC" = c(0, 0, 0),
                           "MESO" = c(0, 0, 0),
                           "OV" = c(0, 128, 128),
                           "PAAD" = c(128, 0, 128),
                           "PCPG" = c(140, 0, 26),
                           "PRAD" = c(173, 216, 230),
                           "READ" = c(0, 0, 139),
                           "SARC" = c(255, 255, 0),
                           "SKCM" = c(0, 0, 0),
                           "STAD" = c(204, 204, 255),
                           "TGCT" = c(186, 85, 211),
                           "THCA" = c(0, 128, 128),
                           "THYM" = c(191, 255, 0),
                           "UCEC" = c(255, 218, 185),
                           "UCS" = c(255, 218, 185),
                           "UVM" = c(0, 0, 0))

cohort_inner_colors <- list("BRCA" = c(255, 192, 203),
                              "COAD" = c(0, 0, 139),
                              "ESCA" = c(204, 204, 255),
                              "HNSC" = c(255, 255, 240),
                              "KICH" = c(255, 165, 0),
                              "KIRC" = c(255, 165, 0),
                              "KIRP" = c(255, 165, 0),
                              "LIHC" = c(80, 200, 120),
                              "LUAD" = c(255, 255, 255),
                              "LUSC" = c(255, 255, 255),
                              "PAAD" = c(128, 0, 128),
                              "READ" = c(0, 0, 139),
                              "STAD" = c(204, 204, 255),
                              "TGCT" = c(186, 85, 211),
                              "THCA" = c(255, 192, 203),
                              "DLBC" = c(191, 255, 0),
                              "SARC" = c(255, 255, 0),
                              "UCEC" = c(255, 218, 185),
                              "CESC" = c(255, 255, 255),
                              "OV" = c(0, 128, 128),
                              "PRAD" = c(173, 216, 230),
                              "SKCM" = c(0, 0, 0),
                              "LGG" = c(128, 128, 128),
                              "GBM" = c(128, 128, 128),
                              "LAML" = c(255, 165, 0),
                              "UCS" = c(255, 218, 185),
                              "UVM" = c(0, 0, 0),
                              "THYM" = c(191, 255, 0),
                              "CHOL" = c(76, 187, 23),
                              "BLCA" = c(234, 162, 33),
                              "ACC" = c(255, 165, 0),
                              "PCPG" = c(255, 255, 240),
                              "MESO" = c(255, 255, 255))

cohort_names <- sort(names(cohort_sizes))
cohort_sizes <- cohort_sizes[cohort_names]
cohort_line_colors <- cohort_line_colors[cohort_names]
cohort_inner_colors <- cohort_inner_colors[cohort_names]

library(rgdal)
library(png)

icon <- readPNG("standing-man.png")
icon <- icon[21:520, 161:380, 1:4]
icon[,,4] <- (icon[,,4] > 0)

inner_region <- matrix(FALSE, nrow(icon), ncol(icon))
for (i in 2:480) {
  for (j in 2:(ncol(icon) - 1)) {
    if (sum(icon[i, 1:(j - 1) , 4]) > 0 & sum(icon[i, (j + 1):ncol(icon), 4]) > 0 & icon[i, j, 4] == 0) {
      inner_region[i, j] <- TRUE
    }
  }
}

aspect_ratio <- nrow(icon) / ncol(icon)

gap_width <- 5
total_gap <- (length(cohort_sizes) - 1) * gap_width
total_width <- sum(sqrt(cohort_sizes)) + total_gap
total_height <- max(sqrt(cohort_sizes)) * aspect_ratio


pdf(file = "cohort_sizes_1_vs_234.pdf", width = 12.0, height = 2.0)
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
plot(0, 0, col = "white", xlim = c(0, total_width), ylim = c(0, total_height * 1.2), axes = FALSE)
for (cohort in cohort_names) {
  icon <- readPNG("standing-man.png")
  icon <- icon[21:520, 161:380, 1:4]
  icon[,,4] <- (icon[,,4] > 0)

  icon[,,1] <- cohort_line_colors[[cohort]][1] / 255
  icon[,,2] <- cohort_line_colors[[cohort]][2] / 255
  icon[,,3] <- cohort_line_colors[[cohort]][3] / 255

  temp1 <- icon[,,1]
  temp1[inner_region] <- cohort_inner_colors[[cohort]][1] / 255
  temp2 <- icon[,,2]
  temp2[inner_region] <- cohort_inner_colors[[cohort]][2] / 255
  temp3 <- icon[,,3]
  temp3[inner_region] <- cohort_inner_colors[[cohort]][3] / 255
  temp4 <- icon[,,4]
  temp4[inner_region] <- 1
  icon[,,1] <- temp1
  icon[,,2] <- temp2
  icon[,,3] <- temp3
  icon[,,4] <- temp4

  c <- which(cohort_names %in% cohort)
  rasterImage(icon, xleft = sum(sqrt(cohort_sizes[1:c])) - sqrt(cohort_sizes[c]) + (c - 1) * gap_width,
                    xright = sum(sqrt(cohort_sizes[1:c])) + (c - 1) * gap_width,
                    ybottom = 0 + max(sqrt(cohort_sizes)) * aspect_ratio * 0.125,
                    ytop = sqrt(cohort_sizes[c]) * aspect_ratio + max(sqrt(cohort_sizes)) * aspect_ratio * 0.125)
  text(x = sum(sqrt(cohort_sizes[1:c])) + (c - 1) * gap_width - 0.5 * sqrt(cohort_sizes[c]), y = 6, labels = cohort, cex = 0.8)
  text(x = sum(sqrt(cohort_sizes[1:c])) + (c - 1) * gap_width - 0.5 * sqrt(cohort_sizes[c]), y = 0, labels = sprintf("(%d)", cohort_sizes[c]), cex = 0.8)
}
dev.off()
system("pdfcrop cohort_sizes_1_vs_234.pdf cohort_sizes_1_vs_234.pdf")
