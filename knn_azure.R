library(AUC)
library(class)

source("classification_helper.R")

cohorts <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH",
             "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
             "TCGA-PAAD", "TCGA-READ", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA")

result <- 0
cohort <- cohorts[12]
data_path <- "./primary_data"
result_path <- "./results_1_vs_234"
pathway <- "none" #replace with "hallmark" if you would like to restrict RF to Hallmark genes only


    #if (file.exists(sprintf("%s/%s/knn_%s_measure_AUROC_replication_%d_result.RData", result_path, cohort, pathway, replication)) == FALSE) {
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
      
      k_set <- seq(from = 1, to = 13, by = 2)
      fold_count <- 4
      train_ratio <- 0.8
      
      if (pathway != "none") {
        pathways <- read_pathways(pathway)
        gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
        X <- X[, which(colnames(X) %in% gene_names)]
      }
      
    
    each_rep <- function(replication){
      set.seed(1606 * replication)
      train_negative_indices <- sample(negative_indices, ceiling(train_ratio * length(negative_indices)))
      train_positive_indices <- sample(positive_indices, ceiling(train_ratio * length(positive_indices)))
      
      auroc_matrix <- matrix(NA, nrow = fold_count, ncol = length(k_set), dimnames = list(1:fold_count, sprintf("%g", k_set)))
      
      negative_allocation <- sample(rep(1:fold_count, ceiling(length(train_negative_indices) / fold_count)), length(train_negative_indices))
      positive_allocation <- sample(rep(1:fold_count, ceiling(length(train_positive_indices) / fold_count)), length(train_positive_indices))
      
      
      for (fold in 1:fold_count) {
        train_indices <- c(train_negative_indices[which(negative_allocation != fold)], train_positive_indices[which(positive_allocation != fold)])
        test_indices <- c(train_negative_indices[which(negative_allocation == fold)], train_positive_indices[which(positive_allocation == fold)])
        
        X_train <- X[train_indices,]
        X_test <- X[test_indices,]
        X_train <- scale(X_train)
        X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
        X_train[is.na(X_train)] <- 0
        X_test[is.na(X_test)] <- 0
        
        X_train[is.infinite(X_train)] <- 0
        X_test[is.infinite(X_test)] <- 0
        
        y_train <- y[train_indices]
        y_test <- y[test_indices]
        
        data <- as.data.frame(X_train)
        #data$label <- as.factor(y_train)
        mytest <- which(is.infinite(X_test), arr.ind = TRUE)
        mytest <- as.data.frame(mytest)
        if( length(mytest) == 0){
          y_test <- y_test[-mytest$row]
        }
        
        X_test <- X_test[is.finite(rowSums(X_test)),]
        X_test <-X_test[complete.cases(X_test),]
        X_test <- na.omit(X_test)
        
        X_train2 <- X_train[is.finite(rowSums(X_train)),]
        X_train2 <-X_train2[complete.cases(X_train2),]
        X_train2 <- na.omit(X_train2)
        
        
        y_train <- y_train[complete.cases(y_train)]
        y_train <- na.omit(y_train)
        
        
        for (k in k_set) {
          print(sprintf("running fold = %d, k = %g", fold, k))
          
          #prediction <- class::knn.cv(data,data$label,k=k,prob=FALSE, use.all = TRUE)
          #prediction <- class:knn(X_train2,X_test,y_train,k = k, l = 0, prob = FALSE, use.all = TRUE)
          prediction <- class::knn(X_train2,as.data.frame(X_test),y_train,k = k, l = 0, prob = FALSE, use.all = TRUE)
          
          auroc_matrix[fold, sprintf("%g", k)] <- AUC::auc(roc(prediction, as.factor(1 * (y_test == +1))))
        }
        
        
      }
      
      knn_star_AUROC <- k_set[max.col(t(colMeans(auroc_matrix, na.rm = TRUE)), ties.method = "last")]
      
      
      train_indices <- c(train_negative_indices, train_positive_indices)
      test_indices <- setdiff(1:length(y), train_indices)
      
      X_train <- X[train_indices,]
      X_test <- X[test_indices,]
      X_train <- scale(X_train)
      X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
      X_train[is.na(X_train)] <- 0
      X_test[is.na(X_test)] <- 0
      invalid_indices <- which(apply(X_train, MARGIN = 2, sd) == 0)
      X_train <- X_train[,-invalid_indices]
      X_test <- X_test[,-invalid_indices]
      
      X_train[is.infinite(X_train)] <- 0
      X_test[is.infinite(X_test)] <- 0
      y_train <- y[train_indices]
      y_test <- y[test_indices]
      
      data <- as.data.frame(X_train)
      #data$label <- as.factor(y_train)
      
     # X_test <- X_test[is.finite(rowSums(X_test)),]
      #X_test <-X_test[complete.cases(X_test),]
      #X_test <- na.omit(X_test)
      
      
      # mytest <- which(is.null(X_test), arr.ind = TRUE)
      # mytest <- as.data.frame(mytest)
      # y_test <- y_test[-mytest$row]
      
      X_train2 <- X_train[is.finite(rowSums(X_train)),]
      X_train2 <-X_train2[complete.cases(X_train2),]
      X_train2 <- na.omit(X_train2)
      
      
      y_train <- y_train[complete.cases(y_train)]
      y_train <- na.omit(y_train)
      #state <- class::knn1(X_train,X_test,y_train)
      prediction <- class::knn(X_train,as.data.frame(X_test),y_train,k=knn_star_AUROC, l = 0, prob = FALSE, use.all =TRUE)
      result <- list()
      result$AUROC <- AUC::auc(roc(prediction, as.factor(1 * (y_test == +1))))
       #sprintf("%s",replication)
      
       #save("result", file = sprintf("%s/%s/knn_%s_measure_AUROC_replication_%d_result.RData", result_path, cohort, pathway, replication))
      return (result)
      }
      start_p <- Sys.time()  
      for (cohort in cohorts) {
        if (dir.exists(sprintf("%s/%s", result_path, cohort)) == FALSE) {
          dir.create(sprintf("%s/%s", result_path, cohort)) 
        }
        
        opt <- list(chunkSize = 10) 
        
        result1 <- foreach (replication = 1:100, .options.azure = opt, .combine='c', .errorhandling = 'remove',.packages = c('class','AUC'))  %dopar% {
          each_rep(replication)
         
        }
       
      }
      end_p <- Sys.time() 
      total_time <- difftime(end_p, start_p, unit = "min")  
      for(replication in 1:100){
        result <- result1[replication]
        save("result", file = sprintf("%s/%s/knn_%s_measure_AUROC_replication_%d_result.RData", result_path, cohort, pathway, replication))
      }
      
#}

