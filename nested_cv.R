library(data.table)
library(glmnet)
library(dplyr)
library(ROCR)
library(cowplot)

# FUNCS -------------------------------------------------------------------

impute <- function(data, replacement) {
  
  # Performing checks of user input
  if (is.data.frame(data) == FALSE) {
    stop("data should be a data frame")
  }
  
  if (replacement != "mean" & replacement != "median") {
    stop("replacement should be either mean or median")
  }
  
  new_data <- data
  # Replacing missing values with the mean
  if (replacement == "mean") {
    for (i in 1:ncol(new_data)) {
      if (is.numeric(new_data[[i]])) {
        new_data[[i]][is.na(new_data[[i]])] <- mean(new_data[[i]], na.rm = TRUE)
      }
    }
  }
  # Replacing missing values with the median
  else if (replacement == "median") {
    for (i in 1:ncol(new_data)) {
      if (is.numeric(new_data[[i]])) {
        new_data[[i]][is.na(new_data[[i]])] <- stats::median(new_data[[i]], na.rm = TRUE)
      }
    }
  }
  return(new_data)
}

# Process rbps
rbps_raw <- fread("~/Desktop/RBP/data/BrainRBPedia.csv", data.table = FALSE)
sfari <- fread("~/Desktop/RBP/data/SFARI-Gene_genes_11-07-2022release_12-14-2022export.csv")
epi <- fread("~/Desktop/RBP/data/Genetic epilepsy syndromes.tsv")
id <- fread("~/Desktop/RBP/data/Intellectual disability.tsv")

rbps_raw$`Autism susceptibility` <- as.numeric(rbps_raw$`Human gene` %in% sfari$`gene-symbol`)
rbps_raw$`Epilepsy susceptibility` <- as.numeric(rbps_raw$`Human gene` %in% epi$`Gene Symbol`)
rbps_raw$`Intellectual disability susceptibility` <- as.numeric(rbps_raw$`Human gene` %in% id$`Gene Symbol`)

rbps <- impute(rbps_raw, "mean")

rbps$extreme_pLI <- as.numeric(rbps$pLI >= 0.99)
rbps$high_pLI <- as.numeric(rbps$pLI >= 0.9 & rbps$pLI < 0.99)
rbps$medium_pLI <- as.numeric(rbps$pLI >= 0.5 & rbps$pLI < 0.9)


# Training ----------------------------------------------------------------

trainCV <- function(data, K = 5, indep_var) {
  
  modelList <- list()
  prList <- list()
  testList <- list()
  coefList <- list()
  aucList <- list()
  aucprList <- list()
  
  switch(indep_var,
         autism = y <- data$`Autism susceptibility`,
         id = y <- data$`Intellectual disability susceptibility`,
         epilepsy = y <- data$`Epilepsy susceptibility`)
  
  # y <- rbps$`Autism susceptibility`
  # y <- rbps$`Intellectual disability susceptibility`
  # y <- rbps$`Epilepsy susceptibility`
  
  # Assign folds evenly using the modulus operator
  set.seed(3019)
  fold0 <- sample.int(sum(y == 0)) %% K
  fold1 <- sample.int(sum(y == 1)) %% K
  foldid <- numeric(length(y))
  foldid[y == 0] <- fold0
  foldid[y == 1] <- fold1
  foldid <- foldid + 1
  
  for (i in 1:K) {
    # Creating training and test data sets
    trainingData <- data[which(foldid != i), ]
    testData <- data[which(foldid == i), ]
    
    # Standardize the training set
    standardizedTrain <- trainingData %>%
      dplyr::mutate_at(c("Neuron enrichment (human)",
                         "Neuron enrichment (mouse)", 
                         "Brain enrichment",
                         "Neuron developmental enrichment (early)",
                         "Neuron developmental enrichment (late)"), ~(scale(.) %>% as.vector))
    
    # Standardize the test set
    standardizedTest <- testData %>%
      dplyr::mutate_at(c("Neuron enrichment (human)",
                         "Neuron enrichment (mouse)", 
                         "Brain enrichment",
                         "Neuron developmental enrichment (early)",
                         "Neuron developmental enrichment (late)"), ~(scale(.) %>% as.vector))
    
    # Adding the test set to testList
    testList[[i]] <- standardizedTest
    
    # Find lambda - size of the penalty
    x <- as.matrix(standardizedTrain[, c(6:11, 17:19)])
    y <- standardizedTrain$`Autism susceptibility`
    # y <- standardizedTrain$`Intellectual disability susceptibility`
    # y <- standardizedTrain$`Epilepsy susceptibility`
    
    fraction_0 <- rep(1 - sum(y == 0) / length(y), sum(y == 0))
    fraction_1 <- rep(1 - sum(y == 1) / length(y), sum(y == 1))
    # assign that value to a "weights" vector
    weights <- numeric(length(y))
    weights[y == 0] <- fraction_0
    weights[y == 1] <- fraction_1
    
    nfold <- 5
    # assign folds evenly using the modulus operator
    fold0_inner <- sample.int(sum(y == 0)) %% nfold
    fold1_inner <- sample.int(sum(y == 1)) %% nfold
    foldid_inner <- numeric(length(y))
    foldid_inner[y == 0] <- fold0_inner
    foldid_inner[y == 1] <- fold1_inner
    foldid_inner <- foldid_inner + 1
    
    set.seed(1028)
    cvfit <- cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = nfold, foldid = foldid_inner,
                       standardize = TRUE, type.measure = "default", weights = weights)
    modelList[[i]] <- cvfit
    
    coefs <- as.matrix(coef(cvfit, s = cvfit$lambda.1se))
    ix <- which(abs(coefs[,1]) > 0)
    
    coefList[[i]] <- coefs[ix, 1, drop = FALSE]
    
    set.seed(2382)
    pred <- predict(cvfit, newx = as.matrix(standardizedTest[, c(6:11, 17:19)]), 
                    type = 'response')[,1]
    
    standardizedTest$predictions <- pred
    predicted <- standardizedTest %>% arrange(desc(predictions))
    
    prList[[i]] <- predicted
    
    pred_ROCR <- prediction(predicted$predictions, predicted$`Autism susceptibility`)
    auc <- performance(pred_ROCR, measure = "auc")
    auc <- auc@y.values[[1]]
    aucList[[i]] <- auc
    
    aucpr <- ROCR::performance(pred_ROCR, "aucpr")
    aucpr <- aucpr@y.values[[1]]
    aucprList[[i]] <- aucpr
    
  }
  
  results <- list(models = modelList,
                  predictions = prList,
                  testSets = testList,
                  betaCoefs = coefList,
                  aucs = aucList,
                  aucprs = aucprList)
  
  class(results) <- "trainCV"
  
  return(results)
}

names(rbps)[17:19] <- c("Extreme pLI", "High pLI", "Medium pLI")
write.table(rbps, file = "~/Desktop/RBP/tables/TableS2_BrainRBPedia.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

autism_res <- trainCV(rbps, K = 5, indep_var = "autism")
id_res <- trainCV(rbps, K = 5, indep_var = "id")
# epi_res <- trainCV(rbps, K = 5, indep_var = "epilepsy")

autism_pred <- do.call("rbind", autism_res$predictions)
id_pred <- do.call("rbind", id_res$predictions)

autism_pred <- autism_pred %>% arrange(-predictions)
