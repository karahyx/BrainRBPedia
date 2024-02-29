# ---
# Title: nested_cv_nonRBPs.R
# Purpose: This script generates the columns in BrainRBPedia from various databases
#          for protein-coding non-RBPs obtained from GENCODE release 44. It also 
#          uses nested cross-validation to construct the machine learning model
#          that predicts these protein-coding non-RBPs' disease susceptibility.
# ---

library(data.table)
library(glmnet)
library(dplyr)
library(ROCR)
library(cowplot)
library(readxl)
library(biomaRt)

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
    for (i in 1:9) {
      if (is.numeric(new_data[[i]])) {
        new_data[[i]][is.na(new_data[[i]])] <- mean(new_data[[i]], na.rm = TRUE)
      }
    }
  }
  # Replacing missing values with the median
  else if (replacement == "median") {
    for (i in 1:9) {
      if (is.numeric(new_data[[i]])) {
        new_data[[i]][is.na(new_data[[i]])] <- stats::median(new_data[[i]], na.rm = TRUE)
      }
    }
  }
  return(new_data)
}

trainCV <- function(data, K = 5, indep_var) {
  
  modelList <- list()
  prList <- list()
  testList <- list()
  coefList <- list()
  aucList <- list()
  aucprList <- list()
  
  switch(indep_var,
         autism = y <- data$`Autism susceptibility`,
         id = y <- data$`Intellectual disability susceptibility`)
  
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
    x <- as.matrix(standardizedTrain[, c(4:9, 11:13)])
    switch(indep_var,
           autism = y <- standardizedTrain$`Autism susceptibility`,
           id = y <- standardizedTrain$`Intellectual disability susceptibility`)
    
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
    pred <- predict(cvfit, newx = as.matrix(standardizedTest[, c(4:9, 11:13)]), 
                    type = 'response')[,1]
    
    standardizedTest$predictions <- pred
    predicted <- standardizedTest %>% arrange(desc(predictions))
    
    prList[[i]] <- predicted
    
    switch(indep_var,
           autism = pred_labels <- predicted$`Autism susceptibility`,
           id = pred_labels <- predicted$`Intellectual disability susceptibility`)
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


# Main --------------------------------------------------------------------

setwd("~/Documents/BrainRBPedia") # CHANGE THIS to the path that contains the BrainRBPedia folder

# Need to generate the gene list and get all columns required for the model:
# Contains canonical RBDs, Neuron enrichment (human), Neuron enrichment (mouse),
# Brain enrichment, Neuron developmental enrichment (early/late), pLI, Disease susceptibility

# The following code is used to generate non_rbps.tsv
# protein-coding genes downloaded from HGNC
pcg <- fread("./data/protein-coding_gene.txt", data.table = F)
non_rbps <- setdiff(pcg$symbol, rbps_raw$`Human gene`)
non_rbps.df <- pcg[which(pcg$symbol %in% non_rbps), c("symbol", "ensembl_gene_id")]

# Contains canonical RBDs
mart <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = mart)
# interpro with short description
domains <- getBM(c("ensembl_gene_id", "interpro", "interpro_short_description", "interpro_start", "interpro_end"), mart = ensembl)

nrbps.domains <- aggregate(interpro_short_description ~ ensembl_gene_id, domains,
                           FUN = toString)
for (i in 1:nrow(nrbps.domains)) {
  nrbps.domains$interpro_short_description[i] <- paste(unique(unlist(strsplit(nrbps.domains$interpro_short_description[i], split = ", "))), collapse = "; ")
}

nrbps.domains$interpro_short_description <- gsub("^; ", "", nrbps.domains$interpro_short_description)
nrbps.domains$interpro_short_description <- gsub("; $", "", nrbps.domains$interpro_short_description)
non_rbps.new <- merge(non_rbps.df, nrbps.domains, by = "ensembl_gene_id", all.x = T)

# # The following code is used to generate canonical_RBDs.txt
# rbds <- read.xlsx("~/Desktop/CisBP-RNA_RBDs.xlsx", sheetIndex = 1)
# rbds2 <- data.frame(RBDs = c("DEAD", "dsrm", "G-patch", "OST-HTH/LOTUS", "ROQ_II", "NGL", "PAZ", "Piwi"),
#                     pfam = c("PF00270", "PF00035", "PF01585", "PF12872", "PF18386", "PF01436", "PF02170", "PF02171"),
#                     interpro = c("IPR011545", "IPR014720", "IPR000467", "IPR025605", "IPR041523", "IPR001258", "IPR003100", "IPR003165"))
# names(rbds) <- c("RBDs", "pfam", "interpro")
# rbds <- rbind(rbds, rbds2)
# write.table(rbds, "./data//canonical_RBDs.txt", row.names = F, col.names = T, quote = F, sep = "\t")
RBD_list <- fread("./data/canonical_RBDs.txt", data.table = F)

# add 'Contains canonical RBDs' column
names(non_rbps.new)[3] <- "protein_domains"
non_rbps.new <- non_rbps.new %>% subset(!is.na(protein_domains))
non_rbps.new$hasCanonicalRBDs <- numeric(nrow(non_rbps.new))

for (i in seq_len(nrow(non_rbps.new))) {
  RBDs <- unlist(strsplit(non_rbps.new$protein_domains[i], "; "))
  interpro_ids <- domains %>% filter_at(vars(interpro_short_description), any_vars(. %in% RBDs))
  if (any(unique(interpro_ids$interpro) %in% RBD_list$interpro)) {
    non_rbps.new$hasCanonicalRBDs[i] <- "1"
  } else {
    non_rbps.new$hasCanonicalRBDs[i] <- "0"
  }
}

# Neuron enrichment (human)
human_neuron_markers <- fread("./data/human_neuron_markers.txt", data.table = F)
non_rbps.hn <- merge(non_rbps.new, human_neuron_markers[, c("avg_log2FC", "gene_symbol")],
                     by.x = "symbol", by.y = "gene_symbol", all.x = T) # 48.0% of non_rbps.new doesn't have an avg_log2FC matched
names(non_rbps.hn)[ncol(non_rbps.hn)] <- "neuronEnrichmentHuman"

# Neuron enrichment (mouse)
mouse_neuron_markers <- fread("./data/mouse_neuron_markers.txt", data.table = F)
mouse_neuron_markers$human_gene_symbol <- toupper(mouse_neuron_markers$gene_symbol)
non_rbps.ms <- merge(non_rbps.hn, mouse_neuron_markers[, c("avg_log2FC", "human_gene_symbol")],
                     by.x = "symbol", by.y = "human_gene_symbol", all.x = T) # 22.9% of non_rbps.ms doesn't have an avg_log2FC matched
names(non_rbps.ms)[ncol(non_rbps.ms)] <- "neuronEnrichmentMouse"

# Brain enrichment
GTEx_data <- fread("./data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", data.table = F)
brain_indices <- grep("Brain", colnames(GTEx_data))

for (i in 1:nrow(GTEx_data)) {
  brain <- GTEx_data[i, brain_indices]
  brain_median <- median(as.numeric(brain))
  non_brain <- GTEx_data[i, -c(1, 2, brain_indices)]
  non_brain_median <- median(as.numeric(non_brain))
  GTEx_data$brainEnrichment[i] <- log2(brain_median + 1) - log2(non_brain_median + 1)
}

non_rbps.bn <- merge(non_rbps.ms, GTEx_data[, c("Description", "brainEnrichment")], by.x = "symbol", by.y = "Description", all.x = TRUE)
# 5.6% of non_rbps.bn doesn't have a brain enrichment value matched
# remove duplicated entries
non_rbps.bn <- non_rbps.bn[-c(485, 1106, 1107, 2495, 3468, 3520, 3725, 4075,
                              4635, 6125, 6430, 7252, 7264, 8622, 8733, 11064,
                              11709, 12087, 12995, 13942, 14176, 14734, 14900,
                              17063, 17459), ]

# Neuron developmental enrichment (early/late)
elife <- read_excel("./data/elife-58124-fig2-data3-v1.xlsx")
elife$Gene_ID <- gsub("__chr(.|..)", "", elife$Gene_ID)
non_rbps.dv <- merge(non_rbps.bn, elife[, c("Gene_ID", "day3/day1_diff", "day7/day1_diff")], by.x = "symbol", by.y = "Gene_ID", all.x = TRUE)

# pLI
pLI <- fread("./data/pLI_original_data.gz", data.table = F)
pLI_no_dup <- pLI[-c(11645, 7027, 4975, 19208, 5970, 2067, 7784, 10148, 7528,
                     2640, 18271, 15770, 19490, 15735, 2990, 14210, 16427, 11384,
                     5939, 13660, 5310, 8785, 18617, 19570, 15149, 11652, 11325,
                     2672, 12947, 8351, 6237, 10921, 16410, 11992, 12743, 17709,
                     15573, 18069, 15244, 7007), ]
non_rbps.pLI <- merge(non_rbps.dv, pLI_no_dup[, c("gene", "pLI")], by.x = "symbol", by.y = "gene", all.x = T)
non_rbps.pLI$extreme_pLI <- as.numeric(non_rbps.pLI$pLI >= 0.99)
non_rbps.pLI$high_pLI <- as.numeric(non_rbps.pLI$pLI >= 0.9 & non_rbps.pLI$pLI < 0.99)
non_rbps.pLI$medium_pLI <- as.numeric(non_rbps.pLI$pLI >= 0.5 & non_rbps.pLI$pLI < 0.9)

# SFARI and Genomics England
sfari <- fread("./data/SFARI-Gene_genes_11-07-2022release_12-14-2022export.csv")
id <- fread("./data/Intellectual disability.tsv")
non_rbps.pLI$`Autism susceptibility` <- as.numeric(non_rbps.pLI$symbol %in% sfari$`gene-symbol`)
non_rbps.pLI$`Intellectual disability susceptibility` <- as.numeric(non_rbps.pLI$symbol %in% id$`Gene Symbol`)

non_rbps.pLI$hasCanonicalRBDs <- as.numeric(non_rbps.pLI$hasCanonicalRBDs)
non_rbps <- impute(non_rbps.pLI, "mean")
non_rbps$extreme_pLI[is.na(non_rbps$extreme_pLI)] <- 0
non_rbps$high_pLI[is.na(non_rbps$high_pLI)] <- 0
non_rbps$medium_pLI[is.na(non_rbps$medium_pLI)] <- 0

names(non_rbps) <- c("Human gene", "Ensembl ID", "Protein domains", "Contains canonical RBDs",
                     "Neuron enrichment (human)", "Neuron enrichment (mouse)", "Brain enrichment",
                     "Neuron developmental enrichment (early)", "Neuron developmental enrichment (late)",
                     "pLI", "Extreme pLI", "High pLI", "Medium pLI", "Autism susceptibility",
                     "Intellectual disability susceptibility")
write.table(non_rbps, "./data/non_rbps.tsv", quote = F,
            row.names = F, col.names = T, sep = "\t")

# After generating non_rbps.tsv, you can just use the following line to read it in
# instead of running the code above again
non_rbps <- fread("/Users/karahan/Documents/RBP/data/non_rbps.tsv", data.table = FALSE)

# positive examples: protein-coding genes that are non-RBPs and related to ASD/ID
# negative examples: protein-coding genes that are non-RBPs and not related to ASD/ID
autism_res2 <- trainCV(non_rbps, K = 5, indep_var = "autism")
id_res2 <- trainCV(non_rbps, K = 5, indep_var = "id")

autism_pred2 <- do.call("rbind", autism_res$predictions)
id_pred2 <- do.call("rbind", id_res$predictions)

autism_pred2 <- autism_pred %>% arrange(-predictions)
id_pred2 <- id_pred %>% arrange(-predictions)


