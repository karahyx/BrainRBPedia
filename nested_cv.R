# ---
# Title: nested_cv.R
# Purpose: This script generates the columns in BrainRBPedia from various databases
#          and uses nested cross-validation to construct the machine learning model
#          that predicts the RBPs' disease susceptibility.
# ---

library(data.table)
library(dplyr)
library(glmnet)
library(ROCR)

setwd("~/Desktop/BrainRBPedia") # CHANGE THIS to the path that contains the BrainRBPedia folder
source("./nested_cv_funcs.R")

rbps <- fread("./data/BrainRBPedia.tsv", data.table = FALSE)

# Train machine learning model
autism_res <- trainCV(rbps, K = 5, indep_var = "autism")
autism_pred <- do.call("rbind", autism_res$predictions)
autism_pred <- autism_pred %>% arrange(-predictions)

id_res <- trainCV(rbps, K = 5, indep_var = "id")
id_pred <- do.call("rbind", id_res$predictions)
id_pred <- id_pred %>% arrange(-predictions)

