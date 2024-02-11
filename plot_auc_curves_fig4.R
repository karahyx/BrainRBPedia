# ---
# Title: plot_auc_curves_fig4.R
# Purpose: This script generates "Figure 4: ASD and ID model performance measured
#          by ROC and PR curves" in the BrainRBPedia manuscript.
# Prerequisite: Run nested_cv.R to obtain objects "autism_res" and "id_res" before
#               running this script.
# ---

library(ggplot2)
library(cowplot)

# ROC
plotROC <- function(results, dependentVarIndex, predictionIndex) {
  
  # Performing checks of user input
  if (class(results) != "trainCV") {
    stop("results should be an S3 object of class trainCV.")
  }
  
  if (is.numeric(dependentVarIndex) == FALSE) {
    stop("dependentVarIndex should be a positive interger.")
  }
  
  if (is.numeric(dependentVarIndex) == TRUE && dependentVarIndex < 1) {
    stop("dependentVarIndex should be a positive interger.")
  }
  
  # Initializing variables
  K <- length(results$models)
  allpf <- data.frame(fpr = numeric(),
                      tpr = numeric(),
                      fold = numeric())
  aucs <- numeric()
  
  # Adding predictions and auc values from each model to allpf
  for (i in 1:K) {
    pred <- ROCR::prediction(results$predictions[[i]][[predictionIndex]],
                             results$predictions[[i]][[dependentVarIndex]])
    perf <- ROCR::performance(pred, "tpr", "fpr")
    
    area <- pROC::auc(results$predictions[[i]][[dependentVarIndex]],
                      results$predictions[[i]][[predictionIndex]])
    aucs <- append(aucs, area)
    
    pf <- data.frame(fpr = perf@x.values,
                     tpr = perf@y.values,
                     fold = i)
    names(pf) <- c("fpr", "tpr", "fold")
    
    allpf <- rbind(allpf, pf)
  }
  
  # Formatting the legend labels
  legendLabels <- character()
  for (i in 1:K) {
    content <- paste("Model ", as.character(i), " (AUC = ",
                     as.character(round(aucs[i], 2)), ")", sep = "")
    legendLabels <- append(legendLabels, content)
  }
  
  # Generating the ROC curves
  rocCurve <-
    ggplot2::ggplot(data = allpf, aes(x = fpr, y = tpr, colour = as.factor(fold))) +
    ggplot2::geom_line() +
    cowplot::theme_cowplot() +
    ggplot2::theme(axis.title = element_text(size = 16, face = "bold"),
                   axis.text.x = element_text(size = 12, face = "bold"),
                   axis.text.y = element_text(size = 12, face = "bold"),
                   axis.line = element_line(size = 1),
                   legend.text = element_text(size = 11, face = "bold"),
                   legend.title = element_blank(),
                   legend.position = c(0.5, 0.3),
                   legend.box.background = element_blank(),
                   legend.key.size = unit(1.0, "cm"),
                   plot.title = element_text(size = 30, hjust = 0.5)) +
    ggsci::scale_color_lancet(labels = legendLabels) +
    ggplot2::labs(x = "False positive rate", y = "True positive rate") +
    ggplot2::ggtitle(label = "ROC Curves")
  
  return(rocCurve)
}

# PR Curves
plotPR <- function(results, dependentVarIndex, predictionIndex) {
  # Performing checks of user input
  if (class(results) != "trainCV") {
    stop("results should be an S3 object of class trainCV.")
  }
  
  if (is.numeric(dependentVarIndex) == FALSE) {
    stop("dependentVarIndex should be a positive interger.")
  }
  
  if (is.numeric(dependentVarIndex) == TRUE && dependentVarIndex < 1) {
    stop("dependentVarIndex should be a positive interger.")
  }
  
  
  # Initializing variables
  K <- length(results$models)
  allpf <- data.frame(fpr = numeric(),
                      tpr = numeric(),
                      fold = numeric())
  aucprs <- numeric()
  
  # Adding predictions and aucpr values from each model to allpf
  for (i in 1:K) {
    pred <- ROCR::prediction(results$predictions[[i]][[predictionIndex]],
                             results$predictions[[i]][[dependentVarIndex]])
    perf <- ROCR::performance(pred, "prec", "rec")
    
    aucprObj <- ROCR::performance(pred, "aucpr")
    aucpr <- aucprObj@y.values[[1]]
    aucprs <- append(aucprs, aucpr)
    
    pf <- data.frame(rec = perf@x.values,
                     prec = perf@y.values,
                     fold = i)
    data.frame(rec = perf@x.values, prec = perf@y.values)
    names(pf) <- c("rec", "prec", "fold")
    
    allpf <- rbind(allpf, pf)
  }
  
  # Formatting the legend labels
  legendLabels <- character()
  for (i in 1:K) {
    content <- paste("Model ", as.character(i), " (AUCPR = ",
                     as.character(round(aucprs[i], 2)), ")", sep = "")
    legendLabels <- append(legendLabels, content)
  }
  
  # Generating the precision-recall curves
  prCurve <-
    ggplot2::ggplot(data = allpf, aes(x = rec, y = prec, colour = as.factor(fold))) +
    ggplot2::geom_line() +
    cowplot::theme_cowplot() +
    ggplot2::theme(axis.title = element_text(size = 16, face = "bold"),
                   axis.text.x = element_text(size = 12, face = "bold"),
                   axis.text.y = element_text(size = 12, face = "bold"),
                   axis.line = element_line(size = 1),
                   legend.text = element_text(size = 11, face = "bold"),
                   legend.title = element_blank(),
                   legend.position = c(0.5, 0.8),
                   legend.box.background = element_blank(),
                   legend.key.size = unit(1.0, "cm"),
                   plot.title = element_text(size = 30, hjust = 0.5)) +
    ggsci::scale_color_lancet(labels = legendLabels) +
    ggplot2::labs(x = "Recall", y = "Precision") +
    ggplot2::ggtitle(label = "Precision-Recall Curve")
  
  return(prCurve)
}

# Run nested_cv.R first to get autism_res and id_res
asd_roc <- plotROC(autism_res, 17, 19)
asd_roc
asd_pr <- plotPR(autism_res, 17, 19)
asd_pr

id_roc <- plotROC(id_res, 18, 19)
id_roc
id_pr <- plotPR(id_res, 18, 19)
id_pr

fig4 <- plot_grid(asd_roc, id_roc, asd_pr, id_pr, labels = "AUTO", align = "hv", label_size = 20)  
fig4

ggsave("./figures/BrainRBPedia_Figure4.jpg", fig4, width = 300, height = 300, units = c("mm"))



