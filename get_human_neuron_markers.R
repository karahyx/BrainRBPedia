# ---
# Title: get_human_neuron_markers.R
# Purpose: Use this script to generate the log2 fold-change of gene expression 
#          in human neurons relative to non-neurons.
# ---

library(Seurat)

setwd("~/Desktop/BrainRBPedia") # CHANGE THIS
output_filename <- "./data/human_neuron_markers.txt"

Seu_AIBS_obj <- readRDS("./data/Seu_AIBS_obj_update_07JUN21.rds")
Idents(Seu_AIBS_obj) <- "NeuN"
human_neuron_markers <- FindMarkers(Seu_AIBS_obj, ident.1 = "Neuronal", ident.2 = "Non-neuronal", 
                                    logfc.threshold = 0.0, min.pct = 0.0, test.use = "MAST")

write.table(human_neuron_markers, output_filename, sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)
