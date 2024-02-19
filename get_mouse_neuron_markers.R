# ---
# Title: get_mouse_neuron_markers.R
# Purpose: Use this script to generate the log2 fold-change of gene expression 
#          in mouse neurons relative to non-neurons.
# ---

library(Seurat)
library(stringr)
library(tidyverse)
library(readxl)
library(data.table)

options(stringsAsFactors = F)
options(echo = TRUE)

# Use the following if running from the script from command line:
# args <- commandArgs(trailingOnly = TRUE)
# visp.dir_path <- args[1]
# visp.matrix_path <- paste0(visp.dir_path, args[2])
# visp.metadata_path <- paste0(visp.dir_path, args[3])
# output_filename <- args[4]

setwd("~/Desktop/BrainRBPedia") # CHANGE THIS
visp.dir_path <- "./data/visp_data/"
visp.matrix_path <- paste0(visp.dir_path, "matrix_VISp.csv")
visp.metadata_path <- paste0(visp.dir_path, "metadata_VISp.csv")
output_filename <- "./data/mouse_neuron_markers.txt"

visp.matrix <- fread(visp.matrix_path, data.table = FALSE)
visp.metadata <- fread(visp.metadata_path, data.table = FALSE)

row.names(visp.metadata) <- visp.metadata$sample_name
row.names(visp.matrix) <- visp.matrix$sample_name
visp.matrix <- visp.matrix[,-c(1,2)]
visp.matrix.t <- t(visp.matrix)

rownames(visp.matrix) <- str_replace_all(rownames(visp.matrix), pattern = "_", replacement = "-")
rownames(visp.metadata) <- str_replace_all(rownames(visp.metadata), pattern = "_", replacement = "-")
colnames(visp.matrix.t) <- str_replace_all(colnames(visp.matrix.t),pattern = "_",replacement = "-")

visp.seurat.t <- CreateSeuratObject(counts = visp.matrix.t, project = "visp", min.cells = 3, min.features = 200, meta.data = visp.metadata)
visp.seurat.t <- NormalizeData(visp.seurat.t, normalization.method = "LogNormalize", scale.factor = 1000000)

Idents(visp.seurat.t) <- 'class_label'
Idents(visp.seurat.t)

visp.seurat.t@meta.data$superclass = factor('Non-neuronal', levels = c('Non-neuronal', "Neuronal", "Others"))
visp.seurat.t@meta.data$superclass[visp.seurat.t@meta.data$class_label == "GABAergic"] = 'Neuronal'
visp.seurat.t@meta.data$superclass[visp.seurat.t@meta.data$class_label == "Glutamatergic"] = 'Neuronal'
visp.seurat.t@meta.data$superclass[visp.seurat.t@meta.data$class_label == ""] = 'Others'
visp.seurat.t@meta.data$superclass[visp.seurat.t@meta.data$class_label == "Non-Neuronal"] = 'Non-neuronal'

Idents(visp.seurat.t) <- "superclass"
Idents(visp.seurat.t)

neuron.markers <- FindMarkers(visp.seurat.t, ident.1 = "Neuronal", ident.2 = "Non-neuronal",
                              logfc.threshold = 0.0, min.pct = 0.0, test.use = "MAST")

write.table(neuron.markers, output_filename, quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)
