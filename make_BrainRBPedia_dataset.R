# ---
# Title: make_BrainRBPedia_dataset.R
# Purpose: This script generates the columns in BrainRBPedia from various databases
#          and uses nested cross-validation to construct the machine learning model
#          that predicts the RBPs' disease susceptibility.
# ---

library(data.table)
library(dplyr)
library(biomaRt)
library(readxl)

setwd("~/Desktop/BrainRBPedia") # CHANGE THIS to the path that contains the BrainRBPedia folder
source("./nested_cv_funcs.R")

rbps.path <- "./data/starting_1072_RBP_list.tsv"
rbps <- fread(rbps.path, data.table = FALSE)
names(rbps) <- gsub(" ", "_", names(rbps))
rbps$Human_gene[rbps$Ensembl_ID == "ENSG00000047597"] <- "XK"
rbps$Mouse_gene[rbps$Ensembl_ID == "ENSG00000047597"] <- "Xk"
rbps$Human_gene[rbps$Human_gene == "C17orf85"] <- "NCBP3"
rbps$Mouse_gene[rbps$Human_gene == "NCBP3"] <- "Ncbp3"

# Need to generate the gene list and get all columns required for the model:
# Contains canonical RBDs, Neuron enrichment (human), Neuron enrichment (mouse),
# Brain enrichment, Neuron developmental enrichment (early/late), pLI, Disease susceptibility

# Add protein domains
mart <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = mart)
# interpro with short description
domains <- getBM(c("ensembl_gene_id", "interpro", "interpro_short_description", "interpro_start", "interpro_end"), mart = ensembl)
rbps.domains <- aggregate(interpro_short_description ~ ensembl_gene_id, domains,
                          FUN = toString)
for (i in 1:nrow(rbps.domains)) {
  rbps.domains$interpro_short_description[i] <- paste(unique(unlist(strsplit(rbps.domains$interpro_short_description[i], split = ", "))), collapse = "; ")
}

rbps.domains$interpro_short_description <- gsub("^; ", "", rbps.domains$interpro_short_description)
rbps.domains$interpro_short_description <- gsub("; $", "", rbps.domains$interpro_short_description)
rbps.domains$interpro_short_description <- gsub(" ;", "", rbps.domains$interpro_short_description)
rbps.new <- merge(rbps, rbps.domains, by.x = "Ensembl_ID", by.y = "ensembl_gene_id", all.x = TRUE)
names(rbps.new)[5] <- "protein_domains"

# Add "Contains canonical RBDs"
RBD_list.path <- "./data/canonical_RBDs.txt"
RBD_list <- fread(RBD_list.path, data.table = FALSE)

rbps.new$hasCanonicalRBDs <- numeric(nrow(rbps.new))
for (i in seq_len(nrow(rbps.new))) {
  RBDs <- unlist(strsplit(rbps.new$protein_domains[i], "; "))
  interpro_ids <- domains %>% filter_at(vars(interpro_short_description), any_vars(. %in% RBDs))
  if (any(unique(interpro_ids$interpro) %in% RBD_list$interpro)) {
    rbps.new$hasCanonicalRBDs[i] <- "1"
  } else {
    rbps.new$hasCanonicalRBDs[i] <- "0"
  }
}
rbps.new$hasCanonicalRBDs <- as.numeric(rbps.new$hasCanonicalRBDs)

# Add "Neuron enrichment (human)"
human_neuron_markers.path <- "./data/human_neuron_markers.txt"
human_neuron_markers <- fread(human_neuron_markers.path, data.table = FALSE)
rbps.new <- merge(rbps.new, human_neuron_markers[, c("avg_log2FC", "gene_symbol")],
                  by.x = "Human_gene", by.y = "gene_symbol", all.x = TRUE)
names(rbps.new)[ncol(rbps.new)] <- "neuronEnrichmentHuman"

# Add "Neuron enrichment (mouse)"
mouse_neuron_markers.path <- "./data/mouse_neuron_markers.txt"
mouse_neuron_markers <- fread(mouse_neuron_markers.path, data.table = FALSE)
mouse_neuron_markers$human_gene_symbol <- toupper(mouse_neuron_markers$gene_symbol)
rbps.new <- merge(rbps.new, mouse_neuron_markers[, c("avg_log2FC", "human_gene_symbol")],
                  by.x = "Human_gene", by.y = "human_gene_symbol", all.x = TRUE) 
names(rbps.new)[ncol(rbps.new)] <- "neuronEnrichmentMouse"

# Add "Brain enrichment"
GTEx_data.path <- "./data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
GTEx_data <- fread(GTEx_data.path, data.table = FALSE)
brain_indices <- grep("Brain", colnames(GTEx_data))

for (i in 1:nrow(GTEx_data)) {
  brain <- GTEx_data[i, brain_indices]
  brain_median <- median(as.numeric(brain))
  non_brain <- GTEx_data[i, -c(1, 2, brain_indices)]
  non_brain_median <- median(as.numeric(non_brain))
  GTEx_data$brainEnrichment[i] <- log2(brain_median + 1) - log2(non_brain_median + 1)
}

rbps.new <- merge(rbps.new, GTEx_data[, c("Description", "brainEnrichment")], by.x = "Human_gene", by.y = "Description", all.x = TRUE)

# Add "Neuron developmental enrichment (early/late)"
elife.path <- "./data/elife-58124-fig2-data3-v1.xlsx"
elife <- read_excel(elife.path)
elife$Gene_ID <- gsub("__chr(.|..)", "", elife$Gene_ID)
rbps.new <- merge(rbps.new, elife[, c("Gene_ID", "day3/day1_diff", "day7/day1_diff")], by.x = "Human_gene", by.y = "Gene_ID", all.x = TRUE)

# Add "Protein expression (caudate/cerebral cortex/hippocampus)"
normal_tissue.path <- "./data/normal_tissue.tsv"
normal_tissue <- fread(normal_tissue.path, data.table = FALSE)
names(normal_tissue) <- gsub(" ", "\\.", names(normal_tissue))
neuronal_types <- normal_tissue %>% filter(Cell.type == "neuronal cells")
aggregated_neuronal <- aggregate(Level ~ Gene, neuronal_types, FUN = toString)
rbps.new <- merge(rbps.new, aggregated_neuronal, by.x = "Ensembl_ID", by.y = "Gene", all.x = TRUE)

# pLI
pLI.path <- "./data/pLI_original_data.gz"
pLI <- fread(pLI.path, data.table = FALSE)
pLI_no_dup <- pLI[-c(11645, 7027, 4975, 19208, 5970, 2067, 7784, 10148, 7528,
                     2640, 18271, 15770, 19490, 15735, 2990, 14210, 16427, 11384,
                     5939, 13660, 5310, 8785, 18617, 19570, 15149, 11652, 11325,
                     2672, 12947, 8351, 6237, 10921, 16410, 11992, 12743, 17709,
                     15573, 18069, 15244, 7007), ] # remove duplicated gene rows
rbps.new <- merge(rbps.new, pLI_no_dup[, c("gene", "pLI")], by.x = "Human_gene", by.y = "gene", all.x = TRUE)

# SFARI and Genomics England
sfari.path <- "./data/SFARI-Gene_genes_11-07-2022release_12-14-2022export.csv"
id.path <- "./data/Intellectual disability.tsv"
sfari <- fread(sfari.path, data.table = FALSE)
id <- fread(id.path, data.table = FALSE)
rbps.new$`Autism susceptibility` <- as.numeric(rbps.new$Human_gene %in% sfari$`gene-symbol`)
rbps.new$`Intellectual disability susceptibility` <- as.numeric(rbps.new$Human_gene %in% id$`Gene Symbol`)

names(rbps.new) <- c("Human gene", "Ensembl ID", "Mouse gene", "Gene description", 
                     "Protein domains", "Contains canonical RBDs", "Neuron enrichment (human)", 
                     "Neuron enrichment (mouse)", "Brain enrichment", "Neuron developmental enrichment (early)", 
                     "Neuron developmental enrichment (late)", "Protein Expression (Caudate/Cerebral Cortex/Hippocampus)",
                     "pLI", "Autism susceptibility", "Intellectual disability susceptibility")
rbps.new <- rbps.new[, c(2, 1, 3:ncol(rbps.new))]

rbps.new <- impute(rbps.new, "mean")
rbps.new$extreme_pLI <- as.numeric(rbps.new$pLI >= 0.99)
rbps.new$high_pLI <- as.numeric(rbps.new$pLI >= 0.9 & rbps.new$pLI < 0.99)
rbps.new$medium_pLI <- as.numeric(rbps.new$pLI >= 0.5 & rbps.new$pLI < 0.9)
rbps.new <- rbps.new[, c(1:13, 16:18, 14:15)]

outfile.path <- "./data/BrainRBPedia_new.tsv" # CHANGE THIS
write.table(rbps.new, outfile.path, quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)

