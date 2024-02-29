# ---
# Title: plot_distribution_fig2.R
# Purpose: This script generates "Figure 2: Distribution of key BrainRBPedia 
#          functional annotations" in the BrainRBPedia manuscript.
# ---

library(readxl)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(data.table)


# Functions ---------------------------------------------------------------

geom_vdensity <- function(data, at, ...) {
  ggplot2::geom_segment(
    data = dplyr::filter(as.data.frame(density(data)[1:2]),
                         seq_along(x) == which.min(abs(x - at))),
    ggplot2::aes(x, 0, xend = x, yend = y), ...)
}


# Main --------------------------------------------------------------------

setwd("~/Documents/BrainRBPedia") # CHANGE THIS to the path that contains the BrainRBPedia folder
rbps <- fread("./data/BrainRBPedia.tsv", data.table = FALSE)
pLI <- read.delim("./data/pLI_original_data.gz")
neuron_markers <- fread("./data/human_neuron_markers.txt", data.table = FALSE)
gtex <- fread("./data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", data.table = FALSE)
dev <- read_excel("./data/elife-58124-fig2-data3-v1.xlsx")

# # Delete rows with duplicated genes: choose rows with higher mouse neuron enrichment value
# n_occur <- data.frame(table(rbps$`Human Gene`))
# n_occur[n_occur$Freq > 1, ]
# duplicates <- rbps[rbps$`Human Gene` %in% n_occur$Var1[n_occur$Freq > 1], ]
# genes_to_keep <- duplicates %>% group_by(`Human Gene`) %>% filter(`Neuron Enrichment (Mouse)` == max(`Neuron Enrichment (Mouse)`))
# genes_to_delete <- duplicates[!duplicates$`Neuron Enrichment (Mouse)` %in% genes_to_keep$`Neuron Enrichment (Mouse)`, ]
# rbps <- rbps[!rbps$`Mouse Gene` %in% genes_to_delete$`Mouse Gene`, ]
# # Check if there are any NA values for HGNC symbol
# na_row <- rbps[is.na(rbps$`Human Gene`) == TRUE, ]
# rbps[rbps$`Ensembl ID` == na_row$`Ensembl ID`, ]$`Human Gene` <- "XK"

# Panel A: comparison of human and mouse neuron enrichment
neuron_enrichment <- data.frame(x = c(rbps$`Neuron enrichment (human)`, 
                                      neuron_markers$avg_log2FC), 
                                type = rep(c("RBPs", "Protein-coding genes"), 
                                           c(length(rbps$`Neuron enrichment (human)`), 
                                             length(neuron_markers$avg_log2FC))))
pA <- ggplot(neuron_enrichment, aes(x = x, fill = type)) + 
      geom_density(alpha = 0.5, size = 0.5) +
      scale_fill_brewer(palette = "PuRd", name = " ") +
      xlab("Neuron enrichment") +
      ylab("Density") +
      theme_cowplot(font_size = 15, 
                    line_size = 1,
                    rel_small = 11/15) +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 9, face = "bold"),
            axis.text.y = element_text(size = 9, face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.text = element_text(size = 9, face = "bold")) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
pA

# Panel B: human brain enrichment
brain_indices <- grep("Brain", colnames(gtex))
for (i in 1:nrow(gtex)) {
  brain <- gtex[i, brain_indices]
  brain_median <- median(as.numeric(brain))
  non_brain <- gtex[i, -c(1, 2, brain_indices)]
  non_brain_median <- median(as.numeric(non_brain))
  gtex$brainEnrichment[i] <- log2(brain_median + 1) - log2(non_brain_median + 1)
}

brain_enrichment <- data.frame(x = c(rbps$`Brain enrichment`, 
                                     gtex$brainEnrichment), 
                               type = rep(c("RBPs", "Protein-coding genes"), 
                                          c(length(rbps$`Brain enrichment`), 
                                            length(gtex$brainEnrichment))))

pB <- ggplot(brain_enrichment, aes(x = x, fill = type)) + 
      geom_density(alpha = 0.5, size = 0.5) +
      scale_fill_brewer(palette = "PuRd", name = " ") +
      xlab("Brain enrichment") +
      ylab("Density") +
      theme_cowplot(font_size = 15, 
                    line_size = 1,
                    rel_small = 11/15) +
      theme(legend.position = c(0.65, 0.9),
            axis.text.x = element_text(size = 9, face = "bold"),
            axis.text.y = element_text(size = 9, face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.text = element_text(size = 7, face = "bold"),
            legend.spacing = unit(0.1, "cm"),
            legend.key.size = unit(0.4, "cm")) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
pB

# Panel C: comparison of early and late neuron development enrichment
develop_enrichment_early <- data.frame(x = c(rbps$`Neuron developmental enrichment (early)`, 
                                             dev$`day3/day1_padj`), 
                                       type = rep(c("RBPs", "Protein-coding genes"), 
                                                  c(length(rbps$`Neuron developmental enrichment (early)`), 
                                                    length(dev$`day3/day1_padj`))))

pC_early <- ggplot(develop_enrichment_early, aes(x = x, fill = type)) + 
            geom_density(alpha = 0.5, size = 0.5) +
            scale_fill_brewer(palette = "PuRd", name = " ") +
            xlab("Neuron development enrichment (early)") +
            ylab("Density") +
            theme_cowplot(font_size = 15, 
                          line_size = 1,
                          rel_small = 11/15) + 
            theme(legend.position = "none",
                  axis.title = element_text(face = "bold"),
                  axis.text.x = element_text(size = 9, face = "bold"),
                  axis.text.y = element_text(size = 9, face = "bold"),
                  legend.text = element_text(size = 9, face = "bold")) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
pC_early

develop_enrichment_late <- data.frame(x = c(rbps$`Neuron developmental enrichment (early)`, 
                                            dev$`day7/day1_padj`), 
                                      type = rep(c("RBPs", "Protein-coding genes"), 
                                                 c(length(rbps$`Neuron developmental enrichment (early)`), 
                                                   length(dev$`day7/day1_padj`))))
pC_late <- ggplot(develop_enrichment_late, aes(x = x, fill = type)) + 
  geom_density(alpha = 0.5, size = 0.5) +
  scale_fill_brewer(palette = "PuRd", name = " ") +
  xlab("Neuron development enrichment (late)") +
  ylab("Density") +
  theme_cowplot(font_size = 15, 
                line_size = 1,
                rel_small = 11/15) + 
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9, face = "bold"),
        axis.text.y = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9, face = "bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
pC_late


# Panel D: histogram for the pLI scores
pLI_no_na <- as.data.frame(na.omit(rbps[, c("pLI")]))
names(pLI_no_na) <- "pLI"

# Version 4: Density
pLI_no_na$breaks <- cut(pLI_no_na$pLI, breaks = c(0, 0.5, 0.9, 0.99, 1))
pD <- ggplot(pLI_no_na, aes(x = pLI)) +
  geom_density(size = 0.7) +
  xlab("pLI (RBPs)") +
  ylab("Density") +
  theme_cowplot(font_size = 15, 
                line_size = 1,
                rel_small = 11/15) + 
  theme(legend.position = c(0.7, 0.9),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
        plot.margin = unit(c(0.2, 1, 0.2, 1), "cm")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) 
pD

# Getting the values of plot
d = ggplot_build(pD)$data[[1]]

# Building shaded area
pD <- pD + 
  geom_area(data = subset(d, x >= 0 & x < 0.5), aes(x = x, y = y), fill = "#9ECAE1", alpha = 1) +
  geom_area(data = subset(d, x >= 0.5 & x < 0.9), aes(x = x, y = y), fill = "#4292C6", alpha = 1) +
  geom_area(data = subset(d, x >= 0.9 & x < 0.99), aes(x = x, y = y), fill = "#08519C", alpha = 1) +
  geom_area(data = subset(d, x >= 0.99), aes(x = x, y = y), fill = "#08306B", alpha = 1)
pD

pE <- ggplot(pLI, aes(x = pLI)) +
  geom_density(size = 0.7) +
  xlab("pLI (protein-coding genes)") +
  ylab("Density") +
  theme_cowplot(font_size = 15, 
                line_size = 1,
                rel_small = 11/15) + 
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 7, face = "bold"),
        legend.spacing = unit(0.1, "cm"), 
        legend.key = element_rect(linewidth = 1.2),
        legend.key.size = unit(0.4, "cm"),
        plot.margin = unit(c(0.2, 1, 0.2, 1), "cm")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) 
pE

# Getting the values of plot
e = ggplot_build(pE)$data[[1]]

# Building shaded area
pE <- pE + 
  geom_area(data = subset(e, x >= 0 & x < 0.5), aes(x = x, y = y, fill = "Low (< 0.5)"), alpha = 1) +
  geom_area(data = subset(e, x >= 0.5 & x < 0.9), aes(x = x, y = y, fill = "Medium (0.5~0.9)"), alpha = 1) +
  geom_area(data = subset(e, x >= 0.9 & x < 0.99), aes(x = x, y = y, fill = "High (0.9~0.99)"), alpha = 1) +
  geom_area(data = subset(e, x >= 0.99), aes(x = x, y = y, fill = "Extreme (> 0.99)"), alpha = 1) +
  theme(legend.position = c(0.65, 0.9)) +
  scale_fill_manual(" ",
                    values = c("Low (< 0.5)" = "#9ECAE1", "Medium (0.5~0.9)" = "#4292C6",
                              "High (0.9~0.99)" = "#08519C", "Extreme (> 0.99)" = "#08306B"),
                    limits = c("Extreme (> 0.99)", "High (0.9~0.99)", "Medium (0.5~0.9)", "Low (< 0.5)"))
pE

fig2 <- plot_grid(pA, pB, pC_early, pC_late, pD, pE,
                  nrow = 3, ncol = 2, align = "hv",
                  labels = "AUTO")
fig2

# Stats
ks.test(rbps$`Neuron enrichment (human)`, neuron_markers$avg_log2FC)
ks.test(rbps$`Brain enrichment`, gtex$brainEnrichment)
ks.test(rbps$`Neuron developmental enrichment (early)`, dev$`day3/day1_padj`)
ks.test(rbps$`Neuron developmental enrichment (late)`, dev$`day7/day1_padj`)
ks.test(pLI_no_na$pLI, pLI$pLI)

ggsave("./figures/BrainRBPedia_Figure2.jpg", fig2, 
       width = 400, height = 250, units = c("mm"))

