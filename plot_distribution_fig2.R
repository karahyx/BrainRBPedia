# ---
# Title: plot_distribution_fig2.R
# Purpose: This script generates "Figure 2: Distribution of key BrainRBPedia 
#          functional annotations" in the BrainRBPedia manuscript.
# ---

library(data.table)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library("wesanderson")

setwd("~/Desktop/BrainRBPedia")
rbps <- fread("./data/BrainRBPedia.tsv")

# # delete rows with duplicated genes: choose rows with higher mouse neuron
# # enrichment value
# n_occur <- data.frame(table(rbps$`Human Gene`))
# n_occur[n_occur$Freq > 1, ]
# duplicates <- rbps[rbps$`Human Gene` %in% n_occur$Var1[n_occur$Freq > 1], ]
# genes_to_keep <- duplicates %>% group_by(`Human Gene`) %>% filter(`Neuron Enrichment (Mouse)` == max(`Neuron Enrichment (Mouse)`))
# genes_to_delete <- duplicates[!duplicates$`Neuron Enrichment (Mouse)` %in% genes_to_keep$`Neuron Enrichment (Mouse)`, ]
# rbps <- rbps[!rbps$`Mouse Gene` %in% genes_to_delete$`Mouse Gene`, ]
# # check if any NA values for HGNC symbol
# na_row <- rbps[is.na(rbps$`Human Gene`) == TRUE, ]
# rbps[rbps$`Ensembl ID` == na_row$`Ensembl ID`, ]$`Human Gene` <- "XK"

# Panel A: comparison of human and mouse neuron enrichment
neuron_enrichment <- data.frame(x = c(rbps$`Neuron enrichment (human)`, 
                                      rbps$`Neuron enrichment (mouse)`), 
                                type = rep(c("Human", "Mouse"), 
                                           c(length(rbps$`Neuron enrichment (human)`), 
                                             length(rbps$`Neuron enrichment (mouse)`))))
pA <- ggplot(neuron_enrichment, aes(x = x, fill = type)) + 
      geom_density(alpha = 0.5, size = 0.9) +
      scale_fill_brewer(palette = "PiYG", name = " ") +
      xlab("log2 neuron/non-neuron expression") +
      ylab("Density") +
      theme_cowplot(font_size = 15, 
                    line_size = 1,
                    rel_small = 11/15) +
      theme(legend.position = c(0.85, 0.9),
            axis.text.x = element_text(size = 14, face = "bold"),
            axis.text.y = element_text(size = 14, face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.text = element_text(size = 14, face = "bold")) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
pA

# Panel B: human brain enrichment
pB <- ggplot(rbps) +
      geom_density(aes(x = `Brain enrichment`), fill = "#FFEDA0", alpha = 0.5, size = 0.9) +
      xlab("log2 brain/non-brain expression") +
      ylab("Density") +
      theme_cowplot(font_size = 15, 
                    line_size = 1) +
      theme(axis.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 14, face = "bold"),
            axis.text.y = element_text(size = 14, face = "bold")) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
pB

# Panel C: comparison of early and late neuron development enrichment
develop_enrichment <- data.frame(x = c(rbps$`Neuron developmental enrichment (early)`, 
                                       rbps$`Neuron developmental enrichment (late)`), 
                                type = rep(c("Early", "Late"), 
                                           c(length(rbps$`Neuron developmental enrichment (early)`), 
                                             length(rbps$`Neuron developmental enrichment (late)`))))
pC <- ggplot(develop_enrichment, aes(x = x, fill = type)) + 
      geom_density(alpha = 0.5, size = 0.9) +
      scale_fill_brewer(palette = "BuGn", name = " ") +
      xlab("log2 neuron developmental expression") +
      ylab("Density") +
      theme_cowplot(font_size = 15, 
                    line_size = 1,
                    rel_small = 11/15) + 
      theme(legend.position = c(0.85, 0.9),
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 14, face = "bold"),
            axis.text.y = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 14, face = "bold")) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
pC

# Panel D: histogram for the pLI scores
pLI_no_na <- na.omit(rbps[, c("pLI")])
pLI_no_na$breaks <- cut(pLI_no_na$pLI, breaks = c(0, 0.5, 0.9, 0.99, 1))
pD <- ggplot(pLI_no_na, aes(x = pLI, fill = breaks)) +
  geom_histogram(bins = 101, boundary = 0, closed = "right", alpha = 2) +
  scale_fill_manual(values = c("#9ECAE1", "#4292C6", "#08519C", "#08306B"), 
                    name = " ", 
                    labels = c("Low (< 0.5)", "Medium (0.5~0.9)", "High (0.9~0.99)", "Extreme (> 0.99)")) +
  xlab("pLI") +
  ylab("Count") +
  theme_cowplot(font_size = 15, 
                line_size = 1,
                rel_small = 11/15) + 
  theme(legend.position = c(0.7, 0.9),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        plot.margin = unit(c(0.2, 1, 0.2, 1), "cm")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) 
pD


fig2 <- plot_grid(pA, pB, pC, pD, labels = "AUTO", align = "hv", label_size = 20)  
fig2

ggsave("./figures/BrainRBPedia_Figure2.jpg", fig2, 
       width = 400, height = 250, units = c("mm"))

