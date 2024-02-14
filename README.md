# :brain: BrainRBPedia

This repository contains the data and scripts required by the BrainRBPedia database.

## :desktop_computer: Instructions
1. Download this repository to your local machine.
2. Files used to generate human and neuron markers are too large to be stored in this repository - please download folder <code>visp_data</code> and <code>Seu_AIBS_obj_update_07JUN21.rds</code> from [Google Drive](https://drive.google.com/drive/folders/15scJPWOfc_yvkqiqA6BZUYtM4y6oiYpI?usp=sharing) and add them to the <code>data</code> folder.

## :green_book: Description

:dna: **make_BrainRBPedia_dataset.R**
* Use this script to generate the BrainRBPedia dataset from various publicly accessible databases.
* Note that you might not be able to generate the current version of BrainRBPedia using this script. This is because features such as “Protein domains” are generated using mart ensembl in library(biomaRt), which is updated periodically. 
* I recommend using get_human_neuron_markers.R and get_mouse_neuron_markers.R to re-generate the neuron markers as I wasn’t able to use the most optimal parameters due to my machine’s computational limit (we can discuss more about this when we meet).

:dna: **nested_cv.R**
* This script uses nested cross-validation to construct the logistic regression model that predicts ASD and ID susceptibility.

:dna: **get_human_neuron_markers.R**
* This script uses the Seurat R package to generate each RBP’s log2 fold-change of gene expression in human neurons relative to non-neurons.

:dna: **get_mouse_neuron_markers.R**
* This script uses the Seurat R package to generate each RBP’s log2 fold-change of gene expression in mouse neurons relative to non-neurons.

:dna: **plot_distribution_fig2.R**
* This script generates “Figure 2: Distribution of key BrainRBPedia functional annotations” in the BrainRBPedia manuscript.

:dna: **plot_beta_coef_fig3.R**
* This script generates “Figure 3: Informative RBP and non-RBP features for ASD and ID susceptibility gene prediction obtained from the penalized multivariate logistic regression models” in the BrainRBPedia manuscript.

:dna: **plot_auc_curves_fig4.R**
* This script generates “Figure 4: ASD and ID model performance measured by ROC and PR curves” in the BrainRBPedia manuscript.


