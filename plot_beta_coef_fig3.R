# ---
# Title: plot_beta_coef_fig3.R
# Purpose: This script generates "Figure 3: Informative RBP and non-RBP features 
#          for ASD and ID susceptibility gene prediction obtained from the penalized 
#          multivariate logistic regression models" in the BrainRBPedia manuscript.
# Prerequisite: Run nested_cv.R and nested_cv_nonRBPs.R to obtain objects "autism_res", 
#               "id_res", "autism_res2", and "id_res2" before running this script.
# ---

library(data.table)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

# Brain enrichment
# Neuron developmental enrichment (early)
# Neuron developmental enrichment (late)
# Extreme pLI
# High pLI
# Medium pLI
# Neuron enrichment (mouse)
# Neuron enrichment (human)

get_beeswarm_tb <- function(res, K) {
  
  df <- as.data.frame(res[["betaCoefs"]][1])
  names(df) <- "fold1"
  df$variable <- rownames(df)
  df <- df[, c(2, 1)]
  
  for (i in 2:K) {
    df2 <- as.data.frame(res[["betaCoefs"]][i])
    names(df2) <- paste0("fold", i)
    df2$variable <- rownames(df2)
    df2 <- df2[, c(2:1)]
    df <- merge(df, df2, 
                by = "variable",
                all.y = T)
  }
  
  df[is.na(df)] <- 0
  
  coefs <- numeric((nrow(df)-1) * 5)
  variables <- character((nrow(df)-1) * 5)
  
  for (i in (2:nrow(df))) {
    ifelse(i == 2, 
           coefs <- as.numeric(df[c(i), c(2:ncol(df))]),
           coefs <- c(coefs, as.numeric(df[c(i), c(2:ncol(df))])))
    
    ifelse(i == 2,
           variables <- rep(df[i, 1], 5),
           variables <- c(variables, rep(df[i, 1], 5)))
  }
  
  tb <- data.frame(coef = coefs,
                   variable = variables)
  
  return(tb)
}

# NOTE: autism_tb is for RBPs; autism_tb2 is for non-RBPs
autism_tb <- get_beeswarm_tb(autism_res, K = 5)
autism_tb2 <- get_beeswarm_tb(autism_res2, K = 5)

autism_tb <- rbind(autism_tb, 
                   data.frame(coef = rep(0, 5),
                              variable = rep("Brain enrichment", 5) ))
autism_tb <- rbind(autism_tb, 
                   data.frame(coef = rep(0, 5),
                              variable = rep("Neuron developmental enrichment (early)", 5) ))
autism_tb <- autism_tb[order(autism_tb$variable),]
autism_tb$group <- "RBP"

autism_tb2 <- autism_tb2[!autism_tb2$variable == "Neuron enrichment (mouse)", ]
autism_tb2$group <- "non-RBP"
autism_all <- rbind(autism_tb, autism_tb2)
autism_all$variable[autism_all$variable == "Neuron developmental enrichment (early)"] <- "Early developmental enrichment"
autism_all$variable[autism_all$variable == "Neuron developmental enrichment (late)"] <- "Late developmental enrichment"
autism_all$odds_ratio <- exp(autism_all$coef)

autism_CIs <- autism_all %>% group_by(variable, group) %>% summarize(Mean = mean(odds_ratio)) %>% ungroup()
autism_all <- autism_all[order(autism_all$variable), ]
groups <- seq(1, 50, by = 5)
CI_low <- numeric(0)
CI_high <- numeric(0)
for (i in groups) {
  members <- autism_all$odds_ratio[i:(i+5-1)]
  if (sum(members) == 5) {
    CI_low <- c(CI_low, 0)
    CI_high <- c(CI_high, 0)
    next
  }
  ci <- as.numeric(t.test(members, conf.level = 0.95)$conf.int)
  CI_low <- c(CI_low, ci[1])
  CI_high <- c(CI_high, ci[2])
}

autism_CIs$CI_low <- CI_low
autism_CIs$CI_high <- CI_high
autism_CIs <- as.data.frame(autism_CIs)
autism_CIs$group <- factor(autism_CIs$group, levels = c("RBP", "non-RBP"))
autism_all$group <- factor(autism_all$group, levels = c("RBP", "non-RBP"))

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p_aut <- ggplot(autism_CIs, aes(x = variable, color = as.factor(group))) +
  geom_beeswarm(data = autism_all, aes(x = variable, y = odds_ratio, fill = as.factor(group)), size = 2, dodge.width = 0.6, shape = 21, color = "black") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), position = position_dodge(0.6), width = 0) +
  geom_point(aes(y = Mean), position = position_dodge(0.6), shape = 3, size = 5) + # to add mean to the error bars
  ylab("Odds ratio") +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  scale_x_discrete(limits = c("Extreme pLI", 
                              "High pLI", 
                              "Late developmental enrichment", 
                              "Brain enrichment",
                              "Early developmental enrichment")) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x.bottom = element_text(face = "bold", size = 9, color = "black"),
        axis.text.y.left = element_text(face = "bold", size = 10, color = "black"),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.82),
        legend.background = element_blank(),
        legend.text = element_text(face = "bold", size = 12)) 
p_aut

id_tb <- get_beeswarm_tb(id_res, K = 5)
id_tb2 <- get_beeswarm_tb(id_res2, K = 5)

id_tb <- rbind(id_tb, 
               data.frame(coef = rep(0, 5),
                          variable = rep("High pLI", 5) ))
id_tb2 <- rbind(id_tb2, 
               data.frame(coef = rep(0, 5),
                          variable = rep("Neuron developmental enrichment (early)", 5) ))
id_tb <- id_tb[order(id_tb$variable), ]
id_tb$group <- "RBP"
id_tb2$group <- "non-RBP"
id_all <- rbind(id_tb, id_tb2)
id_all$odds_ratio <- exp(id_all$coef)

id_CIs <- id_all %>% group_by(variable, group) %>% summarize(Mean = mean(odds_ratio)) %>% ungroup()
id_all <- id_all[order(id_all$variable), ]
id_CIs$group <- factor(id_CIs$group, levels = c("RBP", "non-RBP"))
id_all$group <- factor(id_all$group, levels = c("RBP", "non-RBP"))
groups <- seq(1, 30, by = 5)
CI_low <- numeric(0)
CI_high <- numeric(0)
for (i in groups) {
  members <- id_all$odds_ratio[i:(i+5-1)]
  if (sum(members) == 5) {
    CI_low <- c(CI_low, 0)
    CI_high <- c(CI_high, 0)
    next
  }
  ci <- as.numeric(t.test(members, conf.level = 0.95)$conf.int)
  CI_low <- c(CI_low, ci[1])
  CI_high <- c(CI_high, ci[2])
}

id_CIs$CI_low <- CI_low
id_CIs$CI_high <- CI_high
id_CIs <- as.data.frame(id_CIs)

p_id <- ggplot(id_CIs, aes(x = variable, color = as.factor(group))) +
  geom_beeswarm(data = id_all, aes(x = variable, y = odds_ratio, fill = as.factor(group)), size = 2, dodge.width = 0.5, shape = 21, color = "black") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), position = position_dodge(0.5), width = 0) +
  geom_point(aes(y = Mean), position = position_dodge(0.5), shape = 3, size = 5) + # to add mean to the error bars
  ylab("Odds ratio") +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x.bottom = element_text(face = "bold", size = 9, color = "black"),
        axis.text.y.left = element_text(face = "bold", size = 10, color = "black"),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.82),
        legend.background = element_blank(),
        legend.text = element_text(face = "bold", size = 12)) 
p_id

fig3 <- plot_grid(p_aut, p_id, nrow = 2, ncol = 1, labels = c("A", "B"))
fig3
ggsave("./figures/BrainRBPedia_Figure3.jpg", fig3, 
       width = 300, height = 200, units = c("mm"))
