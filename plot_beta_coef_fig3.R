# ---
# Title: plot_beta_coef_fig3.R
# Purpose: This script generates "Figure 3: Informative RBP and non-RBP features 
#          for ASD and ID susceptibility gene prediction obtained from the penalized 
#          multivariate logistic regression models" in the BrainRBPedia manuscript.
# Prerequisite: Run nested_cv.R to obtain objects "autism_res" and "id_res" before
#               running this script.
# ---

library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(xlsx)

setwd("~/Desktop/BrainRBPedia") # CHANGE THIS


# SETUP -------------------------------------------------------------------

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

# Run nested_cv.R to get autism_res and id_res
autism_tb <- get_beeswarm_tb(autism_res, K = 5)
id_tb <- get_beeswarm_tb(id_res, K = 5)

autism_tb <- rbind(autism_tb, 
                   data.frame(coef = rep(0, 5),
                              variable = rep("Neuron developmental enrichment (early)", 5) ))
id_tb <- rbind(id_tb, 
               data.frame(coef = rep(0, 10),
                          variable = c( rep("high_pLI", 5),
                                        rep("Neuron developmental enrichment (late)", 5) ) ))

autism_tb$disorder <- "ASD"
id_tb$disorder <- "ID"

all_tb <- rbind(autism_tb, id_tb)
all_tb$variable <- gsub("extreme_pLI", "Extreme pLI", all_tb$variable)
all_tb$variable <- gsub("high_pLI", "High pLI", all_tb$variable)
all_tb$variable <- gsub("Neuron developmental enrichment [[:punct:]]early[[:punct:]]", "Early neuron development", all_tb$variable)
all_tb$variable <- gsub("Neuron developmental enrichment [[:punct:]]late[[:punct:]]", "Late neuron development", all_tb$variable)

all_tb$variable <- factor(all_tb$variable, 
                          levels = c("Extreme pLI", 
                                     "High pLI", 
                                     "Late neuron development", 
                                     "Early neuron development"))
all_tb$coef <- exp(all_tb$coef)
names(all_tb)[1] <- "Odds ratio"

CIs <- all_tb %>% group_by(variable, disorder) %>% summarize(Mean = mean(`Odds ratio`)) %>% ungroup()
groups <- seq(1, 40, by = 5)
CI_low <- numeric(0)
CI_high <- numeric(0)
for (i in groups) {
  members <- all_tb$`Odds ratio`[i:(i+5-1)]
  if (sum(members) == 5) {
    CI_low <- c(CI_low, 0)
    CI_high <- c(CI_high, 0)
    next
  }
  ci <- as.numeric(t.test(members, conf.level = 0.95)$conf.int)
  CI_low <- c(CI_low, ci[1])
  CI_high <- c(CI_high, ci[2])
}

CIs$CI_low <- CI_low[c(1,5,2,7,3,8,4,6)]
CIs$CI_high <- CI_high[c(1,5,2,7,3,8,4,6)]

CIs <- as.data.frame(CIs)

write.xlsx(all_tb, "./tables/TableS4_odds_ratios.xlsx", sheetName = "Sheet 1", row.names = FALSE)
write.xlsx(CIs, "./tables/TableS4_odds_ratios.xlsx", sheetName = "Sheet 2", 
           row.names = FALSE, append = TRUE)

# 2-in-1 Beeswarm Plot ----------------------------------------------------

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

level_order <- c("Extreme pLI", "High pLI", "Late neuron development", "Early neuron development")
CIs$variable <- factor(CIs$variable, level = level_order)

fig3 <- ggplot(CIs, aes(x = variable, color = as.factor(disorder))) +
  geom_beeswarm(data = all_tb, aes(x = variable, y = `Odds ratio`, fill = as.factor(disorder)), size = 2, dodge.width = 0.5, shape = 21, color = "black") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), position = position_dodge(0.5), width = 0) +
  geom_point(aes(y = Mean), position = position_dodge(0.5), shape = 3, size = 5) + # to add mean to the error bars
  # geom_point(data = all_tb, aes(x = variable, y = `Odds ratio`, color = as.factor(disorder))) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  scale_x_discrete(limits = level_order) +
  scale_y_continuous(breaks = seq(0, 9, by = 1),
                     limits = c(0.5, 5.5)) +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x.bottom = element_text(face = "bold", size = 14, color = "black"),
        axis.text.y.left = element_text(face = "bold", size = 14, color = "black"),
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
        legend.text = element_text(face = "bold", size = 14)) 
fig3

ggsave("./figures/BrainRBPedia_Figure3.jpg", fig3, width = 300, height = 200, units = c("mm"))
