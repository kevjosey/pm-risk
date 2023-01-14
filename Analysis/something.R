library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(splines)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("all","white", "black"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

load("/Users/kevin/Library/CloudStorage/Dropbox/Projects/ERC-Strata/Output/estimate.RData")

contr <- data.frame()
a.vals <- dat$a.vals
dat$se <- (dat$upper - dat$estimate)/1.96

# contrast indexes

### Create Data

for (i in 1:nrow(scenarios)) {
  
  # QD
  scenario <- scenarios[i,]
  
  dat_tmp <- subset(dat, race == scenario$race & dual == scenario$dual)
  
  idx8 <- which.min(abs(dat_tmp$a.vals - 8))
  idx9 <- which.min(abs(dat_tmp$a.vals - 9))
  idx10 <- which.min(abs(dat_tmp$a.vals - 10))
  idx12 <- which.min(abs(dat_tmp$a.vals - 12))
  
  tmp_1 <- as.numeric(dat_tmp$estimate[idx12]) - as.numeric(dat_tmp$estimate[idx10])
  tmp_2 <- as.numeric(dat_tmp$estimate[idx12]) - as.numeric(dat_tmp$estimate[idx9])
  tmp_3 <- as.numeric(dat_tmp$estimate[idx12]) - as.numeric(dat_tmp$estimate[idx8])
  tmp_4 <- sqrt(as.numeric(dat_tmp$se[idx12])^2 + as.numeric(dat_tmp$se[idx10])^2)
  tmp_5 <- sqrt(as.numeric(dat_tmp$se[idx12])^2 + as.numeric(dat_tmp$se[idx9])^2)
  tmp_6 <- sqrt(as.numeric(dat_tmp$se[idx12])^2 + as.numeric(dat_tmp$se[idx8])^2)
  
  contr_tmp <- data.frame(estimate = c(tmp_1, tmp_2, tmp_3),
                          lower = c(tmp_1 - 1.96*tmp_4, tmp_2 - 1.96*tmp_5, tmp_3 - 1.96*tmp_6),
                          upper = c(tmp_1 + 1.96*tmp_4, tmp_2 + 1.96*tmp_5, tmp_3 + 1.96*tmp_6),
                          pm0 = c(10, 9, 8),
                          pm1 = c(12, 12, 12),
                          race = scenario$race,
                          dual = scenario$dual)
  
  contr_tmp$contrast <- paste0(contr_tmp$pm1, " vs. ", contr_tmp$pm0)
  contr <- rbind(contr, contr_tmp)
  
}

### Contrast Plot

contr <- subset(contr, !(dual %in% c(0,1) & race == "all"))
contr$contrast <- factor(contr$contrast, levels = c("12 vs. 10", "12 vs. 9", "12 vs. 8"))
contr$race_dual <- paste0(str_to_title(contr$race), ifelse(contr$dual == 0, " -\n High SEP", 
                                                           ifelse(contr$dual == 1, " -\n Low SEP",
                                                                  " -\n High + Low SEP")))

contr$race_dual <- ifelse(contr$race_dual == "All -\n High + Low SEP", "All - High\n + Low SEP", contr$race_dual)
contr$race_dual <- ifelse(contr$race_dual == "White -\n High + Low SEP", "White - High\n + Low SEP", contr$race_dual)
contr$race_dual <- ifelse(contr$race_dual == "Black -\n High + Low SEP", "Black - High\n + Low SEP", contr$race_dual)
contr$race_dual <- factor(contr$race_dual, levels = c("All - High\n + Low SEP",
                                                      "White - High\n + Low SEP", 
                                                      "Black - High\n + Low SEP", 
                                                      "White -\n High SEP", "Black -\n High SEP",
                                                      "White -\n Low SEP",  "Black -\n Low SEP"))

contrast_plot <- contr %>% 
  ggplot(aes(x = race_dual, y = 100*estimate, color = contrast)) + 
  geom_pointrange(aes(ymin = 100*lower, ymax = 100*upper), position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  labs(x = "", y = "Risk Difference (%)", title = "Risk Difference Estimates") +
  guides(color = guide_legend(title = ~ PM[2.5]*" Contrasts ("*mu*g*"/"*m^3*")")) +
  theme(legend.position = c(0.17, 0.87),
        legend.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#008080", "#FF00FF","#FFD700")) +
  scale_y_continuous(breaks = round(seq(0, max(100*contr$upper), by = 0.1),1)) +
  grids(linetype = "dashed")

contrast_plot
