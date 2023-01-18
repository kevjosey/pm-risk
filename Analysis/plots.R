.libPaths("~/apps/R_4.0.2")

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
scenarios <- expand.grid(dual = c("both", "low", "high"), race = c("all","white", "black"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(0.00783038, 30.92493, length.out = 101)

# data directories
dir_data <- '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_out <- '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR_All/'

ar <- data.frame()
contrast <- data.frame()
hr <- data.frame()

# contrast indexes
idx8 <- which.min(abs(a.vals - 8))
idx9 <- which.min(abs(a.vals - 9))
idx10 <- which.min(abs(a.vals - 10))
idx11 <- which.min(abs(a.vals - 11))
idx12 <- which.min(abs(a.vals - 12))

### Create Data

for (i in 1:nrow(scenarios)) {
  
  # QD
  scenario <- scenarios[i,]
  load(paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
  
  if (i == 1) {
    
    xiao <- data.frame(a.vals = c(est_data$a.vals, est_data$a.vals), 
                       estimate = c(est_data$estimate.lm, est_data$estimate.cal), 
                       gps_method = c(rep("LM", nrow(est_data)), rep("CAL", nrow(est_data))), 
                       lower = c(est_data[,2] - 1.96*est_data[,3],est_data[,5] - 1.96*est_data[,6]),
                       upper = c(est_data[,2] + 1.96*est_data[,3],est_data[,5] + 1.96*est_data[,6]))
    
    save(xiao, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/xiao.RData')
    
  }
  
  # absolute risks
  ar_tmp <- data.frame(a.vals = c(est_data$a.vals), 
                         estimate = c(est_data$estimate.cal),
                         linear = c(est_data$linear.cal),
                         lower = c(est_data[,2] - 1.96*est_data[,3]),
                         upper = c(est_data[,2] + 1.96*est_data[,3]),
                         race = rep(scenario$race, nrow(est_data)),
                         dual = rep(scenario$dual, nrow(est_data)))
  
  # contrasts
  contr_tmp_1 <- as.numeric(est_data[idx12,2]) - as.numeric(est_data[idx11,2])
  contr_tmp_2 <- as.numeric(est_data[idx12,2]) - as.numeric(est_data[idx10,2])
  contr_tmp_3 <- as.numeric(est_data[idx12,2]) - as.numeric(est_data[idx9,2])
  contr_tmp_4 <- as.numeric(est_data[idx12,2]) - as.numeric(est_data[idx8,2])
  contr_se_1 <- sqrt(as.numeric(est_data[idx12,3])^2 + as.numeric(est_data[idx11,3])^2)
  contr_se_2 <- sqrt(as.numeric(est_data[idx12,3])^2 + as.numeric(est_data[idx10,3])^2)
  contr_se_3 <- sqrt(as.numeric(est_data[idx12,3])^2 + as.numeric(est_data[idx9,3])^2)
  contr_se_4 <- sqrt(as.numeric(est_data[idx12,3])^2 + as.numeric(est_data[idx9,3])^2)
  
  contrast_tmp <- data.frame(estimate = c(contr_tmp_1, contr_tmp_2, contr_tmp_3, contr_tmp_4),
                             lower = c(contr_tmp_1 - 1.96*contr_se_1, 
                                       contr_tmp_2 - 1.96*contr_se_2, 
                                       contr_tmp_3 - 1.96*contr_se_3,
                                       contr_tmp_4 - 1.96*contr_se_4),
                             upper = c(contr_tmp_1 + 1.96*contr_se_1, 
                                       contr_tmp_2 + 1.96*contr_se_2, 
                                       contr_tmp_3 + 1.96*contr_se_3,
                                       contr_tmp_4 + 1.96*contr_se_4),
                             pm0 = c(11,10, 9, 8),
                             pm1 = c(12, 12, 12,12),
                             race = scenario$race,
                             dual = scenario$dual)
  
  contrast_tmp$contrast <- paste0(contrast_tmp$pm1, " vs. ", contrast_tmp$pm0)
  
  # hazard ratio
  hr_tmp_1 <- c(as.numeric(est_data[,2])/as.numeric(est_data[idx8,2]))
  hr_tmp_2 <- c(as.numeric(est_data[,2])/as.numeric(est_data[idx9,2]))
  hr_tmp_3 <- c(as.numeric(est_data[,2])/as.numeric(est_data[idx10,2]))
  hr_tmp_4 <- c(as.numeric(est_data[,2])/as.numeric(est_data[idx11,2]))
  hr_tmp_5 <- c(as.numeric(est_data[,2])/as.numeric(est_data[idx12,2]))
  hr_se_1 <- sqrt(c(est_data[,3]^2)/c(est_data[,2]^2) + 
                    c(est_data[idx8,3]^2)/c(est_data[idx8,2]^2)) 
  hr_se_2 <- sqrt(c(est_data[,3]^2)/c(est_data[,2]^2) + 
                    c(est_data[idx9,3]^2)/c(est_data[idx9,2]^2)) 
  hr_se_3 <- sqrt(c(est_data[,3]^2)/c(est_data[,2]^2) + 
                    c(est_data[idx10,3]^2)/c(est_data[idx10,2]^2)) 
  hr_se_4 <- sqrt(c(est_data[,3]^2)/c(est_data[,2]^2) + 
                    c(est_data[idx11,3]^2)/c(est_data[idx11,2]^2)) 
  hr_se_5 <- sqrt(c(est_data[,3]^2)/c(est_data[,2]^2) + 
                    c(est_data[idx12,3]^2)/c(est_data[idx12,2]^2))
  
  hr_tmp <- data.frame(estimate = c(hr_tmp_1, hr_tmp_2, hr_tmp_3, hr_tmp_4, hr_tmp_5),
                       lower = c(exp(log(hr_tmp_1) - 1.96*hr_se_1), 
                                 exp(log(hr_tmp_2) - 1.96*hr_se_2), 
                                 exp(log(hr_tmp_3) - 1.96*hr_se_3),
                                 exp(log(hr_tmp_4) - 1.96*hr_se_4),
                                 exp(log(hr_tmp_5) - 1.96*hr_se_5)),
                       upper = c(exp(log(hr_tmp_1) + 1.96*hr_se_1), 
                                 exp(log(hr_tmp_2) + 1.96*hr_se_2), 
                                 exp(log(hr_tmp_3) + 1.96*hr_se_3),
                                 exp(log(hr_tmp_4) + 1.96*hr_se_4),
                                 exp(log(hr_tmp_5) + 1.96*hr_se_5)),
                       pm0 = c(rep(8, nrow(est_data)),
                               rep(9, nrow(est_data)),
                               rep(10, nrow(est_data)), 
                               rep(11, nrow(est_data)),
                               rep(12, nrow(est_data))),
                       pm1 = rep(est_data[,1], 5),
                       race = scenario$race,
                       dual = scenario$dual)
  
  ar <- rbind(ar, ar_tmp)
  contrast <- rbind(contrast, contrast_tmp)
  hr <- rbind(hr, hr_tmp)
  
}

save(ar, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/ar.RData')
save(contrast, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/contrast.RData')
save(hr, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/hr.RData')

### Proportion Above

# create proportion data
prop <- data.frame()

for (i in 1:nrow(scenarios)) {

  scenario <- scenarios[i,]

  load(paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))

  prop8 <- weighted.mean(individual_data$pm25 > 8, w = individual_data$time_count)
  prop9 <- weighted.mean(individual_data$pm25 > 9, w = individual_data$time_count)
  prop10 <- weighted.mean(individual_data$pm25 > 10, w = individual_data$time_count)
  prop11 <- weighted.mean(individual_data$pm25 > 11, w = individual_data$time_count)

  prop_tmp <- data.frame(estimate = c(prop8, prop9, prop10, prop11),
                         pm = c(8,9,10,11),
                         race = scenario$race,
                         dual = scenario$dual)

  prop <- rbind(prop, prop_tmp)

}

# make nice labels
prop$race_label <- ifelse(prop$race == "all", "", ifelse(prop$race == "white", " White", " Black"))
prop$dual_label <- ifelse(prop$dual == "both", "All Participants", ifelse(prop$dual == "low", "Low SEP", "High SEP"))

# subset important subgroup
prop$label <- paste(prop$dual_label, prop$race_label, sep = "")
prop_dat <- rbind(subset(prop, race != "all" & dual != "both"),
                  subset(prop, race == "all" & dual == "both"))

prop_plot <- prop_dat %>%
  ggplot(aes(x = label, y = estimate, fill = factor(pm))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Proportion at Risk",
       title = ~ "Proportion of Medicare Recipients with Exceeding Concentrations of "*PM[2.5],
       fill = ~ PM[2.5]*" Concentrations") +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  scale_fill_manual(values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40"),
                    labels = c(~ PM[2.5]*" > 8 "*mu*g*"/"*m^3,
                               ~ PM[2.5]*" > 9 "*mu*g*"/"*m^3,
                               ~ PM[2.5]*" > 10 "*mu*g*"/"*m^3,
                               ~ PM[2.5]*" > 11 "*mu*g*"/"*m^3))

pdf(file = "/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Figures/prop_plot.pdf", width = 10, height = 8)
prop_plot
dev.off()

### Main Plot

scenario <- scenarios[1,]

# histogram and ERF data
load(paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
a_dat <- rep(individual_data$pm25, individual_data$time_count)
hr_tmp <- subset(hr, dual == scenario$dual & race == scenario$race)
hr_tmp <- subset(hr_tmp, pm0 %in% c(8, 9, 10, 11))

# exposure response curve
erf_plot <- hr_tmp %>%
  ggplot(aes(x = pm1, y = estimate, color = factor(pm0))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_segment(x = 8, y = 0.85, xend = 8, yend = 1, linetype = "dashed", color = "#D81B60") +
  geom_segment(x = 9, y = 0.85, xend = 9, yend = 1, linetype = "dashed", color = "#1E88E5") +
  geom_segment(x = 10, y = 0.85, xend = 10, yend = 1, linetype = "dashed", color = "#FFC107") +
  geom_segment(x = 11, y = 0.85, xend = 11, yend = 1, linetype = "dashed", color = "#004D40") +
  theme_bw() +
  coord_cartesian(xlim = c(5,15), ylim = c(0.925,1.05)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Hazard Ratio",
       title = "Hazard Ratios for All Medicare Recipients", color = ~ "Baseline "*PM[2.5]) +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = c(0.925,0.95,0.975,1,1.025,1.05)) +
  scale_x_continuous(breaks = c(5,6,7,8,9,10,11,12,13,14,15)) +
  scale_color_manual(values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40"),
                     labels = c(~ "8 "*mu*g*"/"*m^3, ~ "9 "*mu*g*"/"*m^3,
                                ~ "10 "*mu*g*"/"*m^3, ~ "11 "*mu*g*"/"*m^3))

# histogram
a_hist <- ggplot(data.frame(a = a_dat), mapping = aes(x = a_dat)) +
  geom_density(fill = "grey", alpha = 0.3, adjust = 3) +
  theme_bw() +
  coord_cartesian(xlim = c(5,15), ylim = c(0,0.15)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = ~PM[2.5]*" Density") +
  scale_y_continuous(position = "right", breaks = c(0, 0.05, 0.10, 0.15)) +
  guides(fill = "none") +
  theme_cowplot()

align <- align_plots(a_hist, erf_plot, align = "hv", axis = "tblr")
main_plot <- ggdraw(align[[1]]) + draw_plot(align[[2]])

pdf(file = "/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Figures/erc_plot.pdf", width = 8, height = 8)
main_plot
dev.off()

### ERC by Race

plot_list <- list()
dual.vals <- c(2, 0, 1)

hr$dual_label <- ifelse(hr$dual == "high", "High SEP", ifelse(hr$dual == "low", "Low SEP", "High + Low SEP"))
hr$dual_label <- factor(hr$dual_label, levels = c("High + Low SEP", "High SEP", "Low SEP"))
hr$race <- str_to_title(hr$race)

hr_tmp <- subset(hr, race %in% c("All", "Black", "White") & pm0 == 8)

erf_strata <- hr_tmp %>%
  ggplot(aes(x = pm1, y = estimate, color = factor(race))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_line(size = 1) +
  geom_segment(x = 8, y = 0.8, xend = 8, yend = 1, linetype = "dashed", color = "#D81B60") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ dual_label) +
  coord_cartesian(xlim = c(5.25,14.75), ylim = c(0.85, 1.15)) +
  theme_bw() +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Hazard Ratio",
       color = "Race", title = ~ "Stratified Hazard Ratio Estimates Evaluated with a Baseline "*PM[2.5]*" Concentration of "*8*" "*mu*g*"/"*m^3) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#000000", "#367E18", "#F57328")) +
  scale_y_continuous(breaks = c(0.85,0.9,0.95,1,1.05,1.1,1.15)) +
  scale_x_continuous(breaks = c(5,6,7,8,9,10,11,12,13,14,15))

pdf(file = "/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Figures/erc_strata.pdf", width = 16, height = 8)
erf_strata
dev.off()

### Contrast Plot

contr <- subset(contrast, !(dual %in% c("high","low") & race == "all"))
contr$contrast <- factor(contr$contrast, levels = c("12 vs. 11", "12 vs. 10", "12 vs. 9", "12 vs. 8"))
contr$race_dual <- paste0(str_to_title(contr$race), ifelse(contr$dual == "high", " -\n High SEP",
                                                           ifelse(contr$dual == "low", " -\n Low SEP",
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
  theme_bw() +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = round(seq(0, max(100*contr$upper), by = 0.1),1))

pdf(file = "/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Figures/contrast_plot.pdf", width = 10, height = 8)
contrast_plot
dev.off()
