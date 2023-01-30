library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(cobalt)

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"), race = c("all", "white", "black"))
# scenarios <- expand.grid(dual = c("high", "low"), race = c("asian", "hispanic", "other"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/calibrate.R')
set.seed(42)

### Plot G-Computation, PS Regression, and Calibration Weight Results

# data directories
dir_data <- '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_mod <- '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/GPSReg_All/'
# dir_mod <- '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/GComp_All/'
# dir_mod <- '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/CAL_All/'

ar <- data.frame()
contrast <- data.frame()
hr <- data.frame()

# contrast indexes
idx8 <- which.min(abs(a.vals - 8))
idx9 <- which.min(abs(a.vals - 9))
idx10 <- which.min(abs(a.vals - 10))
idx11 <- which.min(abs(a.vals - 11))
idx12 <- which.min(abs(a.vals - 12))

# Load Data
for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_mod, scenario$dual, "_", scenario$race, ".RData"))
  
  # absolute risks
  ar_tmp <- data.frame(a.vals = c(a.vals), 
                       estimate = mhat.vals,
                       race = rep(scenario$race, length(mhat.vals)),
                       dual = rep(scenario$dual, length(mhat.vals)))
  
  # contrasts
  contr_tmp_11 <- as.numeric(mhat.vals[idx12]) - as.numeric(mhat.vals[idx11])
  contr_tmp_10 <- as.numeric(mhat.vals[idx12]) - as.numeric(mhat.vals[idx10])
  contr_tmp_9 <- as.numeric(mhat.vals[idx12]) - as.numeric(mhat.vals[idx9])
  contr_tmp_8 <- as.numeric(mhat.vals[idx12]) - as.numeric(mhat.vals[idx8])
  
  contrast_tmp <- data.frame(estimate = c(contr_tmp_11, contr_tmp_10, contr_tmp_9, contr_tmp_8),
                             pm0 = c(11,10, 9, 8),
                             pm1 = c(12, 12, 12,12),
                             race = scenario$race,
                             dual = scenario$dual)
  
  contrast_tmp$contrast <- paste0(contrast_tmp$pm1, " vs. ", contrast_tmp$pm0)
  
  # hazard ratio
  hr_tmp_8 <- c(as.numeric(mhat.vals)/as.numeric(mhat.vals[idx8]))
  hr_tmp_9 <- c(as.numeric(mhat.vals)/as.numeric(mhat.vals[idx9]))
  hr_tmp_10 <- c(as.numeric(mhat.vals)/as.numeric(mhat.vals[idx10]))
  hr_tmp_11 <- c(as.numeric(mhat.vals)/as.numeric(mhat.vals[idx11]))
  hr_tmp_12 <- c(as.numeric(mhat.vals)/as.numeric(mhat.vals[idx12]))
  
  hr_tmp <- data.frame(estimate = c(hr_tmp_8, hr_tmp_9, hr_tmp_10, hr_tmp_11, hr_tmp_12),
                       pm0 = c(rep(8, length(mhat.vals)),
                               rep(9, length(mhat.vals)),
                               rep(10, length(mhat.vals)), 
                               rep(11, length(mhat.vals)),
                               rep(12, length(mhat.vals))),
                       pm1 = rep(a.vals, 5),
                       race = scenario$race,
                       dual = scenario$dual)
  
  ar <- rbind(ar, ar_tmp)
  contrast <- rbind(contrast, contrast_tmp)
  hr <- rbind(hr, hr_tmp)
  
}

## Overall ERC

hr_tmp <- subset(hr, dual == "both" & race == "all" & pm0 == 12)

# exposure response curve
hr_plot <- hr_tmp %>%
  ggplot(aes(x = pm1, y = estimate)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_segment(x = 8, y = 0.7, xend = 8, yend = hr_tmp$estimate[idx8], linetype = "dashed", color = "#D81B60") +
  geom_segment(x = 9, y = 0.7, xend = 9, yend = hr_tmp$estimate[idx9], linetype = "dashed", color = "#1E88E5") +
  geom_segment(x = 10, y = 0.7, xend = 10, yend = hr_tmp$estimate[idx10], linetype = "dashed", color = "#FFC107") +
  geom_segment(x = 11, y = 0.7, xend = 11, yend = hr_tmp$estimate[idx11], linetype = "dashed", color = "#004D40") +
  geom_segment(x = 12, y = 0.7, xend = 12, yend = hr_tmp$estimate[idx12], linetype = "dotted") +
  geom_segment(x = 3, y = hr_tmp$estimate[idx8], xend = 8, yend = hr_tmp$estimate[idx8], linetype = "dashed", color = "#D81B60") +
  geom_segment(x = 3, y = hr_tmp$estimate[idx9], xend = 9, yend = hr_tmp$estimate[idx9], linetype = "dashed", color = "#1E88E5") +
  geom_segment(x = 3, y = hr_tmp$estimate[idx10], xend = 10, yend = hr_tmp$estimate[idx10], linetype = "dashed", color = "#FFC107") +
  geom_segment(x = 3, y = hr_tmp$estimate[idx11], xend = 11, yend = hr_tmp$estimate[idx11], linetype = "dashed", color = "#004D40") +
  theme_bw() +
  coord_cartesian(xlim = c(5,13), 
                  ylim = c(0.88, 1.02)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Hazard Ratio",
       title = "GPS Regression Adjustment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = c(0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,
                                0.96,0.97,0.98,0.99,1,1.01,1.02)) +
  scale_x_continuous(breaks = c(5,6,7,8,9,10,11,12,13))

pdf(file = "~/Figures/hr_plot_gpsreg.pdf", width = 10, height = 8)
hr_plot
dev.off()

## Stratified by Race x SEP

hr$dual_label <- ifelse(hr$dual == "high", "High SEP", ifelse(hr$dual == "low", "Low SEP", "All"))
hr$dual_label <- factor(hr$dual_label, levels = c("All", "High SEP", "Low SEP"))
hr$race_title <- str_to_title(hr$race)

hr_tmp <- subset(hr, race %in% c("black", "white") & pm0 == 12)

hr_strata <- hr_tmp %>%
  ggplot(aes(x = pm1, y = estimate, color = factor(race_title))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_segment(x = 12, y = 0.7, xend = 12, yend = hr_tmp$estimate[idx12], linetype = "dotted") +
  facet_wrap(~ dual_label) +
  coord_cartesian(xlim = c(5,13), 
                  ylim = c(0.85, 1.02)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Hazard Ratio",
       color = "Race",
       title = "GPS Regression Adjustment") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#367E18", "#F57328")) +
  scale_y_continuous(breaks = c(0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,
                                0.94,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02)) +
  scale_x_continuous(breaks = c(5,6,7,8,9,10,11,12,13))

pdf(file = "~/Figures/hr_strata_gpsreg.pdf", width = 16, height = 8)
hr_strata
dev.off()

## Compare with Xiao's Plot

DR100 <- read.csv("~/Data/DR100_table.csv")[,-1]
hr_tmp <- subset(hr, dual == "both" & race == "all" & pm0 == 12)

Kevin <- data.frame(a.vals = hr_tmp$pm1, 
                    estimate = hr_tmp$estimate, 
                    gps_method = rep("G-Computation", nrow(hr_tmp)), 
                    lower = rep(NA, nrow(hr_tmp)),
                    upper = rep(NA, nrow(hr_tmp)))

colnames(Kevin) = colnames(DR100)

DR100_table <- rbind(DR100, Kevin)
DR100_table <- subset(DR100_table, !(Methods %in% c("Doubly Robust")))

# standard errors
DR100_table$SE <- (DR100_table$ERC - DR100_table$lower)/1.96

# hazard ratios

plot_dat <- data.frame()

for (j in unique(DR100_table$Methods)) {
  
  tmp <- subset(DR100_table, Methods == j)
  idx0 <- which.min(abs(tmp$a.vals - 12))
  tmp$hr <- c(as.numeric(tmp$ERC)/as.numeric(tmp$ERC[idx0]))
  tmp$log.hr.se <- sqrt(c(tmp$SE^2)/c(tmp$ERC^2) + 
                          c(tmp$SE[idx0]^2)/c(tmp$ERC[idx0]^2))
  tmp$hr.lower <- exp(log(tmp$hr) - 1.96*tmp$log.hr.se)
  tmp$hr.upper <- exp(log(tmp$hr) + 1.96*tmp$log.hr.se)
  plot_dat <- rbind(plot_dat, tmp)
  
}

# Hazard Ratio
# exposure response curve
hr_compare <- ggplot(data=plot_dat, aes(x=a.vals, y = hr, color = Methods)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = hr.lower, ymax = hr.upper), alpha = 0.2, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = c(0.7, 0.2),
        legend.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Hazard Ratio",
       title = "GPS Regression Adjustment") +
  coord_cartesian(xlim = c(5,13), 
                  ylim = c(0.88, 1.02)) +
  scale_color_manual(values=c("blue", "red"),
                     labels = c("G-Computation", "Matching"),
                     breaks = c("G-Computation", "Matching")) +
  scale_y_continuous(breaks = c(0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02)) +
  scale_x_continuous(breaks = c(5,6,7,8,9,10,11,12,13))

pdf(file = "~/Figures/hr_compare_gpsreg.pdf", width = 10, height = 8)
hr_compare
dev.off()