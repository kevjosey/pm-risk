library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)

# scenarios
scenarios <- expand.grid(dual = c("both", "low", "high"), race = c("all","white", "black"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)

### Plot AR and HR Exposure Responses

# data directories
dir_data <- '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_out <- '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR_All/'

pm <- data.frame()
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
  
  scenario <- scenarios[i,]
  load(paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
  load(paste0(dir_out, scenario$dual, "_", scenario$race, "_boot.RData"))
  
  # Replace with Bootstrap SEs
  u.zip <- unique(individual_data$zip)
  m <- 1/log(length(u.zip)) # for m out of n bootstrap
  est_data$se.cal_trunc <- apply(boot_mat, 2, sd, na.rm = T)*sqrt(m)
  perm_mat <- boot_mat[sample(1:nrow(boot_mat), size = nrow(boot_mat), replace = FALSE), ]
  
  if (i == 1) {
    
    xiao <- data.frame(a.vals = rep(est_data$a.vals, 3), 
                       estimate = c(est_data$estimate.lm, est_data$estimate.cal,
                                    est_data$estimate.cal_trunc), 
                       gps_method = c(rep("LM", nrow(est_data)), 
                                      rep("CAL", nrow(est_data)),
                                      rep("CAL_TRUNC", nrow(est_data))), 
                       lower = c(est_data[,2] - 1.96*est_data[,3],
                                 est_data[,4] - 1.96*est_data[,5],
                                 est_data[,6] - 1.96*est_data[,7]),
                       upper = c(est_data[,2] + 1.96*est_data[,3],
                                 est_data[,4] + 1.96*est_data[,5],
                                 est_data[,6] + 1.96*est_data[,7]))
    
    save(xiao, file = '~/Data/xiao.RData')
    
  }
  
  # average PM x Year
  year <- individual_data$year
  list.pm <- split(data.frame(pm = individual_data$pm25, wts = individual_data$time_count), year)
  pm_tmp <- data.frame(pm = do.call(c, lapply(list.pm, function(df) sum(df$pm*df$wts)/sum(df$wts))),
                       year = names(list.pm), dual = scenario$dual, race = scenario$race)
  
  # absolute risks
  ar_tmp <- data.frame(a.vals = c(est_data$a.vals), 
                       estimate = c(est_data[,6]),
                       lower = c(est_data[,6] - 1.96*est_data[,7]),
                       upper = c(est_data[,6] + 1.96*est_data[,7]),
                       race = rep(scenario$race, nrow(est_data)),
                       dual = rep(scenario$dual, nrow(est_data)))
  
  # hazard ratio
  hr_tmp_8 <- c(as.numeric(est_data[,6])/as.numeric(est_data[idx8,6]))
  hr_tmp_9 <- c(as.numeric(est_data[,6])/as.numeric(est_data[idx9,6]))
  hr_tmp_10 <- c(as.numeric(est_data[,6])/as.numeric(est_data[idx10,6]))
  hr_tmp_11 <- c(as.numeric(est_data[,6])/as.numeric(est_data[idx11,6]))
  hr_tmp_12 <- c(as.numeric(est_data[,6])/as.numeric(est_data[idx12,6]))
  
  # bootstrap standard errors
  log_hr_se_8 <- apply(log(boot_mat/perm_mat[,idx8]), 2, sd, na.rm = TRUE)*sqrt(m)
  log_hr_se_9 <- apply(log(boot_mat/perm_mat[,idx9]), 2, sd, na.rm = TRUE)*sqrt(m)
  log_hr_se_10 <- apply(log(boot_mat/perm_mat[,idx10]), 2, sd, na.rm = TRUE)*sqrt(m)
  log_hr_se_11 <- apply(log(boot_mat/perm_mat[,idx11]), 2, sd, na.rm = TRUE)*sqrt(m)
  log_hr_se_12 <- apply(log(boot_mat/perm_mat[,idx12]), 2, sd, na.rm = TRUE)*sqrt(m)
  
  # delta method standard errors
  # log_hr_se_8 <- sqrt(c(est_data[,7]^2)/c(est_data[,6]^2) +
  #                       c(est_data[idx8,7]^2)/c(est_data[idx8,6]^2))
  # log_hr_se_9 <- sqrt(c(est_data[,7]^2)/c(est_data[,6]^2) +
  #                       c(est_data[idx9,7]^2)/c(est_data[idx9,6]^2))
  # log_hr_se_10 <- sqrt(c(est_data[,7]^2)/c(est_data[,6]^2) +
  #                        c(est_data[idx10,7]^2)/c(est_data[idx10,6]^2))
  # log_hr_se_11 <- sqrt(c(est_data[,7]^2)/c(est_data[,6]^2) +
  #                        c(est_data[idx11,7]^2)/c(est_data[idx11,6]^2))
  # log_hr_se_12 <- sqrt(c(est_data[,7]^2)/c(est_data[,6]^2) +
  #                        c(est_data[idx12,7]^2)/c(est_data[idx12,6]^2))
  
  hr_tmp <- data.frame(estimate = c(hr_tmp_8, hr_tmp_9, hr_tmp_10, hr_tmp_11, hr_tmp_12),
                       lower = c(exp(log(hr_tmp_8) - 1.96*log_hr_se_8), 
                                 exp(log(hr_tmp_9) - 1.96*log_hr_se_9), 
                                 exp(log(hr_tmp_10) - 1.96*log_hr_se_10),
                                 exp(log(hr_tmp_11) - 1.96*log_hr_se_11),
                                 exp(log(hr_tmp_12) - 1.96*log_hr_se_12)),
                       upper = c(exp(log(hr_tmp_8) + 1.96*log_hr_se_8), 
                                 exp(log(hr_tmp_9) + 1.96*log_hr_se_9), 
                                 exp(log(hr_tmp_10) + 1.96*log_hr_se_10),
                                 exp(log(hr_tmp_11) + 1.96*log_hr_se_11),
                                 exp(log(hr_tmp_12) + 1.96*log_hr_se_12)),
                       pm0 = c(rep(8, nrow(est_data)),
                               rep(9, nrow(est_data)),
                               rep(10, nrow(est_data)), 
                               rep(11, nrow(est_data)),
                               rep(12, nrow(est_data))),
                       pm1 = rep(est_data[,1], 5),
                       race = scenario$race,
                       dual = scenario$dual)
  
  # hazard contrasts
  contr_tmp_11 <- as.numeric(est_data[idx11,6])/as.numeric(est_data[idx12,6])
  contr_tmp_10 <- as.numeric(est_data[idx10,6])/as.numeric(est_data[idx12,6])
  contr_tmp_9 <- as.numeric(est_data[idx9,6])/as.numeric(est_data[idx12,6])
  contr_tmp_8 <- as.numeric(est_data[idx8,6])/as.numeric(est_data[idx12,6]) 
  
  # contrast standard error - bootstrap
  contr_se_8 <- sd(log(boot_mat[,idx8]/perm_mat[,idx12]))*sqrt(m)
  contr_se_9 <- sd(log(boot_mat[,idx9]/perm_mat[,idx12]))*sqrt(m)
  contr_se_10 <- sd(log(boot_mat[,idx10]/perm_mat[,idx12]))*sqrt(m)
  contr_se_11 <- sd(log(boot_mat[,idx11]/perm_mat[,idx12]))*sqrt(m)
  
  # contrast standard error - delta method
  # contr_se_8 <- sqrt(c(est_data[idx8,7]^2)/c(est_data[idx8,6]^2) +
  #                       c(est_data[idx12,7]^2)/c(est_data[idx12,6]^2))
  # contr_se_9 <- sqrt(c(est_data[idx9,7]^2)/c(est_data[idx9,6]^2) +
  #                      c(est_data[idx12,7]^2)/c(est_data[idx12,6]^2))
  # contr_se_10 <- sqrt(c(est_data[idx10,7]^2)/c(est_data[idx10,6]^2) +
  #                      c(est_data[idx12,7]^2)/c(est_data[idx12,6]^2))
  # contr_se_11 <- sqrt(c(est_data[idx11,7]^2)/c(est_data[idx11,6]^2) +
  #                      c(est_data[idx12,7]^2)/c(est_data[idx12,6]^2))
  
  contrast_tmp <- data.frame(estimate = c(contr_tmp_11, contr_tmp_10, contr_tmp_9, contr_tmp_8),
                             lower = c(exp(log(contr_tmp_11) - 1.96*contr_se_11), 
                                       exp(log(contr_tmp_10) - 1.96*contr_se_10), 
                                       exp(log(contr_tmp_9) - 1.96*contr_se_9),
                                       exp(log(contr_tmp_8) - 1.96*contr_se_8)),
                             upper = c(exp(log(contr_tmp_11) + 1.96*contr_se_11), 
                                       exp(log(contr_tmp_10) + 1.96*contr_se_10), 
                                       exp(log(contr_tmp_9) + 1.96*contr_se_9),
                                       exp(log(contr_tmp_8) + 1.96*contr_se_8)),
                             pm0 = c(11,10, 9, 8),
                             pm1 = c(12, 12, 12,12),
                             race = scenario$race,
                             dual = scenario$dual)
  
  contrast_tmp$contrast <- paste0(contrast_tmp$pm1, " vs. ", contrast_tmp$pm0)
  
  pm <- rbind(pm, pm_tmp)
  ar <- rbind(ar, ar_tmp)
  hr <- rbind(hr, hr_tmp)
  contrast <- rbind(contrast, contrast_tmp)
  
}

save(pm, file = '~/Data/pm.RData')
save(ar, file = '~/Data/ar.RData')
save(contrast, file = '~/Data/contrast.RData')
save(hr, file = '~/Data/hr.RData')

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
prop$dual_label <- ifelse(prop$dual == "both", "All", ifelse(prop$dual == "low", "Lower SEP", "Higher SEP"))

# subset important subgroup
prop$label <- paste(prop$dual_label, prop$race_label, sep = "")
prop_dat <- subset(prop, race != "all" & dual != "both")

prop_plot <- prop_dat %>%
  ggplot(aes(x = label, y = estimate, fill = factor(pm))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "solid") +
  labs(x = "",
       y = "Proportion of Participants at Risk",
       fill = ~ PM[2.5]*" Concentration Thresholds") +
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

pdf(file = "~/Figures/prop_plot.pdf", width = 10, height = 8)
prop_plot
dev.off()

### PM Decrease

# nice labels
pm$race_label <- ifelse(pm$race == "all", "", ifelse(pm$race == "white", " White", " Black"))
pm$dual_label <- ifelse(pm$dual == "both", "All", ifelse(pm$dual == "low", "Lower SEP", "Higher SEP"))
pm_dat <- subset(pm, race != "all" & dual != "both")

exposure_plot <- pm_dat %>%
  ggplot(aes(x = as.numeric(year), y = pm, color = race_label, linetype = dual_label)) +
  geom_line(size = 1) +
  labs(x = "Year",
       y = ~"Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")",
       color = "Race", 
       linetype = "SEP") +
  theme_bw() +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = c(5,6,7,8,9,10,11,12,13,14,15)) +
  scale_color_manual(values = c("#367E18", "#F57328"))

pdf(file = "~/Figures/exposure_plot.pdf", width = 10, height = 8)
exposure_plot
dev.off()

### Main Plot

## Hazard Ratio

hr_tmp <- subset(hr, dual == "both" & race == "all" & pm0 == 12)

# exposure response curve
hr_plot <- hr_tmp %>%
  ggplot(aes(x = pm1, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  geom_segment(x = 8, y = 0.7, xend = 8, yend = hr_tmp$estimate[idx8], linetype = "dashed", color = "#D81B60") +
  geom_segment(x = 9, y = 0.7, xend = 9, yend = hr_tmp$estimate[idx9], linetype = "dashed", color = "#1E88E5") +
  geom_segment(x = 10, y = 0.7, xend = 10, yend = hr_tmp$estimate[idx10], linetype = "dashed", color = "#FFC107") +
  geom_segment(x = 11, y = 0.7, xend = 11, yend = hr_tmp$estimate[idx11], linetype = "dashed", color = "#004D40") +
  geom_segment(x = 12, y = 0.7, xend = 12, yend = hr_tmp$estimate[idx12], linetype = "dotted", color = "black") +
  geom_segment(x = 3, y = hr_tmp$estimate[idx8], xend = 8, yend = hr_tmp$estimate[idx8], linetype = "dashed", color = "#D81B60") +
  geom_segment(x = 3, y = hr_tmp$estimate[idx9], xend = 9, yend = hr_tmp$estimate[idx9], linetype = "dashed", color = "#1E88E5") +
  geom_segment(x = 3, y = hr_tmp$estimate[idx10], xend = 10, yend = hr_tmp$estimate[idx10], linetype = "dashed", color = "#FFC107") +
  geom_segment(x = 3, y = hr_tmp$estimate[idx11], xend = 11, yend = hr_tmp$estimate[idx11], linetype = "dashed", color = "#004D40") +
  theme_bw() +
  coord_cartesian(xlim = c(6,12), 
                  ylim = c(0.92, 1.01)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Hazard Ratio") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = seq(0.92, 1.01, by = 0.01)) +
  scale_x_continuous(breaks = c(6,7,8,9,10,11,12))

pdf(file = "~/Figures/hr_plot.pdf", width = 10, height = 8)
hr_plot
dev.off()

## Mortality Rate

ar_tmp <- subset(ar, dual == "both" & race == "all")

ar_plot <- ar_tmp %>%
  ggplot(aes(x = a.vals, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_line(size = 1) +
  geom_segment(x = 8, y = -0.1, xend = 8, yend = ar_tmp$estimate[idx8], linetype = "dashed", color = "#D81B60") +
  geom_segment(x = 9, y = -0.1, xend = 9, yend = ar_tmp$estimate[idx9], linetype = "dashed", color = "#1E88E5") +
  geom_segment(x = 10, y = -0.1, xend = 10, yend = ar_tmp$estimate[idx10], linetype = "dashed", color = "#FFC107") +
  geom_segment(x = 11, y = -0.1, xend = 11, yend = ar_tmp$estimate[idx11], linetype = "dashed", color = "#004D40") +
  geom_segment(x = 3, y = ar_tmp$estimate[idx8], xend = 8, yend = ar_tmp$estimate[idx8], linetype = "dashed", color = "#D81B60") +
  geom_segment(x = 3, y = ar_tmp$estimate[idx9], xend = 9, yend = ar_tmp$estimate[idx9], linetype = "dashed", color = "#1E88E5") +
  geom_segment(x = 3, y = ar_tmp$estimate[idx10], xend = 10, yend = ar_tmp$estimate[idx10], linetype = "dashed", color = "#FFC107") +
  geom_segment(x = 3, y = ar_tmp$estimate[idx11], xend = 11, yend = ar_tmp$estimate[idx11], linetype = "dashed", color = "#004D40") +
  theme_bw() +
  coord_cartesian(xlim = c(5,15), 
                  ylim = c(0.044,0.051)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")",
       y = "Absolute Mortality Rate") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = seq(0.044, 0.051, by = 0.001)) +
  scale_x_continuous(breaks = seq(5,15,by = 1))

pdf(file = "~/Figures/ar_plot.pdf", width = 10, height = 8)
ar_plot
dev.off()

### ERC by Race

## Hazard Ratio

hr$dual_label <- ifelse(hr$dual == "high", "Higher SEP", ifelse(hr$dual == "low", "Lower SEP", "All"))
hr$dual_label <- factor(hr$dual_label, levels = c("All", "Higher SEP", "Lower SEP"))
hr$race_label <- str_to_title(hr$race)

hr_tmp <- subset(hr, race %in% c("black", "white") & pm0 == 12)

hr_strata <- hr_tmp %>%
  ggplot(aes(x = pm1, y = estimate, color = factor(race_label))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  geom_segment(x = 12, y = 0.7, xend = 12, yend = hr_tmp$estimate[idx12], linetype = "dotted", color = "black") +
  facet_wrap(~ dual_label) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Hazard Ratio",
       color = "Race") +
  theme_bw() +
  coord_cartesian(xlim = c(6,12), 
                  ylim = c(0.86, 1.02)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#367E18", "#F57328")) +
  scale_y_continuous(breaks = seq(0.86, 1.02, by = 0.02)) +
  scale_x_continuous(breaks = c(6,7,8,9,10,11,12))

pdf(file = "~/Figures/hr_strata.pdf", width = 16, height = 8)
hr_strata
dev.off()

## Mortality Rate

ar$dual_label <- ifelse(ar$dual == "high", "Higher SEP", ifelse(ar$dual == "low", "Lower SEP", "All"))
ar$dual_label <- factor(ar$dual_label, levels = c("All", "Higher SEP", "Lower SEP"))
ar$race_label <- str_to_title(ar$race)

ar_tmp <- subset(ar, race %in% c("black", "white"))

ar_strata <- subset(ar_tmp, a.vals < 16 & a.vals > 4) %>%
  ggplot(aes(x = a.vals, y = estimate, color = factor(race_label))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  facet_wrap(~ dual_label) +
  coord_cartesian(xlim = c(5,15), 
                  ylim = c(0.02, 0.11)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")",
       y = "Absolute Mortality Rate", 
       color = "Race") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#367E18", "#F57328")) +
  scale_x_continuous(breaks = c(5,6,7,8,9,10,11,12,13,14,15))

pdf(file = "~/Figures/ar_strata.pdf", width = 16, height = 8)
ar_strata
dev.off()

### Contrast Plot

contr <- subset(contrast, race %in% c("black", "white"))
contr$contrast <- factor(contr$contrast, levels = c("12 vs. 11", "12 vs. 10", "12 vs. 9", "12 vs. 8"))
contr$dual_label <- ifelse(contr$dual == "high", "Higher SEP",
                           ifelse(contr$dual == "low", "Lower SEP", "All"))
contr$race_label <- str_to_title(contr$race)

contrast_plot <- contr %>%
  ggplot(aes(x = contrast, y = 1 - estimate, color = race_label)) +
  geom_pointrange(aes(ymin = 1 - lower, ymax = 1 - upper), 
                  position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0) +
  facet_wrap(~ dual_label) +
  theme_bw() +
  labs(x = ~ PM[2.5]*" Contrasts", y = "1 - (Hazard Ratio)",
       color = "Race", linetype = "Model") +
  theme_bw() +
  coord_cartesian(ylim = c(-0.01, 0.09)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#367E18", "#F57328")) +
  scale_y_continuous(breaks = seq(-0.01, 0.09, by = 0.01)) +
  scale_x_discrete(labels = c(~"12 vs. 11 "*mu*g*"/"*m^3,
                              ~"12 vs. 10 "*mu*g*"/"*m^3,
                              ~"12 vs. 9 "*mu*g*"/"*m^3,
                              ~"12 vs. 8 "*mu*g*"/"*m^3))

pdf(file = "~/Figures/contrast_plot.pdf", width = 16, height = 8)
contrast_plot
dev.off()
