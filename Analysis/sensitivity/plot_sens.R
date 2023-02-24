library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"), race = c("all", "white", "black"))
# scenarios <- expand.grid(dual = c("high", "low"), race = c("asian", "hispanic", "other"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/calibrate.R')
set.seed(42)

### Plot G-Computation and GPS as Regressor Results

# data directories
dir_data <- '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_mod <- c('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR_All/',
             '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/GPSReg_All/',
             '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/GComp_All/')

mod <- c("Doubly-Robust", "GPS as a Regressor", "G-Computation")

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
  
  for (j in 1:3) {
    
    load(paste0(dir_mod[j], scenario$dual, "_", scenario$race, ".RData"))
    load(paste0(dir_mod[j], scenario$dual, "_", scenario$race, "_boot.RData"))
    
    if (j == 1)
      mhat.vals <- est_data[,6]
    
    # Replace with Bootstrap SEs
    u.zip <- unique(individual_data$zip)
    m <- 1/log(length(u.zip)) # for m out of n bootstrap
    mhat.se <- apply(boot_mat, 2, sd, na.rm = T)*sqrt(m)
    perm_mat <- boot_mat[sample(1:nrow(boot_mat), size = nrow(boot_mat), replace = FALSE), ]
    
    # absolute risks
    ar_tmp <- data.frame(a.vals = c(a.vals), 
                         estimate = c(mhat.vals),
                         lower = c(mhat.vals - 1.96*mhat.se),
                         upper = c(mhat.vals + 1.96*mhat.se),
                         race = rep(scenario$race, length(mhat.vals)),
                         dual = rep(scenario$dual, length(mhat.vals)),
                         mod = rep(mod[j], length(mhat.vals)))
    
    # hazard ratio
    hr_tmp_8 <- c(as.numeric(mhat.vals)/as.numeric(mhat.vals[idx8]))
    hr_tmp_9 <- c(as.numeric(mhat.vals)/as.numeric(mhat.vals[idx9]))
    hr_tmp_10 <- c(as.numeric(mhat.vals)/as.numeric(mhat.vals[idx10]))
    hr_tmp_11 <- c(as.numeric(mhat.vals)/as.numeric(mhat.vals[idx11]))
    hr_tmp_12 <- c(as.numeric(mhat.vals)/as.numeric(mhat.vals[idx12]))
    
    # bootstrap standard errors
    log_hr_se_8 <- apply(log(boot_mat/perm_mat[,idx8]), 2, sd, na.rm = TRUE)*sqrt(m)
    log_hr_se_9 <- apply(log(boot_mat/perm_mat[,idx9]), 2, sd, na.rm = TRUE)*sqrt(m)
    log_hr_se_10 <- apply(log(boot_mat/perm_mat[,idx10]), 2, sd, na.rm = TRUE)*sqrt(m)
    log_hr_se_11 <- apply(log(boot_mat/perm_mat[,idx11]), 2, sd, na.rm = TRUE)*sqrt(m)
    log_hr_se_12 <- apply(log(boot_mat/perm_mat[,idx12]), 2, sd, na.rm = TRUE)*sqrt(m)
    
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
                         pm0 = c(rep(8, length(mhat.vals)),
                                 rep(9, length(mhat.vals)),
                                 rep(10, length(mhat.vals)), 
                                 rep(11, length(mhat.vals)),
                                 rep(12, length(mhat.vals))),
                         pm1 = rep(a.vals, 5),
                         race = scenario$race,
                         dual = scenario$dual,
                         mod = mod[j])
    
    # hazard contrasts
    contr_tmp_11 <- as.numeric(mhat.vals[idx11])/as.numeric(mhat.vals[idx12])
    contr_tmp_10 <- as.numeric(mhat.vals[idx10])/as.numeric(mhat.vals[idx12])
    contr_tmp_9 <- as.numeric(mhat.vals[idx9])/as.numeric(mhat.vals[idx12])
    contr_tmp_8 <- as.numeric(mhat.vals[idx8])/as.numeric(mhat.vals[idx12]) 
    
    # contrast standard error - bootstrap
    contr_se_8 <- sd(log(boot_mat[,idx8]/perm_mat[,idx12]))*sqrt(m)
    contr_se_9 <- sd(log(boot_mat[,idx9]/perm_mat[,idx12]))*sqrt(m)
    contr_se_10 <- sd(log(boot_mat[,idx10]/perm_mat[,idx12]))*sqrt(m)
    contr_se_11 <- sd(log(boot_mat[,idx11]/perm_mat[,idx12]))*sqrt(m)
    
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
                               dual = scenario$dual,
                               mod = mod[j])
    
    contrast_tmp$contrast <- paste0(contrast_tmp$pm1, " vs. ", contrast_tmp$pm0)
    
    ar <- rbind(ar, ar_tmp)
    hr <- rbind(hr, hr_tmp)
    contrast <- rbind(contrast, contrast_tmp)
    
  }
  
}

save(ar, file = '~/Data/ar_sensitivity.RData')
save(hr, file = '~/Data/hr_sensitivity.RData')
save(contrast, file = '~/Data/contrast_sensitivity.RData')

load(file = '~/Data/ar_sensitivity.RData')
load(file = '~/Data/hr_sensitivity.RData')
load(file = '~/Data/contrast_sensitivity.RData')

## Overall ERC

# get matching estimators
DR100 <- read.csv("~/Data/DR100_table.csv")[,-1]
DR100 <- subset(DR100, !(Methods %in% c("Doubly Robust")))

ar_tmp <- subset(ar, dual == "both" & race == "all")

Kevin <- data.frame(a.vals = ar_tmp$a.vals, 
                    estimate = ar_tmp$estimate, 
                    gps_method = ar_tmp$mod, 
                    lower = ar_tmp$lower,
                    upper = ar_tmp$upper)

colnames(Kevin) <- colnames(DR100)
DR100_table <- rbind(DR100, Kevin)
DR100_table$SE <- (DR100_table$ERC - DR100_table$lower)/1.96

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

# plot_dat$hr.lower[plot_dat$Methods != c("Doubly-Robust")] <- NA
# plot_dat$hr.upper[plot_dat$Methods != c("Doubly-Robust")] <- NA

plot_dat$Methods <- factor(plot_dat$Methods, levels = c("Doubly-Robust", "Matching", 
                                                        "G-Computation", "GPS as a Regressor"),
                           labels = c("Doubly-Robust (Main Analysis)", "Matching (Wu et al.)", 
                                      "G-Computation (Sensitivity)", "GPS as a Regressor (Sensitivity)"))

# Hazard Ratio
# exposure response curve
hr_compare <- plot_dat %>%
  ggplot(aes(x = a.vals, y = hr, color = Methods, linetype = Methods)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = hr.lower, ymax = hr.upper), alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw() +
  coord_cartesian(xlim = c(6,12), 
                  ylim = c(0.88, 1.04)) +
  theme(legend.position = c(0.7, 0.2),
        legend.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Hazard Ratio") +
  scale_color_manual(values = c("blue", "blue", "blue", "red"),
                     breaks = c("Doubly-Robust (Main Analysis)", "G-Computation (Sensitivity)", 
                                "GPS as a Regressor (Sensitivity)", "Matching (Wu et al.)")) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid"),
                        breaks = c("Doubly-Robust (Main Analysis)", "G-Computation (Sensitivity)", 
                                   "GPS as a Regressor (Sensitivity)", "Matching (Wu et al.)")) +
  scale_y_continuous(breaks = seq(0.88, 1.04, by = 0.02)) +
  scale_x_continuous(breaks = c(6,7,8,9,10,11,12))


pdf(file = "~/Figures/hr_plot_sensitivity.pdf", width = 10, height = 8)
hr_compare
dev.off()

## Stratified by Race x SEP

hr$dual_label <- ifelse(hr$dual == "high", "Panel B: Higher SEP",
                        ifelse(hr$dual == "low", "Panel C: Lower SEP", "Panel A: All"))
hr$race_label <- str_to_title(hr$race)

hr_tmp1 <- hr_tmp2 <- subset(hr, race %in% c("black", "white") & pm0 == 12)
hr_tmp1$CI <- "With Confidence Intervals"
hr_tmp2$CI <- "Without Confidence Intervals"
hr_tmp2$lower[hr_tmp2$mod != "Doubly-Robust"] <- NA
hr_tmp2$upper[hr_tmp2$mod != "Doubly-Robust"] <- NA
hr_tmp <- rbind(hr_tmp1, hr_tmp2)

hr_strata <- hr_tmp %>%
  ggplot(aes(x = pm1, y = estimate, color = factor(race_label), linetype = mod)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_segment(x = 12, y = 0.7, xend = 12, yend = hr_tmp$estimate[idx12], linetype = "dotted") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  facet_grid(CI ~ dual_label) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Hazard Ratio of Mortality",
       color = "Racial Identity",
       linetype = "Model") +
  theme_bw() +
  coord_cartesian(xlim = c(6,12), 
                  ylim = c(0.85, 1.01)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#367E18", "#F57328")) +
  scale_y_continuous(breaks = seq(0.85, 1.01, by = 0.02)) +
  scale_x_continuous(breaks = c(6,7,8,9,10,11,12))

pdf(file = "~/Figures/hr_strata_sensitivity.pdf", width = 16, height = 12)
hr_strata
dev.off()

### Contrast Plot

contr <- subset(contrast, race %in% c("black", "white"))
contr$contrast <- factor(contr$contrast, levels = c("12 vs. 11", "12 vs. 10", "12 vs. 9", "12 vs. 8"))
contr$dual_label <- ifelse(contr$dual == "high", "Panel B: Higher SEP",
                           ifelse(contr$dual == "low", "Panel C: Lower SEP", "Panel A: All"))
contr$race_label <- str_to_title(contr$race)

contrast_plot <- contr %>%
  ggplot(aes(x = contrast, y = 1 - estimate, 
             color = race_label, linetype = mod, shape = mod)) +
  geom_pointrange(aes(ymin = 1 - lower, ymax = 1 - upper), 
                  position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0) +
  facet_wrap(~ dual_label) +
  theme_bw() +
  labs(x = ~ PM[2.5]*" Contrasts", y = "1 - (Hazard Ratio of Mortality)",
       color = "Racial Identity", linetype = "Model", shape = "Model") +
  theme_bw() +
  coord_cartesian(ylim = c(-0.01, 0.1)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#367E18", "#F57328")) +
  scale_y_continuous(breaks = seq(-0.01, 0.1, by = 0.01)) +
  scale_x_discrete(labels = c(~"12 vs. 11 "*mu*g*"/"*m^3,
                              ~"12 vs. 10 "*mu*g*"/"*m^3,
                              ~"12 vs. 9 "*mu*g*"/"*m^3,
                              ~"12 vs. 8 "*mu*g*"/"*m^3))

pdf(file = "~/Figures/contrast_plot_sensitivity.pdf", width = 16, height = 8)
contrast_plot
dev.off()
