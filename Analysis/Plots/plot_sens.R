library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(splines)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(cobalt)

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("all","white", "black"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(4, 16, length.out = 121)

# data directories
dir_data_qd <- '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_data_rm <- '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm/'
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/Other/DR_KWLS/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/Other/DR_rm/'

### Exposure Assessment Sensitivity Plot

scenario <- scenarios[3,]

# QD
load(paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
dat_qd <- data.frame(a.vals = c(est_data$a.vals), 
                     estimate = c(est_data[,6]),
                     lower = c(est_data[,6] - 1.96*est_data[,7]),
                     upper = c(est_data[,6] + 1.96*est_data[,7]),
                     exposure = rep("Di et al. (2019)", nrow(est_data)),
                     race = rep(scenario$race, nrow(est_data)),
                     dual = rep(scenario$dual, nrow(est_data)))

a_dat_sens <- data.frame(a = rep(individual_data$pm25, individual_data$time_count), 
                         exposure = "Di et al. (2019)")

# RM
load(paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
dat_rm <- data.frame(a.vals = c(est_data$a.vals),
                     estimate = c(est_data[,6]),
                     lower = c(est_data[,6] - 1.96*est_data[,7]),
                     upper = c(est_data[,6] + 1.96*est_data[,7]),
                     exposure = rep("van Donkelaar et al. (2016)", nrow(est_data)),
                     race = rep(scenario$race, nrow(est_data)),
                     dual = rep(scenario$dual, nrow(est_data)))

a_dat_sens <- rbind(a_dat_sens, data.frame(a = rep(individual_data$pm25, individual_data$time_count),
                                           exposure = "van Donkelaar et al. (2016)"))

# combine
dat_qd_rm <- rbind(dat_qd, dat_rm)

erf_plot_sens <- dat_qd_rm %>% 
  ggplot(aes(x = a.vals, y = estimate, color = exposure)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(5,15), ylim = c(0.044,0.049)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate",
       color = "Exposure Assessment", title = "Exposure Response Curve under\n Different Exposure Assessments") + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = c(0.044,0.045,0.046,0.047,0.048,0.049)) +
  grids(linetype = "dashed")

leg_plot <- dat_qd_rm %>% 
  ggplot(aes(x = a.vals, y = estimate, color = exposure)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(5,15), ylim = c(0.044,0.049)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate",
       color = "Exposure Assessment") +
  theme_bw() +
  theme(legend.background = element_rect(colour = "black"))

leg <- gtable_filter(ggplot_gtable(ggplot_build(leg_plot)), "guide-box")

a_hist_sens <- ggplot(a_dat_sens, mapping = aes(x = a, fill = exposure)) + 
  geom_density(alpha = 0.3, adjust = 3)+
  coord_cartesian(xlim = c(5,15), ylim = c(0,0.15)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Exposure Density") + 
  theme(panel.grid = element_blank()) +
  scale_y_continuous(position = "right", breaks = c(0, 0.05, 0.10, 0.15)) +
  guides(fill="none") +
  theme_cowplot()

align <- align_plots(a_hist_sens, erf_plot_sens +
                       annotation_custom(leg, xmin = 6, xmax = 7.5, ymin = 0.048, ymax = 0.049), 
                     align = "hv", axis = "tblr")
sens_plot <- ggdraw(align[[1]]) + draw_plot(align[[2]])

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/sens_plot.pdf", width = 8, height = 8)
sens_plot
dev.off()

### Covariate Balance Plot

bal_dat <- function(a, x, weights){

  val <- bal.tab(x, treat = a, weights = weights, method = "weighting", continuous = "raw", s.d.denom = "pooled")
  bal_df <- val$Balance
  labs <- rep(rownames(bal_df), 2)
  vals_tmp <- cbind(bal_df$Corr.Un, bal_df$Corr.Adj)
  vals_year <- c(mean(abs(vals_tmp[1:17,1])), mean(abs(vals_tmp[1:17,2])))
  vals_region <- c(mean(abs(vals_tmp[32:34,1])), mean(abs(vals_tmp[32:34,2])))
  vals_tmp2 <- rbind(cbind(abs(vals_tmp[-c(1:17, 32:34),1]),
                           abs(vals_tmp[-c(1:17, 32:34),2])),
                     vals_year, vals_region)
  rownames(vals_tmp2) <- c("Mean BMI", "Smoking Rate", "% Hispanic", "% Black",
                           "Median Household Income", "Median House Value", "% Below Poverty Level",
                           "% Below High School Education", "Population Density", "% Owner-Occupied Housing",
                           "Summer Temperature","Winter Temperature", "Summer Humidity", "Winter Humidity",
                           "Calendar Year","Census Region")
  vals_tmp2 <- vals_tmp2[order(vals_tmp2[,1], decreasing = TRUE),]

  vals <- c(vals_tmp2[,1], vals_tmp2[,2])
  adjust <- rep(c("Unadjusted", "LM"), each = nrow(vals_tmp2))
  labs <- rep(rownames(vals_tmp2), times = 2)
  df <- data.frame(labs = labs, vals = vals, adjust = adjust)
  df$labs <- factor(df$labs, levels = rev(rownames(vals_tmp2)))

  return(df)

}

load("/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/Other/DR_SuperLearner/2_all_qd.RData")
bdat_1 <- bal_dat(a = zip_data$pm25, x = zip_data[,c(2,4:20)], weights = zip_data$weights.lm)
bdat_2 <- bal_dat(a = zip_data$pm25, x = zip_data[,c(2,4:20)], weights = zip_data$weights.sl)
bdat_3 <- bal_dat(a = zip_data$pm25, x = zip_data[,c(2,4:20)], weights = zip_data$weights.cal)
bdat_tmp_1 <- subset(bdat_2, adjust == "LM")
bdat_tmp_2 <- subset(bdat_3, adjust == "LM")
bdat_tmp_1$adjust <- "XGBoost"
bdat_tmp_2$adjust <- "Calibration"

df <- rbind(bdat_1, bdat_tmp_1, bdat_tmp_2)
df$adjust <- factor(df$adjust, levels = c("Unadjusted", "LM", "XGBoost", "Calibration"))

balance_plot <- ggplot(data = df, aes(x = labs, y = vals, color = adjust)) +
  geom_point(pch = 21, size = 2) +
  geom_line(aes(group = adjust)) +
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(yintercept = 0.1, lty = 3, colour = "black") +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Covariates") + ylab("Absolute Correlation") +
  ggtitle("Covariate Balance by different GPS Implementations") +
  ylim(0, 0.35) +
  guides(color = guide_legend(title = "Implementation")) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 30, hjust = 1),
        legend.position = c(0.85, 0.15),
        legend.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values=c("#F77452","#82F739", "#6867AA", "#575757")) +
  grids(linetype = "dashed")

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/balance_plot.pdf", width = 8, height = 8)
balance_plot
dev.off()
