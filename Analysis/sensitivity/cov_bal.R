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

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/calibrate.R')
set.seed(42)

### Covariate Balance Plots

# data directories
dir_data <- '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_mod <- '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_Data/'

bal_dat <- function(a, x, weights){
  
  val <- bal.tab(x, treat = a, weights = weights, method = "weighting", continuous = "std", s.d.denom = "all")
  bal_df <- val$Balance
  labs <- rep(rownames(bal_df), 2)
  vals_tmp <- bal_df$Corr.Adj
  vals_year <- mean(abs(vals_tmp[1:17]))
  vals_region <- mean(abs(vals_tmp[32:34]))
  vals_tmp2 <- c(abs(vals_tmp[-c(1:17, 32:34)]),
                 vals_year, vals_region)
  names(vals_tmp2) <- c("Mean BMI", "Smoking Rate", "% Hispanic", "% Black",
                        "Median Household Income", "Median House Value", "% Below Poverty Level",
                        "% Below High School Education", "Population Density", "% Owner-Occupied Housing",
                        "Summer Temperature","Winter Temperature", "Summer Humidity", "Winter Humidity",
                        "Calendar Year","Census Region")
  vals <- vals_tmp2[order(vals_tmp2, decreasing = TRUE)]
  adjust <- rep("Unadjusted", each = length(vals_tmp2))
  labs <- names(vals)
  df <- data.frame(labs = labs, vals = vals, adjust = adjust)
  df$labs <- factor(df$labs, levels = rev(names(vals)))
  
  return(df)
  
}

# Standalone Plots + ESS

ess <- NULL

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
  bdat_1 <- bal_dat(a = new_data$x$pm25, x = new_data$x[,c(2,4:20)], weights = rep(1, nrow(new_data$x)))
  bdat_2 <- bal_dat(a = new_data$x$pm25, x = new_data$x[,c(2,4:20)], weights = new_data$x$ipw)
  bdat_3 <- bal_dat(a = new_data$x$pm25, x = new_data$x[,c(2,4:20)], weights = new_data$x$cal)
  bdat_4 <- bal_dat(a = new_data$x$pm25, x = new_data$x[,c(2,4:20)], weights = new_data$x$cal_trunc)
  
  bdat_2$adjust <- "LM"
  bdat_3$adjust <- "Calibration"
  bdat_4$adjust <- "Truncated Calibration"
  
  df <- rbind(bdat_1, bdat_4)
  df$adjust <- factor(df$adjust, levels = c("Unadjusted", "Truncated Calibration"))
  
  # Effective Sample Size
  ess_2 <- sum(new_data$x$ipw)^2/sum(new_data$x$ipw^2)
  ess_3 <- sum(new_data$x$cal)^2/sum(new_data$x$cal^2)
  ess_4 <- sum(new_data$x$cal_trunc)^2/sum(new_data$x$cal_trunc^2)
  ess <- rbind(ess, data.frame(ipw = ess_2, cal = ess_3, cal_trunc = ess_4, race = scenario$race, dual = scenario$dual))
  
  ifelse()
  
  balance_plot <- ggplot(data = df, aes(x = labs, y = vals, color = adjust)) +
    geom_point(pch = 21, size = 2) +
    geom_line(aes(group = adjust)) +
    geom_hline(yintercept = 0, lty = 1) +
    geom_hline(yintercept = 0.1, lty = 3, colour = "black") +
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Covariates") + ylab("Absolute Correlation") +
    ylim(0, 0.35) +
    guides(color = guide_legend(title = "Implementation")) +
    theme_bw() +
    theme(axis.text.y = element_text(angle = 30, hjust = 1),
          legend.position = c(0.85, 0.15),
          legend.background = element_rect(colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(i)
  
  pdf(file = paste0("~/Figures/balance_plot_", scenario$dual, "_", scenario$race,".pdf"), width = 10, height = 8)
  balance_plot
  dev.off()
  
}

save(ess, file = "~/Data/ess.RData")

# Panel Plot

df <- NULL

for (i in c(4:9)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
  race <- paste(str_to_title(scenario$race), "Participants")
  
  if (i %in% c(4,7))
    sep <- "All"
  else if (i %in% c(5,8))
    sep <- "Higher SEP"
  else
    sep <- "Lower SEP"
  
  bdat_1 <- bal_dat(a = new_data$x$pm25, x = new_data$x[,c(2,4:20)], weights = rep(1, nrow(new_data$x)))
  bdat_2 <- bal_dat(a = new_data$x$pm25, x = new_data$x[,c(2,4:20)], weights = new_data$x$ipw)
  bdat_3 <- bal_dat(a = new_data$x$pm25, x = new_data$x[,c(2,4:20)], weights = new_data$x$cal)
  bdat_4 <- bal_dat(a = new_data$x$pm25, x = new_data$x[,c(2,4:20)], weights = new_data$x$cal_trunc)
  
  bdat_2$adjust <- "LM"
  bdat_3$adjust <- "Calibration"
  bdat_4$adjust <- "Truncated Calibration"
  
  df.tmp <- rbind(bdat_1, bdat_4)
  df <- rbind(df, data.frame(df.tmp, race = race, sep = sep))
  
}

df$adjust <- factor(df$adjust, levels = c("Unadjusted", "Truncated Calibration"))
df$sep <- factor(df$sep, levels = c("All", "Higher SEP", "Lower SEP"))
df$race <- factor(df$race, levels = c("Black Participants", "White Participants"))

balance_plot <- ggplot(data = df, aes(x = labs, y = vals, color = adjust)) +
  geom_point(pch = 21, size = 2) +
  facet_grid(race ~ sep) +
  geom_line(aes(group = adjust)) +
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(yintercept = 0.1, lty = 3, colour = "black") +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Covariates") + ylab("Absolute Correlation") +
  ylim(0, 0.35) +
  guides(color = guide_legend(title = "Implementation")) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

pdf(file = "~/Figures/panel_balance_plot.pdf", width = 16, height = 12)
balance_plot
dev.off()
