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

# data directories
dir_data <- '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_mod <- '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_Data/'

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/R/calibrate.R')
set.seed(42)

## Bandwidth Selection

scenario <- scenarios[1,]
load(paste0(dir_mod, scenario$dual, "_", scenario$race, ".RData"))

# fit exposure response curves
target_1 <- count_erf_lm(resid.lm = model_data$resid.lm, 
                         muhat.mat = model_data$muhat.mat, log.pop = model_data$log.pop,
                         w.id = individual_data$id, a = zip_data$pm25, x.id = zip_data$id, 
                         bw = 1, a.vals = a.vals, phat.vals = phat.vals, se.fit = FALSE)

target_xiao <- count_erf_lm(resid.lm = model_data$resid.lm,
                            muhat.mat = model_data$muhat.mat, log.pop = model_data$log.pop,
                            w.id = individual_data$id, a = zip_data$pm25, x.id = zip_data$id, 
                            bw = 1.56147, a.vals = a.vals, phat.vals = phat.vals, se.fit = FALSE)

target_loo <- count_erf_lm(resid.lm = model_data$resid.lm, 
                           muhat.mat = model_data$muhat.mat, log.pop = model_data$log.pop,
                           w.id = individual_data$id, a = zip_data$pm25, x.id = zip_data$id, 
                           bw = 2.8, a.vals = a.vals, phat.vals = phat.vals, se.fit = FALSE)

target_kfold <- count_erf_lm(resid.lm = model_data$resid.lm, 
                             muhat.mat = model_data$muhat.mat, log.pop = model_data$log.pop,
                             w.id = individual_data$id, a = zip_data$pm25, x.id = zip_data$id, 
                             bw = 3.2, a.vals = a.vals, phat.vals = phat.vals, se.fit = FALSE)

est_data <- data.frame(a.vals = a.vals,
                       estimate.loo = target_loo$estimate.lm,
                       estimate.1 = target_1$estimate.lm,
                       estimate.xiao = target_xiao$estimate.lm,
                       estimate.kfold = target_kfold$estimate.lm)

save(est_data,file = paste0(dir_out,"bw_test.RData"))

rm(individual_data, zip_data, model_data, target); gc()

# contrast indexes
idx8 <- which.min(abs(a.vals - 8))
idx9 <- which.min(abs(a.vals - 9))
idx10 <- which.min(abs(a.vals - 10))
idx11 <- which.min(abs(a.vals - 11))
idx12 <- which.min(abs(a.vals - 12))

load(paste0(dir_out, "bw_test.RData"))
idx0 <- which.min(abs(est_data$a.vals - 12))

cv_bw <- rbind(data.frame(a.vals = est_data$a.vals, ERC = est_data$estimate.1, 
                          hr = c(as.numeric(est_data$estimate.1)/as.numeric(est_data$estimate.1[idx0])), bw = "1"),
               data.frame(a.vals = est_data$a.vals, ERC = est_data$estimate.xiao,
                          hr = c(as.numeric(est_data$estimate.xiao)/as.numeric(est_data$estimate.xiao[idx0])), bw = "Xiao"),
               data.frame(a.vals = est_data$a.vals, ERC = est_data$estimate.loo,
                          hr = c(as.numeric(est_data$estimate.loo)/as.numeric(est_data$estimate.loo[idx0])), bw = "Leave-One-Out"),
               data.frame(a.vals = est_data$a.vals, ERC = est_data$estimate.kfold,
                          hr = c(as.numeric(est_data$estimate.kfold)/as.numeric(est_data$estimate.kfold[idx0])), bw = "k-Fold"))


# Hazard Ratio
# exposure response curve
hr_bw <- ggplot(data = cv_bw, aes(x=a.vals, y = hr, color = bw)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = c(0.7, 0.2),
        legend.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Hazard Ratio",
       color = "Bandwidth",
       title = "Badwidth Choices") +
  coord_cartesian(xlim = c(5,13), 
                  ylim = c(0.91, 1.02)) +
  scale_y_continuous(breaks = c(0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02)) +
  scale_x_continuous(breaks = c(5,6,7,8,9,10,11,12,13))

pdf(file = "~/Figures/bw_compare.pdf", width = 10, height = 8)
hr_bw
dev.off()

### Covariate Balance Plot

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

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
  if (scenario$dual == "high") {
    dual0 <- 0
  } else if (scenario$dual == "low") {
    dual0 <- 1
  } else if (scenario$dual == "both") {
    dual0 <- c(0,1)
  }
  
  if (scenario$race == "white") {
    race0 <- 1
  } else if (scenario$race == "black") {
    race0 <- 2
  } else if (scenario$race == "asian") {
    race0 <- 4
  } else if (scenario$race == "hispanic") {
    race0 <- 5
  } else if (scenario$race == "other") {
    race0 <- 3
  } else if (scenario$race == "all") {
    race0 <- c(1,2,3,4,5)
  }
  
  i.zip <- paste(new_data$w$zip, new_data$w$year, sep = "-")
  u.zip <- unique(i.zip[new_data$w$race %in% race0 & new_data$w$dual %in% dual0])
  new_data$x <- subset(new_data$x, id %in% u.zip)
  
  bdat_1 <- bal_dat(a = new_data$x$pm25, x = new_data$x[,c(2,4:20)], weights = rep(1, nrow(new_data$x)))
  bdat_2 <- bal_dat(a = new_data$x$pm25, x = new_data$x[,c(2,4:20)], weights = new_data$x$ipw)
  bdat_3 <- bal_dat(a = new_data$x$pm25, x = new_data$x[,c(2,4:20)], weights = new_data$x$cal)
  bdat_4 <- bal_dat(a = new_data$x$pm25, x = new_data$x[,c(2,4:20)], weights = new_data$x$cal_trunc)
  
  bdat_2$adjust <- "LM"
  bdat_3$adjust <- "Calibration"
  bdat_4$adjust <- "Truncated Calibration"
  
  df <- rbind(bdat_1, bdat_2, bdat_3, bdat_4)
  df$adjust <- factor(df$adjust, levels = c("Unadjusted", "LM", "Calibration", "Truncated Calibration"))
  
  title <- paste(str_to_upper(scenario$dual), str_to_upper(scenario$race), "SUBPOPULATION GPS")
  
  balance_plot <- ggplot(data = df, aes(x = labs, y = vals, color = adjust)) +
    geom_point(pch = 21, size = 2) +
    geom_line(aes(group = adjust)) +
    geom_hline(yintercept = 0, lty = 1) +
    geom_hline(yintercept = 0.1, lty = 3, colour = "black") +
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Covariates") + ylab("Absolute Correlation") +
    ggtitle(title) +
    ylim(0, 0.35) +
    guides(color = guide_legend(title = "Implementation")) +
    theme_bw() +
    theme(axis.text.y = element_text(angle = 30, hjust = 1),
          legend.position = c(0.85, 0.15),
          legend.background = element_rect(colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(i)
  
  pdf(file = paste0("~/Figures/balance_plot_", scenario$dual, "_", scenario$race,"_sub.pdf"), width = 10, height = 8)
  balance_plot
  dev.off()
  
}
