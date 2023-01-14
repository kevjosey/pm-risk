library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(splines)
library(ggplot2)
library(ggpubr)
library(cowplot)

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("all","white", "black"),
                         age = c("65-75", "75-85", "85-95"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios$age <- as.character(scenarios$age)
a.vals <- seq(4, 16, length.out = 121)
n.boot <- 1000

dat <- data.frame()

# race by dual by age plots
for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_', scenario$age, "/",
              scenario$dual, "_", scenario$race, "_qd.RData"))
  
  dat_tmp <- data.frame(a.vals = c(est_data$a.vals), 
                        estimate = c(est_data$estimate.cal),
                        linear = c(est_data$linear.cal),
                        lower = c(est_data[,4] - 1.96*est_data[,5]),
                        upper = c(est_data[,4] + 1.96*est_data[,5]),
                        exposure = rep("Di et al. (2019)", nrow(est_data)),
                        race = rep(scenario$race, nrow(est_data)),
                        sample_size = rep(paste0(str_to_title(scenario$race), " = ", 
                                                 formatC(sum(individual_data$time_count),
                                                         format = "d", big.mark = ",")), 
                                          nrow(est_data)),
                        dual = rep(scenario$dual, nrow(est_data)),
                        age = rep(scenario$age, nrow(est_data)))
  
  dat <- rbind(dat, dat_tmp)
  
}

plot_list <- list()
situations <- expand.grid(dual = c(2, 0, 1), age = c("65-75", "75-85", "85-95"))
situations$dual <- as.numeric(situations$dual)
situations$age <- as.character(situations$age)

for (i in 1:nrow(situations)){
  
  situation <- situations[i,]
  
  if (situation$dual == 1)
    main <- "Low SEP"
  else if (situation$dual == 0)
    main <- "High SEP"
  else
    main <- "High + Low SEP"
  
  # factor race
  dat_tmp <- subset(dat, dual == as.numeric(situation$dual) & race != "all" & age == situation$age)
  dat_tmp$race <- str_to_title(dat_tmp$race)
  dat_tmp$sample_size <- factor(as.character(dat_tmp$sample_size))
  
  # graph breaks
  ylim <- c(min(dat_tmp$lower[dat_tmp$a.vals >= 5 & dat_tmp$a.vals <= 15]),
            max(dat_tmp$upper[dat_tmp$a.vals >= 5 & dat_tmp$a.vals <= 15]))
  breaks <- round(seq(ylim[1], ylim[2], length.out = 6), 4)
    
  # dual ineligible + eligible
  erf_strata_plot <- dat_tmp %>% 
    ggplot(aes(x = a.vals, y = estimate, color = sample_size)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
    geom_line(size = 1) +
    coord_cartesian(xlim = c(5,15), ylim = ylim) +
    labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate", 
         color = "Race", title = main) + 
    theme_bw() +
    guides(color = guide_legend(title = "Race")) + 
    scale_color_manual(values = c("#367E18", "#F57328")) +
    scale_y_continuous(breaks = breaks) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    grids(linetype = "dashed")
  
  leg_plot_tmp <- dat_tmp %>% 
    ggplot(aes(x = a.vals, y = estimate, color = sample_size)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
    geom_line(size = 1) +
    coord_cartesian(xlim = c(5,15), ylim = ylim) +
    labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate", 
         color = "Sample Sizes", title = main) + 
    theme_bw() +
    guides(color = guide_legend(title = "Person-Years at Risk")) +
    theme(legend.background = element_rect(colour = "black")) +
    scale_color_manual(values = c("#367E18", "#F57328"))
  
  leg <- gtable_filter(ggplot_gtable(ggplot_build(leg_plot_tmp)), "guide-box")
  
  if (i %in% c(1,2,4))
    plot_list[[i]] <- erf_strata_plot + annotation_custom(leg, xmin = 5.3, xmax = 8, ymin = breaks[5], ymax = breaks[6])
  else
    plot_list[[i]] <- erf_strata_plot + annotation_custom(leg, xmin = 12, xmax = 14.7, ymin = breaks[1], ymax = breaks[2])
  
    
}

strata_plot_tmp1 <- ggarrange(plotlist = plot_list[1:3], ncol = 3, nrow = 1, legend = "none", common.legend = TRUE)
strata_plot_tmp2 <- ggarrange(plotlist = plot_list[4:6], ncol = 3, nrow = 1, legend = "none", common.legend = TRUE)
strata_plot_tmp3 <- ggarrange(plotlist = plot_list[7:9], ncol = 3, nrow = 1, legend = "none", common.legend = TRUE)

strata_plot1 <- annotate_figure(strata_plot_tmp1, top = text_grob("65 < Age < 75", face = "bold", size = 14))
strata_plot2 <- annotate_figure(strata_plot_tmp2, top = text_grob("75 < Age < 85", face = "bold", size = 14))
strata_plot3 <- annotate_figure(strata_plot_tmp3, top = text_grob("85 < Age < 95", face = "bold", size = 14))

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/age_strata_plot.pdf", width = 16, height = 16)
ggarrange(strata_plot1, strata_plot2, strata_plot3, nrow = 3, ncol = 1, common.legend = FALSE)
dev.off()
