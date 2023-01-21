library(ggplot2)
library(grid) 
library(pBrackets) 
library(gridExtra)
library(fst)
library(data.table)

DR100 <- read.csv("~/Data/DR100_table.csv")[,-1]
load("~/Data/xiao.RData"); Kevin <- xiao
colnames(Kevin) = colnames(DR100)

DR100_table <- rbind(DR100, Kevin)
DR100_table <- subset(DR100_table, !(Methods %in% c("LM", "Doubly Robust")))

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
  coord_cartesian(xlim = c(5,13), ylim = c(0.88, 1.02)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Hazard Ratio") +
  scale_color_manual(values=c("blue", "red"),
                     labels = c("Doubly Robust", "Matching"),
                     breaks = c("CAL", "Matching")) +
  scale_y_continuous(breaks = c(0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02)) +
  scale_x_continuous(breaks = c(5,6,7,8,9,10,11,12,13))

pdf(file = "~/Figures/hr_compare.pdf", width = 8, height = 8)
hr_compare
dev.off()

# Mortality Rate
ar_compare <- ggplot(data = plot_dat, aes(x = a.vals, y = ERC, color = Methods)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = c(0.7, 0.2),
        legend.background = element_rect(colour = "black"),) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Mortality Rate") +
  coord_cartesian(xlim = c(5,15), ylim = c(0.044,0.052)) +
  scale_color_manual(values=c("blue", "red"),
                     labels = c("Doubly Robust", "Matching"),
                     breaks = c("CAL", "Matching")) +
  scale_y_continuous(breaks = c(0.044,0.045,0.046,0.047,0.048,0.049,0.05,0.051,0.052)) +
  scale_x_continuous(breaks = c(5,6,7,8,9,10,11,12,13,14,15))


pdf(file = "~/Figures/ar_compare.pdf", width = 8, height = 8)
ar_compare
dev.off()