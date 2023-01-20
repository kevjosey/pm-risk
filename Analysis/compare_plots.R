library(ggplot2)
library(grid) 
library(pBrackets) 
library(gridExtra)
library(fst)
library(data.table)

DR100 <- read.csv("~/DR100_table.csv")[,-1]
load("~/xiao.RData"); Kevin <- xiao
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
hr_compare <- ggplot(data=plot_dat, aes(x=a.vals, y = hr, color = Methods)) +
  geom_line(aes(x = a.vals, y = hr), lwd = 1.2) +
  geom_ribbon(aes(ymin = hr.lower, ymax = hr.upper), linetype = "dashed", alpha = 0.0, lwd = 1.2) + 
  scale_color_manual(values=c("blue", "red"),
                     labels = c("Doubly Robust", "Matching"),
                     breaks = c("CAL", "Matching")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.7, 0.2)) +
  xlab(expression(paste("PM"[2.5]," (",mu,g/m^3,")"))) +
  ylab("Hazard Ratio") +
  ylim(0.75, 1.1) +
  xlim(2.76,16)

pdf(file = "~/Figures/hr_compare.pdf", width = 10, height = 8)
hr_compare
dev.off()

# Mortality Rate
ar_compare <- ggplot(data=plot_dat, aes(x=a.vals, y = hr, color = Methods)) +
  geom_line(aes(x = a.vals, y = ERC), lwd = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), linetype = "dashed", alpha = 0.0, lwd = 1.2) + 
  scale_color_manual(values=c("blue", "red"),
                     labels = c("Doubly Robust", "Matching"),
                     breaks = c("CAL", "Matching")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.7, 0.2)) +
  xlab(expression(paste("PM"[2.5]," (",mu,g/m^3,")"))) +
  ylab("Absolute Mortality Rate") +
  ylim(0.04, 0.055) +
  xlim(2.76,16)

pdf(file = "~/Figures/ar_compare.pdf", width = 10, height = 8)
ar_compare
dev.off()