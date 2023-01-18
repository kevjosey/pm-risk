library("ggplot2")
library(grid) 
library(pBrackets) 
library(gridExtra)
library(fst)
library(data.table)

DR100 <- read.csv("~/Documents/DR100_table.csv")[,-1]
load("~/Documents/xiao.RData"); Kevin <- xiao
colnames(Kevin) = colnames(DR100)

DR100_table <- rbind(DR100, Kevin)
DR100_table <- subset(DR100_table, Methods != "CAL")

#pdf("/Users/apple/Dropbox/continuous GPS/JASA resubmission/MatchingvsDR.pdf",width = 11.5/1.6, height = 8.5/1.6)
ggplot(data=DR100_table, aes(x=a.vals, y = ERC, color = Methods)) +
  #geom_point(data = subset(DR_all_nosmooth,band=="DR100"),aes(x=a.vals.nosmooth,y=dose_matching.response.nonsmooth), lwd=0.8) +
  geom_line(aes(x= a.vals, y= ERC), lwd = 1.2) +
  geom_ribbon(aes(ymin=lower, ymax=upper), linetype = "dashed", alpha=0.0, lwd = 1.2) + 
  #geom_line(aes(x= DR100_a.vals_boots[1,], y= ERC[1:100]+1.96*DR100_sd),linetype="dashed" ,color="red", lwd=1.2) +
  #geom_line(aes(x= DR100_a.vals_boots[1,], y= ERC[1:100]-1.96*DR100_sd),linetype="dashed" ,color="red" , lwd=1.2) +
  #geom_line(aes(x= DR100_a.vals_boots[1,], y= ERC[1:100]-1.96*DR100_sd),linetype="dashed" ,color="red" , lwd=1.2) +
  #geom_line(aes(x= DR100_a.vals_boots[1,], y= DR_kennedy) ,color="blue", lwd=1.2) +
  #scale_colour_discrete(values=c("#999999", "#E69F00"),
  #  labels=c("Matching", "Doubly Robust")) +
  scale_color_manual(values=c("green", "blue", "red"),
                     labels = c("Kevin - Doubly Robust", "Xiao - Doubly Robust", "Matching"),
                     breaks = c("CAL", "Doubly Robust", "Matching")) +
  theme_bw() +
  #guides(colour = FALSE) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.7, 0.2)) +
  #guides( linetype = FALSE) +
  ggtitle(expression(paste("Causal Exposure-response Curve: PM"[2.5]," v.s. Mortality"))) + 
  xlab(expression(paste("PM"[2.5]," (",mu,g/m^3,")"))) +
  ylab("Mortality Rate (per 100 Medicare Enrollees") +
  ylim(0.027, 0.056) +
  xlim(2.76,17.16)
