#plot for flow chart:
library(ggplot2)
library(dplyr)
library(deSolve)
library(data.table)

source("util_functions")
#gen_immune_func: Function to generate immune proportion from estimated parameters
#calcRinf: Function to calculate R infinity from R0 and R initial
#LL_mcmc_sero_function: function to calculate LL

##outbreak des plots:
Sinit = 1; R0 = 1.5
Rinf1 = calcRinf(R0, Sinit, Iinit)
Rinf2 = calcRinf(R0, Sinit - 0.75*Rinf1, Iinit)

pos_vector = c(numeric(25), rep(Rinf2, 25), rep(Rinf1, 50))
data.plot = data.frame(Age = 0:99, Immune = pos_vector)
sero_data = data.frame(Age = seq(2.5, 100, 10), Immune = pos_vector[seq(5, 100, 10)], 
                       lower = pmax(0, pos_vector[seq(5, 100, 10)] - 0.05),
                       upper =  pos_vector[seq(5, 100, 10)] + 0.05)

#
arrow.pos.both <- data.frame(x = c(2, 27, 37.5, 75), xend = c(23, 48, 37.5, 75),
                  y = c(0.03, Rinf2 + 0.03, Rinf2 - 0.03, Rinf1 - 0.03), 
                  yend = c(0.03, Rinf2 + 0.03, 0.03, 0.03))
arrow.pos.one <- data.frame(x = c(25, 50), xend = c(25, 50),
                             y = c(Rinf2 + 0.23, Rinf1 + 0.23), yend = c(Rinf2 + 0.03, Rinf1 + 0.03))
annot.pot <- data.frame(x = c(37.5, 12.5, 37.5 + 9, 84, 50, 25), 
                        y = c(Rinf2 + 0.08, 0.08, Rinf2/2, Rinf1/2, Rinf1 + 0.26, Rinf2 + 0.26),
                        labels = c(paste0("outbreak~time~", 1:2), 
                                   'list(R[0], S[2](0))', 'list(R[0], S[1](0))', "1^st~outbreak", "2^nd~outbreak"))

intepret_plot_2ob <- 
  ggplot()+
  geom_pointrange(data = sero_data, aes(x = Age, y = Immune, ymin = lower, ymax = upper), size = 0.25)+
  geom_segment(data = arrow.pos.both, aes(x = x, y = y , xend = xend, yend = yend), 
               arrow = arrow(length = unit(0.05, "inches"), ends = "both", type = "closed"), color = "#b2182b", size = 1) +
  geom_segment(data = arrow.pos.one, aes(x = x, y = y , xend = xend, yend = yend), 
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"), color = "#b2182b", size = 1) +
  geom_text(data = annot.pot, aes(x = x, y = y, label = labels), parse = T, size = 3) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        rect = element_rect(fill = "transparent",colour = NA),
        text = element_text(size = 16))
ggsave(filename = "plots/flowchart/intepret_plot_2ob.png", 
       plot = intepret_plot_2ob, units = "mm", bg = "transparent", width = 100, height = 77)
#reconstruct:
reconstruct_plot_2ob <- 
  ggplot(data.plot, aes(Age, Immune))+
  geom_line(size = 1)+
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        rect = element_rect(fill = "transparent",colour = NA),
        text = element_text(size = 16))
ggsave(filename = "plots/flowchart/reconstruct_plot_2ob.png", 
       plot = reconstruct_plot_2ob, units = "mm", bg = "transparent", width = 100, height = 77)

data.plot$Immune <- c(numeric(25), head(data.plot$Immune, - 25))
AR = 0.15
S0_a_over_S0 = (1 - data.plot$Immune)/mean(1 - data.plot$Immune)
R_infty_a = 1 - S0_a_over_S0*(1 - mean(data.plot$Immune) - AR)

project_plot_data = data.frame(Age = data.plot$Age, Immune = c(data.plot$Immune, R_infty_a), 
                               Immunity = rep(c("Before outbreak", "After outbreak"), each = 100))

projected_plot_2ob <-
  ggplot(project_plot_data, aes(Age, Immune, color = Immunity))+
  geom_line(size = 1)+
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("#223559", "#6d8c9c"))+
  theme_bw()+
  labs(color = "")+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        rect = element_rect(fill = "transparent",colour = NA),
        text = element_text(size = 16),
        legend.position="bottom")  

ggsave(filename = "plots/flowchart/projected_plot_2ob.png", 
       plot = projected_plot_2ob, units = "mm", bg = "transparent", width = 100, height = 90)
  