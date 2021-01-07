library(dplyr)
library(ggplot2)
library(directlabels)

source("util_functions.R")
#gen_immune_func: Function to generate immune proportion from estimated parameters
#calcRinf: Function to calculate R infinity from R0 and R initial
#LL_mcmc_sero_function: function to calculate LL

#Pop data:
pop_data <- read.csv("data/pop_data.csv")
pop_data <- pop_data %>%
  mutate(X99 = X99 + X100) %>%
  select(-X100)
#sero data:
CHIKV_sero <- read.csv("data/CHIKV_sero_data.csv")
Study_ID = as.character(unique(CHIKV_sero$Study_id)[17])

#select study:
study_data_raw = filter(CHIKV_sero, Study_id == Study_ID, total != 0)
study_data_raw <- Reduce(rbind, rep(list(study_data_raw), 4))

#Data:
#Getting the demographic data:
country_sel = unique(study_data_raw$country) %>% as.character
country_prop_pop <- filter(pop_data, country == country_sel) %>% 
  select(starts_with("X")) %>%
  apply(1, function(x) x/sum(x)) %>% data.frame
colnames(country_prop_pop) = paste0("X", unique(pop_data$year))

l_study_time = study_data_raw$time %>% unique
study_data_id <- study_data_raw$time %>% table 
l_study_data_id <- list(1:study_data_id[1])
if(length(l_study_time) >= 2 ){
  for(i in 2:length(l_study_time)){
    l_study_data_id <- c(l_study_data_id, list((cumsum(study_data_id)[i - 1] + 1):cumsum(study_data_id)[i]))
  }
}
names(l_study_data_id) = l_study_time
study_time = l_study_time

age_u_lim = max(study_data_raw$time) - min(study_data_raw$time - study_data_raw$age_u)
age_l_lim = max(study_data_raw$time) - min(study_data_raw$time)

#Calculate R infinity:
calcRinf = function(R0, Rinit){
  optimize(
    f = function(Rinf){
      abs(1 - Rinf - exp(log(1 - Rinit) + (-R0 * (Rinf - (Rinit)))))
    },
    interval = c(Rinit + 0.01, 1))$minimum
}

#Gen prob
first_Rinit = 0;

sim_AR_func <- function(ob_year, study_data, R0, total_sample = 2000){
  ob_year = c(sort(ob_year, decreasing = T), 0)
  n_outbreak = length(ob_year)
  study_data$age_l = round(seq(0, 99, 100/nrow(study_data)))
  study_data$age_u = c(tail(study_data$age_l, -1) - 1, 99)
  study_data$age_mid = (study_data$age_l + study_data$age_u )/2
  study_data$total = rep(total_sample/nrow(study_data), nrow(study_data))
  sample_every_age = rep(total_sample/100, 100)
  
  Rinit = first_Rinit
  pos_every_age = numeric(100)
  #Conditional prob:
  ob_dist = c(head(ob_year, -1) - ob_year[-1], tail(ob_year, 1))
  ob_floor = floor(ob_year)
  for(i in 1:length(ob_dist)){
    Rinf = calcRinf(R0, Rinit)#calcRinf_large_init(R0, Rinit, k) #Recovery cov after outbreak
    prev_Rinit = Rinit
    pos_every_age = 1 - (1 - Rinf)*(1 - pos_every_age)/(1 - Rinit)
    #project
    whole_year = floor(ob_dist)[i]; part_year = ob_year[i] - floor(ob_year[i]) 
    pos_every_age <- c(numeric(length = whole_year), (1 - part_year)*pos_every_age[1], pos_every_age[c(1:(99 - whole_year))])
    pos_prop_pop = pos_every_age*country_prop_pop[,paste0("X", ifelse(study_time - ob_floor[i] - 1 < 1950, 1950, study_time - ob_floor[i] - 1))]
    Rinit = sum(pos_prop_pop)
  }
  
  return(Rinf - prev_Rinit)
}

#R0 range:
l_R0 = seq(1, 3, 0.1)
ob_year_1ob = seq(0, 99, 10)

###1ob plot:
#All R0 and number of outbreak:
comb_1obs = expand.grid(R0 = l_R0, ob =  seq(0, 99, 5))

AR_1obs = sapply(1:nrow(comb_1obs), function(x){
  ob_sel = comb_1obs$ob[x]
  R0_sel = comb_1obs$R0[x]
  sim_AR_func(ob_sel, study_data_raw, R0_sel)
})

sim_1_AR_data = cbind(comb_1obs, AR = AR_1obs)

sim_1_contour <- ggplot(sim_1_AR_data, aes(x = ob, y = R0, z = AR))+
  geom_contour_filled(binwidth = 0.05 ) +
  scale_y_continuous(breaks = seq(1, 3, 0.2),
                     label = seq(1, 3, 0.2), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(10, 90, 10), expand = c(0, 0)) +
  labs(fill = "IAR", x = "Time since the last outbreak", y = expression(R[0])) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90),
        text = element_text(size = 18))

ggsave(filename = "plots/simulation/simulation_1ob.png", sim_1_contour, width = 7.5, height = 6, bg = "transparent")

###2ob plot:
#All R0 and number of outbreak:
comb_years = combn(ob_year_1ob, 2) %>% t %>% data.frame
comb_years = comb_years[sapply(1:nrow(comb_years), function(x) {
  ob_year = rev(as.numeric(comb_years[x,]))
  all(c(head(ob_year, -1) - ob_year[-1], tail(ob_year, 1)) > 5)
  }), ]
# comb_years = combn(ob_year_1ob, 2) %>% apply(2, rev) %>% t %>% data.frame 
# comb_years = comb_years[order(comb_years$X1),]
# comb_years = comb_years[sapply(1:nrow(comb_years), function(x) {
#   ob_year = as.numeric(comb_years[x,])
#   all(c(head(ob_year, -1) - ob_year[-1], tail(ob_year, 1)) > 5)
# }), ]
comb_2obs = cbind(R0 = rep(l_R0, each = nrow(comb_years)), comb_years %>% list %>% rep(length(l_R0)) %>% Reduce(rbind, .))

AR_2obs = c()
for(x in 1:nrow(comb_2obs)){
  ob_sel = as.numeric(select(comb_2obs, starts_with("X"))[x,])
  R0_sel = as.numeric(comb_2obs$R0[x])
  AR_2obs = c(AR_2obs, sim_AR_func(ob_sel, study_data_raw, R0_sel))
}

sim_2_AR_data = cbind(comb_2obs, AR = AR_2obs) %>%
  mutate(ob = paste0(X1, "-", X2))

###3ob plot:
#All R0 and number of outbreak:
comb_years = combn(ob_year_1ob, 3) %>% apply(2, rev) %>% t %>% data.frame 
comb_years = comb_years[order(comb_years$X1),]
comb_years = comb_years[sapply(1:nrow(comb_years), function(x) {
  ob_year = as.numeric(comb_years[x,])
  all(c(head(ob_year, -1) - ob_year[-1], tail(ob_year, 1)) > 5)
}), ]
comb_3obs = cbind(R0 = rep(l_R0, each = nrow(comb_years)), comb_years %>% list %>% rep(length(l_R0)) %>% Reduce(rbind, .))

AR_3obs = c()
for(x in 1:nrow(comb_3obs)){
  ob_sel = as.numeric(select(comb_3obs, starts_with("X"))[x,])
  R0_sel = as.numeric(comb_3obs$R0[x])
  AR_3obs = c(AR_3obs, sim_AR_func(ob_sel, study_data_raw, R0_sel))
}

sim_3_AR_data = cbind(comb_3obs, AR = AR_3obs) %>%
  mutate(ob = paste0(X1, "-", X2, "-", X3))

###4ob plot:
#All R0 and number of outbreak:
comb_years = data.frame(t(combn(ob_year_1ob, 4))) 
comb_years = comb_years[sapply(1:nrow(comb_years), function(x) {
  ob_year = rev(as.numeric(comb_years[x,]))
  all(c(head(ob_year, -1) - ob_year[-1], tail(ob_year, 1)) > 5)
}), ]
comb_4obs = cbind(R0 = rep(l_R0, each = nrow(comb_years)), comb_years %>% list %>% rep(length(l_R0)) %>% Reduce(rbind, .))

AR_4obs = c()
for(x in 1:nrow(comb_4obs)){
  ob_sel = as.numeric(select(comb_4obs, starts_with("X"))[x,])
  R0_sel = as.numeric(comb_4obs$R0[x])
  AR_4obs = c(AR_4obs, sim_AR_func(ob_sel, study_data_raw, R0_sel))
}

sim_4_AR_data = cbind(comb_4obs, AR = AR_4obs) %>%
  mutate(ob = paste0(X1, "-", X2, "-", X3, "-", X4))

#focus on pattern:
#2ob:
sim_2_AR_data_mod = rbind(sim_2_AR_data %>% filter(X2 == 90), sim_2_AR_data %>% filter(X1 == 10))
sim_2_AR_data_mod$pattern = rep(c("First outbreak constant", "Second outbreak constant"), each = nrow(sim_2_AR_data_mod)/2)

sim_2_AR_contour = sim_2_AR_data_mod
sim_2_AR_contour$ob_cont = c(rep(seq(10, 80, 10), 21), rep(seq(20, 90, 10), 21))
sim_2_contour <- ggplot(sim_2_AR_contour, aes(x = ob_cont, y = R0, z = AR, subgroup = pattern))+
  geom_contour_filled(binwidth = 0.05 ) +
  scale_y_continuous(breaks = seq(1, 3, 0.2),
                     label = seq(1, 3, 0.2), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(10, 90, 10), expand = c(0, 0)) +
  facet_wrap(.~pattern, scale = "free_x") +
  labs(x = "Timing of previous outbreaks", y = expression(R[0]), fill = "IAR")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90),
        text = element_text(size = 18),
        panel.spacing = unit(2, "lines"))

ggsave(filename = "plots/simulation/simulation_2ob.png", sim_2_contour, width = 10, height = 6, bg = "transparent")

#3ob:
sim_3_AR_data_mod = rbind(sim_3_AR_data %>% filter(X1 == 90, X2 == 80), sim_3_AR_data %>% filter(X3 == 10, X2 == 20), sim_3_AR_data %>% filter(X3 == 10, X1 == 90))
sim_3_AR_data_mod$pattern = mapply(rep, c("First and second outbreaks constant", "Second and third outbreaks constant", "First and third outbreaks constant"), c(147, 147, 147)) %>% as.vector

sim_3_AR_contour = sim_3_AR_data_mod
sim_3_AR_contour$ob_cont = c(rep(seq(10, 70, 10), 21), rep(seq(30, 90, 10), 21), rep(seq(20, 80, 10), 21))
sim_3_contour <- ggplot(sim_3_AR_contour, aes(x = ob_cont, y = R0, z = AR, subgroup = pattern))+
  geom_contour_filled(binwidth = 0.05 ) +
  scale_y_continuous(breaks = seq(1, 3, 0.2),
                     label = seq(1, 3, 0.2), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(10, 90, 10), expand = c(0, 0)) +
  facet_wrap(.~pattern, scale = "free_x", ncol = 2) +
  labs(fill = "IAR", x = "Timing of previous outbreaks", y = expression(R[0])) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90),
        text = element_text(size = 18),
        panel.spacing = unit(2, "lines"))

ggsave(filename = "plots/simulation/simulation_3ob.png", sim_3_contour, width = 10, height = 7, bg = "transparent")

#4ob:
sim_4_AR_data_mod = rbind(sim_4_AR_data %>% filter(X2 == 70, X3 == 80, X4 == 90), sim_4_AR_data %>% filter(X3 == 30, X2 == 20, X1 == 10),
                          sim_4_AR_data %>% filter(X4 == 90, X3 == 80, X1 == 10), sim_4_AR_data %>% filter(X4 == 90, X2 == 20, X1 == 10))
sim_4_AR_data_mod$pattern = mapply(rep, c("First, second, and third outbreaks constant", 
                                          "Second, third, and forth outbreaks constant", 
                                          "First, second, and forth outbreaks constant", 
                                          "First, third, and forth outbreaks constant"), rep(126, 4)) %>% as.vector() # c(84, 84, 105, 105) 

sim_4_AR_contour = sim_4_AR_data_mod
sim_4_AR_contour$ob_cont = c(rep(seq(10, 60, 10), 21), rep(seq(40, 90, 10), 21), rep(seq(20, 70, 10), 21), rep(seq(30, 80, 10), 21))

sim_4_contour <- ggplot(sim_4_AR_contour, aes(x = ob_cont, y = R0, z = AR, subgroup = pattern))+
  geom_contour_filled(binwidth = 0.05 ) +
  scale_y_continuous(breaks = seq(1, 3, 0.2),
                     label = seq(1, 3, 0.2), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(10, 90, 10), expand = c(0, 0)) +
  facet_wrap(.~pattern, nrow = 2, scale = "free_x") +
  labs(fill = "IAR", x = "Timing of previous outbreaks", y = expression(R[0])) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90),
        text = element_text(size = 18),
        panel.spacing = unit(2, "lines"))

ggsave(filename = "plots/simulation/simulation_4ob.png", sim_4_contour, width = 11, height = 7, bg = "transparent")
