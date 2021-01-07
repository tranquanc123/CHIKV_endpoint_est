#Plot for paper
library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
# clear memory
rm(list=ls())
###################
###Summary results:
###################

###Endpoint events:
#loading data:
Rinit_files_name = list.files("results/real_sero_study/immunity/immune_by_age_in_pop/")
CHIKV_sero <- read.csv("data/CHIKV_sero_data.csv")

#Summarise funtion:
quantile_func = function(x) c(quantile(x, c(0.025, 0.5, 0.975)))
expit = function(x) x#exp(x)/(exp(x) + 1)

#file name replace:
file_name_dic <- unique(CHIKV_sero[, c("Study_id", "Study_name")]) 
file_name_dic$file_names = paste0(file_name_dic$Study_id, ".rds")
mod_study_names = file_name_dic$Study_name[sapply(Rinit_files_name, function(x) which(file_name_dic$file_names %in% x))]
study_index = c(paste0("I", 1:8), "D", "O", "C", "M", "G", "P", "Q", "J", "K", "A", "E", "N", "B", "H1", "H2", paste0("L", 1:5), "F", paste0("R", 1:6))

#Overall Endpoints:
Endpoints_l = c()
R0_l = c()
pop2020_l = c()
for(file in 1:length(Rinit_files_name)){
  file_names = str_sub(Rinit_files_name[file], 1, -5)
  endpoints = readRDS(paste0("results/real_sero_study/number_endpoints/" , file_names, ".rds"))
  model_prop = endpoints$model_prop
  
  #get Endpoints:
  Endpoints_l = c(Endpoints_l, list(endpoints$endpoints))
  
  #get R0:
  en_R0 = lapply(1:length(model_prop), function(x) endpoints$R0[[x]]*model_prop[x]) %>% Reduce("+", .)
  R0_l = c(R0_l, list(en_R0))
  
  #get pop:
  pop2020_l = c(pop2020_l, list(endpoints$pop2020))
}
names(Endpoints_l) = mod_study_names

#Summary Endpoints overall and R0:
Overall_Endpoints_R0 = data.frame(Study_ID = mod_study_names, ID = study_index,
                      sapply(Endpoints_l, function(x) quantile_func(rowSums(x))) %>% t, 
                      sapply(R0_l, quantile_func) %>% t)
colnames(Overall_Endpoints_R0) = c("Study_ID", "ID", paste0("Endpoints_", c(0.025, "0.500", 0.975)), paste0("R0_", c(0.025, "0.500", 0.975)))
Overall_Endpoints_R0$Endpoints_sel = ifelse(Overall_Endpoints_R0$Endpoints_0.025 > 0.01, " > 10/1000", " < 10/1000")

#Divide to adult and children:
Children_Endpoints_RO = data.frame(Study_ID = mod_study_names, ID = study_index,
                            sapply(1:length(Endpoints_l), function(x) quantile_func(rowSums(Endpoints_l[[x]][,1:15])/sum(pop2020_l[[x]][1:15]))) %>% t, 
                            sapply(R0_l, quantile_func) %>% t)
Children_Endpoints_RO$Age_group = "Children"

Adult_Endpoints_RO = data.frame(Study_ID = mod_study_names, ID = study_index,
                            sapply(1:length(Endpoints_l), function(x) quantile_func(rowSums(Endpoints_l[[x]][,16:100])/sum(pop2020_l[[x]][16:100]) )) %>% t, 
                            sapply(R0_l, quantile_func) %>% t)
Adult_Endpoints_RO$Age_group = "Adult"
Broad_age_Endpoints_R0 = rbind(Children_Endpoints_RO, Adult_Endpoints_RO)
colnames(Broad_age_Endpoints_R0) = c("Study_ID", "ID",paste0("Endpoints_", c(0.025, "0.500", 0.975)), paste0("R0_", c(0.025, "0.500", 0.975)), "Age_group")
Broad_age_Endpoints_R0$Endpoints_sel = ifelse(Broad_age_Endpoints_R0$Endpoints_0.025 > 0.01, " > 10/1000", " < 10/1000")

#Endpoints for every age group:
age_group_5 = seq(0, 99, 5)
age_group_5_Endpoints_R0 <-
  lapply(1:length(age_group_5), function(y) {
  age = age_group_5[y]
  temp_data = data.frame(Study_ID = mod_study_names, ID = study_index,
             sapply(1:length(Endpoints_l), function(x) quantile_func(rowSums(Endpoints_l[[x]][,(age + 1):(age + 5)]))) %>% t, 
             sapply(R0_l, quantile_func) %>% t)
  colnames(temp_data) =  c("Study_ID", "ID", paste0("Endpoints_", c(0.025, "0.500", 0.975)), paste0("R0_", c(0.025, "0.500", 0.975)))
  temp_data$age_group = paste0(formatC(age, width = 2, flag = 0), "-", formatC(age + 4, width = 2, flag = 0))
  temp_data$pop2020 = sapply(pop2020_l, function(x) sum(x[(age + 1):(age + 5)]))
  temp_data = temp_data[order(temp_data$Endpoints_0.500),]
  return(temp_data)
  }
) %>% Reduce(rbind, .)

age_group_5_Endpoints_R0$Study_ID = factor(age_group_5_Endpoints_R0$Study_ID, levels = unique(age_group_5_Endpoints_R0$Study_ID))

###Sus at 2020:
sus_pop_2020 = c()
immune_pop_age_group_pre2020 = c()
immune_pop_age_group_post2020 = c()
for(file in 1:length(Rinit_files_name)){
  file_names = str_sub(Rinit_files_name[file], 1, -5)
  immunity_pre_ob = readRDS(paste0("results/real_sero_study/immunity/immune_projected/", file_names, ".rds"))
  immunity_post_ob = readRDS(paste0("results/real_sero_study/immunity/immune_post_ob/", file_names, ".rds"))
  model_prop = immunity_pre_ob$model_prop
  
  immune_pop_age_group = lapply(1:length(model_prop), function(x) immunity_pre_ob$projected_immunity_in_pop[[x]]*model_prop[x]) %>% Reduce("+", .)
  immune_pop_age_group_pre2020 = c(immune_pop_age_group_pre2020, list(immune_pop_age_group))
  immune_pop_age_group_post2020 = c(immune_pop_age_group_post2020, list(lapply(1:length(model_prop), function(x) immunity_post_ob$immunity_post_ob[[x]]*model_prop[x]) %>% Reduce("+", .) ))
  
  immune_pop = rowSums(immune_pop_age_group)
  sus_pop_2020 = c(sus_pop_2020, list(1 - immune_pop))
}
sus_pop_2020_sum = sus_pop_2020 %>% sapply(quantile_func) %>% t %>% data.frame
colnames(sus_pop_2020_sum) = paste0("S0_", c(0.025, "0.500", 0.975))
Overall_Endpoints_R0_S0 = cbind(Overall_Endpoints_R0, sus_pop_2020_sum)

#Immunity and S(0) for every age group
age_group = 10
age_group_10 = seq(0, 99, age_group)
age_group_10_Sinit_Rinit_R0 <-
  lapply(1:length(age_group_10), function(y) {
    age = age_group_10[y]
    temp_data = data.frame(Study_ID = mod_study_names, ID = study_index,
                           sapply(1:length(immune_pop_age_group_pre2020), function(x) quantile_func(rowSums(immune_pop_age_group_pre2020[[x]][,(age + 1):(age + age_group)]))) %>% t,
                           sapply(1:length(immune_pop_age_group_post2020), function(x) quantile_func(rowSums(immune_pop_age_group_post2020[[x]][,(age + 1):(age + age_group)]))) %>% t, 
                           sapply(R0_l, quantile_func) %>% t)
    colnames(temp_data) =  c("Study_ID", "ID", paste0("Rinit_", c(0.025, "0.500", 0.975)),paste0("Rinf_", c(0.025, "0.500", 0.975)),  paste0("R0_", c(0.025, "0.500", 0.975)))
    temp_data$age_group = paste0(formatC(age, width = 2, flag = 0), "-", formatC(age + 4, width = 2, flag = 0))
    temp_data$pop2020 = sapply(pop2020_l, function(x) sum(x[(age + 1):(age + age_group)]))
    temp_data = temp_data[order(temp_data$Rinit_0.500),]
    return(temp_data)
  }
  ) %>% Reduce(rbind, .)

age_group_10_Sinit_Rinit_R0$Study_ID = factor(age_group_10_Sinit_Rinit_R0$Study_ID, levels = unique(age_group_5_Endpoints_R0$Study_ID))
###Time since last outbreak:
en_ob_year2020_min = c()
year_proj = 2020
for(file in 1:length(Rinit_files_name)){
  file_names = str_sub(Rinit_files_name[file], 1, -5)
  para.list = readRDS(paste0("results/real_sero_study/RJMCMC/", file_names, ".rds"))
  ob_year = lapply(para.list$l_output, function(x) x[,3:ncol(x)] %>% tail(1000) %>% matrix(nrow = 1000))
  ob_year_min = lapply(ob_year, function(x) apply(x, 1, min))
  
  #get model prob:
  model_prop = readRDS(paste0("results/real_sero_study/immunity/immune_projected/", file_names, ".rds"))$model_prop
  en_ob_year_min = lapply(1:length(model_prop), function(x) ob_year_min[[x]]*model_prop[x]) %>% Reduce("+", .)
  
  #project to 2020:
  study_year = filter(CHIKV_sero, Study_id == file_names)$time %>% max 
  en_ob_year2020_min = c(en_ob_year2020_min, list(en_ob_year_min + year_proj - study_year))
}
en_ob_year2020_min_sum = en_ob_year2020_min %>% sapply(quantile_func) %>% t %>% data.frame
colnames(en_ob_year2020_min_sum) = paste0("ob_", c(0.025, "0.500", 0.975))
Overall_Endpoints_R0_ob = cbind(Overall_Endpoints_R0, en_ob_year2020_min_sum)

#divided by children and adult:
Children_Endpoints_RO_S0 = cbind(Children_Endpoints_RO,
                          sapply(1:length(immune_pop_age_group_pre2020), function(x) {
                            1 - quantile_func(rowSums(immune_pop_age_group_pre2020[[x]][,1:15])/sum(pop2020_l[[x]][1:15]))
                            }) %>% t)

Adult_Endpoints_RO_S0 = cbind(Adult_Endpoints_RO,
                       sapply(1:length(immune_pop_age_group_pre2020), function(x) {
                         1 - quantile_func(rowSums(immune_pop_age_group_pre2020[[x]][,16:100])/sum(pop2020_l[[x]][16:100]))
                       }) %>% t)

Broad_age_Endpoints_R0_S0 = rbind(Children_Endpoints_RO_S0, Adult_Endpoints_RO_S0)
colnames(Broad_age_Endpoints_R0_S0) = c("Study_ID", "ID", paste0("Endpoints_", c(0.025, "0.500", 0.975)), paste0("R0_", c(0.025, "0.500", 0.975)), "Age_group",
                                 paste0("S0_", c(0.025, "0.500", 0.975)))
Broad_age_Endpoints_R0_S0$Endpoints_sel = ifelse(Broad_age_Endpoints_R0_S0$Endpoints_0.025 > 0.01, " > 10/1000", " < 10/1000")
##################
###PLot results
##################

# #Overall Endpoints, R0 plot:
# ggplot(Overall_Endpoints_R0, aes(x = Endpoints_0.500*1000, y = R0_0.500, color = Endpoints_sel))+
#   geom_pointrange(aes(ymin = R0_0.025, ymax = R0_0.975))+
#   geom_errorbar(aes(xmin = Endpoints_0.025*1000, xmax = Endpoints_0.975*1000)) +
#   scale_color_manual(values = c("#4575b4", "#d73027")) +
#   #scale_x_continuous(limits = c(0, 1)) +
#   labs(x = "Number of endpoint events (per thousand)", y = "RO") +
#   theme_bw()+
#   theme(legend.title = element_blank())

# #divided by age group Endpoints, R0 plot:
# ggplot(Broad_age_Endpoints_R0, aes(x = Endpoints_0.500*1000, y = R0_0.500, color = Endpoints_sel))+
#   geom_pointrange(aes(ymin = R0_0.025, ymax = R0_0.975))+
#   geom_errorbar(aes(xmin = Endpoints_0.025*1000, xmax = Endpoints_0.975*1000)) +
#   facet_grid(.~Age_group) + 
#   scale_color_manual(values = c("#4575b4", "#d73027")) +
#   #scale_x_continuous(limits = c(0, 1)) +
#   labs(x = "Number of endpoint events (per thousand)", y = "RO") +
#   theme_bw()+
#   theme(legend.title = element_blank())

#Time from the last outbreak and Endpoints:
Timing_endpoints <- 
  ggplot(Overall_Endpoints_R0_ob, aes(x = ob_0.500, y = Endpoints_0.500*1000, color = R0_0.500))+
  geom_pointrange(aes(xmin = ob_0.025, xmax = ob_0.975))+
  geom_errorbar(aes(ymin = Endpoints_0.025*1000, ymax = Endpoints_0.975*1000, linetype = Endpoints_sel)) +
  #geom_hline(yintercept = 10, color = "#878787", linetype = 4, size = 1.5) +
  geom_text(data = filter(Overall_Endpoints_R0_ob, Endpoints_0.025 > 0.01), aes(x = ob_0.500, y = Endpoints_0.500*1000, label = ID), nudge_x = 1, nudge_y = .1) +
  scale_linetype_manual(values = c(2,1)) +
  scale_colour_gradientn(colors = c('#E6D460','#C7421A', '#1a3a3a'))+
  scale_y_continuous(trans = "log10", breaks = c(10, seq(100, 600, 100))) +
  scale_x_continuous(breaks = seq(0, 60, 10)) +
  labs(y = expression(paste("Number of endpoint events per thousand (", log[10], " scale)")), 
       x = "Time since the last outbreak", color = expression(R[0]), linetype = " ") +
  theme_bw()+
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        text = element_text(size = 20))
ggsave(filename = "CHIKV_endpoints_event/plots/Timing_endpoints.png", Timing_endpoints, width = 10, height = 7.5, bg = "transparent")

  #S0, R0, and Endpoints:
S0_R0_ep <- ggplot(Overall_Endpoints_R0_S0, aes(y = Endpoints_0.500*1000, x = R0_0.500, color = S0_0.500))+
  #geom_hline(yintercept = 10, color = "#878787", linetype = 4, size = 1.5) +
  geom_pointrange(aes(xmin = R0_0.025, xmax = R0_0.975))+
  geom_errorbar(aes(ymin = Endpoints_0.025*1000, ymax = Endpoints_0.975*1000, linetype = Endpoints_sel)) +
  geom_text(data = filter(Overall_Endpoints_R0_S0, Endpoints_0.025 > 0.01), aes(x = R0_0.500, y = Endpoints_0.500*1000, label = ID), nudge_x = .07, nudge_y = .07, color = "black") +
  scale_linetype_manual(values = c(2,1)) +
  scale_colour_gradient(low = "#E6D460", high = "#1a3a3a") +
  #scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(trans = "log10", breaks = c(10, seq(100, 600, 100))) +
  labs(y = expression(paste("Number of endpoint events per thousand (", log[10], " scale)")), 
       x = expression(R[0]), color = "Susceptible proportion \n before the outbreak", linetype = " ") +
  theme_bw()+
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title = element_text(size = 20))
ggsave(filename = "CHIKV_endpoints_event/plots/S0_R0_endpoints.png", S0_R0_ep, width = 10, height = 7.5, bg = "transparent")

# #divided by age group S0, Endpoints, R0 plot:
# ggplot(Broad_age_Endpoints_R0_S0, aes(x = Endpoints_0.500*1000, y = R0_0.500, color = S0_0.500))+
#   geom_pointrange(aes(ymin = R0_0.025, ymax = R0_0.975))+
#   geom_errorbar(aes(xmin = Endpoints_0.025*1000, xmax = Endpoints_0.975*1000, linetype = Endpoints_sel)) +
#   facet_grid(.~Age_group) + 
#   scale_colour_gradient(low = "#2166ac", high = "#b2182b") +
#   #scale_x_continuous(limits = c(0, 1)) +
#   labs(x = "Number of endpoint events (per thousand)", y = "RO",  color = "Susceptible proportion \n before the outbreak", linetype = " ") +
#   theme_bw()+
#   theme(axis.title = element_text(size = 14))

#Children vs Adult:
Children_Adult_comp_Endpoints_R0_S0 = merge(Broad_age_Endpoints_R0_S0 %>% filter(Age_group == "Children"), Broad_age_Endpoints_R0_S0 %>% filter(Age_group == "Adult"),
                                     by = "Study_ID")
Children_adult_ep <- ggplot(Children_Adult_comp_Endpoints_R0_S0, aes(x = Endpoints_0.500.x*1000, y = Endpoints_0.500.y*1000))+
  geom_abline(slope = 1, color = "#878787", linetype = 4, size = 1.5) +
  geom_point(aes(size = R0_0.500.x), shape = 1) +
  geom_errorbar(aes(xmin = Endpoints_0.025.x*1000, xmax = Endpoints_0.975.x*1000, linetype = Endpoints_sel.x, color = S0_0.500.x)) +
  geom_linerange(aes(ymin = Endpoints_0.025.y*1000, ymax = Endpoints_0.975.y*1000, linetype = Endpoints_sel.y, color = S0_0.500.y)) +
  geom_text(data = filter(Children_Adult_comp_Endpoints_R0_S0, Endpoints_0.025.x > 0.01|Endpoints_0.025.y > 0.01), 
            aes(label = ID.x), nudge_x = .1, nudge_y = -.1, color = "black") +
  scale_linetype_manual(values = c(2,1)) +
  scale_colour_gradient(low = "#E6D460", high = "#1a3a3a") +
  #scale_x_continuous(limits = c(0, 1)) +
  labs(x = expression(paste("Number of endpoints per thousand in children (", log[10], " scale)")), 
       y = expression(paste("Number of endpoints per thousand in adults (", log[10], " scale)")),
       color = "Susceptible proportion \n before the outbreak", linetype = " ", size = expression(R[0])) +
  scale_size_continuous(breaks = seq(1, 3, 0.25)) +
  #sacle_color
  scale_y_continuous(trans = "log10", breaks = c(10, seq(100, 600, 100)), limits = c(0.20, 600)) +
  scale_x_continuous(trans = "log10", breaks = c(10, seq(100, 600, 100)), limits = c(0.20, 600)) +
  theme_bw()+
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = -90))
ggsave(filename = "CHIKV_endpoints_event/plots/Children_adult_endpoints.png", Children_adult_ep, width = 10, height = 7.5, bg = "transparent")

# # Endpoints for every age group:
# ggplot(data = age_group_5_Endpoints_R0, mapping = aes(x = age_group)) +
#   geom_bar(aes(y = Endpoints_0.500), stat = "identity", color = "#e31a1c", fill = "#e31a1c")+
#   geom_linerange(aes(ymin = Endpoints_0.025, ymax = Endpoints_0.975), color = "#2166ac") +
#   facet_wrap(.~Study_ID, ncol = 7)+#, scales = "free_y") +
#   labs(x = "Age group", y = "Proportion") +
#   theme_bw()+
#   theme(axis.text.x = element_blank())

#immunity profile for every study:
# ggplot(data = age_group_10_Sinit_Rinit_R0, mapping = aes(x = age_group)) +
#   geom_bar(aes(y = pop2020), stat = "identity", color = "#fb9a99", fill = "#fb9a99") +
#   geom_bar(aes(y = Rinit_0.500), stat = "identity", color = "#e31a1c", fill = "#e31a1c")+
#   geom_linerange(aes(ymin = Rinit_0.025, ymax = Rinit_0.975), color = "#2166ac") +
#   facet_wrap(.~Study_ID, ncol = 7) +
#   labs(x = "Age group", y = "Proportion") +
#   theme_bw()+
#   theme(axis.text.x = element_blank())

age_group_10_Sinit_Rinit_plot_data <- age_group_10_Sinit_Rinit_R0 %>% select(Study_ID, ID, age_group) %>% list() %>% rep(3) %>% Reduce(rbind, .)
age_group_10_Sinit_Rinit_plot_data$Immunity_value = c(age_group_10_Sinit_Rinit_R0$Rinit_0.500, age_group_10_Sinit_Rinit_R0$Rinf_0.500 - 
                                                        age_group_10_Sinit_Rinit_R0$Rinit_0.500, 
                                                      age_group_10_Sinit_Rinit_R0$pop2020 - age_group_10_Sinit_Rinit_R0$Rinf_0.500)
age_group_10_Sinit_Rinit_plot_data$Immunity = factor(rep(c("Before outbreak", "After outbreak", "Susceptible"), each = nrow(age_group_10_Sinit_Rinit_R0)),
                                                     levels = c("Susceptible", "After outbreak", "Before outbreak"))
  
color_pal <- c("#b4d08a","#223559", "#6d8c9c")
ggplot(data = age_group_10_Sinit_Rinit_plot_data, mapping = aes(x = age_group)) +
  #geom_bar(aes(y = pop2020), stat = "identity", color = "#99A776", fill = "#99A776") +
  geom_bar(aes(y = Immunity_value, color = Immunity, fill = Immunity), stat = "identity", position = "stack") +
  geom_text(data = age_group_10_Sinit_Rinit_plot_data %>% select(ID, Study_ID) %>% unique, 
            aes(x = 9.5, y = 0.25, label = ID), color = "black") +
  #geom_bar(aes(y = Rinit_0.500), stat = "identity", color = "#9ebcda", fill = "#9ebcda")+
  scale_color_manual(values = color_pal) +
  scale_fill_manual(values = color_pal) +
  facet_wrap(.~Study_ID, ncol = 7) +
  labs(x = "Age group", y = "Proportion") +
  theme_bw()+
  theme(axis.text.x = element_blank())      

#########Summary model prop:  
# l_model_prop_sum = c()
# for(file in 1:length(Rinit_files_name)){
#   file_names = str_sub(Rinit_files_name[file], 1, -5)
#   country_names = filter(CHIKV_sero, Study_id == str_sub(Rinit_files_name[file], 1, -5))$country %>% unique %>% as.character()
#   print(file_names)
#   
#   #model prop
#   seeds_output = readRDS(paste0("estimate_modelprop/crc/RJMCMC_no_prior_R0/", file_names, ".rds"))
#   model_prop = matrix(as.numeric(seeds_output$seeds$result$`Posterior Model Probabilities`), nrow = 1)
#   colnames(model_prop) = paste0("ob", 1:5)
#   
#   one_row = data.frame(study = file_names, country = country_names, model_prop)
#   l_model_prop_sum  = c(l_model_prop_sum, list(one_row))
# }
# model_prop_sum = Reduce(rbind, l_model_prop_sum)
# 
# m.model_prop_sum = melt(model_prop_sum, id.vars = c("study", "country"), variable.name = "n_outbreak", value.name = "proportion")
# 
#   model_prop_plot <- ggplot(m.model_prop_sum, aes(x = study, y = proportion, group = country, fill = n_outbreak))+
#   geom_bar(stat = "identity")+
#   facet_wrap(.~country, scale = "free") +
#   theme_bw() +
#   scale_y_continuous(limits = c(0, NA))+
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank())#,
# ggsave("chikv_sample_cal/results/summary_model_prop_plot.png", model_prop_plot, width = 10)

#######summary parameters:
l_model_par_sum = c()
l_model_prop_sum = c()
facet_names = c("R[0]", "bar(S)[a]", "~~Outbreak~~time", "~~Model~~Probability")
for(file in 1:length(Rinit_files_name)){
  file_names = str_sub(Rinit_files_name[file], 1, -5)
  country_names = filter(CHIKV_sero, Study_id == str_sub(Rinit_files_name[file], 1, -5))$ISO %>% unique %>% as.character()
  #country_names = str_replace_all(country_names, " ", "~~")
  country_study_years = filter(CHIKV_sero, Study_id == str_sub(Rinit_files_name[file], 1, -5))$time %>% max 
  print(file_names)
  
  seeds_output = readRDS(paste0("results/real_sero_study/RJMCMC/", file_names, ".rds"))
  #extract parameter:
  model_prop = as.numeric(seeds_output$seeds$result$`Posterior Model Probabilities`)
  R0_sum = lapply(1:length(model_prop), function(x) seeds_output$l_output[[x]][,1]*model_prop[x]) %>% Reduce("+", .) %>% tail(1000) %>% quantile_func()
  barSa_sum = lapply(1:length(model_prop), function(x) seeds_output$l_output[[x]][,2]*model_prop[x]) %>% Reduce("+", .) %>% tail(1000) %>% expit %>% quantile_func()
  n_ob_choosen = which.max(model_prop)
  ob_sum = seeds_output$l_output[[n_ob_choosen]][,3:(2 + n_ob_choosen)] %>% matrix(ncol = n_ob_choosen) %>% tail(1000) %>% 
    apply(1, sort, decreasing = T) %>% matrix(nrow = n_ob_choosen) %>% apply(1, quantile_func) %>% t
  ob_year_sum = country_study_years - ob_sum
  
  #format to plot:
  par_sum_data = rbind(R0_sum, barSa_sum, ob_year_sum) %>% data.frame %>% 
    mutate(Study_ID = study_index[file], country = country_names, 
           par_names = c("~~R[0]", "~~bar(S)[a]", rep("Outbreak~~year", n_ob_choosen)), 
           n_ob_id = c("Ensemble", "Ensemble",paste0( 1:n_ob_choosen, " outbreak(s)")))
  
  prop_sum_data = data.frame(model_prop = model_prop, Study_ID = study_index[file], country = country_names, 
                             n_ob_id = paste0(1:length(model_prop), " outbreak(s)"), par_names = "Model~~Probability")
  
  l_model_par_sum = c(l_model_par_sum, list(par_sum_data)); l_model_prop_sum = c(l_model_prop_sum, list(prop_sum_data))
}
model_par_sum = Reduce(rbind, l_model_par_sum); model_prop_sum = Reduce(rbind, l_model_prop_sum)

sum_data_plot <- ggplot(data = NULL, aes(y = Study_ID, color = n_ob_id, fill = n_ob_id))+
  geom_pointrange(data = model_par_sum, aes(x = X50., xmin = X2.5., xmax = X97.5.), position = position_dodge(0.2)) +
  geom_bar(data = model_prop_sum, aes(x = model_prop), stat = "identity") +
  facet_grid(country ~ par_names, scales = "free", space = "free_y", labeller = labeller(par_names = label_parsed))+
  scale_color_manual(values = c("#E6D460","#99A776","#C7421A","#66101f","#1a3a3a", '#304C89')) +
  scale_fill_manual(values = c("#E6D460","#99A776","#C7421A","#66101f","#1a3a3a", '#304C89')) +
  labs(fill = " ", color = " ", x =  " ", y = "Study") +
  theme_bw() +
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(angle = -45),
        panel.spacing.y=unit(0, "lines"),
        panel.spacing.x=unit(0.4, "lines"),
        strip.text.y = element_text(angle = 0))
ggsave(filename = "CHIKV_endpoints_event/plots/Sum_par.png", sum_data_plot, width = 10, height = 7.5, bg = "transparent")
