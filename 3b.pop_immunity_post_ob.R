#OUtbreak in 2020:
library(dplyr)
library(stringr)
# clear memory
rm(list=ls())

Rinit_files_name = list.files("results/real_sero_study/immunity/immune_by_age_in_pop/")

source("util_functions.R")
#gen_immune_func: Function to generate immune proportion from estimated parameters
#calcRinf: Function to calculate R infinity from R0 and R initial
#LL_mcmc_sero_function: function to calculate LL

year_proj = 2020
for(file in 1:length(Rinit_files_name)){
  file_names = str_sub(Rinit_files_name[file], 1, -5)
  print(file_names)
  #Get projected data: 
  projected_immunity = readRDS(paste0("results/real_sero_study/immunity/immune_projected/", file_names, ".rds"))
  pop2020 = projected_immunity$pop2020
  model_prop = projected_immunity$model_prop
  
  #One outbreak in 2020:
  immunity_post_ob = c()
  for(i in 1:length(projected_immunity$R0)){
    immunity_post_ob = c(immunity_post_ob, list(
      sapply(1:length(projected_immunity$R0[[i]]), function(x){
        projected_immunity_in_pop = projected_immunity$projected_immunity_in_pop[[i]][x,]
        pos_every_age = projected_immunity_in_pop/pop2020
        R0 = projected_immunity$R0[[i]][x]
        Rinit = sum(projected_immunity_in_pop)
        Rinf = calcRinf(R0, Rinit) #Recovery cov after outbreak
        pos_every_age_post_ob = 1 - (1 - Rinf)*(1 - pos_every_age)/(1 - Rinit)#(1 - pos_every_age)*(Rinf - Rinit) + pos_every_age
        return(pos_every_age_post_ob*pop2020)
      }) %>% t %>% data.frame()))
  }
  en_immunity_post_ob = lapply(1:length(immunity_post_ob), function(x) immunity_post_ob[[x]]*model_prop[x]) %>% Reduce("+", .)
  print(summary(rowSums(en_immunity_post_ob)))
  #Save file:
  saveRDS(list(immunity_post_ob = immunity_post_ob, pop2020 = pop2020), 
          paste0("results/real_sero_study/immunity/immune_post_ob/", file_names, ".rds"))
}

#Reported AR, by including symptomatic rate:
Symp_rate = readRDS("data/symp_rate_aggregate.rds")
for(file in 1:length(Rinit_files_name)){
  file_names = str_sub(Rinit_files_name[file], 1, -5)
  print(file_names)
  pre_ob_pop_immune = readRDS(paste0("results/real_sero_study/immunity/immune_projected/", file_names, ".rds"))
  post_ob_pop_immune = readRDS(paste0("results/real_sero_study/immunity/immune_post_ob/", file_names, ".rds"))
  endpoint  = list(endpoints = lapply(1:length(pre_ob_pop_immune$model_prop), function(x){
    model_prop = pre_ob_pop_immune$model_prop[x]
    ep = (post_ob_pop_immune$immunity_post_ob[[x]] - pre_ob_pop_immune$projected_immunity_in_pop[[x]])*model_prop*
      matrix(Symp_rate, nrow = 1000, ncol = 100, byrow = F)
    return(ep)
  }) %>% Reduce("+",.),
   pop2020 =  pre_ob_pop_immune$pop2020, R0 = pre_ob_pop_immune$R0, model_prop = pre_ob_pop_immune$model_prop
                   )
  #Save file:
  saveRDS(endpoint, paste0("results/real_sero_study/number_endpoints/" , file_names, ".rds"))
}
    