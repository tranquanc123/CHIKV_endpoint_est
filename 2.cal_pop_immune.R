#current immunnity at the study year.
# clear memory
rm(list=ls())

# loading library:
library(dplyr)
library(stringr)
library(deSolve)
library(data.table)
library(ggplot2)
library(BayesianTools)

#Loading data:
#sero data:
CHIKV_sero_data <- read.csv("data/CHIKV_sero_data.csv", 1)
  
#Pop data:
pop_data <- read.csv("data/pop_data.csv")
pop_data <- pop_data %>%
  mutate(X99 = X99 + X100) %>%
  select(-X100)

source("util_functions.R")
#gen_immune_func: Function to generate immune proportion from estimated parameters
#calcRinf: Function to calculate R infinity from R0 and R initial
#LL_mcmc_sero_function: function to calculate LL

  prob_gen_func <-function(par, study_data){
    R0 = par[1]; Rinit = par[2]; year_epi_start = par[3:length(par)]
    year_epi_start = year_epi_start[year_epi_start > 0]
    
    #Logit Transform Rinit:
    t_Rinit <- expit(Rinit)
    study_time = unique(study_data$time)
    year_epi_start_study = sort(year_epi_start, decreasing = T)
    #if(length(year_dist) != 0) {
    pos_every_age = numeric(100)
    year_dist = c(head(year_epi_start_study, -1) - year_epi_start_study[-1], tail(year_epi_start_study, 1))
    
    #Since S(0) is the susceptible pop before the first study, the S(0) of the first ob of the first study is
    #(assuming susceptible proportion distributed evenly between age groups):
    if(y_time == min(l_study_time)){
      first_epi_year_whole = floor(year_dist[1])
      first_epi_year_part = year_dist[1] - first_epi_year_whole
      if(first_epi_year_whole == 0){
        t_Rinit = (1 - first_epi_year_part)*sum(t_Rinit*country_prop_pop[100:(100 - first_epi_year_whole),paste0("X", ifelse(study_time - first_epi_year_whole < 1950, 1950, study_time - first_epi_year_whole))])
      } else {
        t_Rinit = sum(t_Rinit*country_prop_pop[100:(100 - first_epi_year_whole),paste0("X", ifelse(study_time - first_epi_year_whole - 1 < 1950, 1950, study_time - first_epi_year_whole - 1))]) +
          (1 - first_epi_year_part)*sum(t_Rinit*country_prop_pop[100:(100 - first_epi_year_whole),paste0("X", ifelse(study_time - first_epi_year_whole < 1950, 1950, study_time - first_epi_year_whole))])
      }
    }
    
    for(i in 1:length(year_dist)){
      Rinf = calcRinf(R0, t_Rinit) #Recovery cov after outbreak
      pos_every_age = 1 - (1 - Rinf)*(1 - pos_every_age)/(1 - t_Rinit)#(1 - pos_every_age)*(Rinf - t_Rinit) + pos_every_age
      #project
      whole_year = floor(year_dist)[i]; part_year = year_epi_start_study[i] - floor(year_epi_start_study[i]) 
      pos_every_age <- c(numeric(length = whole_year), (1 - part_year)*pos_every_age[1], pos_every_age[c(1:(99 - whole_year))])
      pos_prop_pop <- pos_every_age*country_prop_pop[,paste0("X", ifelse(study_time - whole_year - 1 < 1950, 1950, study_time - whole_year - 1))]
      t_Rinit <- sum(pos_prop_pop)
    }
    
    return(pos_prop_pop)
  }

all_studies_id = str_sub(list.files("results/real_sero_study/RJMCMC/"), end = -5)
for(which.scenario in all_studies_id){
  
  tryCatch({
    study_data = filter(CHIKV_sero_data, Study_id == which.scenario, total != 0)
    l_study_time = study_data$time %>% unique
    study_data_id <- study_data$time %>% table 
    l_study_data_id <- list(1:study_data_id[1])
    if(length(l_study_time) >= 2 ){
      for(i in 2:length(l_study_time)){
        l_study_data_id <- c(l_study_data_id, list((cumsum(study_data_id)[i - 1] + 1):cumsum(study_data_id)[i]))
      }
    }
    names(l_study_data_id) = l_study_time
    #Data:
    #Getting the demographic data:
    country_sel = unique(study_data$country) %>% as.character
    country_prop_pop <- filter(pop_data, country == country_sel) %>% 
      select(starts_with("X")) %>%
      apply(1, function(x) x/sum(x)) %>% data.frame
    colnames(country_prop_pop) = paste0("X", unique(pop_data$year))
    
    seeds_output = readRDS(paste0("results/real_sero_study/RJMCMC/", which.scenario, ".rds"))
    model_prop = as.numeric(seeds_output$seeds$result$`Posterior Model Probabilities`)
    sim_sample = lapply(seeds_output$l_output, tail, n = 1000)
    
    study_year = unique(study_data$time)
    study_year_dist = max(study_year) - study_year
    for(y_time in study_year){
      sel_study_data <- filter(study_data, time == y_time)
      mod_sim_sample <- 
        lapply(sim_sample, function(x){
          year_est = x[,3:ncol(x)]
          mod_year = year_est - (max(study_year) - y_time)
          #mod_year = mod_year[mod_year > 0]
          mod_year_data = x
          mod_year_data[, 3:ncol(mod_year_data)] = mod_year
          return(mod_year_data)
        })
      gen_pos_prop = lapply(mod_sim_sample, function(y) sapply(1:nrow(y), function(x) {
        par = as.numeric(y[x,]); 
        R0 = par[1]; year_epi_start = par[2:(length(par) - 1)]; Rinit = tail(par, 1)
        prob_gen = gen_immune_func(R0, year_epi_start, Rinit, sel_study_data, y_time, country_prop_pop)
        return(prob_gen)
        }))
      en_gen_pos_prop = lapply(1:length(model_prop), function(x) gen_pos_prop[[x]]*model_prop[x])
      en_immune_by_age_in_pop = t(Reduce("+", en_gen_pos_prop))
      en_immune_by_age = rowSums(en_immune_by_age_in_pop)
  }
  #write data:
  print(which.scenario)
  print(all(en_immune_by_age_in_pop >= 0))
  study_time = max(l_study_time)
  write.csv(en_immune_by_age, paste0("results/real_sero_study/immunity/immune_by_age/", model_chosen, which.scenario, ".csv"), row.names=F, quote=F)
  saveRDS(list(gen_pos_prop = gen_pos_prop, country_pop = country_prop_pop[,paste0("X", study_time)], study_time = study_time, model_prop = model_prop), 
          paste0("results/real_sero_study/immunity/immune_by_age_in_pop/", model_chosen, which.scenario, ".rds"))
    saveRDS(sim_sample, paste0("results/real_sero_study/immunity/parameter_samples/", model_chosen, which.scenario, ".rds"))
  
  }, error = function(e) print("Error"))
}