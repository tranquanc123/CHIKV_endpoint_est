#Project to 2020
# clear memory
rm(list=ls())
library(dplyr)
library(stringr)

#loading data:
Rinit_files_name = list.files("results/real_sero_study/immunity/immune_by_age_in_pop/")
#Pop data:
pop_data <- read.csv("data/pop_data.csv")
pop_data <- pop_data %>%
  mutate(X99 = X99 + X100) %>%
  select(-X100)
#sero data:
CHIKV_sero <- read.csv("data/CHIKV_sero_data.csv")

year_proj = 2020

for(file in 1:length(Rinit_files_name)){
  file_names = str_sub(Rinit_files_name[file], 1, -5)
  country_names = filter(CHIKV_sero, Study_id == str_sub(Rinit_files_name[file], 1, -5))$country %>% unique %>% as.character()
  print(file_names)
  
  #get sim results file => getting R0 and Rinit:
  param.list = readRDS(paste0("results/real_sero_study/parameter_samples/", file_names, ".rds"))
  R0_l = param.list
  R0_l %>% lapply(tail, n = 1000) %>% sapply(quantile, probs = c(0.025, 0.5, 0.975)) %>% print
  
  #get the Rinit value:
  study_inmunity_info <- readRDS(paste0("results/real_sero_study/immunity/immune_by_age_in_pop/", file_names, ".rds"))
  #project to 2020:
  popstudy_matrix = matrix(study_inmunity_info$country_pop, byrow = T, ncol = 100, nrow = 1000)
  study_inmunity_level = lapply(1:length(study_inmunity_info$gen_pos_prop), function(x) t(study_inmunity_info$gen_pos_prop[[x]])/popstudy_matrix) #inmune prob within age
  study_year = study_inmunity_info$study_time
  project_immunity_level = lapply(study_inmunity_level, function(x) cbind(matrix(0, nrow = 1000, ncol = year_proj - study_year), x[, (year_proj - study_year + 1):100]))
  pop2020 = filter(pop_data, country == country_names, year == 2020) %>% select(starts_with("X")) %>% as.numeric()
  pop2020 = pop2020/sum(pop2020)
  project_immunity_in_pop = lapply(project_immunity_level, function(x) x*matrix(pop2020, byrow = T, ncol = 100, nrow = 1000))
    
  #save files:
  saveRDS(list(projected_immunity_in_pop = project_immunity_in_pop, pop2020 = pop2020, R0 = R0_l, model_prop = study_inmunity_info$model_prop), 
          paste0("results/real_sero_study/immunity/immune_projected/", file_names, ".rds"))
}
