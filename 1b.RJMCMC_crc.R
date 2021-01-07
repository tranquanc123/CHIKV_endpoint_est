library(dplyr)
library(rjmcmc)
library(BayesianTools)

#Run the analysis parallel in crc
args = commandArgs(trailingOnly=TRUE)
which.scenario = as.numeric(args[1])

#Pop data:    
pop_data <- read.csv("data/pop_data.csv")
pop_data <- pop_data %>%
  mutate(X99 = X99 + X100) %>%
  select(-X100)
#sero data:
CHIKV_sero <- read.csv("data/CHIKV_sero_data.csv")
all_studies_id = unique(CHIKV_sero$Study_id)[-c(19)]
Study_ID = as.character(all_studies_id[which.scenario])

source("util_functions.R")
#gen_immune_func: Function to generate immune proportion from estimated parameters
#calcRinf: Function to calculate R infinity from R0 and R initial
#LL_mcmc_sero_function: function to calculate LL

#################################
###MCMC
#################################
#select study:
study_data_raw = filter(CHIKV_sero, Study_id == Study_ID, total != 0)

#Data:
#Getting the demographic data:
country_sel = unique(study_data_raw$country) %>% as.character
country_prop_pop <- filter(pop_data, country == country_sel) %>% 
  select(starts_with("X")) %>%
  apply(1, function(x) x/sum(x)) %>% data.frame
colnames(country_prop_pop) = paste0("X", unique(pop_data$year))

#Study time:  
l_study_time = study_data_raw$time %>% unique
study_data_id <- study_data_raw$time %>% table 
l_study_data_id <- list(1:study_data_id[1])
if(length(l_study_time) >= 2 ){
  for(i in 2:length(l_study_time)){
    l_study_data_id <- c(l_study_data_id, list((cumsum(study_data_id)[i - 1] + 1):cumsum(study_data_id)[i]))
  }
}
names(l_study_data_id) = l_study_time
study_data = study_data_raw

##########################
###RJMCMC
##########################
y_id_max = tail(l_study_time, 1) - l_study_time[1]

###Log-likelihood of 5 models
L1 = function(theta) {
  R0 = theta[1]; Rinit = theta[2]; year_epi_start = theta[3] 
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}
L2 = function(theta) {
  R0 = theta[1]; Rinit = theta[2]; year_epi_start = theta[3:4]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}
L3 = function(theta) {
  R0 = theta[1]; Rinit = theta[2]; year_epi_start = theta[3:5]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}
L4 = function(theta) {
  R0 = theta[1]; Rinit = theta[2]; year_epi_start = theta[3:6]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}
L5 = function(theta) {
  R0 = theta[1]; Rinit = theta[2]; year_epi_start = theta[3:7]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}

###Prior of 5 models:
prior1 = function(theta){
  dnorm(theta[2],  0, 10, log = T) +
  sum(dunif(theta[3],  0,  99 + y_id_max, log = T))
}

prior2 = function(theta){
  dnorm(theta[2],  0, 10, log = T) +
  sum(dunif(theta[3:4],  0,  99 + y_id_max, log = T))
}

prior3 = function(theta){
  dnorm(theta[2],  0, 10, log = T) +
  sum(dunif(theta[3:5],  0,  99 + y_id_max, log = T))
}

prior4 = function(theta){
   dnorm(theta[2],  0, 10, log = T) +
  sum(dunif(theta[3:6],  0,  99 + y_id_max, log = T))
}

prior5 = function(theta){
  dnorm(theta[2],  0, 10, log = T) +
  sum(dunif(theta[3:7],  0,  99 + y_id_max, log = T)) 
}

###g and ginv
g1 = g2 = g3 = g4 = g5 = function(psi) {psi}
ginv1 = ginv2 = ginv3 = ginv4 = ginv5 = function(theta) {theta}

###posterior draw from 5 models:
###model runs:
l_n_outbreak = 1:5
l_output = lapply(l_n_outbreak, function(x) {
  one_output = readRDS(paste0("results/real_sero_study/MCMC/", Study_ID, "_", x, "ob.rds"))
  sample_output = getSample(one_output, thin = 100)
  return(sample_output)
})

#reorder par:
for(n_outbreak in l_n_outbreak){
  l_output[[n_outbreak]] <- l_output[[n_outbreak]][,c(1, n_outbreak + 2, 2:(n_outbreak + 1))]
}
### attach posterior for auxiliary var to matrices of posterior draw:
lc = nrow(l_output[[1]])
l_output_att = l_output
for(n_outbreak in l_n_outbreak){
  u_i = matrix(runif(lc*(5 - n_outbreak), 0, 99 + y_id_max), lc, 5 - n_outbreak, byrow = T)
  l_output_att[[n_outbreak]] = cbind(l_output_att[[n_outbreak]], u_i)
}

###function to randomly sample from posterior
getsampler(l_output_att[[1]], "draw1"); getsampler(l_output_att[[2]], "draw2"); getsampler(l_output_att[[3]], "draw3"); 
getsampler(l_output_att[[4]], "draw4"); getsampler(l_output_att[[5]], "draw5")

sel_model = 1:5
seeds = rjmcmcpost(post.draw = mget(paste0("draw", sel_model)),
                   likelihood = mget(paste0("L", sel_model)),
                   g = mget(paste0("g", sel_model)),
                   ginv = mget(paste0("ginv", sel_model)),
                   param.prior = mget(paste0("prior", sel_model)),
                   model.prior = rep(1/length(sel_model), length(sel_model)),
                   chainlength = 1e4)

saveRDS(list(seeds = seeds, l_output = l_output), paste0("results/real_sero_study/RJMCMC/", Study_ID, ".rds"))