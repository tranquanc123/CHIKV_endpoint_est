library(dplyr)
library(data.table)
library(BayesianTools)

#running parallel in crc, which.scenario is the index of 35 study sites:
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

#All sites and number of outbreak:  
sites_obs = expand.grid(sites = all_studies_id, n_ob = 1:5)
  
source("util_functions.R")
#gen_immune_func: Function to generate immune proportion from estimated parameters
#calcRinf: Function to calculate R infinity from R0 and R initial
#LL_mcmc_sero_function: function to calculate LL

#################################
###MCMC run
#################################

n_outbreak = sites_obs$n_ob[which.scenario]
Study_ID = sites_obs$sites[which.scenario] %>% as.character
#select study:
study_data = filter(CHIKV_sero, Study_id == Study_ID, total != 0)
    
#Data:
#Getting the demographic data:
country_sel = unique(study_data$country) %>% as.character
country_prop_pop <- filter(pop_data, country == country_sel) %>% 
  select(starts_with("X")) %>%
  apply(1, function(x) x/sum(x)) %>% data.frame
colnames(country_prop_pop) = paste0("X", unique(pop_data$year))

#Study time:  
l_study_time = study_data$time %>% unique
study_data_id <- study_data$time %>% table 
l_study_data_id <- list(1:study_data_id[1])
if(length(l_study_time) >= 2 ){
  for(i in 2:length(l_study_time)){
    l_study_data_id <- c(l_study_data_id, list((cumsum(study_data_id)[i - 1] + 1):cumsum(study_data_id)[i]))
  }
}
names(l_study_data_id) = l_study_time
l_study_time = sort(l_study_time)

###model runs:
y_id_max = tail(l_study_time, 1) - l_study_time[1] + 99

par_prior_dens = function(par){
  return(dunif(par[1], 1, 10, log = T)+ 
           sum(dunif(par[2:(length(par) - 1)],  0,  98, log = T)) +
           dunif(tail(par, 1), log = T)
  )
}

init_values = function(n = 1){
  return(cbind(runif(n, 1, 10) ,
               matrix(runif(n*n_outbreak,  0,  98), ncol = n_outbreak, nrow = n) ,
               runif(n, 0, 1)
  )
  )
}

par_prior = createPrior(density = par_prior_dens, sampler = init_values, 
                        lower = c(1, rep(0, n_outbreak), -10), upper = c(10, rep(98, n_outbreak), 10),
                        best = c(5, rep(50, n_outbreak), 0))

bayesanSetup <- createBayesianSetup(likelihood = LL_mcmc_sero_function, prior = par_prior)
iter = 1e6
settings <- list(iterations = iter)
output <- runMCMC(bayesianSetup = bayesanSetup, sampler = "DREAMzs", settings = settings)

saveRDS(output, paste0("results/real_sero_study/MCMC/", Study_ID, "_", n_outbreak, "ob.rds"))
