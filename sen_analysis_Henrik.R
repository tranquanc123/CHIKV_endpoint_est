#Simulate Henrik study in 1973 -> project it to 2012, 2013
library(dplyr)
library(data.table)
library(BayesianTools)
library(rjmcmc)
library(Hmisc)

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

#Summarise funtion:
quantile_func = function(x) c(quantile(x, c(0.025, 0.5, 0.975), na.rm = T))

#################################
###MCMC
#################################
Study_ID = "Henrik"
#select study:
study_data_sel = filter(CHIKV_sero, Study_id == Study_ID, total != 0)
study_data = study_data_sel[1:6,]
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

##########################
###RJMCMC
##########################

###Log-likelihood of 5 models
L1 = function(theta) {
  R0 = theta[1]; Rinit = theta[2];# k = theta[3]; 
  year_epi_start = theta[3]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}
L2 = function(theta) {
  R0 = theta[1]; Rinit = theta[2]; #k = theta[3]; 
  year_epi_start =  theta[3:4]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}
L3 = function(theta) {
  R0 = theta[1]; Rinit = theta[2]; #k = theta[3]; 
  year_epi_start =  theta[3:5]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}
L4 = function(theta) {
  R0 = theta[1]; Rinit = theta[2]; #k = theta[3];
  year_epi_start =  theta[3:6]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}
L5 = function(theta) {
  R0 = theta[1]; Rinit = theta[2];# k = theta[3]; 
  year_epi_start =  theta[3:7]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}

###Prior of 5 models:
prior1 = function(theta){
  dnorm(theta[2],  0, 10, log = T) +
  sum(dunif(theta[3],  0,  98, log = T))
}

prior2 = function(theta){
  dnorm(theta[2],  0, 10, log = T) +
  sum(dunif(theta[3:4],  0,  98, log = T)) 
}

prior3 = function(theta){
  dnorm(theta[2],  0, 10, log = T) +
  sum(dunif(theta[3:5],  0,  98, log = T)) 
}

prior4 = function(theta){
  dnorm(theta[2],  0, 10, log = T) +
  sum(dunif(theta[3:6],  0,  98, log = T)) 
}

prior5 = function(theta){
  dnorm(theta[2],  0, 10, log = T) +
  sum(dunif(theta[3:7],  0,  98, log = T)) 
}

###g and ginv
g1 = g2 = g3 = g4 = g5 = function(psi) {theta = psi}
ginv1 = ginv2 = ginv3 = ginv4 = ginv5 = function(theta) {psi = theta}

###posterior draw from 5 models:
###model runs:
n_ob_year = 5
l_n_outbreak = 1:n_ob_year
l_output = c()
for(n_outbreak in l_n_outbreak){
  par_prior_dens = function(par){
      sum(dunif(par[2:(length(par) - 1)],  0,  98, log = T)) +
      dunif(tail(par, 1), log = T)
  }
  init_values = function(n = 1){
    return(cbind(runif(n, 1, 10) ,
                 matrix(runif(n*n_outbreak,  0,  98), ncol = n_outbreak, nrow = n) ,
                 runif(n))
    )
  }
  par_prior = createPrior(density = par_prior_dens, sampler = init_values, 
                          lower = c(1, rep(0, n_outbreak), 0), upper = c(10, rep(98, n_outbreak), 1),
                          best = c(5, rep(50, n_outbreak), 0.5))
  
  bayesanSetup <- createBayesianSetup(likelihood = LL_mcmc_sero_function, prior = par_prior)
  iter = 1e5
  settings <- list(iterations = iter)
  a <- runMCMC(bayesianSetup = bayesanSetup, sampler = "DREAMzs", settings = settings)
  l_output = c(l_output, list(getSample(a, thin = 10)))
}

#reorder par:
for(n_outbreak in 1:n_ob_year){
  l_output[[n_outbreak]] <- l_output[[n_outbreak]][,c(1, n_outbreak + 2, 2:(n_outbreak + 1))]
}
### attach posterior for auxiliary var to matrices of posterior draw:
lc = nrow(l_output[[1]])
l_output_att = l_output
for(n_outbreak in 1:n_ob_year){
  u_i = matrix(runif(lc*(n_ob_year - n_outbreak), 0, 98), lc, n_ob_year - n_outbreak, byrow = T)
  l_output_att[[n_outbreak]] = cbind(l_output_att[[n_outbreak]], u_i)
}

###function to randomly sample from posterior
getsampler(l_output_att[[1]], "draw1"); getsampler(l_output_att[[2]], "draw2"); getsampler(l_output_att[[3]], "draw3"); 
getsampler(l_output_att[[4]], "draw4"); getsampler(l_output_att[[5]], "draw5")

model_list = 1:n_ob_year
seeds = rjmcmcpost(post.draw = mget(paste0("draw", model_list)),
                   likelihood = mget(paste0("L", model_list)),
                   g = mget(paste0("g", model_list)), ginv = mget(paste0("ginv", model_list)),
                   param.prior = mget(paste0("prior", model_list)),
                   model.prior = rep(1/n_ob_year, n_ob_year),
                   chainlength = 1e4)

saveRDS(list(seeds = seeds, l_output = l_output), paste0("results/case_study/", Study_ID, "_1973sim.rds"))

####using this data to get datafit of 2012 and 2013:
Henrik_1973_par = readRDS(paste0("results/case_study/", Study_ID, "_1973sim.rds"))
Henrik_all_par = readRDS(paste0("results/real_sero_study/RJMCMC/Henrik.rds"))

#Simulation function:
sim_function = function(par, study_data, data_fit_time){
  R0 = par[1]; year_epi_start = par[2:(length(par) - 1)]; Rinit = tail(par, 1)
  study_time =  tail(unique(study_data$time), 1)
  
  pos_every_age <- gen_immune_func(R0 = R0, year_epi_start = year_epi_start, Rinit = Rinit, 
    study_data = study_data, l_study_time = study_time, country_prop_pop = country_prop_pop) %>% unlist
  
  #Obtain study-specific age group proportion:
  study_data_fit =  filter(study_data, time == data_fit_time)
  age_l = study_data_fit$age_l + 1
  age_u = study_data_fit$age_u + 1
  pos_prop_study <- sapply(1:length(age_l), function(x) sum(pos_every_age[age_l[x]:age_u[x]]*country_prop_pop[age_l[x]:age_u[x],paste0("X", study_time)])/
                             sum(country_prop_pop[age_l[x]:age_u[x],paste0("X", study_time)]))
  pos_prop_study[pos_prop_study == 1] = 1 - 99e-15
  pos_prop_study[pos_prop_study == 0] = 1e-99
  #derive the true pos:
  ob_pos_prop <- pos_prop_study
  
  return(ob_pos_prop)
}

year_study_fit = study_data_sel$time %>% unique
data_fit_1973_year_sum = c()
for(i in 1:length(year_study_fit)){
  study_year_data = filter(study_data_sel, time %in% year_study_fit[1:i])
  ob_year_1973 = Henrik_1973_par$l_output[[2]][,3:ncol(Henrik_1973_par$l_output[[2]])] %>% 
      tail(1000) %>% matrix(nrow = 1000) %>% apply(1, sort) %>% t
  ob_year_1973_proj = year_study_fit[i] - 1973 + ob_year_1973
  ob_year_2013 = Henrik_all_par$l_output[[5]][,3:ncol(Henrik_all_par$l_output[[5]])] %>% 
      tail(1000) %>% matrix(nrow = 1000) %>% apply(1, sort) %>%  t
  ob_year_2013_proj = year_study_fit[i] - 2013 + ob_year_2013
  
  if(year_study_fit[i] == 1973){
    ob_year_i = ob_year_1973_proj
  }
  if(year_study_fit[i] == 2012){
    ob_year_i = cbind(ob_year_2013_proj[,2:3], ob_year_1973_proj)
  }
  if(year_study_fit[i] == 2013){
    ob_year_i = cbind(ob_year_2013_proj[,1:3], ob_year_1973_proj)
  }
  
  #RO and Rinit from 1973:
  R0_Rinit = Henrik_1973_par$l_output[[2]][,1:2] %>% tail(1000)
  data_fit_i = sapply(1:nrow(R0_Rinit), function(x) sim_function(c(R0_Rinit[x,1], ob_year_i[x,], R0_Rinit[x, 2]), study_year_data, year_study_fit[i])) %>% t
  data_fit_1973_year_sum = c(data_fit_1973_year_sum, list(apply(data_fit_i, 2, quantile_func) %>% t %>% data.frame()))
}

data_fit_all_year_sum = c()
for(i in 1:length(year_study_fit)){
  study_year_data = filter(study_data_sel, time %in% year_study_fit[1:i])
  ob_year_2013 = Henrik_all_par$l_output[[5]][,3:ncol(Henrik_all_par$l_output[[5]])] %>% 
    tail(1000) %>% matrix(nrow = 1000) %>% apply(1, sort) %>%  t
  ob_year_2013_proj = year_study_fit[i] - 2013 + ob_year_2013
  
  if(year_study_fit[i] == 1973){
    ob_year_i = ob_year_2013_proj[,4:5]
  }
  if(year_study_fit[i] == 2012){
    ob_year_i = ob_year_2013_proj[,2:5]
  }
  if(year_study_fit[i] == 2013){
    ob_year_i = ob_year_2013_proj[,1:3]
  }
  
  #RO and Rinit from all years:
  R0_Rinit = Henrik_all_par$l_output[[5]][,1:2] %>% tail(1000)
  data_fit_i = sapply(1:nrow(R0_Rinit), function(x) sim_function(c(R0_Rinit[x,1], ob_year_i[x,], R0_Rinit[x, 2]), study_year_data, year_study_fit[i])) %>% t
  data_fit_all_year_sum = c(data_fit_all_year_sum, list(apply(data_fit_i, 2, quantile_func) %>% t %>% data.frame()))
}


plot_fit = study_data_sel %>% select(time, age_l, age_u, age_mid, total, pos, neg) %>% 
  mutate(age_mid = (age_l + age_u)/2) %>%
  mutate(age_l = formatC(age_l, width = 2, flag = 0)) %>%
  mutate(age_u = formatC(age_u, width = 2, flag = 0)) %>%
  mutate(x_label = paste0(age_l, "-", age_u)) %>% 
  cbind(., binconf(.$pos, .$total)) %>% 
  list %>% rep(2) %>% Reduce(rbind, .) %>%
  cbind(Reduce(rbind, c(data_fit_1973_year_sum, data_fit_all_year_sum))) %>% 
  mutate(Analysis = rep(c("1973 data only", "All years data"), each = nrow(study_data_sel))) %>%
  data.frame()

dat_text <- data.frame(
  label = c("A", "B", "C", "D", "E", "F"),
  time = rep(c(1973, 2012, 2013), 2),
  Analysis = rep(c("1973 data only", "All years data"), each = 3),
  x_label = 1, PointEst = 0.90
)

fit_plot <-
  ggplot(data = plot_fit, aes(x = x_label, y = PointEst, group = 1))+
  geom_pointrange(aes(ymin = Lower, ymax = Upper)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#1a3a3a", alpha = 0.4) +
  geom_line(aes(y = X50.), color = "#C7421A") +
  geom_text(data = dat_text, aes(x = x_label, y = PointEst, label = label), size = 8)+
  facet_grid(Analysis~time) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Age Group", y = "Immune proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        rect = element_rect(fill = "transparent",colour = NA),
        text = element_text(size = 18))

ggsave(filename = "plots/Sen_Henrik_datafit.png", fit_plot, width = 10, height = 7.5, bg = "transparent")

###comparing AR in 2020 when we only have the 1973 data and hold out 2012 and 2013:
Henrik_1973_par = readRDS(paste0("results/case_study/", Study_ID, "_1973sim.rds"))
Henrik_all_par = readRDS(paste0("results/real_sero_study/RJMCMC/Henrik.rds"))
l_study_time = c(1973, 2012, 2013)

#get immune profile in 2013 from estimates of 1973:
ob_year_1973 = Henrik_1973_par$l_output[[2]][,3:ncol(Henrik_1973_par$l_output[[2]])] %>% 
  tail(1000) %>% matrix(nrow = 1000)
ob_year_1973_proj = 2013 - 1973 + ob_year_1973
ob_year_2013 = Henrik_all_par$l_output[[5]][,3:ncol(Henrik_all_par$l_output[[5]])] %>% 
  tail(1000) %>% matrix(nrow = 1000) %>% apply(1, sort) %>%  t
ob_year_2013_proj = 2013 - 2013 + ob_year_2013

ob_year_i = cbind(ob_year_2013_proj[,1:3], ob_year_1973_proj)

#RO and Rinit from 1973:
R0_Rinit = Henrik_1973_par$l_output[[2]][,1:2] %>% tail(1000)

#get pop prop from  2013 all years:
para_sim_all_year = Henrik_all_par$l_output[[5]] %>% tail(1000)

#getting the AR of subsequence outbreaks after 1973 based on AR of 1973 data alone:
par_1973_2013 = cbind(R0_Rinit, ob_year_2013[,1:3], ob_year_1973_proj)
para_sim_all_year
l_sim_run = list(par_1973_2013, para_sim_all_year)

AR_l = c()
for(sim_run in 1:length(l_sim_run)){
  
  #load data:
  par_sim = l_sim_run[[sim_run]]
  
  AR_mat = sapply(1:nrow(par_sim), function(x){
    par = par_sim[x,]
    R0 = par[1]; Rinit = par[2]; year_epi_start = par[3:length(par)]
    #Logit Transform Rinit:
    t_Rinit <- Rinit#exp(Rinit)/(exp(Rinit) + 1)
    
    #which ob in which study:
    year_epi_start = sort(year_epi_start)
    study_year_dist = rev(tail(l_study_time, 1) - rev(l_study_time))
    l_ob_in_study = c()
    for(y in study_year_dist){
      ob_dist = year_epi_start - y; ob_dist = ob_dist[ob_dist > 0]
      l_ob_in_study = c(l_ob_in_study, list(ob_dist))
      year_epi_start = head(year_epi_start, -length(ob_dist))
    }
    
    pos_every_age = rep(0, 100)
    AR_5 = c()
    #for each study time:
    for(time in 1:length(l_study_time)){
      study_time = l_study_time[time]
      year_epi_start_study = sort(l_ob_in_study[[time]], decreasing = T)
      year_dist = c(head(year_epi_start_study, -1) - year_epi_start_study[-1], tail(year_epi_start_study, 1))
      year_floor = floor(year_epi_start_study)
      
      #Since S(0) is the susceptible pop before the first study, the S(0) of the first ob of the first study is
      #(assuming susceptible proportion distributed evenly between age groups):
      if(time == 1){
        first_epi_year_whole = year_floor[1]
        first_epi_year_part = year_epi_start_study[1] - first_epi_year_whole
        if(first_epi_year_whole == 0){
          t_Rinit = (1 - first_epi_year_part)*sum(t_Rinit*country_prop_pop[100:(100 - first_epi_year_whole),paste0("X", ifelse(study_time - first_epi_year_whole < 1950, 1950, study_time - first_epi_year_whole))])
        } else {
          t_Rinit = sum(t_Rinit*country_prop_pop[100:(100 - first_epi_year_whole), paste0("X", ifelse(study_time - first_epi_year_whole - 1 < 1950, 1950, study_time - first_epi_year_whole - 1))]) +
            (1 - first_epi_year_part)*sum(t_Rinit*country_prop_pop[100:(100 - first_epi_year_whole),paste0("X", ifelse(study_time - first_epi_year_whole < 1950, 1950, study_time - first_epi_year_whole))])
        }
      }
      for(i in 1:length(year_dist)){
        Rinf = calcRinf(R0, t_Rinit) #calcRinf_large_init(R0, t_Rinit, k)##Recovery cov after outbreak
        AR_5 = c(AR_5, Rinf - t_Rinit)
        pos_every_age = 1 - (1 - Rinf)*(1 - pos_every_age)/(1 - t_Rinit)#(1 - pos_every_age)*(Rinf - t_Rinit) + pos_every_age
        #project to the current study year:
        whole_year = floor(year_dist)[i]; part_year = year_epi_start_study[i] - floor(year_epi_start_study[i]) 
        pos_every_age <- c(numeric(length = whole_year), (1 - part_year)*pos_every_age[1], pos_every_age[c(1:(99 - whole_year))])
        pos_prop_pop = pos_every_age*country_prop_pop[,paste0("X", ifelse(study_time - year_floor[i] - 1 < 1950, 1950, study_time - year_floor[i] - 1))]
        t_Rinit = sum(pos_prop_pop)
      }
    }
    return(AR_5)
  }) %>% t
  AR_l = c(AR_l, list(AR_mat))
}

#AR in 2020
dataset_names = mapply(rep, c("1973 data only", "1973 data with \n outbreak time estimated from \n all years data", "All years data"), times = c(2, 3, 5)) %>% unlist
AR_plot_data = lapply(1:length(AR_l), function(x){
  mod_data = AR_l[[x]] %>% apply(2, quantile, prob = c(0.025, 0.25, 0.5, 0.75, 0.975)) %>% round(6) %>% t %>% data.frame 
  mod_data$Outbreak = paste0(1:5, " outbreak")
  return(mod_data)
}) %>% Reduce(rbind, .)
AR_plot_data$Dataset = dataset_names

AR_comp_plot  <- ggplot(AR_plot_data , aes(x = Outbreak, fill = Dataset))+
  geom_boxplot(aes(ymin = X2.5., lower = X25., middle = X50., upper = X75., ymax = X97.5.),
               stat = "identity", color = "grey10")+
  scale_fill_manual(values = c("#E6D460","#99A776","#C7421A")) +
  theme_bw()+
  labs(y = "IAR")

ggsave(filename = "plots/Sen_Henrik_AR.png", AR_comp_plot, width = 7.5, height = 5, bg = "transparent")
