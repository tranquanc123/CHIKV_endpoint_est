library(dplyr)
library(data.table)
library(BayesianTools)
library(reshape2)
library(rjmcmc)
library(ggplot2)
library(Hmisc)
library(stringr)

#Pop data:
pop_data <- read.csv("data/pop_data.csv")
pop_data <- pop_data %>%
  mutate(X99 = X99 + X100) %>%
  select(-X100)

#sero data:
CHIKV_sero <- read.csv("data/CHIKV_sero_data.csv")
Study_ID = unique(CHIKV_sero$Study_id)[17] %>% as.character

source("util_functions.R")
#gen_immune_func: Function to generate immune proportion from estimated parameters
#calcRinf: Function to calculate R infinity from R0 and R initial
#LL_mcmc_sero_function: function to calculate LL

#################################
###Simulation data
#################################
l_scen = c("original", "larger_agegroup", "smaller_samplesize") #Choode which scenario, the last 2 are for the sensitivity analyses;
scen = l_scen[1]

R0 = 1.5; first_Rinit = 0
ob_year = c(85, 45, 5) 
n_ob_year = 5

#select study:
study_data_raw = filter(CHIKV_sero, Study_id == Study_ID, total != 0)
study_data_raw <- Reduce(rbind, rep(list(study_data_raw), 
                                    ifelse(scen %in% c("original", "smaller_samplesize"), 4, 1))) #<=== change this for dif age group

#Data:
#Getting the demographic data:
country_sel = unique(study_data_raw$country) %>% as.character
country_prop_pop <- filter(pop_data, country == country_sel) %>% 
  select(starts_with("X")) %>%
  apply(1, function(x) x/sum(x)) %>% data.frame
colnames(country_prop_pop) = paste0("X", unique(pop_data$year))

#Getting list of study time:
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

#Getting limit of age:
age_u_lim = max(study_data_raw$time) - min(study_data_raw$time - study_data_raw$age_u)
age_l_lim = max(study_data_raw$time) - min(study_data_raw$time)

sim_data_func <- function(ob_year, study_data, total_sample){ #<==+ change this for dif sample size
  #Modified study data
  n_outbreak = length(ob_year)
  study_data$age_l = round(seq(0, 99, 100/nrow(study_data)))
  study_data$age_u = c(tail(study_data$age_l, -1) - 1, 99)
  study_data$age_mid = (study_data$age_l + study_data$age_u )/2
  study_data$total = rep(total_sample/nrow(study_data), nrow(study_data))
  sample_every_age = rep(total_sample/100, 100)
  
  #Obtain immunity profile:
  pos_every_age <- gen_immune_func(R0, year_epi_start = ob_year, Rinit = first_Rinit, 
                                   study_data, l_study_time, country_prop_pop) %>% unlist
    
  real_pos = sapply(1:nrow(study_data), function(x) {
    one_row_data = study_data[x,]
    n_pos_every_age = pos_every_age*sample_every_age
    sum(n_pos_every_age[(one_row_data$age_l + 1):(one_row_data$age_u + 1)])
  })
  
  study_data$pos = round(real_pos)
  study_data$neg = study_data$total - study_data$pos
  
  return(study_data)
}
study_data = sim_data_func(ob_year, study_data_raw, 
                           total_sample = ifelse(scen %in% c("original", "larger_agegroup"), 2000, 200))

#################################
###MCMC
#################################

###model run:
l_n_outbreak = 1:n_ob_year
l_output = c()
for(n_outbreak in l_n_outbreak){
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
  iter = 1e5
  settings <- list(iterations = iter)
  a <- runMCMC(bayesianSetup = bayesanSetup, sampler = "DREAMzs", settings = settings)
  l_output = c(l_output, list(getSample(a, thin = 10)))
}

#reorder par:
for(n_outbreak in 1:n_ob_year){
  l_output[[n_outbreak]] <- l_output[[n_outbreak]][,c(1, n_outbreak + 2, 2:(n_outbreak + 1))]
}

##########################
###RJMCMC
##########################
###Log-likelihood of 5 models
L1 = function(theta) {
  R0 = theta[1]; Rinit = theta[2];
  year_epi_start = theta[3]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}
L2 = function(theta) {
  R0 = theta[1]; Rinit = theta[2]; 
  year_epi_start =  theta[3:4]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}
L3 = function(theta) {
  R0 = theta[1]; Rinit = theta[2]; 
  year_epi_start =  theta[3:5]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}
L4 = function(theta) {
  R0 = theta[1]; Rinit = theta[2];
  year_epi_start =  theta[3:6]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}
L5 = function(theta) {
  R0 = theta[1]; Rinit = theta[2];
  year_epi_start =  theta[3:7]
  return(LL_mcmc_sero_function(c(R0, year_epi_start, Rinit)))
}

###Prior of 5 models:
prior1 = function(theta){
  dunif(theta[2], 0, 1, log = T) +
  sum(dunif(theta[3],  0,  98, log = T))
}
prior2 = function(theta){
  dunif(theta[2], 0, 1, log = T) +  
  sum(dunif(theta[3:4],  0,  98, log = T)) 
}
prior3 = function(theta){
  dunif(theta[2], 0, 1, log = T) +
  sum(dunif(theta[3:5],  0,  98, log = T)) 
}
prior4 = function(theta){
  dnorm(theta[2],  0, 10, log = T) +
  sum(dunif(theta[3:6],  0,  98, log = T)) 
}
prior5 = function(theta){
  dunif(theta[2], 0, 1, log = T) +  
  sum(dunif(theta[3:7],  0,  98, log = T)) 
}

###g and ginv
g1 = g2 = g3 = g4 = g5 = function(psi) {theta = psi}
ginv1 = ginv2 = ginv3 = ginv4 = ginv5 = function(theta) {psi = theta}

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
                   chainlength = 1e5)

saveRDS(list(seeds = seeds, l_output = l_output, study_data  = study_data), paste0("results/case_study/sim_", scen, 
                                                         "3obs.rds"))

############################
###Results from the model
############################
###Data fit:
datafit_function <-function(par, study_data){
  par = as.numeric(par)
  R0 = par[1];  Rinit = par[2]; year_epi_start = tail(par, length(par) - 2)
  
  pos_every_age <- gen_immune_func(R0, year_epi_start, Rinit, study_data, l_study_time, country_prop_pop)
  
  ob_pos_prop = c()
  for(time in 1:length(l_study_time)){
    study_time = l_study_time[time]
    current_study_year_dist = tail(l_study_time, 1) - study_time
    study_data_sel = study_data[l_study_data_id[[as.character(study_time)]],]
    
    #Imunity for every age group:
    if(pos_every_age[[time]] == -1e99) return(-1e99)
    age_l = study_data_sel$age_l + 1
    age_u = study_data_sel$age_u + 1
    pos_prop_study <- sapply(1:length(age_l), function(x) sum(pos_every_age[[time]][age_l[x]:age_u[x]]*country_prop_pop[age_l[x]:age_u[x],paste0("X", study_time)])/
                               sum(country_prop_pop[age_l[x]:age_u[x],paste0("X", study_time)]))
    pos_prop_study[pos_prop_study == 1] = 1 - 99e-15
    pos_prop_study[pos_prop_study == 0] = 1e-99
    ob_pos_prop = c(ob_pos_prop, pos_prop_study)
  }
  return(ob_pos_prop)
}

#Choose which scenario
l_sim_run = c("original", "larger_agegroup", "smaller_samplesize")
sim_run = l_sim_run[3]

###plot data fit
seeds_output = readRDS(paste0("results/case_study/sim_", sim_run, "3obs.rds"))
study_data = seeds_output$study_data
model_prop = Re(seeds_output$seeds$result$`Posterior Model Probabilities`)
sim_sample = lapply(seeds_output$l_output, tail, n = 1000)

#generate positive proportion for each age group:
gen_pos_prop = lapply(sim_sample, function(y) sapply(1:nrow(y), function(x) datafit_function(y[x,], study_data)))
gen_pos_prop = lapply(1:length(model_prop), function(x) gen_pos_prop[[x]]*model_prop[x])
gen_pos_prop = Reduce("+", gen_pos_prop)
gen_pos_prop_sum = apply(gen_pos_prop, 1, function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975)))) 
  
plot_fit = study_data %>% select(age_l, age_u, age_mid, total, pos, neg) %>% 
    mutate(age_mid = (age_l + age_u)/2) %>%
    mutate(age_l = formatC(age_l, width = 2, flag = 0)) %>%
    mutate(age_u = formatC(age_u, width = 2, flag = 0)) %>%
    mutate(x_label = paste0(age_l, "-", age_u)) %>% 
    cbind(., binconf(.$pos, .$total)) %>%
    cbind(t(gen_pos_prop_sum)) %>% data.frame()

fit_plot <-
    ggplot(data = plot_fit, aes(x = x_label, y = PointEst, group = 1))+
    geom_pointrange(aes(ymin = Lower, ymax = Upper)) +
    geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), color = "grey", alpha = 0.4) +
    geom_line(aes(y = mean), color = "red") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Age Group", y = "Immune proportion") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -90),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA), 
          rect = element_rect(fill = "transparent",colour = NA),
          text = element_text(size = 14))
      
ggsave(paste0("datafit_", sim_run, ".jpeg"), fit_plot, device = "jpeg", width = 8,
       path = "plots/case_study/")

###Marginal plot:
#Choose which scenario
l_sim_run = c("original", "larger_agegroup", "smaller_samplesize")
sim_run = l_sim_run[1]

seeds_output = readRDS(paste0("results/case_study/sim_", sim_run, "3obs.rds"))
study_data = seeds_output$study_data
model_prop = Re(seeds_output$seeds$result$`Posterior Model Probabilities`)
sim_sample = lapply(seeds_output$l_output, function(x) tail(x, 500))
  
facet_names = c("R[0]", "bar(S)[a]", "~~Outbreak~~time", "~~Model~~Probability")
#get samples:
margin_data = lapply(1:length(sim_sample), function(x) {
    mcmc_re <- sim_sample[[x]] %>% data.frame %>% mutate_at(paste0("par.", x + 2), function(x) 1 - x)
    colnames(mcmc_re) = c(facet_names[1:2], paste0("outbreak time ", 1:x))
    mcmc_re[,3:(3 + x - 1)] = select(mcmc_re, starts_with("outbreak time")) %>%
      apply(1, sort, decreasing = T) %>% matrix(ncol = x, byrow = T)
    mcmc_re$n_ob = x
    return(mcmc_re)
    })
    
#get ensemble data:
mcmc_ensem = lapply(1:length(margin_data), function(x) margin_data[[x]][,facet_names[1:2]]*model_prop[x])
mcmc_ensem = Reduce("+", mcmc_ensem)

#modifiy data to plot
m.par_data = margin_data %>% lapply(melt, id.vars = "n_ob", variable.name = "parameters") %>% Reduce(rbind, .)
m.mcmc_ensem = melt(mcmc_ensem, variable.name = "parameters") %>% mutate(n_ob = "ens")
m.par_true = data.frame(parameters = rep(c(facet_names[1:2], paste0("outbreak time ", 1:5)), 50),
                        value = rep(c(R0, 1 - first_Rinit, ob_year[1:5]), 50), n_ob = "true")
m.par_data = m.par_data %>% rbind(m.mcmc_ensem) %>% rbind(m.par_true)
ob_string = which(str_length(m.par_data$parameters) > 10)
m.par_data$facet = m.par_data$parameters %>% as.character()
m.par_data[ob_string, "facet"] <- facet_names[3]
m.par_data$facet = factor(m.par_data$facet, levels = facet_names)
    
model_prop_data = data.frame(x = 1, value = model_prop, n_ob = paste0(1:5), facet = facet_names[4])
    
margin_plot <- ggplot(data = NULL, aes(y = value, fill = n_ob))+
      geom_violin(data = m.par_data %>% filter(facet == facet_names[3]), aes(x = n_ob, color = parameters), scale = "width") + 
      geom_violin(data = m.par_data %>% filter(parameters %in% c(facet_names[1:2])), aes(x = n_ob), scale = "width")+
      geom_bar(data = model_prop_data, aes(x = x), stat = "identity")+
      scale_color_manual(values = c(rep("black", 6))) +
      scale_fill_manual(values = c("#E6D460","#99A776","#C7421A","#66101f","#1a3a3a",'#304C89', "white"),
                          labels = c("1 outbreak", paste0(2:5, " outbreaks"), "Ensemble", "True")) +
      facet_wrap(.~facet, scales = "free", labeller = "label_parsed", nrow = 2) +
      labs(x = NULL, fill = NULL) +
      guides(color = FALSE) +
      theme_bw() +
      theme(panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            rect = element_rect(fill = "transparent",colour = NA),
            text = element_text(size = 18))

ggsave(paste0("marginplot_", sim_run, ".jpeg"), margin_plot, device = "jpeg", 
       path = "plots/case_study/", width = 12, height = 8)

###Final AR
###Compare 
l_sim_run = c("original", "larger_agegroup", "smaller_samplesize")

AR_ense_l = c()
for(sim_run in l_sim_run){
  #load data and ensemble:
  model_output = readRDS(paste0("results/case_study/sim_", sim_run,"3obs.rds"))
  
  model_prop = Re(model_output$seeds$result$`Posterior Model Probabilities`)
  sim_output = lapply(model_output$l_output, function(x) tail(x, 500))
  study_data = model_output$study_data
    
  ###project and final epi size:
  prob_output = lapply(sim_output, function(x) {
    sapply(1:nrow(x), function(y) gen_immune_func(R0 = x[y,1], Rinit = x[y, 2], year_epi_start = tail(x[y, ], -2),
                                                  l_study_time = l_study_time, country_prop_pop = country_prop_pop,
                                                  study_data = study_data)) %>% Reduce(rbind, .)
    })
  
  #project to 2020:
  whole_year = 20
  proj2020 <- lapply(prob_output, function(x) {cbind(matrix(0 ,ncol = whole_year, nrow = nrow(x)),
                    x[,c(1:(100 - whole_year))])
    })
  
  R_2020_distr = lapply(proj2020, function(x) (x*matrix(country_prop_pop$X2020, byrow = T, nrow = nrow(proj2020[[1]]), ncol =  ncol(proj2020[[1]]))) %>% rowSums %>% as.numeric) 
  
  #Rinf after an outbreak in 2020:
  AR_est = lapply(1:length(R_2020_distr), function(x){
    R0_sel = sim_output[[x]][,"par 1"] %>% as.numeric()
    Rinf_post_2020 = sapply(1:length(R0_sel), function(y) calcRinf(R0_sel[y], R_2020_distr[[x]][y]) - R_2020_distr[[x]][y])
  })
  AR_ense = lapply(1:length(AR_est), function(x) AR_est[[x]]*model_prop[x]) %>% Reduce("+", .) %>% as.numeric()
  AR_ense_l = c(AR_ense_l, list(AR_ense))
}

#true value:
R0; first_Rinit
ob = ob_year
whole_year

prop = gen_immune_func(R0 = R0, Rinit = 0, year_epi_start = ob,
                        l_study_time = l_study_time, country_prop_pop = country_prop_pop,
                        study_data = study_data) %>% unlist
prop_proj2020 = c(numeric(whole_year), prop[1:(100 - whole_year)])
S_2020 = 1 - sum(prop_proj2020*country_prop_pop[,"X2020"])

Rinf = calcRinf(R0, 1 - S_2020)
  AR_true = Rinf- (1 - S_2020)

#plot:
AR_plot_data = Reduce(rbind, AR_ense_l) %>% t %>% data.frame
colnames(AR_plot_data) = l_sim_run
AR_sum_data = AR_plot_data %>% apply(2, quantile, prob = c(0.025, 0.25, 0.5, 0.75, 0.975)) %>% round(4) %>% t %>% data.frame
AR_sum_data$datatype = rownames(AR_sum_data)
AR_plot_data.m = melt(AR_plot_data)
x_labels = c("Original", "Larger age bin", "Smaller samples")
AR_sum_data$datatype = factor(x_labels, levels = x_labels)

  AR_comp_plot  <- ggplot(AR_sum_data, aes(x = datatype))+
  geom_boxplot(aes(ymin = X2.5., lower = X25., middle = X50., upper = X75., ymax = X97.5.),
               stat = "identity")+
  geom_hline(yintercept = AR_true, linetype = 2) +
  theme_bw()+
  labs(x = "Different dataset", y = "IAR")

ggsave(paste0("plots/case_study/AR_", sim_run, ".jpeg"), AR_comp_plot)  
