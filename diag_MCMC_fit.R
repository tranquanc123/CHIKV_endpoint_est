library(dplyr)
library(Hmisc)
library(coda)
library(ggplot2)
library(data.table)
library(gridExtra)
library(stringr)
library(BayesianTools)

#Pop data:
pop_data <- read.csv("data/pop_data.csv")
pop_data <- pop_data %>%
  mutate(X99 = X99 + X100) %>%
  select(-X100)
#sero data:
CHIKV_sero <- read.csv("data/CHIKV_sero_data.csv")

#Study_ID vs Study_name:
ID_name = unique(CHIKV_sero[,c("Study_id", "Study_name")])

source("util_functions.R")
#gen_immune_func: Function to generate immune proportion from estimated parameters
#calcRinf: Function to calculate R infinity from R0 and R initial
#LL_mcmc_sero_function: function to calculate LL

##################
####Diagnosis plot
##################
all_studies_id = list.files(paste0("results/real_sero_study/MCMC/"))

for(file_names in all_studies_id){
  tryCatch({
  print(file_names)
    mcmc_model = readRDS(paste0("results/real_sero_study/MCMC/", file_names))
    Study_ID = str_sub(file_names, 1, -9)
    Study_name = ID_name$Study_name[ID_name$Study_id == Study_ID]
    mcmc_chain = getSample(mcmc_model, start = 5e5, thin = 1e3) %>% data.frame
    n_outbreak = ncol(mcmc_chain) - 2
    if(n_outbreak != 1) mcmc_chain[,2:(ncol(mcmc_chain) - 1)] = t(apply(mcmc_chain[,2:(ncol(mcmc_chain) - 1)], 1, sort, decreasing = T))
    colnames(mcmc_chain) = c("~~R[0]", paste0("Outbreak.", 1:n_outbreak), "~~bar(S)[a]" )
    mcmc_chain =  mcmc_chain %>% mutate(chain = as.character(rep(1:3, nrow(mcmc_chain)/3)), iter = rep(1:(nrow(mcmc_chain)/3), each = 3))
    m.mcmc_chain = melt(mcmc_chain, id.vars = c("chain", "iter"), variable.name = "data_type")
    trace_plot <- 
      ggplot(m.mcmc_chain, aes(x = iter, y = value, group = chain, color = chain)) +
      geom_line() +
      facet_wrap(~data_type, scales = "free", labeller = label_parsed) +
      theme_bw() +
      ggtitle(Study_name)
    
    density_plot <- 
      ggplot(m.mcmc_chain, aes(x = value, group = chain, fill = chain, color = chain)) +
      geom_density(alpha = 0.2) +
      facet_wrap(~data_type, scales = "free", labeller = label_parsed) +
      theme_bw()
    
    diag_plot <- gridExtra::grid.arrange(trace_plot, density_plot)
    ggsave(paste0("plots/diag_plot/" , str_sub(file_names, 1, -5), ".png"), diag_plot, width = 10, height = 10)
  }, error = function(e) print("Error"))
}

#################
####Data fit
#################
#ensemble:
###Fit data:
datafit_function <-function(par, study_data, y_time, l_study_time){
  R0 = par[1]; Rinit = par[2]; year_epi_start = par[3:length(par)]
  year_epi_start = year_epi_start[year_epi_start > 0]
  
  #Logit Transform Rinit:
  t_Rinit <- Rinit#exp(Rinit)/(exp(Rinit) + 1)
  study_time = unique(study_data$time)
    year_epi_start_study = sort(year_epi_start, decreasing = T)
    
    #if(length(year_dist) != 0) {
    pos_every_age = numeric(100)
    year_dist = c(head(year_epi_start_study, -1) - year_epi_start_study[-1], tail(year_epi_start_study, 1))
    year_floor = floor(year_epi_start_study)
    #Since S(0) is the susceptible pop before the first study, the S(0) of the first ob of the first study is
    #(assuming susceptible proportion distributed evenly between age groups):
    
    if(y_time == min(l_study_time)){
      first_epi_year_whole = year_floor[1]
      first_epi_year_part = year_epi_start_study[1] - first_epi_year_whole
      if(first_epi_year_whole == 0){
        t_Rinit = (1 - first_epi_year_part)*sum(t_Rinit*country_prop_pop[100:(100 - first_epi_year_whole),paste0("X", ifelse(study_time - first_epi_year_whole < 1950, 1950, study_time - first_epi_year_whole))])
      } else {
        t_Rinit = sum(t_Rinit*country_prop_pop[100:(100 - first_epi_year_whole),paste0("X", ifelse(study_time - first_epi_year_whole - 1 < 1950, 1950, study_time - first_epi_year_whole - 1))]) +
          (1 - first_epi_year_part)*sum(t_Rinit*country_prop_pop[100:(100 - first_epi_year_whole),paste0("X", ifelse(study_time - first_epi_year_whole < 1950, 1950, study_time - first_epi_year_whole))])
      }
    }
    
    for(i in 1:length(year_dist)){
      Rinf = calcRinf(R0, t_Rinit) #Recovery cov after outbreak
      pos_every_age = 1 - (1 - Rinf)*(1 - pos_every_age)/(1 - t_Rinit) #(1 - pos_every_age)*(Rinf - t_Rinit) + pos_every_age
      #project
      whole_year = floor(year_dist)[i]; part_year = year_epi_start_study[i] - floor(year_epi_start_study[i]) 
      pos_every_age <- c(numeric(length = whole_year), (1 - part_year)*pos_every_age[1], pos_every_age[c(1:(99 - whole_year))])
      pos_prop_pop <- pos_every_age*country_prop_pop[,paste0("X", ifelse(study_time - whole_year - 1 < 1950, 1950, study_time - whole_year - 1))]
      t_Rinit <- sum(pos_prop_pop)
    }
    #Obtain study-specific age group proportion:
    age_l = study_data$age_l + 1
    age_u = study_data$age_u + 1
    pos_prop_study <- sapply(1:length(age_l), function(x) sum(pos_every_age[age_l[x]:age_u[x]]*country_prop_pop[age_l[x]:age_u[x],paste0("X", study_time)])/
                               sum(country_prop_pop[age_l[x]:age_u[x],paste0("X", study_time)]))
    pos_prop_study[pos_prop_study == 1] = 1 - 99e-15
    pos_prop_study[pos_prop_study == 0] = 1e-99
    
    #derive the true pos:
    ob_pos_prop <- pos_prop_study

  return(ob_pos_prop)
}

l_datafit_plot = c()
all_sites_id = list.files(paste0("results/real_sero_study/RJMCMC/"))

for(file_names in all_sites_id){
  print(file_names)
  Study_ID = str_sub(file_names, 1, -5)
  study_data = filter(CHIKV_sero, Study_id == Study_ID, total != 0)
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
  #sens_specs_data = filter(sens_specs_CHIKV, Study_id == Study_ID)

  seeds_output = readRDS(paste0("results/real_sero_study/RJMCMC/", Study_ID, ".rds"))
  model_prop = as.numeric(seeds_output$seeds$result$`Posterior Model Probabilities`)
  sim_sample = lapply(seeds_output$l_output, function(x) tail(x, 1000))
  
  study_year = unique(study_data$time)
  study_year_dist = max(study_year) - study_year
  for(y_time in study_year){
    sel_study_data <- filter(study_data, time == y_time)
    mod_sim_sample <- 
      lapply(sim_sample, function(x){
      #mod year:
      year_est = x[,3:ncol(x)]
      mod_year = year_est - (max(study_year) - y_time)
      #mod_year = mod_year[mod_year > 0]
      mod_year_data = x
      mod_year_data[, 3:ncol(mod_year_data)] = mod_year
      return(mod_year_data)
    })
    gen_pos_prop = lapply(mod_sim_sample, function(y) sapply(1:nrow(y), function(x) datafit_function(y[x,], sel_study_data, y_time, study_year)))
    gen_pos_prop = lapply(1:length(model_prop), function(x) gen_pos_prop[[x]]*model_prop[x])
    gen_pos_prop = Reduce("+", gen_pos_prop)
    gen_pos_prop_sum = apply(gen_pos_prop, 1, function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975)))) 
    
    plot_fit = sel_study_data %>% select(Study_name, time, age_l, age_u, age_mid, total, pos, neg) %>% 
      mutate(plot_tittle = paste0(Study_name, " ", time)) %>%
      mutate(age_mid = (age_l + age_u)/2) %>%
      mutate(age_l = formatC(age_l, width = 2, flag = 0)) %>%
      mutate(age_u = formatC(age_u, width = 2, flag = 0)) %>%
      mutate(x_label = paste0(age_l, "-", age_u)) %>% 
      cbind(., binconf(.$pos, .$total)) %>%
      cbind(t(gen_pos_prop_sum)) %>% data.frame()
    l_datafit_plot = c(l_datafit_plot, list(plot_fit))
  }
}
datafit_plot = Reduce(rbind, l_datafit_plot)

fit_plot <-
  ggplot(data = datafit_plot, aes(x = x_label, y = PointEst, group = 1))+
  geom_pointrange(aes(ymin = Lower, ymax = Upper)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), color = "grey", alpha = 0.4) +
  geom_line(aes(y = mean), color = "red") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Age Group", y = "Immune proportion") +
  facet_wrap(~plot_tittle, scale = "free", ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        rect = element_rect(fill = "transparent",colour = NA),
        text = element_text(size = 12),
        axis.text.x.bottom = element_text(size = 7))
ggsave("data_fit_all.png", fit_plot, device = "png", 
           path = "plots/datafit/", width = 12, height = 15)
