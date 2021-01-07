#function to simulate immunity profile from R0, Rinit, outbreak years:
gen_immune_func <-function(R0, year_epi_start, Rinit, study_data, l_study_time, country_prop_pop){
  
  l_pos = c()
  pos_every_age = rep(0, 100)
  
  #which ob in which study:
  year_epi_start = sort(year_epi_start)
  study_year_dist = rev(tail(l_study_time, 1) - rev(l_study_time))
  l_ob_in_study = c()
  for(y in study_year_dist){
    ob_dist = year_epi_start - y; ob_dist = ob_dist[ob_dist > 0]
    l_ob_in_study = c(l_ob_in_study, list(ob_dist))
    if(length(ob_dist != 0)){
      year_epi_start = head(year_epi_start, -length(ob_dist))
    }
  }
  
  #for each study time:
  for(time in 1:length(l_study_time)){
    study_time = l_study_time[time]
    current_study_year_dist = tail(l_study_time, 1) - study_time
    study_data_sel = study_data[l_study_data_id[[as.character(study_time)]],]
    year_epi_start_study = sort(l_ob_in_study[[time]], decreasing = T)
    
    #project pop immunity from the previous study time to the current study time 
    if(current_study_year_dist < max(l_study_time) - min(l_study_time)){
      year_dist_to_the_last_study = l_study_time[time] - l_study_time[time - 1]
      pos_every_age <- c(numeric(year_dist_to_the_last_study), head(pos_every_age, - year_dist_to_the_last_study))
      pos_prop_pop = pos_every_age*country_prop_pop[,paste0("X", ifelse(study_time < 1950, 1950, study_time))]
      Rinit = sum(pos_prop_pop)
    }
    
    year_dist = c(head(year_epi_start_study, -1) - year_epi_start_study[-1], tail(year_epi_start_study, 1))
    year_floor = floor(year_epi_start_study)
    
    #Since S(0) is the susceptible pop before the first study, the S(0) of the first ob of the first study is
    #(assuming susceptible proportion distributed evenly between age groups):
    if(time == 1){
      if(length(year_dist) == 0){
        Rinit = 0 #if the first study has no outbreak, Rinit of the 1st outbreak of the next study will be 0
      } else {
        first_epi_year_whole = floor(year_dist[1])
        first_epi_year_part = year_dist[1] - first_epi_year_whole
        if(first_epi_year_whole == 0){
          Rinit = (1 - first_epi_year_part)*sum(Rinit*country_prop_pop[100:(100 - first_epi_year_whole),paste0("X", ifelse(study_time - first_epi_year_whole < 1950, 1950, study_time - first_epi_year_whole))])
        } else {
          Rinit = sum(Rinit*country_prop_pop[100:(100 - first_epi_year_whole),paste0("X", ifelse(study_time - first_epi_year_whole - 1 < 1950, 1950, study_time - first_epi_year_whole - 1))]) +
            (1 - first_epi_year_part)*sum(Rinit*country_prop_pop[100:(100 - first_epi_year_whole),paste0("X", ifelse(study_time - first_epi_year_whole < 1950, 1950, study_time - first_epi_year_whole))])
        }
      }
    }
    
    if(length(year_dist) != 0){
      for(i in 1:length(year_dist)){
        Rinf = calcRinf(R0, Rinit)
        #if(Rinf - Rinit < 0.01) return(-1e99)
        pos_every_age = 1 - (1 - Rinf)*(1 - pos_every_age)/(1 - Rinit) 
        #project to the current study year:
        whole_year = floor(year_dist)[i]; part_year = year_dist[i] - floor(year_dist[i]) 
        pos_every_age <- c(numeric(length = whole_year), (1 - part_year)*pos_every_age[1], pos_every_age[c(1:(99 - whole_year))])
        pos_prop_pop = pos_every_age*country_prop_pop[,paste0("X", ifelse(study_time - year_floor[i] - 1 < 1950, 1950, study_time - year_floor[i] - 1))]
        Rinit = sum(pos_prop_pop)
      }
    }
    
    l_pos = c(l_pos, list(pos_every_age)) 
  }
  return(l_pos)
}

#Function to calculate R infinity from R0 and R initial
calcRinf = function(R0, Rinit){
  optimize(
    f = function(Rinf){
      abs(1 - Rinf - exp(log(1 - Rinit) + (-R0 * (Rinf - (Rinit)))))
    },
    interval = c(Rinit + 0.001, 1))$minimum
}

#Function to calculate LL:
LL_mcmc_sero_function <-function(par){
  par = as.numeric(par)
  R0 = par[1]; year_epi_start = par[2:(length(par) - 1)]; Rinit = tail(par, 1)
  
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
  
  return(sum(dbinom(study_data_sel$pos, study_data_sel$total, ob_pos_prop, log = T)))
}