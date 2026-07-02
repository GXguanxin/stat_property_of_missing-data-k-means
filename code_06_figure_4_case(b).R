#############################################
########  converge to true centers    #######
#############################################
# case (b)

rho_starstar_list= c(1,2,3,4, seq(5,7,by=0.2), 8,9,10)

### MCAR 0.1 ------ 

ExperimentResults_5_k2_rho_haszero_MCAR_10pct <- list(
  allsettings=list(
    rho_starstar_list= c(1,2,3,4, seq(5,7,by=0.2), 8,9,10) ,p=2,k=2,
    variance_vec=rep(1,2), beta_starstar=3,
    missmechanism="MCAR",MCARmissprob=rep(0.1,2)
  ),
  allresults=list(MSEvsrho=array(NA,dim=c(18,100,2)),
                  CERvsrho=array(NA,dim=c(18,100,2)) )
)
for (rho_id in 1:18){
  rho_starstar = rho_starstar_list[rho_id]
  for (rept in 1:100){
    mustarstar_k2 <- rbind( c(0,0), c(rho_starstar,0)  )
    data <- generate_truncated_data_2016(n=10000, p=2, k=2, 
                                         centers=mustarstar_k2,
                                         variance_vec=rep(1,2), beta_starstar=3,
                                         missmechanism="MCAR",
                                         MCARmissprob=rep(0.1,2))
    res_once <- kmmd_2026(x=data$Missing,k=2,
                          initial_centers_array = get_initial_centers_2026(x=data$Missing,k=2,nstart = 30,init_method = "random") )
    res_km <- kmeans(data$Orig,centers=2,iter.max = 100,nstart = 30)
    ExperimentResults_5_k2_rho_haszero_MCAR_10pct$allresults$MSEvsrho[rho_id,rept,] <- c(
      MSEcenters(mustarstar_k2,res_km$centers),
      MSEcenters(mustarstar_k2,res_once$centers)   ) 
    ExperimentResults_5_k2_rho_haszero_MCAR_10pct$allresults$CERvsrho[rho_id,rept,] <- c(
      classError(classification = res_km$cluster, class = data$Origlabel)$errorRate,
      classError(classification = res_once$cluster, class = data$Misslabel)$errorRate   ) 
  }
  print(rho_starstar)
}


### MCAR 0.3 ------ 

ExperimentResults_5_k2_rho_haszero_MCAR_30pct <- list(
  allsettings=list(
    rho_starstar_list= c(1,2,3,4, seq(5,7,by=0.2), 8,9,10) ,p=2,k=2,
    variance_vec=rep(1,2), beta_starstar=3,
    missmechanism="MCAR",MCARmissprob=rep(0.3,2)
  ),
  allresults=list(MSEvsrho=array(NA,dim=c(18,100,2)),
                  CERvsrho=array(NA,dim=c(18,100,2)) )
)
for (rho_id in 1:18){
  rho_starstar = rho_starstar_list[rho_id]
  for (rept in 1:100){
    mustarstar_k2 <- rbind( c(0,0), c(rho_starstar,0)  )
    data <- generate_truncated_data_2016(n=10000, p=2, k=2, 
                                         centers=mustarstar_k2,
                                         variance_vec=rep(1,2), beta_starstar=3,
                                         missmechanism="MCAR",
                                         MCARmissprob=rep(0.3,2))
    res_once <- kmmd_2026(x=data$Missing,k=2,
                          initial_centers_array = get_initial_centers_2026(x=data$Missing,k=2,nstart = 30,init_method = "random") )
    res_km <- kmeans(data$Orig,centers=2,iter.max = 100,nstart = 30)
    ExperimentResults_5_k2_rho_haszero_MCAR_30pct$allresults$MSEvsrho[rho_id,rept,] <- c(
      MSEcenters(mustarstar_k2,res_km$centers),
      MSEcenters(mustarstar_k2,res_once$centers)   ) 
    ExperimentResults_5_k2_rho_haszero_MCAR_30pct$allresults$CERvsrho[rho_id,rept,] <- c(
      classError(classification = res_km$cluster, class = data$Origlabel)$errorRate,
      classError(classification = res_once$cluster, class = data$Misslabel)$errorRate   ) 
  }
  print(rho_starstar)
}


### MCAR 0.5 ------ 

ExperimentResults_5_k2_rho_haszero_MCAR_50pct <- list(
  allsettings=list(
    rho_starstar_list= c(1,2,3,4, seq(5,7,by=0.2), 8,9,10) ,p=2,k=2,
    variance_vec=rep(1,2), beta_starstar=3,
    missmechanism="MCAR",MCARmissprob=rep(0.5,2)
  ),
  allresults=list(MSEvsrho=array(NA,dim=c(18,100,2)),
                  CERvsrho=array(NA,dim=c(18,100,2)) )
)
for (rho_id in 1:18){
  rho_starstar = rho_starstar_list[rho_id]
  for (rept in 1:100){
    mustarstar_k2 <- rbind( c(0,0), c(rho_starstar,0)  )
    data <- generate_truncated_data_2016(n=10000, p=2, k=2, 
                                         centers=mustarstar_k2,
                                         variance_vec=rep(1,2), beta_starstar=3,
                                         missmechanism="MCAR",
                                         MCARmissprob=rep(0.5,2))
    res_once <- kmmd_2026(x=data$Missing,k=2,
                          initial_centers_array = get_initial_centers_2026(x=data$Missing,k=2,nstart = 30,init_method = "random") )
    res_km <- kmeans(data$Orig,centers=2,iter.max = 100,nstart = 30)
    ExperimentResults_5_k2_rho_haszero_MCAR_50pct$allresults$MSEvsrho[rho_id,rept,] <- c(
      MSEcenters(mustarstar_k2,res_km$centers),
      MSEcenters(mustarstar_k2,res_once$centers)   ) 
    ExperimentResults_5_k2_rho_haszero_MCAR_50pct$allresults$CERvsrho[rho_id,rept,] <- c(
      classError(classification = res_km$cluster, class = data$Origlabel)$errorRate,
      classError(classification = res_once$cluster, class = data$Misslabel)$errorRate   ) 
  }
  print(rho_starstar)
}


### draw figures ------

# summarize data
ExperimentResults_5_k2_rho_haszero_MCAR_allpct <- list(
  summarized = list(), figures = list()  )
ExperimentResults_5_k2_rho_haszero_MCAR_allpct$summarized <- list(
  mse = data.frame(
    rho_starstar = rep(rho_starstar_list, 4),
    mse_mean = c(
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_10pct$allresults$MSEvsrho[,,1], 1, mean),
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_10pct$allresults$MSEvsrho[,,2], 1, mean),
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_30pct$allresults$MSEvsrho[,,2], 1, mean),
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_50pct$allresults$MSEvsrho[,,2], 1, mean)
    ),
    mse_sd = c(
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_10pct$allresults$MSEvsrho[,,1], 1, sd),
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_10pct$allresults$MSEvsrho[,,2], 1, sd),
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_30pct$allresults$MSEvsrho[,,2], 1, sd),
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_50pct$allresults$MSEvsrho[,,2], 1, sd)
    ),
    method= factor( rep(c("fully observed data", "10% missing data","30% missing data","50% missing data"), each = 18) , 
                    levels =  c("fully observed data", "10% missing data","30% missing data","50% missing data") )
  ),
  cer = data.frame(
    rho_starstar = rep(rho_starstar_list, 4),
    cer_mean = c(
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_10pct$allresults$CERvsrho[,,1], 1, mean),
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_10pct$allresults$CERvsrho[,,2], 1, mean),
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_30pct$allresults$CERvsrho[,,2], 1, mean),
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_50pct$allresults$CERvsrho[,,2], 1, mean)
    ),
    cer_sd = c(
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_10pct$allresults$CERvsrho[,,1], 1, sd),
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_10pct$allresults$CERvsrho[,,2], 1, sd),
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_30pct$allresults$CERvsrho[,,2], 1, sd),
      apply(ExperimentResults_5_k2_rho_haszero_MCAR_50pct$allresults$CERvsrho[,,2], 1, sd)
    ),
    method= factor( rep(c("fully observed data", "10% missing data","30% missing data","50% missing data"), each = 18) , 
                    levels =  c("fully observed data", "10% missing data","30% missing data","50% missing data") )
  )
)

# mcar 0.1
ggplot(
  subset(ExperimentResults_5_k2_rho_haszero_MCAR_allpct$summarized$mse,
         method %in% c("fully observed data", "10% missing data") &
           rho_starstar %in% c(1:10) ),
  aes(x = rho_starstar, y = mse_mean, color = method, linetype=method, shape=method, group=method)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) + 
  scale_linetype_manual(values = c("solid", "dashed"))+
  scale_color_manual(values = c("black", "blue"))+
  scale_shape_manual(values = c(16, 17))+
  scale_x_continuous(breaks = seq(1, 10, by = 1))+
  labs(x = expression(rho), y = "MSE",color = NULL, shape=NULL, linetype = NULL) + 
  theme_bw() + coord_cartesian(ylim = c(0, 0.2)) +
  theme( legend.position = c(1, 1), legend.justification = c("right", "top"),
         legend.background = element_rect(fill = "white", color = "black", size = 0.3))

# mcar 0.3 
ggplot(
  subset(ExperimentResults_5_k2_rho_haszero_MCAR_allpct$summarized$mse,
         method %in% c("fully observed data", "30% missing data") &
           rho_starstar %in% c(1:10) ),
  aes(x = rho_starstar, y = mse_mean, color = method, linetype=method, shape=method, group=method)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) + 
  scale_linetype_manual(values = c("solid", "dashed"))+
  scale_color_manual(values = c("black", "blue"))+
  scale_shape_manual(values = c(16, 17))+
  scale_x_continuous(breaks = seq(1, 10, by = 1))+
  labs(x = expression(rho), y = "MSE",color = NULL, shape=NULL, linetype = NULL) + 
  theme_bw() + coord_cartesian(ylim = c(0, 0.2)) +
  theme( legend.position = c(1, 1), legend.justification = c("right", "top"),
         legend.background = element_rect(fill = "white", color = "black", size = 0.3))

# mcar 0.5 
ggplot(
  subset(ExperimentResults_5_k2_rho_haszero_MCAR_allpct$summarized$mse,
         method %in% c("fully observed data", "50% missing data") &
           rho_starstar %in% c(1:10) ),
  aes(x = rho_starstar, y = mse_mean, color = method, linetype=method, shape=method, group=method)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) + 
  scale_linetype_manual(values = c("solid", "dashed"))+
  scale_color_manual(values = c("black", "blue"))+
  scale_shape_manual(values = c(16, 17))+
  scale_x_continuous(breaks = seq(1, 10, by = 1))+
  labs(x = expression(rho), y = "MSE",color = NULL, shape=NULL, linetype = NULL) + 
  theme_bw() + coord_cartesian(ylim = c(0, 0.3)) +
  theme( legend.position = c(1, 1), legend.justification = c("right", "top"),
         legend.background = element_rect(fill = "white", color = "black", size = 0.3))

