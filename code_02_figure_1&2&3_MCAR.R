#################################################################
########  consistency & convergence rate & excess risk    #######
#################################################################
# MCAR

mustarstar=rbind( c(-sqrt(6)/2, -sqrt(6)/2) , 
                  c((sqrt(6)+3*sqrt(2))/4, (sqrt(6)-3*sqrt(2))/4), 
                  c((sqrt(6)-3*sqrt(2))/4, (sqrt(6)+3*sqrt(2))/4)  )

### MCAR: 0.1 ------

Lfun_minvalue_vec <- rep(0,100)
mustar_array <- array(0,dim=c(3,2,100))
for (rept in 1:100){
  data <- generate_data_2026(n=10000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MCAR", MCARmissprob = rep(0.1,2))
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random") )
  Lfun_minvalue_vec[rept] <- res_once$obj_val/(10000*3)
  reorderindex <- order(res_once$centers[,1])
  mustar_array[,,rept] <- res_once$centers[reorderindex,]
  if (rept%%10==0){print(rept)}
}
Lfun_minvalue <- mean(Lfun_minvalue_vec)
mustar <- apply(mustar_array,1:2,mean)


ExperimentResults_123_MCAR_10pct <- list( 
  Lfun_minvalue_vec=Lfun_minvalue_vec, mustar_array=mustar_array,
  Lfun_minvalue=Lfun_minvalue, mustar=mustar,
  lossdiff=matrix(0,100,8),
  mse=matrix(0,100,8) )
each_n <- c(100,300,500, 1000, 2000, 3000, 4000, 5000)
for (ni in 1:8){
  for (rept in 1:100){
    data <- generate_data_2026(n=each_n[ni]*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                               missmechanism = "MCAR", MCARmissprob = rep(0.1,2))
    res_once <- kmmd_2026(x=data$Missing,k=3,
                          initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random")  )
    ExperimentResults_123_MCAR_10pct$mse[rept,ni] <- MSEcenters(mustar,res_once$centers)
    ExperimentResults_123_MCAR_10pct$lossdiff[rept,ni] <- res_once$obj_val/(each_n[ni]*3) - Lfun_minvalue
    if (rept%%10==0){print(rept)}
  }
  print(ni)
}


ExperimentResults_123_MCAR_10pct$summarized <- data.frame(
  n=each_n*3, 
  mse_mean= apply(ExperimentResults_123_MCAR_10pct$mse,2,mean),
  mse_sd= apply(ExperimentResults_123_MCAR_10pct$mse,2,sd),
  lossdiff_mean= apply(abs(ExperimentResults_123_MCAR_10pct$lossdiff),2,mean),
  lossdiff_sd= apply(abs(ExperimentResults_123_MCAR_10pct$lossdiff),2,sd)
)



#loss convergence 
ggplot(ExperimentResults_123_MCAR_10pct$summarized, aes(x = n, y = lossdiff_mean)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_errorbar(aes(ymin = lossdiff_mean - lossdiff_sd/sqrt(100), ymax = lossdiff_mean + lossdiff_sd/sqrt(100)), width = 50, color = "black") +
  labs(x = "Sample Size (n)",
       y = "Excess Risk") + ylim(0,0.07) +
  theme_bw()

#sqrt(n)-excess risk bound
ols_fit <- lm(log10(ExperimentResults_123_MCAR_10pct$summarized$lossdiff_mean) ~ 1, 
              offset = -0.5*log10(ExperimentResults_123_MCAR_10pct$summarized$n))
ggplot(ExperimentResults_123_MCAR_10pct$summarized, aes(x = n, y = lossdiff_mean)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_abline(slope = -0.5, 
              intercept = as.numeric(coef(ols_fit)),
              linetype = "dashed", color = "red") +
  labs(x = "n",
       y = TeX("$|\\hat{m}_n - m^*|$")) +
  theme_bw()+
  scale_x_log10() + scale_y_log10() 

#consistency
ggplot(ExperimentResults_123_MCAR_10pct$summarized, aes(x = n, y = mse_mean)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_errorbar(aes(ymin = mse_mean - mse_sd/sqrt(100), ymax = mse_mean + mse_sd/sqrt(100)), width = 200, color = "black") +
  labs(x = "n",
       y = TeX("$\\|\\hat{M}_{n,t} - M^*\\|_F^2$")) +
  theme_bw()

#convergence rate: log(MSE^2) = intercept - log(n) 
ggplot(ExperimentResults_123_MCAR_10pct$summarized, aes(x = n, y = mse_mean)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_abline(slope = -1, 
              intercept = mean(log10(ExperimentResults_123_MCAR_10pct$summarized$mse_mean) + log10(ExperimentResults_123_MCAR_10pct$summarized$n)), 
              linetype = "dashed", color = "red") +
  labs(x = "n",
       y = TeX("$Avg(\\|\\hat{M}_{n,t} - M^*\\|_F^2)$")) +
  theme_bw()+
  scale_x_log10() + scale_y_log10() 


### MCAR: 0.3 ------

Lfun_minvalue_vec <- rep(0,100)
mustar_array <- array(0,dim=c(3,2,100))
for (rept in 1:100){
  data <- generate_data_2026(n=10000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MCAR", MCARmissprob = rep(0.3,2))
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random") )
  Lfun_minvalue_vec[rept] <- res_once$obj_val/(10000*3)
  reorderindex <- order(res_once$centers[,1])
  mustar_array[,,rept] <- res_once$centers[reorderindex,]
  print(rept)
}
Lfun_minvalue <- mean(Lfun_minvalue_vec)
mustar <- apply(mustar_array,1:2,mean)


ExperimentResults_123_MCAR_30pct <- list( 
  Lfun_minvalue_vec=Lfun_minvalue_vec, mustar_array=mustar_array,
  Lfun_minvalue=Lfun_minvalue, mustar=mustar,
  lossdiff=matrix(0,100,8),
  mse=matrix(0,100,8) )
each_n <- c(100,300,500, 1000, 2000, 3000, 4000, 5000)
for (ni in 1:8){
  for (rept in 1:100){
    data <- generate_data_2026(n=each_n[ni]*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                               missmechanism = "MCAR", MCARmissprob = rep(0.3,2))
    res_once <- kmmd_2026(x=data$Missing,k=3,
                          initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random")  )
    ExperimentResults_123_MCAR_30pct$mse[rept,ni] <- MSEcenters(mustar,res_once$centers)
    ExperimentResults_123_MCAR_30pct$lossdiff[rept,ni] <- res_once$obj_val/(each_n[ni]*3) - Lfun_minvalue
    if (rept%%10==0){print(rept)}
  }
  print(ni)
}


## draw figures
ExperimentResults_123_MCAR_30pct$summarized <- data.frame(
  n=each_n*3, 
  mse_mean= apply(ExperimentResults_123_MCAR_30pct$mse,2,mean),
  mse_sd= apply(ExperimentResults_123_MCAR_30pct$mse,2,sd),
  lossdiff_mean= apply(abs(ExperimentResults_123_MCAR_30pct$lossdiff),2,mean),
  lossdiff_sd= apply(abs(ExperimentResults_123_MCAR_30pct$lossdiff),2,sd)
)



#loss convergence 
ggplot(ExperimentResults_123_MCAR_30pct$summarized, aes(x = n, y = lossdiff_mean)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_errorbar(aes(ymin = lossdiff_mean - lossdiff_sd/sqrt(100), ymax = lossdiff_mean + lossdiff_sd/sqrt(100)), width = 50, color = "black") +
  labs(x = "Sample Size (n)",
       y = "Excess Risk") + ylim(0,0.07) +
  theme_bw()

#sqrt(n)-excess risk bound
ols_fit <- lm(log10(ExperimentResults_123_MCAR_30pct$summarized$lossdiff_mean) ~ 1, 
              offset = -0.5*log10(ExperimentResults_123_MCAR_30pct$summarized$n))
ggplot(ExperimentResults_123_MCAR_30pct$summarized, aes(x = n, y = lossdiff_mean)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_abline(slope = -0.5, 
              intercept = as.numeric(coef(ols_fit)),
              linetype = "dashed", color = "red") +
  labs(x = "n",
       y = TeX("$|\\hat{m}_n - m^*|$")) +
  theme_bw()+
  scale_x_log10() + scale_y_log10() 


#consistency
ggplot(ExperimentResults_123_MCAR_30pct$summarized, aes(x = n, y = mse_mean)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_errorbar(aes(ymin = mse_mean - mse_sd/sqrt(100), ymax = mse_mean + mse_sd/sqrt(100)), width = 200, color = "black") +
  labs(x = "n",
       y = TeX("$\\|\\hat{M}_{n,t} - M^*\\|_F^2$")) +
  theme_bw()

#convergence rate: log(MSE^2) = intercept - log(n) 
ggplot(ExperimentResults_123_MCAR_30pct$summarized, aes(x = n, y = mse_mean)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_abline(slope = -1, 
              intercept = mean(log10(ExperimentResults_123_MCAR_30pct$summarized$mse_mean) + log10(ExperimentResults_123_MCAR_30pct$summarized$n)), 
              linetype = "dashed", color = "red") +
  labs(x = "n",
       y = TeX("$Avg(\\|\\hat{M}_{n,t} - M^*\\|_F^2)$")) +
  theme_bw()+
  scale_x_log10() + scale_y_log10() 



### MCAR: 0.5 ------

Lfun_minvalue_vec <- rep(0,100)
mustar_array <- array(0,dim=c(3,2,100))
for (rept in 1:100){
  data <- generate_data_2026(n=10000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MCAR", MCARmissprob = rep(0.5,2))
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random") )
  Lfun_minvalue_vec[rept] <- res_once$obj_val/(nrow(data$Missing))
  reorderindex <- order(res_once$centers[,1])
  mustar_array[,,rept] <- res_once$centers[reorderindex,]
  if (rept%%10==0){print(rept)}
}
Lfun_minvalue <- mean(Lfun_minvalue_vec)
mustar <- apply(mustar_array,1:2,mean)


ExperimentResults_123_MCAR_50pct <- list( 
  Lfun_minvalue_vec=Lfun_minvalue_vec, mustar_array=mustar_array,
  Lfun_minvalue=Lfun_minvalue, mustar=mustar,
  lossdiff=matrix(0,100,8),
  mse=matrix(0,100,8) )
each_n <- c(100,300,500, 1000, 2000, 3000, 4000, 5000)
for (ni in 1:8){
  for (rept in 1:100){
    data <- generate_data_2026(n=each_n[ni]*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                               missmechanism = "MCAR", MCARmissprob = rep(0.5,2))
    res_once <- kmmd_2026(x=data$Missing,k=3,
                          initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random")  )
    ExperimentResults_123_MCAR_50pct$mse[rept,ni] <- MSEcenters(mustar,res_once$centers)
    ExperimentResults_123_MCAR_50pct$lossdiff[rept,ni] <- res_once$obj_val/(nrow(data$Missing)) - Lfun_minvalue
    if (rept%%10==0){print(rept)}
  }
  print(ni)
}



## draw figures
ExperimentResults_123_MCAR_50pct$summarized <- data.frame(
  n=each_n*3, 
  mse_mean= apply(ExperimentResults_123_MCAR_50pct$mse,2,mean),
  mse_sd= apply(ExperimentResults_123_MCAR_50pct$mse,2,sd),
  lossdiff_mean= apply(abs(ExperimentResults_123_MCAR_50pct$lossdiff),2,mean),
  lossdiff_sd= apply(abs(ExperimentResults_123_MCAR_50pct$lossdiff),2,sd)
)


#loss convergence 
ggplot(ExperimentResults_123_MCAR_50pct$summarized, aes(x = n, y = lossdiff_mean)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_errorbar(aes(ymin = lossdiff_mean - lossdiff_sd/sqrt(100), ymax = lossdiff_mean + lossdiff_sd/sqrt(100)), width = 50, color = "black") +
  labs(x = "Sample Size (n)",
       y = "Excess Risk") + ylim(0,0.08) +
  theme_bw()

#sqrt(n)-excess risk bound
ols_fit <- lm(log10(ExperimentResults_123_MCAR_50pct$summarized$lossdiff_mean) ~ 1, 
              offset = -0.5*log10(ExperimentResults_123_MCAR_50pct$summarized$n))
ggplot(ExperimentResults_123_MCAR_50pct$summarized, aes(x = n, y = lossdiff_mean)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_abline(slope = -0.5, 
              intercept = as.numeric(coef(ols_fit)),
              linetype = "dashed", color = "red") +
  labs(x = "n",
       y = TeX("$|\\hat{m}_n - m^*|$")) +
  theme_bw()+
  scale_x_log10() + scale_y_log10() 


#consistency
ggplot(ExperimentResults_123_MCAR_50pct$summarized, aes(x = n, y = mse_mean)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_errorbar(aes(ymin = mse_mean - mse_sd/sqrt(100), ymax = mse_mean + mse_sd/sqrt(100)), width = 200, color = "black") +
  labs(x = "n",
       y = TeX("$\\|\\hat{M}_{n,t} - M^*\\|_F^2$")) +
  theme_bw()

#convergence rate: log(MSE^2) = intercept - log(n) 
ggplot(ExperimentResults_123_MCAR_50pct$summarized, aes(x = n, y = mse_mean)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_abline(slope = -1, 
              intercept = mean(log10(ExperimentResults_123_MCAR_50pct$summarized$mse_mean) + log10(ExperimentResults_123_MCAR_50pct$summarized$n)), 
              linetype = "dashed", color = "red") +
  labs(x = "n",
       y = TeX("$Avg(\\|\\hat{M}_{n,t} - M^*\\|_F^2)$")) +
  theme_bw()+
  scale_x_log10() + scale_y_log10() 



