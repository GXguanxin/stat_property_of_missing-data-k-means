#########################################
########  asymptotic normality    #######
#########################################
# MCAR & MNAR

### MCAR: 0.1 ------

#calculate mustar
Lfun_minvalue_vec <- rep(0,100)
mustar_array <- array(0,dim=c(3,2,100))
for (rept in 1:100){
  data <- generate_data_2026(n=10000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MCAR",MCARmissprob = rep(0.1,2))
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random") )
  Lfun_minvalue_vec[rept] <- res_once$obj_val/(nrow(data$Missing))
  reorderindex <- order(res_once$centers[,1])
  mustar_array[,,rept] <- res_once$centers[reorderindex,]
  if (rept%%10==0){print(rept)}
}
Lfun_minvalue <- mean(Lfun_minvalue_vec)
mustar <- apply(mustar_array,1:2,mean)


est_mu_array <- array(0,dim=c(3,2,1000))
thetahat <- matrix(0,3*2,1000)
for (rept in 1:1000){
  data <- generate_data_2026(n=2000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MCAR",MCARmissprob = rep(0.1,2))
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random")  )
  est_mu_array[,,rept] <- align_centers(res_once$centers,mustar)
  thetahat[,rept] <- sqrt(nrow(data$Missing))*( as.vector(t(est_mu_array[,,rept] - mustar)) )
  if (rept%%100==0){print(rept)}
}
#Malahabios distance^2
Smat_inv <- solve(cov(t(thetahat)))
Mdist2 <- apply(thetahat, 2, function(z) (t(z))%*%Smat_inv%*%z )

ExperimentResults_4_MCAR_10pct <- list(
  Lfun_minvalue_vec=Lfun_minvalue_vec, mustar_array=mustar_array,
  Lfun_minvalue=Lfun_minvalue, mustar=mustar,
  est_mu_array=est_mu_array, thetahat=thetahat,
  Smat_inv=Smat_inv, Mdist2=Mdist2
)

#multi-normal test of vec(cmat)
print(mvnorm.etest(t(ExperimentResults_4_MCAR_10pct$thetahat), R = 199))


### MCAR: 0.3 ------

#calculate mustar
Lfun_minvalue_vec <- rep(0,100)
mustar_array <- array(0,dim=c(3,2,100))
for (rept in 1:100){
  data <- generate_data_2026(n=10000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MCAR",MCARmissprob = rep(0.3,2))
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random") )
  Lfun_minvalue_vec[rept] <- res_once$obj_val/(nrow(data$Missing))
  reorderindex <- order(res_once$centers[,1])
  mustar_array[,,rept] <- res_once$centers[reorderindex,]
  if (rept%%10==0){print(rept)}
}
Lfun_minvalue <- mean(Lfun_minvalue_vec)
mustar <- apply(mustar_array,1:2,mean)


est_mu_array <- array(0,dim=c(3,2,1000))
thetahat <- matrix(0,3*2,1000)
for (rept in 1:1000){
  data <- generate_data_2026(n=2000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MCAR",MCARmissprob = rep(0.3,2))
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random")  )
  est_mu_array[,,rept] <- align_centers(res_once$centers,mustar)
  thetahat[,rept] <- sqrt(nrow(data$Missing))*( as.vector(t(est_mu_array[,,rept] - mustar)) )
  if (rept%%100==0){print(rept)}
}
#Malahabios distance^2
Smat_inv <- solve(cov(t(thetahat)))
Mdist2 <- apply(thetahat, 2, function(z) (t(z))%*%Smat_inv%*%z )

ExperimentResults_4_MCAR_30pct <- list(
  Lfun_minvalue_vec=Lfun_minvalue_vec, mustar_array=mustar_array,
  Lfun_minvalue=Lfun_minvalue, mustar=mustar,
  est_mu_array=est_mu_array, thetahat=thetahat,
  Smat_inv=Smat_inv, Mdist2=Mdist2
)

#multi-normal test of vec(cmat)
print(mvnorm.etest(t(ExperimentResults_4_MCAR_30pct$thetahat), R = 199))


### MCAR: 0.5 ------

#calculate mustar
Lfun_minvalue_vec <- rep(0,100)
mustar_array <- array(0,dim=c(3,2,100))
for (rept in 1:100){
  data <- generate_data_2026(n=10000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MCAR",MCARmissprob = rep(0.5,2))
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random") )
  Lfun_minvalue_vec[rept] <- res_once$obj_val/(nrow(data$Missing))
  reorderindex <- order(res_once$centers[,1])
  mustar_array[,,rept] <- res_once$centers[reorderindex,]
  if (rept%%10==0){print(rept)}
}
Lfun_minvalue <- mean(Lfun_minvalue_vec)
mustar <- apply(mustar_array,1:2,mean)


est_mu_array <- array(0,dim=c(3,2,1000))
thetahat <- matrix(0,3*2,1000)
for (rept in 1:1000){
  data <- generate_data_2026(n=2000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MCAR",MCARmissprob = rep(0.5,2))
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random")  )
  est_mu_array[,,rept] <- align_centers(res_once$centers,mustar)
  thetahat[,rept] <- sqrt(nrow(data$Missing))*( as.vector(t(est_mu_array[,,rept] - mustar)) )
  if (rept%%100==0){print(rept)}
}
#Malahabios distance^2
Smat_inv <- solve(cov(t(thetahat)))
Mdist2 <- apply(thetahat, 2, function(z) (t(z))%*%Smat_inv%*%z )

ExperimentResults_4_MCAR_50pct <- list(
  Lfun_minvalue_vec=Lfun_minvalue_vec, mustar_array=mustar_array,
  Lfun_minvalue=Lfun_minvalue, mustar=mustar,
  est_mu_array=est_mu_array, thetahat=thetahat,
  Smat_inv=Smat_inv, Mdist2=Mdist2
)


#multi-normal test of vec(cmat)
print(mvnorm.etest(t(ExperimentResults_4_MCAR_50pct$thetahat), R = 199))


### MNAR0: 0.1 ------

#calculate mustar
Lfun_minvalue_vec <- rep(0,100)
mustar_array <- array(0,dim=c(3,2,100))
for (rept in 1:100){
  data <- generate_data_2026(n=10000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MNAR0",MNAR0lambda_star = 10)
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random") )
  Lfun_minvalue_vec[rept] <- res_once$obj_val/(nrow(data$Missing))
  reorderindex <- order(res_once$centers[,1])
  mustar_array[,,rept] <- res_once$centers[reorderindex,]
  if (rept%%10==0){print(rept)}
}
Lfun_minvalue <- mean(Lfun_minvalue_vec)
mustar <- apply(mustar_array,1:2,mean)


est_mu_array <- array(0,dim=c(3,2,1000))
thetahat <- matrix(0,3*2,1000)
for (rept in 1:1000){
  data <- generate_data_2026(n=2000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MNAR0",MNAR0lambda_star = 10)
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random")  )
  est_mu_array[,,rept] <- align_centers(res_once$centers,mustar)
  thetahat[,rept] <- sqrt(nrow(data$Missing))*( as.vector(t(est_mu_array[,,rept] - mustar)) )
  if (rept%%100==0){print(rept)}
}
#Malahabios distance^2
Smat_inv <- solve(cov(t(thetahat)))
Mdist2 <- apply(thetahat, 2, function(z) (t(z))%*%Smat_inv%*%z )

ExperimentResults_4_MNAR0_10pct <- list(
  Lfun_minvalue_vec=Lfun_minvalue_vec, mustar_array=mustar_array,
  Lfun_minvalue=Lfun_minvalue, mustar=mustar,
  est_mu_array=est_mu_array, thetahat=thetahat,
  Smat_inv=Smat_inv, Mdist2=Mdist2
)

#multi-normal test of vec(cmat)
print(mvnorm.etest(t(ExperimentResults_4_MNAR0_10pct$thetahat), R = 199))


### MNAR0: 0.3 ------

#calculate mustar
Lfun_minvalue_vec <- rep(0,100)
mustar_array <- array(0,dim=c(3,2,100))
for (rept in 1:100){
  data <- generate_data_2026(n=10000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MNAR0",MNAR0lambda_star = 0.8)
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random") )
  Lfun_minvalue_vec[rept] <- res_once$obj_val/(nrow(data$Missing))
  reorderindex <- order(res_once$centers[,1])
  mustar_array[,,rept] <- res_once$centers[reorderindex,]
  if (rept%%10==0){print(rept)}
}
Lfun_minvalue <- mean(Lfun_minvalue_vec)
mustar <- apply(mustar_array,1:2,mean)


est_mu_array <- array(0,dim=c(3,2,1000))
thetahat <- matrix(0,3*2,1000)
for (rept in 1:1000){
  data <- generate_data_2026(n=2000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MNAR0",MNAR0lambda_star = 0.8)
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random")  )
  est_mu_array[,,rept] <- align_centers(res_once$centers,mustar)
  thetahat[,rept] <- sqrt(nrow(data$Missing))*( as.vector(t(est_mu_array[,,rept] - mustar)) )
  if (rept%%100==0){print(rept)}
}
#Malahabios distance^2
Smat_inv <- solve(cov(t(thetahat)))
Mdist2 <- apply(thetahat, 2, function(z) (t(z))%*%Smat_inv%*%z )

ExperimentResults_4_MNAR0_30pct <- list(
  Lfun_minvalue_vec=Lfun_minvalue_vec, mustar_array=mustar_array,
  Lfun_minvalue=Lfun_minvalue, mustar=mustar,
  est_mu_array=est_mu_array, thetahat=thetahat,
  Smat_inv=Smat_inv, Mdist2=Mdist2
)

#multi-normal test of vec(cmat)
print(mvnorm.etest(t(ExperimentResults_4_MNAR0_30pct$thetahat), R = 199))


### MNAR0: 0.5 ------

#calculate mustar
Lfun_minvalue_vec <- rep(0,100)
mustar_array <- array(0,dim=c(3,2,100))
for (rept in 1:100){
  data <- generate_data_2026(n=10000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MNAR0",MNAR0lambda_star = 0.3)
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random") )
  Lfun_minvalue_vec[rept] <- res_once$obj_val/(nrow(data$Missing))
  reorderindex <- order(res_once$centers[,1])
  mustar_array[,,rept] <- res_once$centers[reorderindex,]
  if (rept%%10==0){print(rept)}
}
Lfun_minvalue <- mean(Lfun_minvalue_vec)
mustar <- apply(mustar_array,1:2,mean)


est_mu_array <- array(0,dim=c(3,2,1000))
thetahat <- matrix(0,3*2,1000)
for (rept in 1:1000){
  data <- generate_data_2026(n=2000*3,p=2,k=3,centers=mustarstar,variance_vec = rep(1,2),
                             missmechanism = "MNAR0",MNAR0lambda_star = 0.3)
  res_once <- kmmd_2026(x=data$Missing,k=3,
                        initial_centers_array = get_initial_centers_2026(x=data$Missing,k=3,nstart = 30,init_method = "random")  )
  est_mu_array[,,rept] <- align_centers(res_once$centers,mustar)
  thetahat[,rept] <- sqrt(nrow(data$Missing))*( as.vector(t(est_mu_array[,,rept] - mustar)) )
  if (rept%%100==0){print(rept)}
}
#Malahabios distance^2
Smat_inv <- solve(cov(t(thetahat)))
Mdist2 <- apply(thetahat, 2, function(z) (t(z))%*%Smat_inv%*%z )

ExperimentResults_4_MNAR0_50pct <- list(
  Lfun_minvalue_vec=Lfun_minvalue_vec, mustar_array=mustar_array,
  Lfun_minvalue=Lfun_minvalue, mustar=mustar,
  est_mu_array=est_mu_array, thetahat=thetahat,
  Smat_inv=Smat_inv, Mdist2=Mdist2
)

#multi-normal test of vec(cmat)
print(mvnorm.etest(t(ExperimentResults_4_MNAR0_50pct$thetahat[,-668]), R = 199))


###### normality test for each component of cmat ----

for (i in 1:6){  print(shapiro.test(ExperimentResults_4_MCAR_10pct$thetahat[i,])$p.value)  }
for (i in 1:6){  print(shapiro.test(ExperimentResults_4_MCAR_30pct$thetahat[i,])$p.value)  }
for (i in 1:6){  print(shapiro.test(ExperimentResults_4_MCAR_50pct$thetahat[i,])$p.value)  }
for (i in 1:6){  print(shapiro.test(ExperimentResults_4_MNAR0_10pct$thetahat[i,])$p.value) }
for (i in 1:6){  print(shapiro.test(ExperimentResults_4_MNAR0_30pct$thetahat[i,])$p.value) }
for (i in 1:6){  print(shapiro.test(ExperimentResults_4_MNAR0_50pct$thetahat[i,])$p.value) }



