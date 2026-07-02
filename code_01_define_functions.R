###################################
########  main functions    #######
###################################

#packages--------------
library(kpodclustr)
library(MASS)
library(Rfast)
library(Rcpp)
library(RcppArmadillo) 
library(abind)
library(RSKC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(clue)
library(energy)
library(mclust)


sourceCpp("cpp_initialcenters_kmmd_once.cpp")
# define functions ------------------

get_initial_centers_2026 <- function(x,k,nstart=10,
                                     init_method=c("random","kmpp")){
  n=dim(x)[1]
  p=dim(x)[2]
  
  initial_centers_array=array(0,dim = c(k,p,nstart) )
  
  if (init_method[1]=="random"){
    #mean imputation
    x_imputed <- x
    mu_mat <- matrix(colMeans(x,na.rm = TRUE), nrow = n, ncol = p, byrow = TRUE)
    x_imputed[is.na(x)==1] <- mu_mat[is.na(x)==1] 
    #check NA columns in x_imputed 
    NAcols <- which( apply(x_imputed, 2, function(z) sum(is.na(z)) ) == n )
    if (length(NAcols)>0){x_imputed[,NAcols] <- rep(0,n)}
    #generate initial centers
    for (rept in 1:nstart){
      row_index <- sample(1:n, k, replace = FALSE)
      cmat <- x_imputed[row_index,]
      if(any(duplicated(cmat)) == TRUE){ cmat <- cmat + matrix(rnorm(k*p,0,0.01),k,p) }
      initial_centers_array[,,rept] <- cmat
    }
  }
  
  if (init_method[1]=="kmpp"){
    #mean imputation
    x_imputed <- x
    mu_mat <- matrix(colMeans(x,na.rm = TRUE), nrow = n, ncol = p, byrow = TRUE)
    x_imputed[is.na(x)==1] <- mu_mat[is.na(x)==1] 
    #check NA columns in x_imputed 
    NAcols <- which( apply(x_imputed, 2, function(z) sum(is.na(z)) ) == n )
    if (length(NAcols)>0){x_imputed[,NAcols] <- rep(0,n)}
    #generate initial centers
    for (rept in 1:nstart){
      initial_centers_array[,,rept] <- kmpp_init_centers_cpp(x_imputed,k)$centers
    }
  }
  
  return(initial_centers_array=initial_centers_array)
}

kmmd_2026 <- function(x,k,initial_centers_array=NULL,itermax=100,eps=1e-4){
  n=dim(x)[1]
  p=dim(x)[2]
  
  if (is.null(initial_centers_array)){
    initial_centers_array <- get_initial_centers_2026(x,k,nstart=10,init_method = "kmpp")
  }
  nstart=dim(initial_centers_array)[3]
  
  obj_old=Inf
  obj_nstart=rep(0,nstart)
  for (start in 1:nstart){
    #run for each initialization
    initial_centers=initial_centers_array[,,start]
    tmp_res=kmmd_once_2026cpp(x=x,k=k,init_centers=initial_centers, itermax = itermax,eps=eps)
    #record
    obj=tmp_res$obj_val
    obj_nstart[start]=obj
    #compare
    if(is.finite(obj) && (obj<obj_old)){
      obj_old=obj
      result=tmp_res
    }
  }
  result$obj_each_initial=obj_nstart
  return(result)
}




generate_data_2026 <- function(n,p,k,centers=matrix(0,k,p),variance_vec=rep(1,p),
                               missmechanism=c("MCAR","MAR","MNAR1","MNAR2","MNAR0"),
                               MCARmissprob=rep(0,p),
                               MARpsi=c(0.1,0.1),
                               MNAR1phi=c(3,2),
                               MNAR2percent=rep(0.1,p),
                               MNAR0lambda_star=0.1,
                               nomissing=FALSE                    ){
  ##original data
  x <- matrix(0, nrow = n, ncol = p)
  clusterlabel <- sample(1:k, n, replace = TRUE)
  for (kk in 1:k){
    kk_cluster_size <- sum(clusterlabel==kk)
    x_kk_cluster <- MASS::mvrnorm(kk_cluster_size, mu=centers[kk,], Sigma=diag(variance_vec) ) 
    x[which(clusterlabel==kk),] <- x_kk_cluster
  }
  
  if(nomissing==TRUE){
    return(list( Orig=x, Missing=NULL, CompleteCase=x, 
                 Origlabel=clusterlabel, Misslabel=NULL, Complabel=clusterlabel  ))
  }
  
  ##missingness
  
  if (missmechanism[1] == "MCAR"){
    miss_mat <- matrix(0,n,p)
    for (j in 1:p){
      miss_mat[,j] <- rbinom(n, 1, prob=MCARmissprob[j])
    }
  }
  
  if (missmechanism[1] == "MAR"){
    missprob=rep(0,n)
    miss_mat=matrix(0,n,p)
    for (i in 1:n){
      missprob[i]=1/(1+exp(- (MARpsi[1])*( x[i,1] - (MARpsi[2]) ) ) )
      for (j in 2:p){
        miss_mat[i,j]=rbinom(n=1,size=1,prob=missprob[i])
      }
    }
  }
  
  if (missmechanism[1] == "MNAR1"){
    missprob=1/(1+exp(- (MNAR1phi[1])*( x - (MNAR1phi[2]) ) ) )
    miss_mat=matrix(0,n,p)
    for (j in 1:p){
      for (i in 1:n){
        miss_mat[i,j]=rbinom(n=1,size=1,prob=missprob[i,j])
      }
    }
  }
  
  if (missmechanism[1] == "MNAR2"){
    miss_mat=matrix(0,n,p)
    for (j in 1:p){
      threshold=quantile(abs(x[,j]),probs=MNAR2percent[j])
      miss_mat[,j]=( abs(x[,j])< threshold )
    }
  }
  
  if (missmechanism[1] == "MNAR0"){
    missprob=exp(-MNAR0lambda_star*(x^2) )
    miss_mat=matrix(0,n,p)
    for (j in 1:p){
      for (i in 1:n){
        miss_mat[i,j]=rbinom(n=1,size=1,prob=missprob[i,j])
      }
    }
  }
  
  ##missing data and complete data
  x_tmp <- x
  x_tmp[miss_mat==1] <- NA
  nullrows <- which( (Rfast::rowsums(miss_mat)) == p )
  comprows <- which( (Rfast::rowsums(miss_mat)) == 0 )
  is_null_row <- ( (Rfast::rowsums(miss_mat)) == p )
  
  if(length(nullrows)==0){
    x_miss <- x_tmp
    clusterlabel_miss <- clusterlabel
  }else{
    x_miss <- x_tmp[-nullrows,]
    clusterlabel_miss <- clusterlabel[-nullrows]
  }
  
  if(length(comprows)==0){
    x_comp <- NULL
    clusterlabel_comp <- NULL
  }else{
    x_comp <- x_tmp[comprows,]
    clusterlabel_comp <- clusterlabel[comprows]
  }
  
  return(list(Orig=x, Missing=x_miss, CompleteCase=x_comp, 
              Origlabel=clusterlabel, Misslabel=clusterlabel_miss, Complabel=clusterlabel_comp,
              is_null_row=is_null_row  ))
}



MSEcenters<-function(codebook0,codebook){
  k=dim(codebook0)[1]
  p=dim(codebook0)[2]
  
  mse=rep(0,k)
  for (kk in 1:k){
    dist=apply(codebook0,1,function(z) sum((codebook[kk,]-z)^2)  )
    mse[k]=min(dist)
  }
  return(sum(mse))
}


align_centers <- function(est_centers, true_centers){
  k <- nrow(true_centers)
  dist_mat <- matrix(0, k, k)
  for (kk in 1:k){
    for (kkk in 1:k){
      dist_mat[kk,kkk] <- sqrt( sum( (est_centers[kk,] - true_centers[kkk,])^2 ) )
    }
  }
  mapping <- solve_LSAP(dist_mat)   
  new_order <- order(mapping)       
  return(est_centers[new_order,])
}


generate_truncated_data_2016 <- function(n=10000, p=2, k=2, centers=matrix(0,k,p),
                                         variance_vec=rep(1,p), beta_starstar=3,
                                         missmechanism=c("MCAR","MAR","MNAR1","MNAR2","MNAR0"),
                                         MCARmissprob=rep(0,p),
                                         MARpsi=c(0.1,0.1),
                                         MNAR1phi=c(3,2),
                                         MNAR2percent=rep(0.1,p),
                                         MNAR0lambda_star=0.1,
                                         nomissing=FALSE            ){
  
  ##original data
  x <- matrix(0, nrow = n, ncol = p)
  clusterlabel <- sample(1:k, n, replace = TRUE)
  for (kk in 1:k){
    kk_cluster_size <- sum(clusterlabel==kk)
    x_kk_cluster_tmp <- MASS::mvrnorm(n, mu=centers[kk,], Sigma=diag(variance_vec) )
    tmp_dist_l2 <- apply(x_kk_cluster_tmp, 1, function(z) sqrt(sum((z - centers[kk,])^2)) )
    x_kk_cluster_keep <- x_kk_cluster_tmp[tmp_dist_l2 <= beta_starstar, ]
    if (nrow(x_kk_cluster_keep) >= kk_cluster_size){
      x[which(clusterlabel==kk),] <- x_kk_cluster_keep[sample(1:nrow(x_kk_cluster_keep), kk_cluster_size, replace=FALSE),]
    }else{
      x[which(clusterlabel==kk),] <- rbind(x_kk_cluster_keep,
                                           matrix(NA, kk_cluster_size-nrow(x_kk_cluster_keep), p))
    }
  }
  zero_rows <- which( is.na(rowSums(x) ) )
  if (length(zero_rows)>0){
    x <- x[-zero_rows,]
  }
  
  if(nomissing==TRUE){
    return(list( Orig=x, Missing=NULL, CompleteCase=x, 
                 Origlabel=clusterlabel, Misslabel=NULL, Complabel=clusterlabel  ))
  }
  
  ##missingness
  
  if (missmechanism[1] == "MCAR"){
    miss_mat <- matrix(0,n,p)
    for (j in 1:p){
      miss_mat[,j] <- rbinom(n, 1, prob=MCARmissprob[j])
    }
  }
  
  if (missmechanism[1] == "MAR"){
    missprob=rep(0,n)
    miss_mat=matrix(0,n,p)
    for (i in 1:n){
      missprob[i]=1/(1+exp(- (MARpsi[1])*( x[i,1] - (MARpsi[2]) ) ) )
      for (j in 2:p){
        miss_mat[i,j]=rbinom(n=1,size=1,prob=missprob[i])
      }
    }
  }
  
  if (missmechanism[1] == "MNAR1"){
    missprob=1/(1+exp(- (MNAR1phi[1])*( x - (MNAR1phi[2]) ) ) )
    miss_mat=matrix(0,n,p)
    for (j in 1:p){
      for (i in 1:n){
        miss_mat[i,j]=rbinom(n=1,size=1,prob=missprob[i,j])
      }
    }
  }
  
  if (missmechanism[1] == "MNAR2"){
    miss_mat=matrix(0,n,p)
    for (j in 1:p){
      threshold=quantile(abs(x[,j]),probs=MNAR2percent[j])
      miss_mat[,j]=( abs(x[,j])< threshold )
    }
  }
  
  if (missmechanism[1] == "MNAR0"){
    missprob=exp(-MNAR0lambda_star*(x^2) )
    miss_mat=matrix(0,n,p)
    for (j in 1:p){
      for (i in 1:n){
        miss_mat[i,j]=rbinom(n=1,size=1,prob=missprob[i,j])
      }
    }
  }
  
  ##missing data and complete data
  x_tmp <- x
  x_tmp[miss_mat==1] <- NA
  nullrows <- which( (Rfast::rowsums(miss_mat)) == p )
  comprows <- which( (Rfast::rowsums(miss_mat)) == 0 )
  is_null_row <- ( (Rfast::rowsums(miss_mat)) == p )
  
  if(length(nullrows)==0){
    x_miss <- x_tmp
    clusterlabel_miss <- clusterlabel
  }else{
    x_miss <- x_tmp[-nullrows,]
    clusterlabel_miss <- clusterlabel[-nullrows]
  }
  
  if(length(comprows)==0){
    x_comp <- NULL
    clusterlabel_comp <- NULL
  }else{
    x_comp <- x_tmp[comprows,]
    clusterlabel_comp <- clusterlabel[comprows]
  }
  
  return(list(Orig=x, Missing=x_miss, CompleteCase=x_comp, 
              Origlabel=clusterlabel, Misslabel=clusterlabel_miss, Complabel=clusterlabel_comp,
              is_null_row=is_null_row  ))
}
