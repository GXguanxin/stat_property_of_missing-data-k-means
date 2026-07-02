// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// [[Rcpp::export]]
List kmpp_init_centers_cpp(const arma::mat xmat, const int k){
  int n = xmat.n_rows, p = xmat.n_cols;
  double tmp_dist = 0;
  arma::mat init_cmat = arma::zeros(k,p);
  arma::vec all_index = arma::regspace(0,n-1), tmp_dist_vec = arma::ones(n), init_index(1);
  //arma::vec init_centers_index_vec(k);
  for(int kk=0; kk<k; ++kk){
    //sample kk-th center with probability=tmp_dist_vec_copy
    init_index = RcppArmadillo::sample(all_index, 1, false, tmp_dist_vec);
    init_cmat(kk,arma::span::all) = xmat(init_index(0),arma::span::all);
    
    for(int i = 0; i < n; ++i){//update weights for sampling
      tmp_dist = sum( pow(xmat(i,arma::span::all) - init_cmat(kk,arma::span::all), 2) );
      if(kk==0){
        tmp_dist_vec(i) = tmp_dist;
      }else if(tmp_dist < tmp_dist_vec(i)){
        tmp_dist_vec(i) = tmp_dist;
      }
    }
  }
  
  List res;
  res["centers"] = init_cmat;
  
  return res;
}


// [[Rcpp::export]]
List kmmd_once_2026cpp(const arma::mat x, const int k,
                        const arma::mat init_centers, const int itermax, const double eps){
  int n = x.n_rows, p = x.n_cols;
  arma::mat missing = arma::conv_to<arma::mat>::from(x != x); // (NAN != NAN) = TRUE，(NAN == NAN) = FALSE
  arma::mat observed = 1.0 - missing;
  
  double obj_old = 1e-8, obj_new = 0, term1 = 0;
  int iter_end = 0;
  arma::uvec cluster(n);
  arma::rowvec masked_center(p), numerator(p), denominator(p);
  arma::vec obj_list(itermax); 
  arma::mat cmat= init_centers, cmat_old= init_centers, dist_mat = arma::zeros(n,k), ab_mat = arma::zeros(n,p);
  arma::umat label_list(n, itermax, arma::fill::zeros);
  arma::cube cmat_list = arma::zeros(k,p,itermax);
  
  arma::mat Pomega_x = x;
  Pomega_x.replace(arma::datum::nan, 0);
  
  for (int iter = 0; iter < itermax; ++iter){
    
    //given tmp_cmat, estimate tmp_cluster
    dist_mat = arma::zeros(n,k);
    for (int i = 0; i < n; ++i){
      for (int kk = 0; kk < k; ++kk){
        //partial distance between x_i and mu_kk on observed dimensions
        masked_center = cmat(kk,arma::span::all)%observed(i,arma::span::all);
        dist_mat(i,kk) = sum( pow( Pomega_x(i,arma::span::all) - masked_center , 2) ); 
      }
    }
    cluster = arma::index_min(dist_mat,1);
    
    //update cmat
    cmat = arma::zeros(k,p);
    for (int kk = 0; kk < k; ++kk){
      arma::uvec rowindex = arma::find(cluster == kk); 
      if (rowindex.n_elem > 0){
        numerator = sum( Pomega_x.rows( rowindex ) , 0);
        denominator = sum( observed.rows( rowindex ) , 0);
        for(int j=0; j<p; ++j) {
          if(denominator(j) > 0){ cmat(kk, j) = numerator(j)/denominator(j); }
          else{ cmat(kk, j) = cmat_old(kk, j); }
        }
      }else{ cmat.row(kk) = cmat_old.row(kk); }
    }
    
    //record
    ab_mat = cmat.rows(cluster);
    term1 = arma::accu( pow( (Pomega_x - ab_mat)%observed ,2) );
    obj_new = term1;
    
    //stopping
    if(iter > 0){
      if( abs(obj_new-obj_old)/abs(obj_old) < eps ){
        obj_list(iter)=obj_new;
        label_list.col(iter)=cluster;
        cmat_list.slice(iter)=cmat;
        iter_end=iter;
        break;}
    }
    //update
    obj_old=obj_new;
    iter_end=iter;
    cmat_old = cmat;
    obj_list(iter)=obj_new;
    label_list.col(iter)=cluster;
    cmat_list.slice(iter)=cmat;
  }
  
  
  List result;
  result["cluster"] = cluster+1;
  result["centers"] = cmat;
  result["obj_val"] = obj_new;
  result["term1"] = term1;
  result["obj_list"] = obj_list.subvec(0,iter_end);
  result["label_list"] = label_list.cols(0,iter_end)+1;
  result["cmat_list"] = cmat_list.slices(0,iter_end);
  return(result);
}


