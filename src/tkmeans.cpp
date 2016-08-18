//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//'@importFrom Rcpp sourceCpp
//'@useDynLib tkmeans



// constructing priority queues
#include <iostream>
#include <queue>
#include <vector>
#include <functional>
#include <tuple>
// maths
#include <math.h>


const double log2pi = std::log(2.0 * M_PI);

typedef std::tuple<double, int, int>  di_pair;

class CompareDist
{
public:
  bool operator()(di_pair n1, di_pair n2) {
    return std::get<0>(n1)<std::get<0>(n2);
  }
};




arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov) {
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);
}



arma::vec dmvnorm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool log = false) {
  arma::vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;

  if (log) {
    return(logretval);
  } else {
    return(exp(logretval));
  }
}



arma::vec distCentre(int k, arma::mat ui, arma::mat centres, double lambda, int d){
  //set up sigma
  arma::vec sigmaDiag(d);  sigmaDiag.fill(lambda);
  arma::mat sigma = diagmat(sigmaDiag);

   //store distances
  arma::vec dist =  arma::zeros(k);

  for (int j=0; j<k; j++)
  {
    dist.at(j) = 1.0/k * dmvnorm_arma(ui, centres.row(j), sigma, false)[0];
  }
  //Rcpp::Rcout << "size: " << size(dist) << std::endl;
  return(dist);
}


arma::vec distCentre2(int k, arma::mat ui, arma::mat centres, double lambda, int d){

  //store distances
  arma::vec dist =  arma::zeros(k);

  for (int j=0; j<k; j++)
  {
    dist.at(j) = sqrt(arma::sum(arma::pow(centres.row(j) - ui, 2.0)));
  }


  return(dist);
}


arma::vec distCentre3(int k, arma::mat ui, arma::mat centres, double lambda, int d){
  arma::vec dist = arma::sqrt(arma::sum(arma::pow(centres - repmat(ui,k,1), 2.0),1));
  return(dist);
}



int max_index(arma::vec x){

  //returns index of biggest value in x.
  Rcpp::NumericVector y = NumericVector(x.begin(),x.end());
  Rcpp::NumericVector::iterator it =  std::max_element(y.begin(), y.end());

  //Rcpp::Rcout << "it  " << it <<  " y " << y.begin() <<  " diff " << std::distance(y.begin(), it) << std::endl;
  return (std::distance(y.begin(), it));
}



int min_index(arma::vec x){

  //returns index of biggest value in x.
  Rcpp::NumericVector y = NumericVector(x.begin(),x.end());
  Rcpp::NumericVector::iterator it =  std::min_element(y.begin(), y.end());

  //Rcpp::Rcout << "it  " << it <<  " y " << y.begin() <<  " diff " << std::distance(y.begin(), it) << std::endl;
  return (std::distance(y.begin(), it));
}



arma::mat init_centres(arma::mat M, int k){
  int n = size(M)[0];
  int d = size(M)[1];

  //randomly pick centres
  arma::mat B = arma::randi<arma::mat>(k,1,arma::distr_param(0,n));

  arma::mat cent =  arma::zeros(k,d);

  for(int i=0; i<k; i++){
    //Rcpp::Rcout << "i = " << i << " n = " << B.at(i,0) <<  std::endl;

    cent.row(i) = M.row(B.at(i,0));
  }

  return cent;
}


//'@title Trimmed k-means clustering
//'@description
//'Performs trimmed k-means clustering algorithm on a matrix of data. Each row is an observation.
//'@details
//'k is the number of clusters. Alpha is the proportion of data that will be excluded in the clustering.
//'
//'Algorithm will halt if either maximum number of iterations is reached or the change between iterations drops below tol.
//'
//'When n_starts is greater than 1, the algorithm will run multiple times and the result with the best BIC will be returned.
//'The centres are intialised by picking k observations.
//'
//'The function only returns the k cluster centres. To calculate the nearest cluster centre for each observation use the function \code{nearest_cluster}.
//'
//'@param M matrix (n x m). Rows are observations, columns are predictors.
//'@param k number of clusters
//'@param alpha proportion of data to be trimmed
//'@param nstart number of restarts
//'@param iter maximum number of iterations
//'@param tol criteria for algorithm convergence
//'@param verbose If true will output more information on algorithm progress.
//'@return matrix of cluster means (k x m).
//'@examples
//'iris_mat <- as.matrix(iris[,1:4])
//'iris_centres2 <- tkmeans(iris_mat, 2 , 0.1, 1, 10, 0.001) # 2 clusters
//'@export
// [[Rcpp::export]]
arma::mat tkmeans(arma::mat& M, int k , double alpha, int nstart = 1, int iter = 10, double tol = 0.0001, bool verbose = false){
  int n = size(M)[0];
  int d = size(M)[1];

  if(alpha>=1 | alpha <=0){
    stop("alpha must be in (0,1).");
  }

  if(nstart<1){
    stop("Number of starts must be at least 1.");
  }

  if(iter<5){
    stop("Maximum number of iterations is too low.");
  }

  if(verbose){
    Rcpp::Rcout << "n = " << n << " d = " << d <<  " k = " << k << std::endl;
    Rcpp::Rcout << "n starts = " << nstart << " max iterations " << iter <<  " alpha = " << alpha << std::endl;
  }
  arma::mat centres = init_centres(M, k);
  arma::mat means = centres;

  double lambda = 1;
  //double tol = 0.001;
  double diff = 1;
  int m=0;

  // cluster membership record
  //std::vector<int> cluster_membership;
  //cluster_membership.reserve(n);

  //keep list of smallest using priority queue
  std::priority_queue< di_pair, std::vector<di_pair>, CompareDist  >  smallest;

  int queue_len;
  if(alpha<=0.5){

    if(alpha != 0.0){
      queue_len = floor((n-1)*(alpha))+1;
    }else{
      queue_len = 0;
    }


    //Rcpp::Rcout << "number of outliers = " <<  queue_len   << std::endl;
    //Rcpp::Rcout << "Begin main loop..."  << std::endl;


    while(m<iter & diff > tol){

      //Rcpp::Rcout << m << " of " << iter  << " iterations" << std::endl;

      arma::mat new_centres = arma::zeros(k, d);
      arma::mat centre_members = arma::zeros(k, 1);

      for(int i=0; i<n; i++)
      {
        arma::mat Dc = distCentre2(k, M.row(i), centres, lambda, d);
        //distCentre(k, M.row(i), centres, lambda, d).print();

        int cluster = min_index(Dc);

        di_pair temp = di_pair(Dc.at(cluster), cluster, i);

        new_centres.row(cluster) =  new_centres.row(cluster)+ M.row(i);
        centre_members(cluster,0) =  centre_members(cluster,0) + 1;

        if(queue_len > 0){
          if(smallest.size() < queue_len){
            smallest.push(temp);
          }else{
            if(std::get<0>(smallest.top()) > std::get<0>(temp)){
              //Rcpp::Rcout << "out: "   <<  std::get<0>(smallest.top()) <<  " in: " << std::get<0>(temp) << std::endl;
              smallest.pop();
              smallest.push(temp);
            }
          }
        }
        //Rcpp::Rcout << "size of heap:"  <<  smallest.size() << std::endl;
      }


      //remove all the outliers
      if(queue_len > 0){

        while (!smallest.empty())
        {
          new_centres.row(std::get<1>(smallest.top())) = new_centres.row(std::get<1>(smallest.top())) - M.row(std::get<2>(smallest.top()));

          centre_members(std::get<1>(smallest.top()),0) = centre_members(std::get<1>(smallest.top()),0) - 1;

          //Rcpp::Rcout << " " << std::get<0>(smallest.top()) <<  " " << std::get<1>(smallest.top()) << " " << std::get<2>(smallest.top()) << std::endl;

          smallest.pop();
        }
      }

      if(centre_members.min() > 0) {
        means = new_centres.each_col() / centre_members;
      }else{
        stop("Empty cluster");
      }

      //centre_members.print();
      //Rcpp::Rcout << diff << std::endl;

      diff = arma::accu(arma::abs(centres-means));
      centres = means;
      m++;
    }
  }
  else{
    if(alpha != 0.0){
      queue_len = floor((n-1)*(1-alpha))+1;
    }else{
      queue_len = 0;
    }


    //Rcpp::Rcout << "number of outliers = " <<  queue_len   << std::endl;
    //Rcpp::Rcout << "Begin main loop..."  << std::endl;


    while(m<iter & diff > tol){

      //Rcpp::Rcout << m << " of " << iter  << " iterations" << std::endl;

      arma::mat new_centres = arma::zeros(k, d);
      arma::mat centre_members = arma::zeros(k, 1);

      for(int i=0; i<n; i++)
      {
          arma::mat Dc = distCentre2(k, M.row(i), centres, lambda, d);
          //distCentre(k, M.row(i), centres, lambda, d).print();

          int cluster = min_index(Dc);

          di_pair temp = di_pair(Dc.at(cluster), cluster, i);

          //new_centres.row(cluster) =  new_centres.row(cluster)+ M.row(i);
          //centre_members(cluster,0) =  centre_members(cluster,0) + 1;

          if(queue_len > 0){
            if(smallest.size() < queue_len){
               smallest.push(temp);
            }else{
              if(std::get<0>(smallest.top()) > std::get<0>(temp)){
                //Rcpp::Rcout << "out: "   <<  std::get<0>(smallest.top()) <<  " in: " << std::get<0>(temp) << std::endl;
                smallest.pop();
                smallest.push(temp);
              }
            }
          }
          //Rcpp::Rcout << "size of heap:"  <<  smallest.size() << std::endl;
      }


    //remove all the outliers
      if(queue_len > 0){
          while (!smallest.empty())
          {
            new_centres.row(std::get<1>(smallest.top())) = new_centres.row(std::get<1>(smallest.top())) + M.row(std::get<2>(smallest.top()));

            centre_members(std::get<1>(smallest.top()),0) = centre_members(std::get<1>(smallest.top()),0) + 1;

            //Rcpp::Rcout << " " << std::get<0>(smallest.top()) <<  " " << std::get<1>(smallest.top()) << " " << std::get<2>(smallest.top()) << std::endl;

            smallest.pop();
          }
      }

      if(centre_members.min() > 0) {
          means = new_centres.each_col() / centre_members;
      }else{
         stop("Empty cluster");
      }

      //centre_members.print();
      //Rcpp::Rcout << diff << std::endl;

      diff = arma::accu(arma::abs(centres-means));
      centres = means;
      m++;
    }
  }
  return means;
}

//'@title Rescales a matrix in place.
//'@description
//'Recales matrix so that each col to have mean 0 and sd 1.
//'The original matrix is overwritten in place. The function returns the means and standard deviations of each column used to rescale it.
//'@details
//'The key advantage of this method is that it can be applied to very matrices without having to make a second copy in memory and the orginal can be restored using the saved information.
//'
//'@param M matrix of data (n x m)
//'@return Returns a matrix of size 2 x m. The first row contains the means. The second row  contains the sd. Scales original matrix in place.
//'@examples
//'m = matrix(rnorm(24, 1, 2),4, 6)
//'scale_params = scale_mat_inplace(m)
//'sweep(sweep(m,2,scale_params[2,],'*'),2,scale_params [1,], '+') # orginal matrix restored
//'@export
// [[Rcpp::export]]
arma::mat scale_mat_inplace(arma::mat& M){
  arma::rowvec means = mean(M);
  arma::rowvec sds = stddev(M);

  for(int i=0; i< M.n_cols; i++){
    M.col(i) -= means.at(i);
    M.col(i) /= sds.at(i);
  }

  arma::mat means_sds(2,means.n_elem);
  means_sds.row(0) = means;
  means_sds.row(1) = sds;

  return means_sds;
  }


//'@title Calculates BIC for a given clustering.
//'@description
//'Computes Bayesian information criterion for a given clustering of a data set.
//'@details
//'Bayesian information criterion (BIC) is calculated using the formula, BIC =  -2 * log(L) + k*log(n).
//'k is the number of free parameters, in this case is m*k + k - 1.
//'n is the number of observations (rows of data).
//'L is the liklihood for the given set of cluster centres.
//'
//'@param data a matrix (n x m). Rows are observations, columns are predictors.
//'@param centres matrix of cluster means (k x m), where k is the number of clusters.
//'@return BIC value
//'@examples
//'iris_mat <- as.matrix(iris[,1:4])
//'iris_centres2 <- tkmeans(iris_mat, 2 , 0.1, 1, 10, 0.001) # 2 clusters
//'iris_centres3 <- tkmeans(iris_mat, 3 , 0.1, 1, 10, 0.001) # 3 clusters
//'cluster_BIC(iris_mat, iris_centres2)
//'cluster_BIC(iris_mat, iris_centres3)
//'@export
// [[Rcpp::export]]
double  cluster_BIC(arma::mat& data, arma::mat& centres){
          int x = data.n_rows;
          int m = centres.n_cols;
          int k = centres.n_rows;
          double PI_val = log(1./k);


          if(data.n_cols!=m){
            stop("Cluster centre dimensionality does not match data.");
          }

          arma::mat temp_row(1, k);


          double log_like_accum = 0.;

          for(int i=0;i<x;i++){

            for(int j=0;j<k;j++){

              temp_row.at(0,j) =

                -0.5*(arma::accu(arma::pow(data.row(i),2.))
                + arma::accu(arma::pow(centres.row(j),2.))
                - 2.*dot(data.row(i),centres.row(j)))
                - (m/2.)*log(2.*M_PI) - log(m)/2. + PI_val;

            }

            double temp_max = temp_row.max();


            temp_row = temp_row - temp_max;

            double log_like = temp_max + log(arma::accu(arma::exp(temp_row)));
            //Rcpp::Rcout <<  log_like << std::endl;
            log_like_accum += log_like;
          }

          //Rcpp::Rcout <<  log_like_accum << std::endl;


  return -2*log_like_accum  + log(x)*(m*k + k - 1);
}




//'@title Allocates each rw (observation) in data to the nearest cluster centre.
//'@description
//'For each observation the euclidean distance to each of the cluster centres is calculated and cluster with the smallest distance is return for that observation.
//'@param data a matrix (n x m) to be clustered
//'@param centres matrix of cluster means (k x m), wher k is the number of clusters.
//'@return vector of cluster allocations, n values ranging from 1 to k.
//'@examples
//'iris_mat <- as.matrix(iris[,1:4])
//'centres<- tkmeans(iris_mat, 3 , 0.2, 1, 10, 0.001)
//' nearest_cluster(iris_mat, centres)
//'@export
// [[Rcpp::export]]
arma::uvec nearest_cluster(arma::mat& data, arma::mat& centres){
  int k =  centres.n_rows;
  double lambda = 1.0;
  int d = size(data)[0];

  if(data.n_cols!=centres.n_cols){
    stop("Cluster centre dimensionality does not match data.");
  }

  arma::uvec clusters = arma::zeros<arma::uvec>(d);

  for(int i=0;i<d;i++){
    arma::mat Dc = distCentre2(k, data.row(i), centres, lambda, d);
    int cluster = min_index(Dc);
    clusters.at(i) = cluster+1;//going back to 1 indexed R
  }

  return clusters;
}

