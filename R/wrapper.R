#'@title Trimmed k-means clustering
#'@description
#'Performs trimmed k-means clustering algorithm [1] on a matrix of data. Each row  in the data is an observation, each column is  a variable.
#'For optimal use columns should be scaled to have the same means and variances using \code{scale_mat_inplace}.
#'@details
#'k is the number of clusters. alpha is the proportion of data that will be excluded in the clustering.
#'
#'Algorithm will halt if either maximum number of iterations is reached or the change between iterations drops below tol.
#'
#'When n_starts is greater than 1, the algorithm will run multiple times and the result with the best BIC will be returned.
#'The centres are intialised by picking k observations.
#'
#'The function only returns the k cluster centres. To calculate the nearest cluster centre for each observation use the function \code{nearest_cluster}.
#'
#'@param M matrix (n x m). Rows are observations, columns are predictors.
#'@param k number of clusters
#'@param alpha proportion of data to be trimmed
#'@param weights weightings for variables (columns).
#'@param nstart number of restarts
#'@param iter maximum number of iterations
#'@param tol criteria for algorithm convergence
#'@param verbose If true will output more information on algorithm progress.
#'@return Returns a matrix of cluster means (k x m).
#'@references
#' [1] Cuesta-Albertos, J. A., Gordaliza A., and Matran C. "Trimmed K-Means: An Attempt to Robustify Quantizers." The Annals of Statistics 25.2 (1997): 553-76.
#'@examples
#'iris_mat <- as.matrix(iris[,1:4])
#'scale_params<-scale_mat_inplace(iris_mat)
#'iris_cluster<- tkmeans(iris_mat, 2 , 0.1, c(1,1,1,1), 1, 10, 0.001) # 2 clusters
#'@export
tkmeans <- function(M, k, alpha, weights = rep(1,ncol(M)),  nstart = 1L, iter = 10L, tol = 0.0001, verbose = FALSE) {
  .Call('tkmeans_tkmeans', PACKAGE = 'tkmeans', M, k, alpha, weights, nstart, iter, tol, verbose)
}


