#' Build the single-particle Hamiltonian for impurity + finite bath
#'
#'@description Analyze differential scanning calorimetry (DSC) data, detect thermal transitions and calculate thermodynamic parameters.
#'
#'
#' @param eps_d impurity level (numeric)
#' @param eps_k numeric vector of bath levels
#' @param V_k numeric vector of hybridizations (same length as eps_k)
#' @return single-particle Hamiltonian matrix of size (1 + length(eps_k))
#' @export

build_hamiltonian <- function(eps_d, eps_k, V_k){
  if(length(eps_k) != length(V_k)) stop("eps_k and V_k must have same length")
  n_bath <- length(eps_k)
  N <- 1 + n_bath
  H0 <- matrix(0, nrow = N, ncol = N)
  # index 1 is impurity
  H0[1,1] <- eps_d
  if(n_bath > 0){
    H0[2:(N), 2:(N)] <- diag(eps_k)
    # hybridization
    H0[1, 2:(N)] <- V_k
    H0[2:(N), 1] <- V_k
  }
  return(H0)
}


#' Matrix exponential via eigen decomposition (small matrices)
#' @param A detail
#' @param dt detail
matrix_exp_dt <- function(A, dt){
  ev <- eigen(A)
  vals <- ev$values
  vecs <- ev$vectors
  expD <- diag(exp(-dt * vals))
  # reconstruct
  return(Re(vecs %*% expD %*% solve(vecs)))
}