#' Hirsch-Fye QMC using Rcpp-accelerated core
#' 
#' @description Analyze differential scanning calorimetry (DSC) data, detect thermal transitions and calculate thermodynamic parameters.
#'
#' @param eps_d impurity level
#' @param eps_k bath energy levels
#' @param V_k hybridization strengths
#' @param U on-site Hubbard interaction
#' @param beta inverse temperature
#' @param L number of time slices
#' @param n_steps Monte Carlo steps
#' @param thermal thermalization steps
#' @param measure_interval measurement interval
#' @param seed RNG seed

#' @useDynLib SC25014049, .registration = TRUE
#' 
#' @export hf_step
#' @importFrom stats runif
#' @importFrom Rcpp sourceCpp
#' @export



simulate_hirsch_fye <- function(eps_d = 0, eps_k = c(-1,1), V_k = c(0.5,0.5), U = 2.0,
                                beta = 5.0, L = 40, n_steps = 20000, thermal = 2000,
                                measure_interval = 10, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  H0 <- build_hamiltonian(eps_d, eps_k, V_k)
  s <- sample(c(-1,1), L, replace = TRUE)
  
  
  meas_nup <- numeric(0)
  meas_ndn <- numeric(0)
  
  
  for(step in seq_len(n_steps)){
    l <- sample.int(L, 1)
    s_new <- s; s_new[l] <- -s_new[l]
    
    
    old <- hf_step(H0, s, U, beta)
    new <- hf_step(H0, s_new, U, beta)
    
    
    if(log(runif(1)) < (new$logdet - old$logdet)){
      s <- s_new
    }
    
    
    if (step > thermal && ((step - thermal) %% measure_interval == 0)) {
      
      cur <- hf_step(H0, s, U, beta)
      
      Gup <- solve(diag(nrow(H0)) + cur$Bup)
      Gdn <- solve(diag(nrow(H0)) + cur$Bdn)
      
      n_up <- 1 - Gup[1,1]
      n_dn <- 1 - Gdn[1,1]
      
      meas_nup <- c(meas_nup, Re(n_up))
      meas_ndn <- c(meas_ndn, Re(n_dn))
    }
  }
  
  
  list(
    measurements = list(
      n_up = meas_nup,
      n_dn = meas_ndn
    )
  )
}