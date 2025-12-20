## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=FALSE---------------------------------------------------------------
set.seed(123)
m<-1000
sim <- 10000
alpha <-0.1

null_n <- 950
al_n <- 50


fwer_1 <- numeric(sim)
fwer_2 <- numeric(sim)

fdr_1 <- numeric(sim)
fdr_2 <- numeric(sim)

tpr_1 <- numeric(sim)
tpr_2 <- numeric(sim)

for (i in 1:sim){
  null_p <- runif(null_n)
  al_p <-  rbeta(al_n,0.1,1)
  
  truth_results <- c(rep(0, null_n), rep(1, al_n))
  p_values <- c(null_p, al_p)
  
  #bonfer
  
  p1 <- p.adjust(p_values, method = "bonferroni")
  truth_e0_1 <- 0 
  for (j in 1:950){
    if (p1[j]<=0.1){
      truth_e0_1=truth_e0_1+1
    }
  }
  fwer_1[i]=as.numeric(truth_e0_1>=1)
  reject_1 <- as.numeric(p1 <= alpha)
  fdr_1[i] <- ifelse(sum(reject_1[1:null_n])>0,sum(reject_1[1:null_n]) / sum(reject_1),0)
  tpr_1[i] <- sum(reject_1[(null_n+1):m]) / al_n
  #B-H
  
  p2 <- p.adjust(p_values, method = "BH")
  truth_e0_2 <- 0 
  for (j in 1:950){
    if (p2[j]<=0.1){
      truth_e0_2=truth_e0_2+1
    }
  }
  fwer_2[i]=as.numeric(truth_e0_2>=1)
  reject_2 <- as.numeric(p2 <= alpha)
  fdr_2[i] <- ifelse(sum(reject_2[1:null_n])>0,sum(reject_2[1:null_n]) / sum(reject_2),0)
  tpr_2[i] <- sum(reject_2[(null_n+1):m]) / al_n
}

results <- data.frame(
  FWER = c(mean(fwer_1),mean(fwer_2)),
  FDR = c(mean(fdr_1),mean(fdr_2)),
  TPR = c(mean(tpr_1),mean(tpr_2))
)

row.names(results) <- c("Bonferroni correction", "B-H correction")


knitr::kable(t(results))



## ----echo=FALSE---------------------------------------------------------------
rm(list = ls())

times <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
n <- length(times) 

# 计算λ的MLE (危险率估计)
mean_t <- mean(times) 
lambda_mle <- 1/ mean_t

cat("λ的原始MLE =", round(lambda_mle, 6), "/小时\n")

## ----echo=FALSE---------------------------------------------------------------


set.seed(123) 
B <- 10000   
boot_lambda <- numeric(B) 


for (i in 1:B) {
  boot_sample <- sample(times, size = n, replace = TRUE)
  boot_lambda[i] <- 1/ mean(boot_sample)
}

bias_boot <- mean(boot_lambda) - lambda_mle 
se_boot <- sd(boot_lambda)              

knitr::kable(t(
  data.frame(
    "boostrap" = mean(boot_lambda),
    "Bootstrap_bias " =bias_boot,
    "Bootstrap_sd " = se_boot
  ))
)


## ----echo=FALSE---------------------------------------------------------------
rm(list = ls())
library(bootstrap)
data(scor, package = "bootstrap")


Sigma_hat <- cov(scor)
eigen_values <- eigen(Sigma_hat)$values

theta_hat <- eigen_values[1] / sum(eigen_values)

knitr::kable(t(
  data.frame(
    eigen_values = eigen_values,
    theta_hat = theta_hat
  ))
)

## ----echo=FALSE---------------------------------------------------------------

B <- 10000
n <- nrow(scor)
theta_hat_b <- numeric(B)

for (i in 1:B) {

  boot_indices <- sample(1:n, size = n, replace = TRUE)
  boot_sample <- scor[boot_indices, ]

  Sigma_boot <- cov(boot_sample)
  eigen_boot <- eigen(Sigma_boot)$values
  theta_hat_b[i] <- eigen_boot[1] / sum(eigen_boot)
  
}


boot_mean <- mean(theta_hat_b)
boot_bias <- boot_mean - theta_hat
boot_se <- sd(theta_hat_b)


knitr::kable(t(
  data.frame(
    bootstrap_est = boot_mean,
    bias = boot_bias,
    sd = boot_se
  ))
)


