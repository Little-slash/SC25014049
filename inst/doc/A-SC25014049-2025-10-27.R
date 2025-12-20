## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=FALSE---------------------------------------------------------------
library(bootstrap)
data(scor, package = "bootstrap")

set.seed(123)
Sigma_hat <- cov(scor)
eigen_values <- eigen(Sigma_hat)$values

theta_ <- eigen_values[1] / sum(eigen_values)

knitr::kable(t(
  data.frame(
    eigen_values = eigen_values,
    theta_ = theta_
  ))
)

## ----echo=FALSE---------------------------------------------------------------

n <- nrow(scor)
theta_hat <- numeric(n)


for (i in 1:n) {
  jack_sample <- scor[(1:n)[-i], ]

  Sigma_jack <- cov(jack_sample)
  eigen_jack <- eigen(Sigma_jack)$values
  theta_hat[i] <- eigen_jack[1] / sum(eigen_jack)
  
}


jack_mean <- mean(theta_hat)
jack_bias <- (n-1)*(jack_mean - theta_)
jack_se <- sqrt((n-1)*mean((jack_mean - theta_hat)^2))


knitr::kable(t(
  data.frame(
    jackknife_est = jack_mean,
    bias = jack_bias,
    sd = jack_se
  ))
)



## ----echo=FALSE---------------------------------------------------------------

library(DAAG)

n<-length(ironslag$magnetic)
e1<-e2<-e3<-e4 <- array(numeric(n), dim = c(n, n, 2))

for(i in 1:n){
  for(j in (1:n)[-i]){
    k <- c(i, j)
    y <- ironslag$magnetic[-k]
    x<- ironslag$chemical[-k]
    J1<-lm(y~x)
    yhat11<-J1$coef[1]+ J1$coef[2]* ironslag$chemical[i]
    yhat12<-J1$coef[1]+ J1$coef[2]* ironslag$chemical[j]
    e1[i,j,2]<-ironslag$magnetic[i]+ironslag$magnetic[j]-yhat11-yhat12
    
    J2<-lm(y ~x+I(x^2))
    yhat21<-J2$coef[1]+J2$coef[2]*ironslag$chemical[i]+J2$coef[3]*ironslag$chemical[i]^2
    yhat22<-J2$coef[1]+J2$coef[2]*ironslag$chemical[j]+J2$coef[3]*ironslag$chemical[j]^2
    e2[i,j,2]<-ironslag$magnetic[i]+ironslag$magnetic[j]-yhat21-yhat22
    
    J3<-lm(log(y)~x)
    logyhat31<-J3$coef[1]+J3$coef[2]* ironslag$chemical[i]
    yhat31<-exp(logyhat31)
    logyhat32<-J3$coef[1]+J3$coef[2]* ironslag$chemical[j]
    yhat32<-exp(logyhat32)
    e3[i,j,2]<-ironslag$magnetic[i]+ironslag$magnetic[j]-yhat31-yhat32
    
    J4<-lm(log(y)~log(x))
    logyhat41<-J4$coef[1]+ J4$coef[2]* log(ironslag$chemical[i])
    yhat41<-exp(logyhat41)
    logyhat42<-J4$coef[1]+ J4$coef[2]* log(ironslag$chemical[j])
    yhat42<-exp(logyhat42)
    e4[i,j,2]<-ironslag$magnetic[i]+ironslag$magnetic[j]-yhat41-yhat42
  }
}
knitr::kable(data.frame(
  error_mean1  = mean(e1^2),
  error_mean2  = mean(e2^2),
  error_mean3  = mean(e3^2),
  error_mean4  = mean(e4^2)
))

## ----echo=FALSE---------------------------------------------------------------

set.seed(123)

n <- 30; B <- 1000; M <- 1000
true_mean <- 0; alpha <- 0.05

results <- matrix(0, nrow = M, ncol = 9)
colnames(results) <- c("Norm_cov", "Basic_cov", "Perc_cov",
                      "Norm_left", "Basic_left", "Perc_left", 
                      "Norm_right", "Basic_right", "Perc_right")

for (i in 1:M) {
  data <- rnorm(n, true_mean, 1)
  samp_mean <- mean(data)
  
  boot_means <- replicate(B, mean(sample(data, n, replace = TRUE)))
  boot_se <- sd(boot_means)
  boot_quants <- quantile(boot_means, c(alpha/2, 1-alpha/2))
  
  z <- qnorm(1 - alpha/2)
  ci_norm <- samp_mean + c(-1, 1) * z * boot_se
  ci_basic <- 2 * samp_mean - rev(boot_quants)
  ci_perc <- boot_quants

  results[i, 1:3] <- c(
    ci_norm[1] <= true_mean & true_mean <= ci_norm[2],
    ci_basic[1] <= true_mean & true_mean <= ci_basic[2], 
    ci_perc[1] <= true_mean & true_mean <= ci_perc[2]
  )
  
  results[i, 4:6] <- c(
    true_mean < ci_norm[1],
    true_mean < ci_basic[1],
    true_mean < ci_perc[1]
  )
  
  results[i, 7:9] <- c(
    true_mean > ci_norm[2],
    true_mean > ci_basic[2], 
    true_mean > ci_perc[2]
  )
}

# 汇总结果
final_results <- data.frame(
  Method = c("Normal", "Basic", "Percentile"),
  Coverage = colMeans(results[, 1:3]),
  Miss_Left = colMeans(results[, 4:6]),
  Miss_Right = colMeans(results[, 7:9])
)

knitr::kable(final_results)


## ----echo=FALSE---------------------------------------------------------------
bootstrap_1 <-function(n){
  B<-1000
  u<-runif(n)
  x<-log(1-u)*(-2)
  mean_t <- mean(x) 
  lambda_mle <- 1/ mean_t
  boot_lambda <- numeric(B)
  
  for (i in  1:B){
    boot_sample<-sample(x,size = n,replace = TRUE)
    boot_lambda[i]<-1/mean(boot_sample)
  }
  boot_bias <- mean(boot_lambda)-lambda_mle
  boot_se <- sd(boot_lambda)
  return (c(boot_bias,boot_se))
}
simulation<-function(n){
  m<-1000
  thero_bias <- 2/(n-1)
  thero_sd <- 2*n/((n-1)*sqrt(n-2))
  simulation_n_bias <- numeric(m)
  simulation_n_sd <- numeric(m)
  
  for (i in 1:m){
    s_i = bootstrap_1(n)
    simulation_n_bias[i] <- s_i[1]
    simulation_n_sd[i] <- s_i[2]
  }
  mean_bootstrap_bias <- mean(simulation_n_bias)
  mean_bootstrap_se <- mean(simulation_n_sd)
  
  return (list(
    mean_bootstrap_bias = mean_bootstrap_bias,
    mean_bootstrap_se = mean_bootstrap_se,
    n = n,
    theoretical_bias = thero_bias,
    theoretical_se = thero_sd
  ))
  
  
}

re_5 <- simulation(5)
re_10 <- simulation(10)
re_20 <- simulation(20)

sp <- data.frame(
  sample_size = c(re_5$n, re_10$n, re_20$n),
  mean_bootstrap_bias = c(re_5$mean_bootstrap_bias, re_10$mean_bootstrap_bias, re_20$mean_bootstrap_bias),
  mean_bootstrap_se = c(re_5$mean_bootstrap_se, re_10$mean_bootstrap_se, re_20$mean_bootstrap_se),
  theoretical_bias = c(re_5$theoretical_bias, re_10$theoretical_bias, re_20$theoretical_bias),
  theoretical_se = c(re_5$theoretical_se, re_10$theoretical_se, re_20$theoretical_se)
)


knitr::kable(sp, digits = 4)


