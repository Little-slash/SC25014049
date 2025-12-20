## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,fig.width=7, fig.height=5, out.width="100%")

## ----echo=FALSE,warning=FALSE-------------------------------------------------
library(coda)
set.seed(123)

lap <- function(x) {
  0.5 * exp(-abs(x))
}


lap_m <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0  
  
  for(i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)

    alpha_accept <- lap(y) / lap(x[i-1])
    
    if(u[i] < alpha_accept) {
      x[i] <- y
    } else {
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
  
  acceptance_rate <- 1 - k / (N-1)
  return(list(x = x, acceptance_rate = acceptance_rate))
}

run_until_converged <- function(sigma, x0_values, max_iter = 50000, target_rhat = 1.2) {
  n_chains <- length(x0_values)
  chains <- list()
  acceptance_rates <- numeric(n_chains)
  

  iter <- 2000
  for (j in 1:n_chains) {
    result <- lap_m(sigma, x0_values[j], iter)
    chains[[j]] <- result$x
    acceptance_rates[j] <- result$acceptance_rate
  }
  
  mcmc_chains <- lapply(chains, function(chain) mcmc(chain))
  mcmc_list <- mcmc.list(mcmc_chains)
  rhat <- gelman.diag(mcmc_list)$psrf[1]
  

  while (rhat > target_rhat && iter < max_iter) {
    iter <- iter + 1000
    for (j in 1:n_chains) {

      new_result <- lap_m(sigma, tail(chains[[j]], 1), 1000)
      chains[[j]] <- c(chains[[j]], new_result$x[-1])  # 去掉重复的第一个值
      acceptance_rates[j] <- (acceptance_rates[j] * (iter-1000) + 
                              new_result$acceptance_rate * 1000) / iter
    }
    
    mcmc_chains <- lapply(chains, function(chain) mcmc(chain))
    mcmc_list <- mcmc.list(mcmc_chains)
    rhat <- gelman.diag(mcmc_list)$psrf[1]

  }
  
  return(list(
    chains = chains,
    acceptance_rates = acceptance_rates,
    final_iter = iter,
    rhat = rhat
  ))
}

# 主程序
set.seed(123)
N <- 10000
sigma <- c(0.05, 0.5, 2, 16)
x0_values <- c(-10, 0, 10, 25)

results <- list()

for(i in 1:length(sigma)) {
  result <- run_until_converged(sigma[i], x0_values)
  results[[as.character(sigma[i])]] <- result
}

par(mfrow = c(2, 2))


for(i in 1:length(sigma)) {
  sigma_val <- sigma[i]
  result <- results[[as.character(sigma_val)]]
  
  plot(result$chains[[1]], type = "l", 
       main = paste("sigma =", sigma_val, 
                   "\n平均接受率 =", round(mean(result$acceptance_rates), 3),
                   "\nRhat =", round(result$rhat, 3)),
       xlab = "迭代", ylab = "x值",
       ylim = c(-10, 30))
  

  for(j in 2:length(result$chains)) {
    lines(result$chains[[j]], col = j)
  }
}


## ----echo=FALSE---------------------------------------------------------------

summary_table <- data.frame(
  sigma = sigma,
  avg_acceptance_rate = sapply(sigma, function(s) mean(results[[as.character(s)]]$acceptance_rates)),
  final_iterations = sapply(sigma, function(s) results[[as.character(s)]]$final_iter),
  final_rhat = sapply(sigma, function(s) results[[as.character(s)]]$rhat)
)

knitr::kable(summary_table)


## ----echo=FALSE---------------------------------------------------------------

set.seed(123)
library(coda)


a <- 2
b <- 3
n <- 20
nc <- 4 
max_iter <- 50000 
min_iter <- 5000   
burn_ratio <- 0.3 

# 收敛监控函数
monitor_convergence <- function() {
  cl <- list()  # 链列表
  
  for(cid in 1:nc) {

    sm <- matrix(NA, nrow = max_iter, ncol = 2)
    colnames(sm) <- c("x", "y")
    

    x_cur <- sample(0:n, 1)
    y_cur <- runif(1)

    for(i in 1:max_iter) {
      x_cur <- rbinom(1, n, y_cur)
      y_cur <- rbeta(1, x_cur + a, n - x_cur + b)
      sm[i, ] <- c(x_cur, y_cur)
    }
    cl[[cid]] <- mcmc(sm)
  }
  
  # 合并链
  cc <- mcmc.list(cl)
  return(cc)
}


iter <- min_iter
converged <- FALSE
final_chains <- NULL

while(!converged) {

  
  # 运行采样
  chains <- monitor_convergence()
  
  # 计算燃烧期后样本
  burn_in <- floor(iter * burn_ratio)
  chains_burned <- mcmc.list(
    lapply(chains, function(x) mcmc(x[(burn_in+1):iter, ]))
  )
  
  # 计算R_hat
  gd <- gelman.diag(chains_burned, multivariate = FALSE)
  r_hat <- gd$psrf[, 1]
  # 
  # cat("当前R_hat值 - x:", round(r_hat[1], 4), "y:", round(r_hat[2], 4), "\n")
  
  if(all(r_hat < 1.2)) {
    converged <- TRUE
    final_chains <- chains_burned
    cat("最终R_hat值 - x:", round(r_hat[1], 4), "y:", round(r_hat[2], 4), "\n")
  } 
  else {
    iter <- iter * 2
    if(iter > max_iter) {
      iter <- max_iter
      final_chains <- chains_burned
    }
  }
}

final_gd <- gelman.diag(final_chains)
# print(final_gd)

par(mfrow = c(1, 2))
# 绘制收敛图
gelman.plot(final_chains)

## ----echo=FALSE---------------------------------------------------------------

# 提取合并样本
all_samps <- do.call(rbind, final_chains)
x_samps <- all_samps[, "x"]
y_samps <- all_samps[, "y"]


df <- data.frame(x = x_samps, y = y_samps)
chain1_x <- final_chains[[1]][, "x"]
plot(chain1_x, type = "l", col = "blue", 
     main = "X的轨迹图（第一条链）", xlab = "迭代", ylab = "x值")

# 6. Y的轨迹图（第一条链）
chain1_y <- final_chains[[1]][, "y"]
plot(chain1_y, type = "l", col = "red", 
     main = "Y的轨迹图（第一条链）", xlab = "迭代", ylab = "y值")

## ----echo=FALSE---------------------------------------------------------------
set.seed(12345)

model <-function(N,b1,b2,b3,f0){
  X1 <- rpois(N, 1)    
  X2 <- rexp(N, 1)             
  X3 <- rbinom(N, 1,0.5)
  
  g <- function(alpha) {

    tmp <- alpha + b1 * X1 + b2 * X2 + b3 * X3
    p <- 1 / (1 + exp(tmp))
    
    mean(p) - f0
  }
  alpha <- uniroot(g,c(-20,20))$root
  return (alpha)
  
}



N <- 10^6
b1 <- 1
b2 <- 1
b3 <- -1
f0_values <- c(0.1, 0.01, 0.001, 0.0001)

results <- data.frame(
  f0 = f0_values,
  alpha1 = numeric(length(f0_values))
  
)

for (i in seq_along(f0_values)) {

  alpha_result <- model(N, b1, b2, b3, f0_values[i])
  results$alpha1[i] <- alpha_result
}
knitr::kable(t(results))


