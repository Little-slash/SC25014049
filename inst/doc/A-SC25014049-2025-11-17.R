## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=FALSE---------------------------------------------------------------

p <- 0.3
t <- 5
n <- 10000  

set.seed(123)
T <- rgeom(n, prob = 1 - p)

T_t <- T[T >= t]
S <- T_t - t 

mean_S <- mean(S)
phi_p <- p / (1 - p)


cat("模拟剩余寿命均值:", round(mean_S, 4), "\n")
cat("理论phi(p):", round(phi_p, 4), "\n")
cat("差异:", round(mean_S - phi_p, 4))

## ----echo=FALSE---------------------------------------------------------------

library(boot)
a <- c(-4, -2, -9)
A1 <- rbind(c(2, 1, 1),   
            c(1, -1, 3))

b1 <- c(2, 3)

res <- simplex(a = a, A1 = A1, b1 = b1, maxi = TRUE)

cat("最优解:\n")
cat("x =", res$soln[1], "\n")
cat("y =", res$soln[2], "\n") 
cat("z =", res$soln[3], "\n")
cat("最优目标函数值:", -res$value, "(最小化形式)\n")  # 注意转换回最小化问题

## ----echo=FALSE---------------------------------------------------------------
u <- c(11,8,27,13,16,0,23,10,24,2)
v <- c(12,9,28,14,17,1,24,11,25,3)
n <- length(u)

l0 <- 0.1

# 直接极大化似然函数
dir_mle <- function(l, u, v) {
  ll <- sum(log(exp(-l*u) - exp(-l*v)))
  return(-ll)
}

res <- optimize(dir_mle, c(0.001, 1), u=u, v=v)
l_direct <- res$minimum

# EM算法
em_alg <- function(u, v, l0, tol=1e-8, max_iter=1000) {
  l_curr <- l0
  l_hist <- numeric(max_iter)
  l_hist[1] <- l0
  
  for (i in 2:max_iter) {
    num <- u*exp(-l_curr*u) - v*exp(-l_curr*v)
    den <- exp(-l_curr*u) - exp(-l_curr*v)
    ex <- 1/l_curr + num/den
    

    l_new <- n / sum(ex)
    
    l_hist[i] <- l_new
    
 
    if (abs(l_new - l_curr) < tol) {
      l_hist <- l_hist[1:i]
      break
    }
    l_curr <- l_new
  }
  
  return(list(l_est = l_curr, l_hist = l_hist, iter = i-1))
}

em_res <- em_alg(u, v, l0)
cat("直接法:", l_direct, "\n")
cat("EM算法:", em_res$l_est, "\n")
cat("相对差异:", abs(l_direct - em_res$l_est)/l_direct, "\n")

## ----echo=FALSE---------------------------------------------------------------

l_em <- em_res$l_hist
err_em <- abs(l_em - l_direct)

ratios <- err_em[2:length(err_em)] / err_em[1:(length(err_em)-1)]
avg_ratio <- mean(ratios[!is.na(ratios) & !is.infinite(ratios)])

# 绘制收敛过程
plot(1:length(l_em), l_em, type="b", xlab="迭代次数", ylab="λ估计值", 
     main="EM算法收敛过程", col="blue")
abline(h=l_direct, col="red", lty=2, lwd=2)
legend("topright", legend=c("EM算法", "直接法MLE"), 
       col=c("blue", "red"), lty=c(1,2))

# 绘制误差收敛
plot(1:length(err_em), log(err_em), type="b", 
     xlab="迭代次数", ylab="log(误差)", 
     main="EM算法收敛速度", col="darkgreen")


## ----echo=FALSE---------------------------------------------------------------

# 比较结果
cat("\n方法比较:\n")
cat("EM算法迭代次数:", em_res$iter, "\n")

