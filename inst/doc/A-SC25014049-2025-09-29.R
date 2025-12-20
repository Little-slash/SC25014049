## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=FALSE---------------------------------------------------------------

monte_beta <- function(x, n=100000){
  
  samples <- rbeta(n, 3, 3)
  return(mean(samples <= x))
  
}
x_ <-seq(0.1,0.9,0.1) 

est = rep(0,9)

for(i in 1:9){
  est[i] <- monte_beta(x_[i])
  
}

true_values <- pbeta(x_, 3, 3)

result = cbind(x_,est,true_values)
df1 = data.frame(result)
colnames(df1) <- c("X值", "蒙特卡洛积分估计值","pbeta-真实值")
knitr::kable(df1)


## ----echo=FALSE---------------------------------------------------------------

var_eu <- -1*exp(1)*exp(1)/2-3/2+2*exp(1)
var_sum <- 0.01565

variance_ratio <- (1/4 * var_sum) / var_eu

percent_reduction <- (1 - variance_ratio) * 100

message <- sprintf("方差减少百分比 = %f%%",percent_reduction)
print(message)


## ----echo=FALSE---------------------------------------------------------------

g <- function(x) {
  (x^2 / sqrt(2*pi)) * exp(-x^2/2)
}

f1 <- function(x) {
  dnorm(x) / (1 - pnorm(1))
}


f2 <- function(x) {
  return(exp(-1*(x-1)))
}

x_3 <- seq(1, 5, 0.001)
g_x <- g(x_3)
f1_x <- f1(x_3)
f2_x <- f2(x_3)

plot(x_3, g_x, type = "l", lwd = 2, col = "black", 
     ylab = "Density", main = "函数比较")
lines(x_3, f1_x, col = "red", lwd = 2)
lines(x_3, f2_x, col = "blue", lwd = 2)
legend("topright", legend = c("g(x)", "f1(x): 截断正态", "f2(x): 截断Gamma"),
       col = c("black", "red", "blue"), lwd = 2)

## ----echo=FALSE---------------------------------------------------------------

set.seed(123)
n <- 100000

u1 <- runif(n)
samples_f1 <- qnorm(pnorm(1) + u1 * (1 - pnorm(1)))
weights_f1 <- g(samples_f1) / f1(samples_f1)
est_f1 <- mean(weights_f1)
variance_f1 <- var(weights_f1) / n


sample_f2 <- function(n) {
  rexp(n, rate = 1) + 1
}
samples_f2 <- sample_f2(n)
weights_f2 <- g(samples_f2) / f2(samples_f2)
est_f2 <- mean(weights_f2)
variance_f2 <- var(weights_f2) / n

data3 <- data.frame(t(c(variance_f1,variance_f2)))
colnames(data3) <- c( "使用f1(x)的方差","使用f2(x)的方差")
knitr::kable(data3)

## ----echo=FALSE---------------------------------------------------------------
library(rbenchmark)

n <- c(1e4, 2*1e4, 4 * 1e4, 6 * 1e4, 8 * 1e4)

fast_sort <- function(n){
  sp = sample(1:n)  
  sort(sp)
}

result <- benchmark(
  "a_1"={
    re <- fast_sort(n[1])
  },
  "a_2"={
    re <- fast_sort(n[2])
  },
  "a_3"={
    re <- fast_sort(n[3])
  },
  "a_4"={
    re <- fast_sort(n[4])
  },
  "a_5"={
    re <- fast_sort(n[5])
  },
  replications = 1000,
  columns = c("test", "replications", "elapsed"))

knitr::kable(data.frame(result))


## ----echo=FALSE---------------------------------------------------------------

data <- data.frame(
  n = n,
  t_n = n * log(n),
  a_n = result$elapsed
)

fit <- lm(a_n ~ t_n, data = data)

plot(data$t_n,data$a_n,type="p",main="线性拟合回归",
     xlab="t", ylab="n*log(n)")
abline(fit, col = "red", lwd = 2)

