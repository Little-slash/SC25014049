## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,fig.width=7, fig.height=5, out.width="100%")

## ----echo=FALSE---------------------------------------------------------------
set.seed(123)

r <- runif(1000)
sigma1 <- 1
x <- sigma1*sqrt(-2*log(1-r))
knitr::kable(head(x))

## ----echo=FALSE---------------------------------------------------------------
hist(x,breaks = 30,prob = TRUE, main="Rayleigh distribution")
curve((x/sigma1^2)*exp(-x^2/2/sigma1^2),add = TRUE, col="red",lwd = 2)


## ----echo=FALSE---------------------------------------------------------------

set.seed(123)

r <- runif(1000)
sigma1 <- 2
x <- sigma1*sqrt(-2*log(1-r))
hist(x,breaks = 30,prob = TRUE, main="Rayleigh distribution")
curve((x/sigma1^2)*exp(-x^2/2/sigma1^2),add = TRUE, col="red",lwd = 2)


## ----echo=FALSE---------------------------------------------------------------
my_data <- data.frame(
  x = c(0,1,2,3,4),
  "p(x)" = c(0.1, 0.2, 0.2, 0.2,0.3) 
)

knitr::kable(t(my_data))


## ----echo=FALSE---------------------------------------------------------------

x <- c(-1, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5)
y <- c(0, 0, 0.1, 0.1, 0.3, 0.3, 0.5, 0.5, 0.7, 0.7, 1, 1)

plot(x, y, type = "s", lwd = 2, col = "blue",
     main = "离散随机变量X的累积分布函数",
     xlab = "x", ylab = "F(x)",
     xlim = c(-0.5, 4.5), ylim = c(0, 1.1),
     xaxt = "n", yaxt = "n")
axis(1, at = -1:5)
axis(2, at = seq(0, 1, by = 0.1))
points(4, 1, pch = 16, col = "red", cex = 1.5)
points(c(0,1,2,3), c(0.1,0.3,0.5,0.7), pch = 1, col = "blue", cex = 1.5)


## ----echo=FALSE---------------------------------------------------------------
set.seed(123)
x <- c(0,1,2,3,4)
p <- c(0.1, 0.2, 0.2, 0.2, 0.3)
cp <- cumsum(p)
m <- 1000
U = runif(m)
r <- x[findInterval(U,cp)+1]
knitr::kable(head(r))

## ----echo=FALSE---------------------------------------------------------------
ct1 <- as.vector(table(r))
set.seed(123)
compare2 = sample(0:4,size = 1000,replace = TRUE,prob = c(.1,.2,.2,.2,.3))
ct2<- as.vector(table(compare2))

knitr::kable(t(data.frame("Inverse_Transform"=ct1/1000,
                          "Sample_Function"=ct2/1000,
                          "Theoretical" = c(0.1,0.2,0.2,0.2,0.3)
                          )),caption = "comparison of empirical and theorectical")

## ----echo=FALSE---------------------------------------------------------------

set.seed(123)

count <- 0

y <- numeric(1000)
c <- 4/27
while (count < 1000) {
  
    x_candi <- runif(1)
    
    f_x <- x_candi^2 * (1-x_candi)
    u <- runif(1)
    
    if (u <= f_x / c) {
      count <- count + 1
      y[count] <- x_candi
    }
}
knitr::kable(head(y))


## ----echo=FALSE---------------------------------------------------------------

x_vals <- seq(0, 1, length.out = 1000)
theory_density <- dbeta(x_vals, 3, 2)

hist(y, prob = TRUE, breaks = 30, col = "blue",
     main = "Beta(3,2)分布的接受-拒绝采样结果",
     xlab = "x", ylab = "密度", ylim = c(0, 2.0))
lines(x_vals, theory_density, col = "red", lwd = 2)

legend("topleft", legend = c("sample", "theoretical"), 
       fill = c("blue", "red"))



## ----echo=FALSE---------------------------------------------------------------
mix_normal <- function(num, p1){
  set.seed(123)
  count <- 1
  y <- numeric(num+1)
  
  while(count<=num){
    r <- runif(1)
    if(r<=p1) {
      re  <- rnorm(1,0,1)
      y[count] <- re
      count <- count+1
    }
    else{
      re <- rnorm(1,3,1)
      y[count] <- re
      count <- count+1
    }
  }
  return(y)
}

knitr::kable(head(mix_normal(1000,0.75)))


## ----echo=FALSE---------------------------------------------------------------

paint <- function(p1, n = 1000){
  sample = mix_normal(1000,p1)
  
  mixture_density <- function(x) {
      p1 * dnorm(x, mean = 0, sd = 1) + (1 - p1) * dnorm(x, mean = 3, sd = 1)
  }
  
  hist(sample, breaks = 50, freq = FALSE, main = paste("Mixture with p1 =", p1),xlab = "x", ylab = "Density", col = "blue", border = "white")
  
  x <- seq(min(sample), max(sample), length.out = 1000)
  lines(x, mixture_density(x), col = "red", lwd = 2)

  legend("topright", legend = c("Histogram", "Density Curve"), 
         col = c("blue", "red"), lwd = c(2, 2), bty = "n")
}
paint(0.75)


## ----echo=FALSE---------------------------------------------------------------

par(mfrow = c(3, 3))
p1_values <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

for (p1 in p1_values) {
  paint(p1)
}

par(mfrow = c(1, 1))



## ----echo=FALSE---------------------------------------------------------------

set.seed(123)
r <- 4
theta_scale <- 2
n <- 1000
lambda <- rgamma(n, shape = r, scale = theta_scale)
knitr::kable(head(lambda, 10))


## ----echo=FALSE---------------------------------------------------------------

Y <- rexp(n, rate = lambda)

knitr::kable(head(Y, 10),caption = "\n前10个观测值:\n")

## ----echo=FALSE---------------------------------------------------------------
mean_Y <- mean(Y)
print(paste("生成变量的均值为", round(mean_Y, 4)))


## ----echo=FALSE---------------------------------------------------------------

hist(Y, breaks = 30, main = "Exponential-Gamma Mixture (r=4, scale=2)", 
     xlab = "Y", col = "lightblue")


