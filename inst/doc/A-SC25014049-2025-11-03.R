## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=FALSE---------------------------------------------------------------
cvm <- function(x, y) {
  x1 <- length(x)
  x2 <- length(y)
  r <- rank(c(x, y))
  r_x <- r[1:x1]
  r_y <- r[(x1 + 1):(x1 + x2)]
  U <- x1 * sum((r_x - 1:x1)^2) + x2 * sum((r_y - 1:x2)^2)
  W2 <- U / (x1 * x2 * (x1 + x2)) - (4 * x2 * x1 - 1) / (6 * (x2 + x1))
  return(W2)
}
set.seed(123)

x <- rnorm(100,0,1)
y <- rnorm(100,0,1)
w2 = cvm(x, y)

sample3 <- rnorm(100, 0,1)
sample4 <- rnorm(100,100,10) 
w2_1 <- cvm(sample3, sample4)

knitr::kable(data.frame(
  "相同分布的W^2:"= w2,
  "不同分布的W^2:"= w2_1
))


## ----echo=FALSE---------------------------------------------------------------

centered <- function(x, y) {
  x <- x - mean(x)
  y <- y - mean(y)
  outx <- sum(x > max(y)) + sum(x < min(y))
  outy <- sum(y > max(x)) + sum(y < min(x))
  return(max(c(outx, outy)))
}

p_test <- function(x, y, R = 1000) {
  n1 <- length(x)
  n2 <- length(y)
  combined <- c(x, y)
  T_obs <- centered(x, y)
  
  T_perm <- replicate(R, {
    perm <- sample(combined)
    x_perm <- perm[1:n1]
    y_perm <- perm[(n1 + 1):(n1 + n2)]
    centered(x_perm, y_perm)
  })
  
  p_ <- mean(T_perm >= T_obs)
  return(p_)
}

n1 <- 20
n2 <- 30
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
m <- 10000
R <- 100   

set.seed(123) 
e_rej <- numeric(m) 

for (i in 1:m) {
  x <- rnorm(n1, mu1, sigma1)
  y <- rnorm(n2, mu2, sigma2)
  p_ <- p_test(x, y, R)
  e_rej[i] <- (p_ < 0.05)
}


alphahat <- mean(e_rej)
se <- sqrt(alphahat * (1 - alphahat) / m)

knitr::kable(data.frame(
  "TIe rate-est"= alphahat,
  "stad.error"= se
))

