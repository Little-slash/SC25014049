## ----star, include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=FALSE---------------------------------------------------------------
set.seed(123)

data1 <- rnorm(n = 100, mean = 50, sd = 10) #均值为50，标准差为10
data2 <- data.frame(data1)
#mode(data2)
knitr::kable(head(data2))

## ----echo=FALSE---------------------------------------------------------------

x_mean <- mean(data1)


x_low <- x_mean-1.96*10/10
x_high <- x_mean+1.96*10/10


## ----echo=FALSE---------------------------------------------------------------
# 绘制直方图并绘制理论正太分布曲线，评估置信区间。
hist(data1, breaks = 30, col = "lightblue", main = "正态分布数据与95%置信区间", 
     xlab = "值", 
     ylab = "频率",
     prob = TRUE,  
     ylim = c(0, 0.07))

curve(dnorm(x, mean = 50, sd = 10), 
      col = "red", lwd = 2, lty = 2, add = TRUE)

# 添加置信区间线，并填充颜色
abline(v = x_low, col = "orange", lwd = 2, lty = 2)
abline(v = x_high, col = "orange", lwd = 2, lty = 2)

x_fill <- seq(x_low, x_high, length.out = 100)
y_fill <- dnorm(x_fill, mean = x_mean, sd = sd(data1))
polygon(c(x_fill, rev(x_fill)), c(y_fill, rep(0, length(y_fill))), 
        col = rgb(1, 0.65, 0, 0.3), border = NA)

legend("topright", 
       legend = c("总体分布", "95%置信区间"),
       col = c( "red", "orange"),
       lty = c(1, 1),
       lwd = c(2, 2),
       cex = 0.8)

## ----echo=FALSE---------------------------------------------------------------
library(car)
data(cars)
data_e2 <- cars
knitr::kable(data_e2[1:10,])


## ----echo=FALSE---------------------------------------------------------------
poly_model <- lm(dist ~ poly(speed, 2), data = data_e2)
poly_mode2 <- lm(dist ~ poly(speed, 3), data = data_e2)
co1 <- coef(poly_model)
co2 <- coef(poly_mode2)


## ----echo=FALSE---------------------------------------------------------------
plot(data_e2$speed, data_e2$dist, 
     main = "汽车速度与刹车距离关系",
     xlab = "速度", ylab = "刹车距离")

speed_seq <- seq(min(data_e2$speed), max(data_e2$speed), length.out = 100)
predict1 <- predict(poly_model, data.frame(speed = speed_seq))
predict2 <- predict(poly_mode2, data.frame(speed = speed_seq))
lines(speed_seq, predict1, col = "red", lwd = 2)
lines(speed_seq, predict2, col = "blue", lwd = 2)
legend("topleft", legend = c("二次多项式拟合","三次多项式拟合"), col = c("red","blue"), lwd = c(2,2))

## ----echo=FALSE---------------------------------------------------------------

set.seed(123)
n_sam <- 200
set.seed(123)
n <- 200
x1 <- rnorm(n, mean = 10, sd = 2)
x2 <- rnorm(n, mean = 5, sd = 1.5)
x3 <- rnorm(n, mean = 8, sd = 3)

beta0 <- 2.5  
beta1 <- 1.8  
beta2 <- -0.5 
beta3 <- 1

error <- rnorm(n, mean = 0, sd = 1.5)

y <- beta0 + beta1*x1 + beta2*x2 + beta3*x3 + error

data_e3 <- data.frame(y, x1, x2, x3)
knitr::kable(data_e3[1:10,])


## ----echo=FALSE---------------------------------------------------------------

duoyuan_model <- lm(y ~ x1 + x2 + x3, data = data_e3)
co3 = coef(duoyuan_model)

## ----echo=FALSE---------------------------------------------------------------


residuals <- resid(duoyuan_model)
hist(residuals, breaks = 20, 
     xlab = "残差", ylab = "频率",
     main = "残差分布直方图", col = "lightblue", border = "black")

curve(dnorm(x, mean = 0, sd = 1.5)* 100, 
      from = -4.5, to = 4.5,  # 覆盖±3个标准差的范围
      main = "正态分布曲线",
      col = "blue", lwd = 2, add=TRUE)


