## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=FALSE---------------------------------------------------------------

stratified_importance_sampling <- function(n = 5, M = 10000) {
  n_per <- M / n
  stratum_vars <- numeric(n)
  stratum_estimates <- numeric(n)
  for (j in 0:(n-1)) {
    a <- j / n
    b <- (j + 1) / n
    u <- runif(n_per,a,b)
    # x <- -log(1 - u * (1 - exp(-1)) * (1 - exp(-1/n)) / 
    #             (1 - exp(-1)) * exp(j/n))
    x <- -log(1 - u * (1 - exp(-1)))
    
    hx <- exp(-x) / (1 + x^2)
    stratum_f_x <- n * exp(-x) / (1 - exp(-1))
    
    weights <- hx / stratum_f_x
    stratum_estimates[j+1] <- mean(weights)
    stratum_vars[j+1] <- var(weights) / n_per
  }
  theta_hat <- sum(stratum_estimates)
  se_hat <- sqrt(sum(stratum_vars))
  
  return(list(
    estimate = theta_hat,
    std_error = se_hat,
    stratum_estimates = stratum_estimates
  ))
}

# 运行分层重要性抽样
set.seed(12345)
result_stratified <- stratified_importance_sampling(n = 5, M = 10000)

knitr::kable(data.frame(
  "method"= c("sample", "stratified_sampling"),
  "est"= c(0.5257801,result_stratified$estimate),
  "stad.error"=c(0.0970314,result_stratified$std_error)
))

## ----echo=FALSE---------------------------------------------------------------
rm(list=ls())

library(ggplot2)
plot_power_curve <- function(n){
  m<-1000
  mu0<-500
  sigma<-100
  mu<-c(seq(450,650,10))#alternatives
  M <- length(mu)
  power1 <- numeric(M)
  for(i in 1:M){
    mu1 <- mu[i]
    pvalues<-replicate(m,expr={
      x<-rnorm(n,mean=mu1,sd =sigma)
      ttest <- t.test(x,alternative ="greater",
                      mu = mu0)
      ttest$p.value  })
    power1[i]<-mean(pvalues<=.05)
  }
  se<-sqrt(power1 *(1-power1)/m)
  # return(list(
  #   mu <- mu,
  #   power1 <- power1,
  #   se=se
  # ))
  return(data.frame(
    n = n,
    mu = mu,
    power = power1,
    se = se,
    upper = power1 + 2 * se,
    lower = power1 - 2 * se
  ))
}



sample_size <- c(seq(10,50,10))
all_results <- data.frame()

for (i in sample_size){
  result_n <- plot_power_curve(i)
  all_results <- rbind(all_results, result_n)
}

all_results$n <- as.factor(all_results$n)


ggplot(all_results, aes(x = mu, y = power, color = n, linetype = n)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = c(0, 0.05), linetype = c("solid", "dashed"), color = "gray50") +
  labs(
    title = "power曲线与样本量(n)的关系",
    subtitle = "不同样本量下的power随均值变化情况",
    x = "备择假设均值 (mu)",
    y = "统计(Power)",
    color = "样本量 (n)",
    linetype = "样本量 (n)"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )


## ----echo=FALSE---------------------------------------------------------------
n <- 20      
df_chisq <- 2    
mu_true <- df_chisq 
alpha <- 0.05   
conf_level <- 1 - alpha  
M <- 10000    

cover <- logical(M)

set.seed(123)
for (i in 1:M) {

  x <- rchisq(n, df = df_chisq)
  x_bar <- mean(x)
  s <- sd(x)
  t_quantile <- qt(1 - alpha/2, df = n-1)
  margin_error <- t_quantile * s / sqrt(n)
  ci_lower <- x_bar - margin_error
  ci_upper <- x_bar + margin_error
  

  cover[i] <- (ci_lower <= mu_true) & (mu_true <= ci_upper)
}


coverage_prob <- mean(cover)

# 输出结果
cat("基于", M, "次蒙特卡洛模拟的结果：\n")
cat("样本大小：", n, "\n")
cat("数据分布：卡方，理论均值 =", mu_true, "\n")
cat("95% t区间的实际覆盖概率：", round(coverage_prob, 4), "\n")



## ----echo=FALSE---------------------------------------------------------------

sigma2 <- 2 * df_chisq 
cover_var <- logical(M)  

set.seed(123)
for (i in 1:M) {
  x <- rchisq(n, df = df_chisq)
  s2 <- var(x)
  chi2_critical <- qchisq(alpha, df = n - 1)
  upper_var <- (n - 1) * s2 / chi2_critical
  cover_var[i] <- (sigma2 <= upper_var)
}
coverage_prob_var <- mean(cover_var)

knitr::kable(data.frame(
 "method"=c("t-interval","varience"),
  "coverage"=c(coverage_prob,coverage_prob_var) 
)
)

## ----echo=FALSE---------------------------------------------------------------
rm(list=ls())


set.seed(123)
M <- 10000     
n <- 30        
alpha <- 0.05  
results <- list()


mu0_chisq <- 1

reject_chisq <- numeric(M)

for(i in 1:M) {
  sample_data <- rchisq(n, df = 1)
  
  test_result <- t.test(sample_data, mu = mu0_chisq)
  
  reject_chisq[i] <- test_result$p.value < alpha
}

empirical_alpha_chisq <- mean(reject_chisq)
results$chisq <- empirical_alpha_chisq

cat("x²(1)分布的经验第一类错误率:", round(empirical_alpha_chisq, 4), "\n")


## ----echo=FALSE---------------------------------------------------------------


mu0_unif <- 1
reject_unif <- numeric(M)

for(i in 1:M) {
  sample_data <- runif(n, min = 0, max = 2)
  test_result <- t.test(sample_data, mu = mu0_unif)
  
  reject_unif[i] <- test_result$p.value < alpha
}

empirical_alpha_unif <- mean(reject_unif)
results$uniform <- empirical_alpha_unif

cat("Uniform(0,2)分布的经验第一类错误率:", round(empirical_alpha_unif, 4), "\n")

## ----echo=FALSE---------------------------------------------------------------


mu0_exp <- 1

reject_exp <- numeric(M)

for(i in 1:M) {
  sample_data <- rexp(n, rate = 1)
  test_result <- t.test(sample_data, mu = mu0_exp)
  reject_exp[i] <- test_result$p.value < alpha
}

empirical_alpha_exp <- mean(reject_exp)
results$exponential <- empirical_alpha_exp

cat("Exponential(1)分布的经验第一类错误率:", round(empirical_alpha_exp, 4), "\n")


## ----echo=FALSE---------------------------------------------------------------


# 汇总结果
summary_df <- data.frame(
  Distribution = c("χ²(1)", "Uniform(0,2)", "Exponential(1)"),
  Nominal_Alpha = alpha,
  Empirical_Alpha = c(results$chisq, results$uniform, results$exponential),
  Bias = c(results$chisq - alpha, results$uniform - alpha, results$exponential - alpha)
)
knitr::kable(summary_df)


## ----echo=FALSE---------------------------------------------------------------
# 可视化结果
library(ggplot2)

ggplot(summary_df, aes(x = Distribution, y = Empirical_Alpha)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = alpha, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "t检验在不同分布下的经验第一类错误率",
       subtitle = paste("显著性水平=", alpha, "，样本量 n =", n, "，模拟次数 M =", M),
       y = "经验第一类错误率",
       x = "分布类型") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.1)) +
  geom_text(aes(label = round(Empirical_Alpha, 4)), vjust = -0.5)


