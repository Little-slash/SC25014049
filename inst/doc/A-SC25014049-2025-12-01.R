## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,fig.width=7, fig.height=5, out.width="100%")

## ----echo=FALSE---------------------------------------------------------------
# 设置随机种子确保结果可重现
set.seed(123)

# 加密变量名定义
a1_obs_counts <- c(125, 18, 20, 34)
b2_total_samples <- sum(a1_obs_counts) 
c3_chains <- 4
d4_iter_per_chain <- 4000 
e5_burnin <- 3000 
f6_thin_interval <- 5 
g7_proposal_sd <- 0.08 


h8_calc_log_posterior <- function(i9_param_theta) {
  if (i9_param_theta < 0 || i9_param_theta > 1) {
    return(-Inf)
  }
  
  j10_prob_vector <- c(
    0.5 + i9_param_theta/4,     
    (1 - i9_param_theta)/4,   
    (1 - i9_param_theta)/4,      
    i9_param_theta/4           
  )
  
  if (any(j10_prob_vector < 0) || any(j10_prob_vector > 1)) {
    return(-Inf)
  }

  k11_log_likelihood <- sum(a1_obs_counts * log(j10_prob_vector))
  
  l12_log_prior <- 0
  return(k11_log_likelihood + l12_log_prior)
}


m13_run_mh_chain <- function(n14_start_theta, o15_chain_length, p16_burnin, q17_thin, r18_proposal_sd) {

  s19_samples <- numeric(o15_chain_length)
  t20_current_theta <- n14_start_theta
  u21_current_log_post <- h8_calc_log_posterior(t20_current_theta)
  v22_accept_count <- 0
  

  for (w23_iter in 1:o15_chain_length) {

    x24_candidate <- rnorm(1, mean = t20_current_theta, sd = r18_proposal_sd)
    

    y25_candidate_log_post <- h8_calc_log_posterior(x24_candidate)

    z26_accept_prob <- exp(y25_candidate_log_post - u21_current_log_post)
    z26_accept_prob <- min(1, z26_accept_prob)
    
    # 决定是否接受候选值
    if (runif(1) < z26_accept_prob) {
      t20_current_theta <- x24_candidate
      u21_current_log_post <- y25_candidate_log_post
      v22_accept_count <- v22_accept_count + 1
    }
    

    s19_samples[w23_iter] <- t20_current_theta
  }
  

  aa27_accept_rate <- v22_accept_count / o15_chain_length

  bb28_post_burn_samples <- s19_samples[(p16_burnin + 1):o15_chain_length]
  cc29_thinned_samples <- bb28_post_burn_samples[seq(1, length(bb28_post_burn_samples), by = q17_thin)]
  
  # 返回结果
  return(list(
    chain_samples = cc29_thinned_samples,
    full_chain = s19_samples,
    accept_rate = aa27_accept_rate,
    chain_id = as.character(n14_start_theta)
  ))
}

# 运行多条独立MCMC链

dd30_all_chains <- list()
ee31_initial_values <- seq(0.1, 0.9, length.out = c3_chains)  # 分散的初始值

for (ff32_chain_idx in 1:c3_chains) {
  gg33_chain_result <- m13_run_mh_chain(
    n14_start_theta = ee31_initial_values[ff32_chain_idx],
    o15_chain_length = d4_iter_per_chain,
    p16_burnin = e5_burnin,
    q17_thin = f6_thin_interval,
    r18_proposal_sd = g7_proposal_sd
  )
  
  dd30_all_chains[[ff32_chain_idx]] <- gg33_chain_result

}

# Gelman-Rubin收敛诊断函数
hh34_gelman_rubin_diagnostic <- function(ii35_chain_list) {
  jj36_chain_samples <- lapply(ii35_chain_list, function(x) x$chain_samples)
  
  kk37_chain_lengths <- sapply(jj36_chain_samples, length)
  ll38_chain_means <- sapply(jj36_chain_samples, mean)
  mm39_chain_vars <- sapply(jj36_chain_samples, var)
  

  nn40_overall_mean <- mean(unlist(jj36_chain_samples))
  
  oo41_num_chains <- length(jj36_chain_samples)
  pp42_n <- kk37_chain_lengths[1]
  
  # 计算链内方差(Between-chain variance)
  qq43_b <- pp42_n * sum((ll38_chain_means - nn40_overall_mean)^2) / (oo41_num_chains - 1)
  
  # 计算链间方差(Within-chain variance)
  rr44_w <- mean(mm39_chain_vars)
  
  # 计算潜在尺度缩减因子(Potential Scale Reduction Factor)
  # 方差估计
  ss45_var_hat <- ((pp42_n - 1) / pp42_n) * rr44_w + (1 / pp42_n) * qq43_b
  
  # R-hat统计量
  tt46_r_hat <- sqrt(ss45_var_hat / rr44_w)
  
  # 计算有效样本量(Effective Sample Size)
  # 首先计算自相关
  uu47_combined_samples <- unlist(jj36_chain_samples)
  vv48_auto_corr <- acf(uu47_combined_samples, plot = FALSE, lag.max = 1000)$acf[, 1, 1]
  
  # 找到第一个负自相关的滞后
  ww49_lag_tau <- which(vv48_auto_corr < 0)[1]
  if (is.na(ww49_lag_tau)) ww49_lag_tau <- length(vv48_auto_corr)
  
  # 整合自相关时间
  xx50_integrated_act <- 1 + 2 * sum(vv48_auto_corr[2:min(ww49_lag_tau, length(vv48_auto_corr))])
  
  # 有效样本量
  yy51_ess <- length(uu47_combined_samples) / xx50_integrated_act
  
  return(list(
    r_hat = tt46_r_hat,
    between_var = qq43_b,
    within_var = rr44_w,
    var_hat = ss45_var_hat,
    ess = yy51_ess,
    chain_means = ll38_chain_means,
    chain_vars = mm39_chain_vars
  ))
}


zz52_gr_diag <- hh34_gelman_rubin_diagnostic(dd30_all_chains)

# 合并所有链的样本进行后验分析
aaa53_all_post_samples <- unlist(lapply(dd30_all_chains, function(x) x$chain_samples))

# 后验分布统计量
bbb54_post_mean <- mean(aaa53_all_post_samples)
ccc55_post_median <- median(aaa53_all_post_samples)
ddd56_post_sd <- sd(aaa53_all_post_samples)
eee57_post_ci <- quantile(aaa53_all_post_samples, c(0.025, 0.975))

cat("链数量:", c3_chains, "\n")
cat("每条链迭代次数:", d4_iter_per_chain, "\n")
cat("每条链后验样本数:", length(dd30_all_chains[[1]]$chain_samples), "\n")
cat("总后验样本数:", length(aaa53_all_post_samples), "\n")


cat("Gelman-Rubin R-hat值:", round(zz52_gr_diag$r_hat, 4), "\n")
if (zz52_gr_diag$r_hat < 1.1) {
  cat("结论: R-hat < 1.1，链已收敛 ✓\n")
} else if (zz52_gr_diag$r_hat < 1.2) {
  cat("结论: R-hat < 1.2，链基本收敛\n")
} else {
  cat("结论: R-hat ≥ 1.2，链可能未收敛，建议增加迭代次数\n")
}


cat("\n========== 后验分布统计量 ==========\n")
cat("后验均值:", round(bbb54_post_mean, 4), "\n")
cat("后验中位数:", round(ccc55_post_median, 4), "\n")
cat("后验标准差:", round(ddd56_post_sd, 4), "\n")
cat("95%后验可信区间: (", round(eee57_post_ci[1], 4), ", ", round(eee57_post_ci[2], 4), ")\n")

# 绘制收敛诊断图
par(mfrow = c(2, 2))

# 1. 多条链轨迹图（重叠）
colors <- c("red", "blue", "green", "purple")
plot(1, type = "n", xlim = c(1, d4_iter_per_chain), ylim = c(0, 1),
     main = "多条MCMC链轨迹图", xlab = "迭代次数", ylab = expression(theta))

for (fff58_idx in 1:c3_chains) {
  lines(dd30_all_chains[[fff58_idx]]$full_chain, 
        col = adjustcolor(colors[fff58_idx], alpha.f = 0.6), lwd = 1)
}
abline(v = e5_burnin, col = "black", lty = 2, lwd = 2)
legend("topright", legend = paste("链", 1:c3_chains), 
       col = colors[1:c3_chains], lty = 1, lwd = 2, cex = 0.8)

# 2. 链均值演化图
plot(1, type = "n", xlim = c(1, c3_chains), ylim = c(0, 1),
     main = "各链后验均值比较", xlab = "链编号", ylab = "后验均值")
points(1:c3_chains, zz52_gr_diag$chain_means, 
       col = colors[1:c3_chains], pch = 19, cex = 1.5)
abline(h = bbb54_post_mean, col = "black", lty = 2, lwd = 2)
segments(1:c3_chains, zz52_gr_diag$chain_means - sqrt(zz52_gr_diag$chain_vars),
         1:c3_chains, zz52_gr_diag$chain_means + sqrt(zz52_gr_diag$chain_vars),
         col = colors[1:c3_chains], lwd = 2)

plot(density(dd30_all_chains[[1]]$chain_samples), col = colors[1], lwd = 2,
     main = "各链后验密度对比", xlab = expression(theta), ylab = "密度",
     xlim = c(0, 1), ylim = c(0, max(sapply(dd30_all_chains, 
     function(x) max(density(x$chain_samples)$y))) * 1.2))

for (ggg59_idx in 2:c3_chains) {
  lines(density(dd30_all_chains[[ggg59_idx]]$chain_samples), 
        col = colors[ggg59_idx], lwd = 2)
}
legend("topright", legend = paste("链", 1:c3_chains), 
       col = colors[1:c3_chains], lty = 1, lwd = 2, cex = 0.8)


barplot(c(zz52_gr_diag$r_hat, 1.1), names.arg = c("R-hat", "阈值"),
        col = c(ifelse(zz52_gr_diag$r_hat < 1.1, "green", "orange"), "red"),
        main = "Gelman-Rubin诊断", ylab = "R-hat值", ylim = c(0, max(1.2, zz52_gr_diag$r_hat)))
abline(h = 1.1, col = "red", lty = 2, lwd = 2)




## ----echo=FALSE---------------------------------------------------------------
# 设置随机种子保证结果可重现
set.seed(123)

sigma_value <- 2.0
sample_size <- 10000

# 逆变换
generate_rayleigh_inverse <- function(n, sigma) {
  uniform_random <- runif(n)
  rayleigh_random <- sigma * sqrt(-2 * log(uniform_random))
  return(rayleigh_random)
}

# BOX
generate_rayleigh_boxmuller <- function(n, sigma) {

  if (n %% 2 != 0) {
    n <- n + 1
  }
  
  U1_value <- runif(n/2)
  U2_value <- runif(n/2)
  
  Z1_vector <- sqrt(-2 * log(U1_value)) * cos(2 * pi * U2_value)
  Z2_vector <- sqrt(-2 * log(U1_value)) * sin(2 * pi * U2_value)
  
  R_values <- sqrt(Z1_vector^2 + Z2_vector^2)

  rayleigh_random <- sigma * c(R_values)
  
  # 如果n原本是奇数，去掉最后一个值
  if (length(rayleigh_random) > n) {
    rayleigh_random <- rayleigh_random[1:n]
  }
  
  return(rayleigh_random)
}

# 生成三组随机数
rayleigh_inverse <- generate_rayleigh_inverse(sample_size, sigma_value)
rayleigh_boxmuller <- generate_rayleigh_boxmuller(sample_size, sigma_value)



knitr::kable(data.frame(head(rayleigh_inverse),head(rayleigh_boxmuller)))

## ----echo=FALSE---------------------------------------------------------------
library(ggplot2)
theoretical_mean <- sigma_value * sqrt(pi/2)
theoretical_sd <- sigma_value * sqrt((4 - pi)/2)

summary_stats <- data.frame(
  方法 = c("逆变换法", "Box-Muller法","理论"),
  均值 = c(mean(rayleigh_inverse), mean(rayleigh_boxmuller),theoretical_mean),
  标准差 = c(sd(rayleigh_inverse), sd(rayleigh_boxmuller),theoretical_sd)
)

knitr::kable(summary_stats)

# 4.3 QQ图比较
qrayleigh <- function(p, sigma = 1) {
  ifelse(p < 0 | p > 1, NaN, sigma * sqrt(-2 * log(1 - p)))
}

qq_data_inverse <- data.frame(
  sample_quantiles = quantile(rayleigh_inverse, probs = ppoints(100)),
  theoretical_quantiles = qrayleigh(ppoints(100), sigma_value)
)

ggplot(qq_data_inverse, aes(theoretical_quantiles, sample_quantiles)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("逆变换法 vs 理论瑞利分布")


par(mfrow = c(1, 2))


# 逆变换法 vs Box-Muller法
qqplot(rayleigh_inverse, rayleigh_boxmuller, 
       main = "逆变换法 vs Box-Muller法",
       xlab = "逆变换法分位数",
       ylab = "Box-Muller法分位数",
       col = "purple", pch = 19, cex = 0.5)
abline(0, 1, col = "red", lwd = 2)


plot(density(rayleigh_inverse), col = "blue", lwd = 2,
     main = "密度曲线比较", xlab = "x", ylab = "密度",
     xlim = c(0, max(rayleigh_inverse, rayleigh_boxmuller)))
lines(density(rayleigh_boxmuller), col = "darkgreen", lwd = 2)
legend("topright", legend = c("逆变换法", "Box-Muller法"),
       col = c("blue", "darkgreen"), lwd = 2, lty = c(1, 1))



## ----echo=FALSE---------------------------------------------------------------
library(microbenchmark)

scale_parameter <- 2.0
sample_size <- 10000
timing_results_list <- list()


  
# 使用microbenchmark进行性能测试
benchmark_result <- microbenchmark(
  逆变换法 = generate_rayleigh_inverse(sample_size, scale_parameter),
  BoxMuller法 = generate_rayleigh_boxmuller(sample_size, scale_parameter),
  times = 100, 
  unit = "ms" 
)



# 提取统计摘要
summary_stats <- summary(benchmark_result)


print(summary_stats)

time_ratio <- summary_stats$mean[2] / summary_stats$mean[1]



