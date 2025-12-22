## ----echo=FALSE---------------------------------------------------------------

knitr::include_graphics("./dsc_analyze.png")

## ----echo=TRUE,eval=FALSE-----------------------------------------------------
# 
# temp <- seq(0, 300, length.out = 200) # 从0到300度，200个点
# heat_flow<- 0.5 - 0.3 * dnorm(temp, mean = 150, sd = 15) + 0.8 * dnorm(temp, mean = 250, sd = 12)
# 
# rel <- dsc_analyze(temperature = temp, heat_flow = heat_flow, baseline_method = "linear")
# 
# 
# print(rel)
# 

## ----echo=FALSE---------------------------------------------------------------

knitr::include_graphics("./tga_analyze.png")

## ----echo=TRUE----------------------------------------------------------------
library(SC25014049)
simulate_tga <- function(temp, stages) {
  # stages是一个列表，每个元素包含: loss (总损失%), T0 (中心温度), k (速率), T_start (起始温度)
  total_loss <- 0
  for (s in stages) {
    effective_temp <- pmax(temp - s$T_start, 0)
    stage_loss <- s$loss / (1 + exp(-s$k * (temp - s$T0)))
    total_loss <- total_loss + stage_loss
  }
  mass_remaining <- 100 - total_loss
  return(mass_remaining)
}

stages <- list(
  list(loss = 5, T0 = 100, k = 0.1, T_start = 50),   # 阶段1：低温脱水
  list(loss = 40, T0 = 350, k = 0.07, T_start = 200), # 阶段2：主分解[citation:2]
  list(loss = 30, T0 = 550, k = 0.05, T_start = 400)  # 阶段3：高温碳化[citation:2]
)

temp <- seq(30, 800, length.out = 200)
mass <- simulate_tga(temp, stages)

plot(temp, mass, type = "l", col = "blue", lwd = 2,
     xlab = "Temperature (°C)", ylab = "Mass Remaining (%)",
     main = "Simulated TGA Curve (Three-Stage Decomposition)",
     ylim = c(0, 100))
grid()

rel <- tga_analyze(temperature = temp, mass = mass, n_components = 3)
print(rel)

