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
temp <- seq(30, 800, length.out = 200) # 从30到800度，200个点
mass <- 100 - 0.1*(temp-100)^0.5 - 0.3*(temp-300)^0.7 - 0.4*(temp-500)^0.8

rel <- tga_analyze(temperature = temp, mass = mass, n_components = 3)

print(rel) 

