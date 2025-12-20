#' @title TGA Data Analysis
#'
#' @description Analyze thermogravimetric analysis (TGA) data, decompose multi-component systems, and calculate thermal stability parameters.
#'
#' @param temperature Temperature vector(degree)
#' @param mass mass
#' @param n_components Expected number of components(1-4)
#' 
#' @return A list containing analysis results with the following components
#' 
#' 
#' @export
#' @useDynLib SC25014049, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @examples
#' \dontrun{
#' temp <- seq(30, 800, length=200)
#' mass <- 100 - 0.1*(temp-100)^0.5 - 0.3*(temp-300)^0.7 - 0.4*(temp-500)^0.8
#' result <- tga_analyze(temp, mass, n_components=3)
#' print(result$summary)
#' plot(result)
#' }

tga_analyze <- function(temperature, mass, n_components = 2) {
  
  if(length(temperature) != length(mass)) {
    stop("temperature and mass must length same")
  }
  
  result <- tga_decompose_cpp(temperature, mass, n_components)
  
  n_peaks <- result$n_components
  
  summary <- data.frame(
    Parameter = c(
      "start mass(%)",
      "re_carbon(%)",
      "total weight loss(%)",
      if(n_peaks > 0) paste("component", 1:n_peaks, "weight loss(%)"),
      if(n_peaks > 0) paste("component", 1:n_peaks, "peak temp(degree)"),
      "T5%(degree)",
      "T95%(degree)"
    ),
    Value = c(
      mass[1],
      result$char_yield,
      100 - result$char_yield,
      if(n_peaks > 0) round(result$weight_loss, 2),
      if(n_peaks > 0) round(result$peak_temperatures, 1),
      round(result$T5, 1),
      round(result$T95, 1)
    ),
    stringsAsFactors = FALSE
  )
  
  summary <- summary[!is.na(summary$Value), ]

  structure(
    list(
      data = data.frame(
        temperature = temperature,
        mass = mass,
        dtg = c(result$dtg, NA)
      ),
      decomposition = list(
        n_components = n_peaks,
        weight_loss = result$weight_loss,
        peak_temperatures = result$peak_temperatures
      ),
      stability = list(
        char_yield = result$char_yield,
        T5 = result$T5,
        T95 = result$T95
      ),
      summary = summary,
      settings = list(
        n_components = n_components
      )
    ),
    class = "tga_result"
  )
}

#' Plot TGA results
#'
#' @param x An object of class \code{tga_result}
#' @param type Plot type: "both", "tga", or "dtg"
#' @param ... Additional graphical parameters
#'
#' @export
plot.tga_result <- function(x, type = c("both", "tga", "dtg"), ...) {
  type <- match.arg(type)
  
  if(type %in% c("both", "tga")) {
    plot(x$data$temperature, x$data$mass,
         type = "l", lwd = 2, col = "blue",
         xlab = "temp(degree)", ylab = "weight (%)",
         main = "TGA curve", ...)
    if(x$decomposition$n_components > 0) {
      peaks <- x$decomposition$peak_temperatures
      abline(v = peaks, col = "red", lty = 2, lwd = 1)
      text(peaks, max(x$data$mass), 
           labels = paste(round(x$decomposition$weight_loss, 1), "%"),
           pos = 3, col = "red")
    }
  }
  
  if(type == "both") {
    par(new = TRUE)
    plot(x$data$temperature, x$data$dtg,
         type = "l", lwd = 1.5, col = "darkgreen",
         xlab = "", ylab = "", axes = FALSE)
    axis(4, col = "darkgreen")
    mtext("DTG (%/degree)", side = 4, line = 3, col = "darkgreen")
    
    legend("topright", 
           legend = c("TGA", "DTG"),
           col = c("blue", "darkgreen"),
           lty = 1, lwd = 2)
  } else if(type == "dtg") {
    plot(x$data$temperature, x$data$dtg,
         type = "l", lwd = 2, col = "darkgreen",
         xlab = "temp(degree)", ylab = "DTG (%/degree)",
         main = "DTG curve", ...)
  }
}

#' Print TGA analysis summary
#'
#' @param x An object of class \code{tga_result}
#' @param ... Ignored
#'
#' @export
print.tga_result <- function(x, ...) {
  cat("TGA analysis\n")
  cat("===========\n")
  cat(sprintf("data point: %d\n", nrow(x$data)))
  cat(sprintf("detected component: %d\n", x$decomposition$n_components))
  cat("\nThermal stability parameter:\n")
  cat(sprintf("carbon yield: %.1f %%\n", x$stability$char_yield))
  cat(sprintf("T5%%: %.1f degree\n", x$stability$T5))
  cat(sprintf("T95%%: %.1f degree\n", x$stability$T95))
  
  if(x$decomposition$n_components > 0) {
    cat("\nComponent decomposition:\n")
    for(i in 1:x$decomposition$n_components) {
      cat(sprintf("component%d: %.1f %%weight loss, peak-temp: %.1f degree\n",
                  i, x$decomposition$weight_loss[i],
                  x$decomposition$peak_temperatures[i]))
    }
  }
}
