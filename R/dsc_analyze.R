#' @title DSC Data Analysis
#'
#' @description Analyze differential scanning calorimetry (DSC) data, detect thermal transitions and calculate thermodynamic parameters.
#'
#' @param temperature Temperature vector(degree)
#' @param heat_flow Heat flow signal(mW)
#' @param baseline_method Baseline correction method("linear"and"convex")
#'
#' @return A list containing analysis results with the following components
#' 
#' @export
#' @useDynLib SC25014049, .registration = TRUE

#' @importFrom graphics abline axis legend lines mtext par points text
#' @importFrom stats sd
#' 
#' @examples
#' \dontrun{
#' temp <- seq(0, 300, length=200)
#' heat <- 0.5 - 0.3*dnorm(temp, 150, 15) + 0.8*dnorm(temp, 250, 12)
#' result <- dsc_analyze(temp, heat)
#' print(result$transitions)
#' }

dsc_analyze <- function(temperature, heat_flow, baseline_method = "linear") {
  
  baseline <- baseline_correct_cpp(temperature, heat_flow, baseline_method)
  corrected <- heat_flow - baseline
  
  transitions <- detect_transitions(temperature, corrected)
  
  if(length(transitions$peaks) > 0) {
    params <- calculate_transition_params(temperature, corrected, transitions)
  } else {
    params <- list(
      enthalpies = numeric(0),
      peak_temps = numeric(0),
      onset_temps = numeric(0)
    )
  }
  
  
  
  tg <- detect_glass_transition(temperature, corrected)
  

  structure(
    list(
      data = data.frame(
        temperature = temperature,
        heat_flow = heat_flow,
        baseline = baseline,
        corrected = corrected
      ),
      transitions = list(
        n_peaks = length(transitions$peaks),
        peak_indices = transitions$peaks,
        enthalpies = params$enthalpies,
        peak_temperatures = params$peak_temps,
        onset_temperatures = params$onset_temps,
        glass_transition = tg
      ),
      summary = data.frame(
        Parameter = c(
          "Detected peak",
          "Total enthalpy change",
          "glass transition temperature"
        ),
        Value = c(
          length(transitions$peaks),
          sum(params$enthalpies, na.rm = TRUE),
          if(!is.null(tg)) tg$midpoint else NA
        ),
        Units = c(" ", "MJ", "degree")
      )
    ),
    class = "dsc_result"
  )
}

detect_transitions <- function(temp, signal) {
  n <- length(signal)
  peaks <- c()
  

  for(i in 3:(n-2)) {

    if(signal[i] > signal[i-1] && 
       signal[i] > signal[i+1] &&
       signal[i] > mean(signal) + sd(signal)) {
      peaks <- c(peaks, i)
    }

    else if(signal[i] < signal[i-1] && 
            signal[i] < signal[i+1] &&
            signal[i] < mean(signal) - sd(signal)) {
      peaks <- c(peaks, i)
    }
  }
  
  list(peaks = peaks)
}

calculate_transition_params <- function(temp, signal, transitions) {
  enthalpies <- numeric(0)
  peak_temps <- numeric(0)
  onset_temps <- numeric(0)
  
  for(idx in transitions$peaks) {
    peak_temp <- temp[idx]
    peak_signal <- signal[idx]

    start <- max(1, idx - 20)
    end <- min(length(temp), idx + 20)

    area <- sum(diff(temp[start:end]) * 
                  (signal[start:(end-1)] + signal[(start+1):end]) / 2)
    
    enthalpies <- c(enthalpies, abs(area))
    peak_temps <- c(peak_temps, peak_temp)
    onset_temps <- c(onset_temps, temp[max(1, idx - 15)])
  }
  
  list(
    enthalpies = enthalpies,
    peak_temps = peak_temps,
    onset_temps = onset_temps
  )
}

detect_glass_transition <- function(temp, signal) {
  n <- length(signal)
  

  slopes <- diff(signal) / diff(temp)
  
  slope_changes <- diff(slopes)
  max_change_idx <- which.max(abs(slope_changes))
  
  if(length(max_change_idx) > 0) {
    tg_idx <- max_change_idx
    return(list(
      midpoint = temp[tg_idx],
      onset = temp[max(1, tg_idx - 10)],
      endset = temp[min(n, tg_idx + 10)]
    ))
  }
  
  return(NULL)
}

#' Plot DSC results
#'
#' @param x An object of class \code{dsc_result}
#' @param show_baseline Logical; whether to show the baseline
#' @param ... Additional graphical parameters
#'
#' @export
plot.dsc_result <- function(x, show_baseline = TRUE, ...) {

  plot(x$data$temperature, x$data$heat_flow,
       type = "l", lwd = 2, col = "blue",
       xlab = "temp(degree)", ylab = "HEAT-flow(mW)",
       main = "DSC-curve", ...)
  
  if(show_baseline) {
    lines(x$data$temperature, x$data$baseline,
          col = "gray", lty = 2, lwd = 1)
  }
  

  if(x$transitions$n_peaks > 0) {
    peaks <- x$transitions$peak_indices
    peak_temps <- x$data$temperature[peaks]
    
    abline(v = peak_temps, col = "red", lty = 3, lwd = 1)
    points(peak_temps, x$data$heat_flow[peaks],
           pch = 17, col = "red", cex = 1.5)
    

    text(peak_temps, x$data$heat_flow[peaks],
         labels = paste(round(x$transitions$enthalpies, 1), "mJ"),
         pos = 3, col = "red")
  }

  if(!is.null(x$transitions$glass_transition)) {
    tg <- x$transitions$glass_transition$midpoint
    abline(v = tg, col = "green", lty = 2, lwd = 2)
    text(tg, min(x$data$heat_flow),
         labels = paste("Tg =", round(tg, 1), "degree"),
         pos = 1, col = "green")
  }
  
  legend("topright",
         legend = c("DSC_data", if(show_baseline) "baseline",
                    "heat transition", "Tg"),
         col = c("blue", if(show_baseline) "gray",
                 "red", "green"),
         lty = c(1, if(show_baseline) 2, 3, 2),
         lwd = c(2, if(show_baseline) 1, 1, 2))
}

#' Print DSC analysis summary
#'
#' @param x An object of class \code{dsc_result}
#' @param ... Ignored
#'
#' @export
print.dsc_result <- function(x, ...) {
  cat("DSC analysis\n")
  cat("===========\n")
  cat(sprintf("data point: %d\n", nrow(x$data)))
  cat(sprintf("heat-tranisition: %d\n", x$transitions$n_peaks))
  
  if(!is.null(x$transitions$glass_transition)) {
    tg <- x$transitions$glass_transition
    cat(sprintf("\nglass_transition:\n"))
    cat(sprintf("midum_temp: %.1f degree\n", tg$midpoint))
    cat(sprintf("start_temp: %.1f degree\n", tg$onset))
    cat(sprintf("end_temp: %.1f degree\n", tg$endset))
  }
  
  if(x$transitions$n_peaks > 0) {
    cat("\nheat-tranisition-peak:\n")
    for(i in 1:x$transitions$n_peaks) {
      cat(sprintf("peak%d: %.1f degree, dH = %.2f MJ\n",
                  i, x$transitions$peak_temperatures[i],
                  x$transitions$enthalpies[i]))
    }
    cat(sprintf("\ntotal enthalpy change: %.2f mj\n", sum(x$transitions$enthalpies)))
  }
}
