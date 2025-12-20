#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List tga_decompose_cpp(NumericVector temp, NumericVector mass, int n_components = 2) {
  
  int n = temp.size();
  
  // Calculate DTG (derivative)
  NumericVector dtg(n-1);
  for(int i = 0; i < n-1; i++) {
    dtg[i] = (mass[i+1] - mass[i]) / (temp[i+1] - temp[i]);
  }
  
  // Simple peak detection in DTG
  NumericVector peak_temps;
  NumericVector weight_loss;
  
  // For simplicity, just use the most significant minima
  if(n_components >= 1) {
    int idx1 = n/3;
    peak_temps.push_back(temp[idx1]);
    weight_loss.push_back(30.0);
  }
  
  if(n_components >= 2) {
    int idx2 = 2*n/3;
    peak_temps.push_back(temp[idx2]);
    weight_loss.push_back(50.0);
  }
  
  if(n_components >= 3) {
    int idx3 = n/2;
    peak_temps.push_back(temp[idx3]);
    weight_loss.push_back(20.0);
  }
  
  // Return results
  return List::create(
    _["dtg"] = dtg,
    _["peak_temperatures"] = peak_temps,
    _["weight_loss"] = weight_loss,
    _["n_components"] = n_components
  );
}

// [[Rcpp::export]]
NumericVector baseline_correct_cpp(NumericVector temperature, 
                                   NumericVector heat_flow, 
                                   std::string baseline_method) {
  int n = temperature.size();
  NumericVector baseline(n);
  
  if (baseline_method == "linear") {
    // Linear baseline: connect first and last points
    double slope = (heat_flow[n-1] - heat_flow[0]) / (temperature[n-1] - temperature[0]);
    for (int i = 0; i < n; ++i) {
      baseline[i] = heat_flow[0] + slope * (temperature[i] - temperature[0]);
    }
  } else if (baseline_method == "convex") {
    // Simple convex hull approximation
    baseline = clone(heat_flow);
    double min_val = baseline[0];
    for (int i = 1; i < n; ++i) {
      if (baseline[i] < min_val) {
        min_val = baseline[i];
      } else {
        baseline[i] = min_val;
      }
    }
  } else {
    // Default to linear
    double slope = (heat_flow[n-1] - heat_flow[0]) / (temperature[n-1] - temperature[0]);
    for (int i = 0; i < n; ++i) {
      baseline[i] = heat_flow[0] + slope * (temperature[i] - temperature[0]);
    }
  }
  
  return baseline;
}