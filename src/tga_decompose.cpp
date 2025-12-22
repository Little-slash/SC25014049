#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Rcpp;

// 辅助函数：寻找局部极值点
// [[Rcpp::export]]
std::vector<int> find_local_minima(NumericVector x, int window = 5, double threshold = 0.1) {
  int n = x.size();
  std::vector<int> minima;
  
  for (int i = window; i < n - window; ++i) {
    bool is_minimum = true;
    
    // 检查左侧窗口
    for (int j = 1; j <= window; ++j) {
      if (x[i] > x[i - j] * (1 - threshold)) {
        is_minimum = false;
        break;
      }
    }
    
    // 检查右侧窗口
    if (is_minimum) {
      for (int j = 1; j <= window; ++j) {
        if (x[i] > x[i + j] * (1 - threshold)) {
          is_minimum = false;
          break;
        }
      }
    }
    
    // 确保是最小值且为显著峰（绝对值大于阈值）
    if (is_minimum && std::abs(x[i]) > 0.01) {
      // 检查是否已经记录了附近的极小值
      if (!minima.empty() && (i - minima.back()) < window) {
        // 如果当前值比前一个记录的最小值更小，则替换
        if (x[i] < x[minima.back()]) {
          minima.back() = i;
        }
      } else {
        minima.push_back(i);
      }
    }
  }
  
  return minima;
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

// 辅助函数：计算积分（梯形法）
// [[Rcpp::export]]
double trapezoidal_integral(NumericVector x, NumericVector y, int start, int end) {
  double integral = 0.0;
  for (int i = start; i < end; ++i) {
    integral += (x[i+1] - x[i]) * (y[i] + y[i+1]) / 2.0;
  }
  return integral;
}

// 主函数：TGA分解
// [[Rcpp::export]]
List tga_decompose_cpp(NumericVector temp, NumericVector mass, int n_components = 2) {
  
  int n = temp.size();
  
  // 1. 计算DTG（微分热重曲线）
  NumericVector dtg(n-1);
  for(int i = 0; i < n-1; i++) {
    dtg[i] = (mass[i+1] - mass[i]) / (temp[i+1] - temp[i]);
  }
  
  // 2. 平滑DTG曲线（移动平均滤波）
  int smooth_window = 5;
  NumericVector dtg_smooth(n-1);
  for(int i = 0; i < n-1; i++) {
    double sum = 0.0;
    int count = 0;
    for(int j = std::max(0, i - smooth_window/2); 
        j <= std::min(n-2, i + smooth_window/2); j++) {
      sum += dtg[j];
      count++;
    }
    dtg_smooth[i] = sum / count;
  }
  
  // 3. 寻找DTG曲线的局部极小值（因为质量下降，DTG为负，极小值对应分解最快点）
  std::vector<int> min_indices = find_local_minima(dtg_smooth, 10, 0.15);
  
  // 4. 对检测到的峰按深度（DTG值）排序，选择最显著的峰
  std::vector<std::pair<double, int>> peaks_with_depth;
  for(int idx : min_indices) {
    peaks_with_depth.push_back(std::make_pair(dtg_smooth[idx], idx));
  }
  
  // 按DTG值升序排序（因为DTG为负，更小的值表示更深的峰）
  std::sort(peaks_with_depth.begin(), peaks_with_depth.end());
  
  // 5. 选择指定数量的峰（不超过n_components）
  int n_peaks = std::min((int)peaks_with_depth.size(), n_components);
  
  // 存储峰的相关信息
  NumericVector peak_temps;
  NumericVector weight_loss;
  NumericVector peak_dtg_values;
  
  // 6. 确定每个峰的边界并计算失重
  if(n_peaks > 0) {
    // 先收集峰的位置（按温度排序）
    std::vector<int> selected_indices;
    for(int i = 0; i < n_peaks; ++i) {
      selected_indices.push_back(peaks_with_depth[i].second);
    }
    std::sort(selected_indices.begin(), selected_indices.end());
    
    // 为每个峰确定边界并计算失重
    for(int i = 0; i < n_peaks; ++i) {
      int peak_idx = selected_indices[i];
      peak_dtg_values.push_back(dtg_smooth[peak_idx]);
      peak_temps.push_back(temp[peak_idx]);
      
      // 确定峰的边界
      int left_bound = 0;
      int right_bound = n-2;
      
      if(i == 0) {
        // 第一个峰：左边界是0，右边界是到下一个峰的中点或曲线拐点
        left_bound = 0;
        if(n_peaks > 1) {
          // 寻找两个峰之间的最低点作为分界
          int next_peak = selected_indices[1];
          double min_val = dtg_smooth[peak_idx];
          int min_idx = peak_idx;
          for(int j = peak_idx; j < next_peak; ++j) {
            if(dtg_smooth[j] < min_val) {
              min_val = dtg_smooth[j];
              min_idx = j;
            }
          }
          // 取两个峰之间DTG值最大的点作为分界（实际是分解最慢的点）
          double max_val = dtg_smooth[peak_idx];
          int max_idx = peak_idx;
          for(int j = peak_idx; j < next_peak; ++j) {
            if(dtg_smooth[j] > max_val) {
              max_val = dtg_smooth[j];
              max_idx = j;
            }
          }
          right_bound = max_idx;
        }
      } else if(i == n_peaks - 1) {
        // 最后一个峰：左边界是前一个峰的右边界，右边界是终点
        int prev_peak = selected_indices[i-1];
        double max_val = dtg_smooth[prev_peak];
        int max_idx = prev_peak;
        for(int j = prev_peak; j < peak_idx; ++j) {
          if(dtg_smooth[j] > max_val) {
            max_val = dtg_smooth[j];
            max_idx = j;
          }
        }
        left_bound = max_idx;
        right_bound = n-2;
      } else {
        // 中间的峰：左边界是前一个峰的右边界，右边界是到下一个峰的中点
        int prev_peak = selected_indices[i-1];
        int next_peak = selected_indices[i+1];
        
        // 与前一个峰的边界
        double max_val1 = dtg_smooth[prev_peak];
        int max_idx1 = prev_peak;
        for(int j = prev_peak; j < peak_idx; ++j) {
          if(dtg_smooth[j] > max_val1) {
            max_val1 = dtg_smooth[j];
            max_idx1 = j;
          }
        }
        left_bound = max_idx1;
        
        // 与后一个峰的边界
        double max_val2 = dtg_smooth[peak_idx];
        int max_idx2 = peak_idx;
        for(int j = peak_idx; j < next_peak; ++j) {
          if(dtg_smooth[j] > max_val2) {
            max_val2 = dtg_smooth[j];
            max_idx2 = j;
          }
        }
        right_bound = max_idx2;
      }
      
      // 计算该峰的失重（对DTG的绝对值在边界内积分）
      // 注意：DTG为负，我们取绝对值积分得到正的质量损失
      double area = 0.0;
      for(int j = left_bound; j < right_bound; ++j) {
        area += std::abs(dtg_smooth[j]) * (temp[j+1] - temp[j]);
      }
      weight_loss.push_back(area);
    }
    
    // 7. 计算残炭率和其他热稳定性参数
    // 总失重 = 起始质量 - 最终质量
    double total_weight_loss = mass[0] - mass[n-1];
    
    // 如果检测到峰，重新归一化每个峰的失重百分比
    if(total_weight_loss > 0) {
      double sum_areas = 0.0;
      for(int i = 0; i < n_peaks; ++i) {
        sum_areas += weight_loss[i];
      }
      
      // 归一化：使各峰失重之和等于总失重
      for(int i = 0; i < n_peaks; ++i) {
        weight_loss[i] = (weight_loss[i] / sum_areas) * total_weight_loss;
      }
    }
  }
  
  // 8. 计算热稳定性参数：T5和T95
  double T5 = NA_REAL;
  double T95 = NA_REAL;
  
  // 寻找失重5%和95%时的温度
  double mass_5 = mass[0] * 0.95;  // 失重5%，剩余95%
  double mass_95 = mass[0] * 0.05; // 失重95%，剩余5%
  
  for(int i = 0; i < n-1; ++i) {
    if(std::isnan(T5) && mass[i] <= mass_5) {
      // 线性插值求T5
      double t1 = temp[i-1], t2 = temp[i];
      double m1 = mass[i-1], m2 = mass[i];
      T5 = t1 + (t2 - t1) * (mass_5 - m1) / (m2 - m1);
    }
    if(std::isnan(T95) && mass[i] <= mass_95) {
      // 线性插值求T95
      double t1 = temp[i-1], t2 = temp[i];
      double m1 = mass[i-1], m2 = mass[i];
      T95 = t1 + (t2 - t1) * (mass_95 - m1) / (m2 - m1);
      break;
    }
  }
  
  // 如果没找到，设为最后一个温度
  if(std::isnan(T5)) T5 = temp[n-1];
  if(std::isnan(T95)) T95 = temp[n-1];
  
  // 9. 返回结果
  return List::create(
    _["dtg"] = dtg_smooth,
    _["peak_temperatures"] = peak_temps,
    _["weight_loss"] = weight_loss,
    _["peak_dtg_values"] = peak_dtg_values,
    _["n_components"] = n_peaks,
    _["char_yield"] = mass[n-1],  // 残炭率 = 最终质量
                          _["T5"] = T5,
                          _["T95"] = T95
  );
}