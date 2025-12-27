#include<Rcpp.h>
using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]
NumericMatrix matmul(const NumericMatrix &A, const NumericMatrix &B) {
  int n = A.nrow();
  int m = B.ncol();
  int k = A.ncol();
  NumericMatrix C(n, m);
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      double s = 0.0;
      for(int l = 0; l < k; l++){
        s += A(i,l) * B(l,j);
      }
      C(i,j) = s;
    }
  }
  return C;
}


// identity matrix
NumericMatrix eye(int n){
  NumericMatrix I(n,n);
  for(int i=0;i<n;i++) I(i,i)=1.0;
  return I;
}


// matrix exponential via truncated Taylor series (small matrices only)
NumericMatrix mat_exp(const NumericMatrix &A, double dt, int order = 12){
  int n = A.nrow();
  NumericMatrix X = eye(n);
  NumericMatrix term = eye(n);
  
  
  for(int k = 1; k <= order; k++){
    term = matmul(term, A);
    double coef = pow(-dt, k) / tgamma(k + 1.0);
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        X(i,j) += coef * term(i,j);
      }
    }
  }
  return X;
}


// determinant via LU decomposition (no pivoting, small matrices assumed)
double det_lu(NumericMatrix A){
  int n = A.nrow();
  double det = 1.0;
  for(int k = 0; k < n; k++){
    if(std::abs(A(k,k)) < 1e-12) return 0.0;
    det *= A(k,k);
    for(int i = k+1; i < n; i++){
      double factor = A(i,k) / A(k,k);
      for(int j = k; j < n; j++){
        A(i,j) -= factor * A(k,j);
      }
    }
  }
  return det;
}


//' Hirsch-Fye算法单步计算核心函数
//'
//' 这是Hirsch-Fye量子蒙特卡洛算法的核心计算函数，使用C++实现以获得高性能。
//' 该函数计算给定辅助场配置下的传播子和行列式。
//'
//' @param H0 单粒子哈密顿量矩阵
//' @param s 辅助场配置向量，长度为时间切片数L
//' @param U Hubbard相互作用强度
//' @param beta 逆温度
//'
//' @return 包含以下元素的列表：
//'   \item{Bup}{向上自旋的传播子矩阵}
//'   \item{Bdn}{向下自旋的传播子矩阵}
//'   \item{logdet}{对数行列式值，用于Metropolis接受率计算}
//' @examples
//' \dontrun{
//' # 构建哈密顿量
//' H0 <- build_hamiltonian(0, c(-1, 1), c(0.5, 0.5))
//' 
//' # 随机辅助场配置
//' set.seed(123)
//' s <- sample(c(-1, 1), 20, replace = TRUE)
//' 
//' # 运行单步计算
//' result <- hf_step(H0, s, U = 2.0, beta = 5.0)
//' 
//' # 查看结果
//' cat("对数行列式:", result$logdet, "\n")
//' cat("Bup矩阵维度:", dim(result$Bup), "\n")
//' }
//'
//' @export
// [[Rcpp::export]]
List hf_step(const NumericMatrix &H0,
             const IntegerVector &s,
             double U,
             double beta){

  int L = s.size();
  int N = H0.nrow();
  double dt = beta / L;
  
  
  NumericMatrix B0 = mat_exp(H0, dt);
  
  
  double lambda = acosh(exp(dt * U / 2.0));
  
  
  auto build_prod = [&](int spin){
    NumericMatrix M = eye(N);
    for(int l = 0; l < L; l++){
      NumericMatrix D = eye(N);
      D(0,0) = exp(spin * lambda * s[l]);
      M = matmul(matmul(D, B0), M);
    }
    return M;
  };
  NumericMatrix Bup = build_prod(+1);
  NumericMatrix Bdn = build_prod(-1);
  
  
  NumericMatrix I = eye(N);
  NumericMatrix Mup = matmul(I, I);
  NumericMatrix Mdn = matmul(I, I);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      Mup(i,j) += Bup(i,j);
      Mdn(i,j) += Bdn(i,j);
    }
  }
  
  
  double det_up = det_lu(Mup);
  double det_dn = det_lu(Mdn);
  
  
  double logdet = log(std::abs(det_up)) + log(std::abs(det_dn));
  
  
  return List::create(
    Named("Bup") = Bup,
    Named("Bdn") = Bdn,
    Named("logdet") = logdet
  );
}