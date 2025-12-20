## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=FALSE---------------------------------------------------------------

lapply_2 <- function(FUN, ..., FUN.VALUE, USE.NAMES = TRUE) {

  temp_result <- Map(FUN, ...)
  
  if (is.matrix(FUN.VALUE)) {

    final_output <- vapply(
      temp_result, 
      function(x) if (is.matrix(x)) x else matrix(x, nrow = nrow(FUN.VALUE)),
      FUN.VALUE
    )
    return(final_output)
  } else {

    final_output <- vapply(temp_result, identity, FUN.VALUE, USE.NAMES = USE.NAMES)
    return(final_output)
  }
}

## ----echo=FALSE---------------------------------------------------------------

library(rbenchmark)
chisq_2 <- function(P4k9q, R2m6s) {
  if (length(P4k9q) != length(R2m6s)) {
    stop("Vectors must have the same length")
  }
  
  if (any(is.na(P4k9q)) || any(is.na(R2m6s))) {
    stop("Vectors must not contain missing values")
  }

  J8d1f <- table(P4k9q, R2m6s)
  
  W5t2r <- rowSums(J8d1f)
  C9n4p <- colSums(J8d1f)
  T6b8v <- sum(J8d1f)
  
  E3k7z <- outer(W5t2r, C9n4p) / T6b8v
  
  Q1x5y <- sum((J8d1f - E3k7z)^2 / E3k7z)
  
  return(Q1x5y)
}

# test
set.seed(123)
test_vector1 <- sample(1:3, 100, replace = TRUE)
test_vector2 <- sample(1:2, 100, replace = TRUE)

fast_result <- chisq_2(test_vector1, test_vector2)

standard_result <- chisq.test(test_vector1, test_vector2)$statistic

df1 = data.frame(fast_result,standard_result,abs(fast_result - standard_result))
colnames(df1) <- c("快速版本卡方统计量", "标准版本卡方统计量","结果差异")
knitr::kable(df1)

## ----echo=FALSE---------------------------------------------------------------
F7j2k_results <- benchmark(
  "标准chisq.test" = {
    chisq.test(test_vector1, test_vector2)$statistic
  },
  "快速chisq" = {
    chisq_2(test_vector1, test_vector2)
  },
  replications = 100,
  columns = c("test", "replications", "elapsed", "relative", "user.self", "sys.self")
)

# print(F7j2k_results)
knitr::kable(F7j2k_results)


## ----echo=FALSE---------------------------------------------------------------
library(rbenchmark)

fast_table_int <- function(X9s8d, Y7k9f) {

  if (!is.integer(X9s8d) || !is.integer(Y7k9f)) {
    stop("Vectors must be integer type")
  }
  

  A7s9d <- sort(unique(X9s8d))
  B8j7k <- sort(unique(Y7k9f))
  X_idx <- match(X9s8d, A7s9d)
  Y_idx <- match(Y7k9f, B8j7k)

  R8n7m <- length(A7s9d)
  C9p8o <- length(B8j7k)
  U6h5g <- (X_idx - 1) * C9p8o + Y_idx
  
  S5l4k <- tabulate(U6h5g, nbins = R8n7m * C9p8o)
  dim(S5l4k) <- c(R8n7m, C9p8o)

  rownames(S5l4k) <- as.character(A7s9d)
  colnames(S5l4k) <- as.character(B8j7k)
  
  return(S5l4k)
}

chisq_2 <- function(P4k9q, R2m6s) {
  if (length(P4k9q) != length(R2m6s)) {
    stop("Vectors must have the same length")
  }
  
  if (any(is.na(P4k9q)) || any(is.na(R2m6s))) {
    stop("Vectors must not contain missing values")
  }

  J8d1f <- fast_table_int(P4k9q, R2m6s)
  
  W5t2r <- rowSums(J8d1f)
  C9n4p <- colSums(J8d1f)
  T6b8v <- sum(J8d1f)
  
  E3k7z <- outer(W5t2r, C9n4p) / T6b8v
  
  Q1x5y <- sum((J8d1f - E3k7z)^2 / E3k7z)
  
  return(Q1x5y)
}

set.seed(123)
test_vec1 <- as.integer(sample(1:100, size = 1e6, replace = TRUE))
test_vec2 <- as.integer(sample(1:50, size = 1e6, replace = TRUE))

chisq_2_original <- function(P4k9q, R2m6s) {
  if (length(P4k9q) != length(R2m6s)) {
    stop("Vectors must have the same length")
  }
  
  if (any(is.na(P4k9q)) || any(is.na(R2m6s))) {
    stop("Vectors must not contain missing values")
  }

  J8d1f <- table(P4k9q, R2m6s)
  
  W5t2r <- rowSums(J8d1f)
  C9n4p <- colSums(J8d1f)
  T6b8v <- sum(J8d1f)
  
  E3k7z <- outer(W5t2r, C9n4p) / T6b8v
  
  Q1x5y <- sum((J8d1f - E3k7z)^2 / E3k7z)
  
  return(Q1x5y)
}

bench_res <- benchmark(
  "原始版(table())" = chisq_2_original(test_vec1, test_vec2),
  "加速版(fast_table)" = chisq_2(test_vec1, test_vec2),
  replications = 1,
  columns = c("test", "elapsed"),
  order = "elapsed"
)

knitr::kable(bench_res)

cat("\n结果一致性：", all.equal(
  chisq_2_original(test_vec1, test_vec2),
  chisq_2(test_vec1, test_vec2)
), "\n")

