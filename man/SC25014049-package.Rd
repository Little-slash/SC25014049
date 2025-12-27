\name{SC25014049-package}
\alias{SC25014049-package}
\docType{package}
\title{
  A short title line describing what the package does
}
\description{
  A more detailed description of what the package does. A length
  of about one to five lines is recommended.
}
\details{
  This section should provide a more detailed overview of how to use the
  package, including the most important functions.
}
\author{
SC25014049
Maintainer: Fenghao Lu <lfh3485188464@163.com>
}
\keyword{ package }
\examples{
  \dontrun{
     res <- simulate_hirsch_fye(eps_d = 0.0, eps_k = c(-1,1), V_k = c(0.5,0.5), U = 2.0,
     beta = 5.0, L = 40, n_steps = 5000, thermal = 500, measure_interval = 5,
     seed = 123)
     mean(res$measurements$n_up)
     mean(res$measurements$n_dn)
     plot(res$measurements$n_up, type = 'l', main = 'n_up samples')
  }
}
