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
\references{
  This optional section can contain literature or other references for
  background information.
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
     temp <- seq(30, 800, length.out = 100)
     mass <- 100 - 15*exp(-((temp-200)/50)^2) - 35*exp(-((temp-350)/60)^2)
     result <- tga_analyze(temp, mass, n_components = 2)
     print(result)
     plot(result)
     
     # Simple DSC analysis
     temp <- seq(0, 300, length.out = 100)
     heat <- 10 + 5*sin(temp/50) - 8*dnorm(temp, 150, 20)
     result <- dsc_analyze(temp, heat)
     print(result)
     plot(result)
  }
}
