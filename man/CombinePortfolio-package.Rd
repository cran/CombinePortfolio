\name{CombinePortfolio-package}
\docType{package}
\alias{CombinePortfolio}
\alias{CombinePortfolio-package}
\title{
Estimation of optimal combined portfolios based on an 8-fund rule.
}
\description{This package computes optimal portfolio weights as combination of simple portfolio strategies, like the tangency, GMV or naive (1/N). It is based on an 8-fund rule.
}
\details{
\tabular{ll}{
Package: \tab CombinePortfolio\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-06-01\cr
License: \tab GPL-3\cr
Depends: \tab R (>= 3.0), methods \cr
URL: \tab http://www.cran.r-project.org, http://www.bioconductor.org, http://www.statomics.com \cr
}

}
\author{
Code: Florian Ziel \cr
Documentation: Florian Ziel\cr
Maintainer: Florian Ziel <florian.ziel@uni-due.de>
}
\references{
(list of references)
}

\keyword{ Package }
\keyword{ Portfolio }
\keyword{ Combination rule }
\keyword{ Restricted portfolio rule}


\examples{
 
	ret<- diff(log(EuStockMarkets)) ## sample asset returns
	crule<- combination.rule(ret,detailed.output=TRUE)
	crule$w["1'",] ## Adjusted Kan-Zhou(2007) 2-fund rule
	crule$w["1''2",] ## Adjusted Kan-Zhou(2007) 3-fund rule
	crule$w["124",] ## Combination rule: Tangency+GMV+naive 4-fund rule, plug-in estimator
	crule$delta["124",] ## Combination weights
	crule$V[,c(1,2,4)] ## Combination targets: Tangency, GMV and naive (1/N)

}


