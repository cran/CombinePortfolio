\name{combination.rule}
\alias{combination.rule}
\title{ Function for estimating portfolio weights by the 8fund rule }
\description{
  This function computes optimal portfolio weights based on an 8-fund rule. 
}
\usage{
combination.rule(ret, gamma=1, superset=1:7, subset=NULL, detailed.output=FALSE, 
		RHO.grid.size= 100, Kmax.init= 500, tail.cut.exp= 20)
}
\arguments{
	\item{ret}{ Matrix or data.frame of excess returns }
  \item{gamma}{Relative risk aversion parameter}
  \item{superset}{ Vector of integers from 1,2,...,7. It gives the possible included target rules, \code{1:7} provides all full 8-fund rule solutions.}
  \item{subset}{ Vector of integers of subset. It gives the target rules that must be included in the model, \code{NULL} provides all possible solutions.}
  \item{detailed.output}{ If \code{FALSE} only the estimated portfolio weight vectors of the models are returned. If \code{TRUE} a list of the portfolio weight vectors, the combination weights, and the target rules is provided.}
	\item{RHO.grid.size}{Just for convergence issues, the larger the more time-consuming, but the higher the precision of the results, only relevant if one of 5, 6 or 7 rule is included.}
	\item{Kmax.init}{See description of \code{RHO.grid.size}}
	\item{tail.cut.exp}{See description of \code{RHO.grid.size}}
 }
\details{
 The target vectors are scaled so that their weights sum up to 1. Thus target rules are interpretable, i.e. 1 = tancency, 2 = GMV and 4 = naive (1/N).
The function computes optimal portfolio weights given any combination rule of the riskfree asset and several target rule. These rules are called (and ordered) by and proportional to
\deqn{1 \equiv \widehat{\boldsymbol{\Sigma}}^{-1} \widehat{\boldsymbol{\mu}} }{'1' = Sigma^(-1) mu}
\deqn{2 \equiv \widehat{\boldsymbol{\Sigma}}^{-1} \boldsymbol{1}}{'2' = Sigma^(-1) 1}
\deqn{3 \equiv \widehat{\boldsymbol{\mu}} \ \ \ \ \ \ \ \ }{'3' = mu}
\deqn{4 \equiv \boldsymbol{1} \ \ \ }{'4' = 1}
\deqn{5 \equiv \widehat{\boldsymbol{S}}^{-2} \widehat{\boldsymbol{\mu}}}{'5' = S^(-2) mu}
\deqn{6 \equiv \widehat{\boldsymbol{S}}^{-2} \boldsymbol{1}}{'6' = S^(-2) 1}
\deqn{7 \equiv \widehat{\boldsymbol{S}}^{-1} \boldsymbol{1}}{'7' = S^(-1) 1}
where \eqn{\widehat{\boldsymbol{\mu}} }{mu} and \eqn{\widehat{\boldsymbol{\Sigma}} }{Sigma}
are the Gaussian ML-estimators of the asset mean vector \eqn{\boldsymbol{\mu} }{mu}  and the covariance matrix \eqn{\boldsymbol{\Sigma} }{Sigma}.
Moreover, we use the decomposition \eqn{\widehat{\boldsymbol{\Sigma}} = \widehat{\boldsymbol{S}}\widehat{\boldsymbol{R}}\widehat{\boldsymbol{S}} }{Sigma = SRS} with
\eqn{\widehat{\boldsymbol{R}} }{R} as sample correlation matrix and \eqn{\widehat{\boldsymbol{S}} }{S} as diagonal matrix with the sample standard deviations on the diagonal.
}
\value{
	Returns matrix of estimated weights for possible combination rules. If \code{detailed.output}  is \code{TRUE} \code{TRUE} a list of the portfolio weight vectors, the combination weights, and the target rules is provided. The names of the combination rule are coded by their portfolio that is incorporated. If "'" is contained is the name \eqn{\theta^2 }{theta^2}-adjusted estimation is used,
if   "''" is contained is the name \eqn{\theta^2 }{theta^2}-adjusted estimation is used. Hence e.g. "1'" represents the \eqn{\theta^2}{theta^2}-adjusted 2-fund rule of Kan-Zhou(2007)
and "1''2" represents the \eqn{\psi^2 }{psi^2}-adjusted 3-fund rule of Kan-Zhou(2007).
}

\author{ Florian Ziel
\cr
\email{ziel@europa-uni.de}
}

\seealso{\code{\link{combination.rule}} }

\examples{
	ret<- diff(log(EuStockMarkets))

	combination.rule(ret) ## all 8-fund rule estimates

	crule<- combination.rule(ret,gamma=5,detailed.output=TRUE)
	crule$w["1'",] ## Adjusted Kan-Zhou(2007) 2-fund rule
	crule$w["1''2",] ## Adjusted Kan-Zhou(2007) 3-fund rule
	crule$w["124",] ## Combination rule: Tangency+GMV+naive 4-fund rule, plug-in estimator
	crule$delta["124",] ## Combination weights
	crule$V[,c(1,2,4)] ## Combination targets: Tangency, GMV and naive

	## only models that can contain Tangency, GMV and naive, but must contain GMV
	crule2<- combination.rule(ret, superset=c(1,2,4), subset=2, detailed.output=TRUE) 
	crule2$w # weights
	crule2$delta # combination weights
	crule2$V # target vectors

	## case where T <= N - 4
	ret2<- cbind(ret[1:10,], ret[11:20,], ret[21:30,]) ## (TxN) 10x12-matrix
	combination.rule(ret2) ## only accessible solutions

}

\keyword{ Portfolio }
\keyword{ Combination rule }
