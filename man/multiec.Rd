\name{multiec}
\alias{multiec}
\alias{print.mlteqc}

\title{Multiple Equating Coefficients}
\description{Calculates the equating coefficients between multiple forms.}
\usage{
multiec(mods, base = 1, method = "mean-mean", se = TRUE, nq = 30, start = NULL, 
	iter.max = 100000, obsinf = TRUE, trace = TRUE)
}
\arguments{
  \item{mods}{an object of the class \code{\link[equateIRT]{modIRT}} containing item parameter
  coefficients and their covariance matrix of the forms to be equated.
}
\item{base}{
  integer value indicating the base form.
}
\item{method}{
  the method used to compute the equating coefficients.
  This should be one of "\code{mean-mean}", "\code{mean-gmean}", "\code{irf}",  "\code{trf}" or "\code{lik}" (see details).
}
\item{se}{
  logical; if \code{TRUE} the standard errors of the equating coefficients and
  the item parameters on a common scale are computed.
}
\item{nq}{
  number of quadrature points used for the Gauss-Hermite quadrature for
  methods "\code{irf}" or "\code{trf}".
}
\item{start}{
  initial values. This can be a vector containing the A and B equating coefficients excluding
  the base form, or an object of class \code{mlteqc} returned by function \code{multiec}.
  Used only with methods "\code{irf}", "\code{trf}" and "\code{lik}".
}
\item{iter.max}{
  maximum number of iterations allowed.
  Used only with methods "\code{irf}", "\code{trf}" and "\code{lik}".
}
\item{obsinf}{
  logical; if \code{TRUE} the standard errors are computed using the observed 
  information matrix. Used only with method "\code{lik}",
  see Battauz and  Leoncio (2023).
}
\item{trace}{
  logical; if \code{TRUE} tracing information is produced.
}
}
\details{
  The methods implemented for the computation of the multiple equating coefficients are
  the multiple mean-mean method ("\code{mean-mean}"),
  the multiple mean-geometric mean method ("\code{mean-gmean}"),
  the multiple item response function method ("\code{irf}"), 
  the multiple test response function method ("\code{trf}"),
  and likelihood-based equating ("\code{lik}").

}
\value{
	An object of class \code{mlteqc} with components
  \item{A}{A equating coefficients.}
  \item{B}{B equating coefficients.}
  \item{se.A}{standard errors of A equating coefficients.}
  \item{se.B}{standard errors of B equating coefficients.}
  \item{varAB}{covariance matrix of equating coefficients.}
  \item{as}{discrimination parameters on a common scale \eqn{\hat a_j^*}.}
  \item{bs}{difficulty parameters on a common scale \eqn{\hat b_j^*}.}
  \item{se.as}{standard errors of discrimination parameters on a common scale.}
  \item{se.bs}{standard errors of difficulty parameters on a common scale.}
  \item{tab}{data frame containing item parameter names (\code{Item}), 
	item parameter estimates of all the forms (e.g. \code{T1}, \dots, \code{T3}), 
	and item parameter estimates of all the forms converted 
	in the scale of the base form (e.g. \code{T3.as.T1}).}
  \item{varFull}{list of covariance matrices of the item parameter estimates 
  of every form.}
  \item{partial}{partial derivatives of equating coefficients 
  with respect to the item parameters.}
  \item{itmp}{number of item parameters of the IRT model.}
  \item{method}{the equating method used.}
  \item{basename}{the name of the base form.}
  \item{convergence}{An integer code. 0 indicates successful convergence.
  Returned only with methods "\code{irf}", "\code{trf}" and "\code{lik}".}
}
\references{
	Battauz, M. (2017). Multiple equating of separate IRT calibrations.
	\emph{Psychometrika}, \bold{82}, 610--636.
	
  Battauz, M., Leoncio, W. (2023). 
	A Likelihood Approach to Item Response Theory Equating of Multiple Forms
	\emph{Applied Psychological Measurement}, \bold{47}, 200-220.
	doi: 10.1177/01466216231151702.
	
  Haberman, S. J. (2009). Linking parameter estimates derived from an
  item response model through separate calibrations.
  ETS Research Report Series, 2009, i-9.
}
\author{
  Michela Battauz,
  Waldir Leoncio [ctb]
}

\seealso{
  \code{\link{score.mlteqc}}, \code{\link[equateIRT]{modIRT}}
}
\examples{
data(est2pl)
# prepare the data
mods <- modIRT(coef = est2pl$coef, var = est2pl$var, display = FALSE)
# Estimation of the equating coefficients with the multiple mean-mean method
eqMM <- multiec(mods = mods, base = 1, method = "mean-mean")
summary(eqMM)

# Estimation of the equating coefficients with the 
# multiple mean-geometric mean method (Haberman, 2009)
eqMGM <- multiec(mods = mods, base = 1, method = "mean-gmean")
summary(eqMGM)

# Estimation of the equating coefficients with the multiple item response function method
eqIRF <- multiec(mods = mods, base = 1, method = "irf")
summary(eqIRF)

# Estimation of the equating coefficients with the multiple item response function method
# using as initial values the estimates obtained with the multiple mean-geometric mean method
eqMGM <- multiec(mods = mods, base = 1, method = "mean-gmean", se = FALSE)
eqIRF <- multiec(mods = mods, base = 1, method = "irf", start = eqMGM)
summary(eqIRF)

# Estimation of the equating coefficients with the multiple test response function method
eqTRF <- multiec(mods = mods, base = 1, method = "trf")
summary(eqTRF)

# Estimation of the equating coefficients with the likelihood-based method
eqLIK <- multiec(mods = mods, base = 1, method = "lik")
summary(eqLIK)
}

