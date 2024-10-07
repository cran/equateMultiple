\name{item.common}
\alias{item.common}
\alias{item.common.mlteqc}
\title{
	Estimates the Item Parameters on a Common Scale
}
\description{
	\code{item.common} is a generic function which extracts the estimates of the 
	item parameters on the scale of the base form, which is obtained using 
	the item parameter estimates across all forms.
}
\usage{
item.common(x, \dots)

\method{item.common}{mlteqc}(x, \dots)
}
\arguments{
	\item{x}{object of the class \code{mlteqc} returned by function 
	  \code{\link{multiec}}} 
	\item{\dots}{further arguments passed to or from other methods.}
}
\value{
	A data frame containing item names (\code{Item}),
	estimated value of the item parameter (\code{Estimate})
	and its standard error (\code{StdErr}).
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link[equateMultiple]{itm}}, \code{\link{multiec}}
}
\examples{
data(est2pl)
# prepare the data
mods <- modIRT(coef = est2pl$coef, var = est2pl$var, display = FALSE)
# Estimation of the equating coefficients with the multiple item response function method
eqIRF <- multiec(mods = mods, base = 1, method = "irf")

# estimates of the item parameters on the scale of the base form
item.common(eqIRF)
}
