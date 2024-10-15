\name{plot.mlteqc}
\alias{plot.mlteqc}
\title{
	Plot of Item Parameter Estimates
}
\description{
	\code{plot} method for class \code{mlteqc}, which plots item parameter 
	estimates of one form agains the item parameter estimates of the base form 
	before and after conversion.
}
\usage{
\method{plot}{mlteqc}(x, form = 2, ask = prod(par("mfcol")) < x$itmp*2 && dev.interactive(), ...)
}
\arguments{
	\item{x}{
		an object of the class \code{mlteqc} returned by function 
		\code{\link{multiec}}.
	}
	\item{form}{
	  test form. It can be specified by name or number.
	}
	\item{ask}{
	  logical; if \code{TRUE}, the user is asked before each plot.
	}
	\item{\dots}{
		further arguments passed to or from other methods.
	}
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{multiec}}
}
\examples{
data(est2pl)
# prepare the data
mods <- modIRT(coef = est2pl$coef, var = est2pl$var, display = FALSE)
# Estimation of the equating coefficients with the multiple mean-mean method
eqMM <- multiec(mods = mods, base = 1, method = "mean-mean")
plot(eqMM, form = "T2")
}
