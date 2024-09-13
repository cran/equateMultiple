\name{mathTest}
\alias{mathTest}
\docType{data}
\title{
Math Test Data
}
\description{
List of 5 datasets with the binary responses (correct/wrong) to 5 forms of a math test.
}
\usage{data("mathTest")}
\format{
	A list of length 5, containing 5 data frames. 
	Each dataset contains rows of responses from individuals to various items, 
	with the item labels as the column headers.
}
\examples{
data(mathTest)
mathTest[[1]][1:3,]
}
\keyword{datasets}
