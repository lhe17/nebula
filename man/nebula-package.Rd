\name{nebula-package}
\alias{nebula-package}
\docType{package}
\title{
\packageTitle{nebula}
}
\description{
\packageDescription{nebula}
}
\details{
nebula is an R package for performing association analysis using a fast negative binomial mixed model for multi-subject single-cell data.
}
\author{
\packageAuthor{nebula}

Maintainer: \packageMaintainer{nebula}
}
\references{
He. L, Kulminski, A.M., NEBULA: a fast negative binomial mixed model for differential expression and co-expression analysis of large-scale single-cell data. 2020
}

\keyword{ package }

\examples{
library(nebula)
data(sample_data)
pred = model.matrix(~X1+X2+cc,data=sample_data$pred)
re = nebula(count=sample_data$count,id=sample_data$sid,pred=pred)
}
