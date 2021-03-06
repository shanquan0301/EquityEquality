% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ti.R
\name{ti}
\alias{ti}
\title{Theil Index (TI) and Decomposition}
\usage{
ti(y = y, s = s, group1 = group1, group2 = group2)
}
\arguments{
\item{y}{Outcome varialbe.}

\item{s}{Variable based to rank y.}

\item{group1}{First level of cases.}

\item{group2}{Second level of cases.}
}
\description{
Calculate the theil index, and decompose it between groups.
}
\examples{
library(dineq)
data(mex_inc_2016)
mex_inc_2016 \%$\% ti(y = age,
                 s = income,
                 group1 = hh_number,
                 group2 = hh_structure)
}
\author{
Shanquan CHEN \email{shanquan0301@gmial.com}
}
\keyword{Index}
\keyword{Theil}
