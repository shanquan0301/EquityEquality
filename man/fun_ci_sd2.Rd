% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci.R
\name{fun_ci_sd2}
\alias{fun_ci_sd2}
\title{Concentration Index (CI) and Decomposition}
\usage{
fun_ci_sd2(inc, y, x, fir.st = TRUE, nam_adj = nam_adj)
}
\arguments{
\item{inc}{Variable based to rank y.}

\item{y}{is the dependent variable.}

\item{x}{is the data frame of independent variables.}

\item{fir.st}{if whether fit for first step. Default is TRUE.}

\item{nam_adj}{variables needed to be transfered into dummy variavles. Variables list in nam_adj must be factors.}
}
\description{
Calculate the concentration index, and decompose it among factors.
}
\examples{
res <- dat_2011_rural_imp \%$\%
           fun_ci_sd2(inc = hous_inc,
           y = outp_yes,
           x = data.frame(health_professionals.10.000., te,
                          edu, live_alone, nation, smoke_year, drink_year, phy_act, liv_con, mental_num, fam_n_inc,
                          hous_inc, ncmi, pens, comi,
                          age_dec, gender, self_stat, adl_heal, pain_num, chro_num),
           nam_adj = c("edu", "self_stat"))
res <- dat_2011_rural_imp \%>\% filter(outp_tim > 0) \%$\%
         fun_ci_sd2(inc = hous_inc,
         y = outp_tim,
         x = data.frame(health_professionals.10.000., te,
                        edu, live_alone, nation, smoke_year, drink_year, phy_act, liv_con, mental_num, fam_n_inc,
                        hous_inc, ncmi, pens, comi,
                        age_dec, gender, self_stat, adl_heal, pain_num, chro_num),
            fir.st = FALSE,
            nam_adj = c("edu", "self_stat"))
}
\author{
Shanquan CHEN \email{shanquan0301@gmial.com}
}
\keyword{Concentration}
\keyword{Index}
