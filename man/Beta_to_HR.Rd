\name{Beta_to_HR}
\alias{Beta_to_HR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert score statistic to hazard ratio
}
\description{
Convert score statistic to hazard ratio
}
\usage{
Beta_to_HR(beta, g, resids, cumhaz, frail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
  Score statistics
}
  \item{g}{
  Genotype values
}
  \item{resids}{
  Martingale residuals
}
  \item{cumhaz}{
  Cumulative hazards
}
  \item{frail}{
  Frailty variance of null model
}
}
\details{
Hazard ratios are computed, using the SNP values, frailties and cumulative hazard
}
\value{
Hazard ratios
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (beta, g, resids, cumhaz, frail)
{
    constant = sum(g^2 * cumhaz/(frail * cumhaz + 1))/sum(g^2)
    HR = beta/constant
    return(HR)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
