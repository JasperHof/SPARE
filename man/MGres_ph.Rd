\name{MGres_ph}
\alias{MGres_ph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Compute martingale residuals from a coxph null model
}
\description{
   The cumulative baseline hazards are computed from the null model, based on Efrons formula. Subsequently, martingale residuals are computed.
}
\usage{
MGres_ph(
  fitph,
  data
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitph}{
    A 'coxph' object that is fitted for the null model. This object is obtained from the 'coxph' function of the survival R package.
}
  \item{data}{
    The data frame that is used for fitting the coxph object.
}
}
\details{
The data frame should include the samples under column name 'subjects'.
}
\value{
    A numeric vector containing outcome martingale residuals for all subjects
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Jasper Hof
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(survival)

# Read recurrent event data
datafile = system.file('extdata', 'Recurrencedata', package = 'SPARE')
data <- read.csv2(datafile, stringsAsFactors = FALSE)
data = as.data.frame(data)

fitph = coxph(Surv(tstart, tstop, Status) ~ X1 + X2 + X3 + frailty.gaussian(subject), data = data)
mresids = MGres_ph(fitph, data)

}

