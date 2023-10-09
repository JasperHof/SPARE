\name{MGres}
\alias{MGres}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Compute martingale residuals from a coxme null model
}
\description{
   The cumulative baseline hazards are computed from the null model, based on Breslows formula. Subsequently, martingale residuals are computed.
}
\usage{
MGres(
  fitme,
  data
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitme}{
    A 'coxme' object that is fitted for the null model. This object is obtained from the 'coxme' function.
}
  \item{data}{
    The data frame that is used for fitting the coxme object.
}
}
\details{
The data frame should include the samples under column name 'subjects'.
}
\value{
    An object containing: (i) martingale residuals for all subjects; (ii) cumulative hazards for all subjects, and (iii) the frailty variance of the null model. The latter two are used for computing the approximate hazard ratios in SPARE.
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
library(coxme)

# Read recurrent event data
datafile = system.file('extdata', 'Recurrencedata', package = 'SPARE')
data <- read.csv2(datafile, stringsAsFactors = FALSE)
data = as.data.frame(data)

# Examples for four different recurrent event models:
# Fit a nonstratified calendar time model (AG model)
fitme = coxme(Surv(tstart, tstop, Status) ~ X1 + X2 + X3 + (1|subject), data = data, ties = 'breslow')
mresids = MGres(fitme, data)

# Stratified calendar time model (PWP-CT model)
fitme = coxme(Surv(tstart, tstop, Status) ~ X1 + X2 + X3 + strata(Recurrence) + (1|subject), data = data, ties = 'breslow')
mresids = MGres(fitme, data)

# Create new data set for gap time data
GTdata = data; GTdata$tstop = GTdata$tstop - GTdata$tstart

# Nonstratified gap time model (GT-UR model)
fitme = coxme(Surv(tstop, Status) ~ X1 + X2 + X3 + (1|subject), data = GTdata, ties = 'breslow')
mresids = MGres(fitme, data)

# Stratified gap time model (PWP-GT)
fitme = coxme(Surv(tstop, Status) ~ X1 + X2 + X3 + strata(Recurrence) + (1|subject), data = GTdata, ties = 'breslow')
mresids = MGres(fitme, data)

}

