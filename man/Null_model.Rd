\name{Null_model}
\alias{Null_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute martingale residuals and implement saddle-point approximation
}
\description{
  First computes martingale residuals based on a null model, then computes a saddle-point approximation for the distribution of score statistics in SPARE
}
\usage{
Null_model(
  fitme,
  data,
  IDs = NULL,
  range = c(-20,20),
  length.out = 50000
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitme}{
    A 'coxme' or 'coxph' object that is fitted for the null model. This object is obtained from the 'coxme' function of the coxme R package, or the 'coxph' function of the survival package.
}
  \item{data}{
    The data frame that is used for fitting the coxme object.
}
  \item{IDs}{
    Character vector of subjects in data frame
  }
  \item{range}{
    A two-element numeric vector to specify the domain of the empirical CGF of the martingale residuals (default: (-20     ,20)).
  }
  \item{length.out}{
    A positive integer (default: 50000), specifying the number of points used to compute the empirical CGF of the        martingale residuals
  }
}
\details{
The data frame should include the samples under column name 'subjects'.
}
\value{
    Object of class 'NULL_model' which includes martingale residuals, the derivatives of the CGF and the IDs.
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

# Example code for computing a Null model
library(coxme)

# Read phenotype data
datafile = system.file('extdata', 'Recurrencedata', package = 'SPARE')
data <- read.csv2(datafile, stringsAsFactors = FALSE)
data = as.data.frame(data)

# Fit a null model
fitme = coxme(Surv(tstart, tstop, Status) ~ X1 + X2 + X3 + (1|subject), data = data, ties = 'breslow')

# Compute the null model object, including martingale residuals and saddle-point approximation
obj.null = Null_model(fitme, data, IDs = unique(data$subject))

}
