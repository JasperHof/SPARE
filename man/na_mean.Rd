\name{na_mean}
\alias{na_mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Impute missing genotype values
}
\description{
Replaces missing genotype values with the mean value of the SNP
}
\usage{
na_mean(x, option = "mean", maxgap = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{option}{
%%     ~~Describe \code{option} here~~
}
  \item{maxgap}{
%%     ~~Describe \code{maxgap} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
function (x, option = "mean", maxgap = Inf)
{
    data <- x
    if (!is.null(dim(data)[2]) && dim(data)[2] > 1) {
        for (i in 1:dim(data)[2]) {
            if (!anyNA(data[, i])) {
                next
            }
            tryCatch(data[, i] <- na_mean(data[, i], option,
                maxgap), error = function(cond) {
                warning(paste("imputeTS: No imputation performed for column",
                  i, "because of this", cond), call. = FALSE)
            })
        }
        return(data)
    }
    else {
        missindx <- is.na(data)
        if (!anyNA(data)) {
            return(data)
        }
        if (any(class(data) == "tbl")) {
            data <- as.vector(as.data.frame(data)[, 1])
        }
        if (all(missindx)) {
            stop("Input data has only NAs. Input data needs at least 1 non-NA data point for applying na_mean")
        }
        if (!is.null(dim(data)[2]) && !dim(data)[2] == 1) {
            stop("Wrong input type for parameter x")
        }
        if (!is.null(dim(data)[2])) {
            data <- data[, 1]
        }
        else if (option == "mean") {
            mean <- mean(data, na.rm = TRUE)
            data[missindx] <- mean
        }
        if (is.finite(maxgap) && maxgap >= 0) {
            rlencoding <- rle(is.na(x))
            rlencoding$values[rlencoding$lengths <= maxgap] <- FALSE
            en <- inverse.rle(rlencoding)
            data[en == TRUE] <- NA
        }
        if (!is.null(dim(x)[2])) {
            x[, 1] <- data
            return(x)
        }
        return(data)
    }
  }
}
