\name{MGres_check}
\alias{MGres_check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check the input format of null object
}
\description{
Check the input format of null object
}
\usage{
MGres_check(fitme, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitme}{
%%     ~~Describe \code{fitme} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
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
function (fitme, data)
{
    if (is.null(fitme))
        stop("no coxme object included")
    if (class(fitme) != "coxme")
        stop("object not of class coxme")
    if (is.null(data))
        stop("no data object included")
    if (class(data) != "data.frame")
        stop("data is not a data frame object")
    if (!"subject" \%in\% colnames(data))
        stop("please include individuals as \"subject\" in dataframe")
    check_strat = strsplit(as.character(fitme$formulaList$fixed)[3],
        "strata")[[1]]
    if (length(check_strat) > 2) {
        stop("do not include stratum/covariates with name strata")
    }
    if (length(check_strat) > 1) {
        name = substring(strsplit(check_strat[2], ")")[[1]][1],
            2)
        strats = data[, which(colnames(data) == name)]
        try(if (length(strats) != dim(data)[1])
            stop("please include the strata once in the data frame"))
    }
    else {
        strats = rep(0, dim(data)[1])
    }
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
