\name{check_input}
\alias{check_input}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check the input format of null object
}
\description{
Check the input format of null object
}
\usage{
check_input(data, IDs, mresid, range)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{IDs}{
%%     ~~Describe \code{IDs} here~~
}
  \item{mresid}{
%%     ~~Describe \code{mresid} here~~
}
  \item{range}{
%%     ~~Describe \code{range} here~~
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
function (data, IDs, mresid, range)
{
    if (is.null(IDs))
        stop("Argument 'IDs' is required in case of potential errors. For more information, please refer to 'Details'.")
    IDs = as.character(IDs)
    if (any(!is.element(IDs, data$subject)))
        stop("All elements in IDs should be also in data as 'subject'")
    if (anyDuplicated(IDs) != 0)
        stop("Argument 'IDs' should not have a duplicated element.")
    if (range[2] != -1 * range[1])
        stop("range[2] should be -1*range[1]")
    if (length(mresid) != length(IDs))
        stop("length(mresid)!=length(IDs) where mresid are the martingale residuals.")
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
