\name{check_input_SPARE}
\alias{check_input_SPARE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check the input for the SPARE function
}
\description{
Check the input for the SPARE function
}
\usage{
check_input_SPARE(obj.null, Geno.mtx, par.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.null}{
%%     ~~Describe \code{obj.null} here~~
}
  \item{Geno.mtx}{
%%     ~~Describe \code{Geno.mtx} here~~
}
  \item{par.list}{
%%     ~~Describe \code{par.list} here~~
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
function (obj.null, Geno.mtx, par.list)
{
    if (class(obj.null) != "NULL_Model")
        stop("obj.null should be a returned outcome from 'Null_model'")
    if (is.null(rownames(Geno.mtx)))
        stop("Row names of 'Geno.mtx' should be given.")
    if (is.null(colnames(Geno.mtx)))
        stop("Column names of 'Geno.mtx' should be given.")
    if (!is.numeric(Geno.mtx) | !is.matrix(Geno.mtx))
        stop("Input 'Geno.mtx' should be a numeric matrix.")
    if (length(intersect(obj.null$IDs, rownames(Geno.mtx))) ==
        0)
        stop("None of 'IDs' are included in rownames(Geno.mtx).")
    print(paste0("In total, ", length(intersect(obj.null$IDs,
        rownames(Geno.mtx))), " samples with phenotype and genotype information"))
    if (!is.numeric(par.list$min.maf) | par.list$min.maf < 0 |
        par.list$min.maf > 0.5)
        stop("Argument 'min.maf' should be a numeric value >= 0 and <= 0.5.")
    if (!is.numeric(par.list$missing.cutoff) | par.list$missing.cutoff <
        0 | par.list$missing.cutoff > 1)
        stop("Argument 'missing.cutoff' should be a numeric value between 0 and 1.")
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
