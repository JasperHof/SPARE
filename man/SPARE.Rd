\name{SPARE}
\alias{SPARE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Saddlepoint approximation for recurrent events
}
\description{
Saddlepoint approximation for recurrent events based on input null model and genotype matrix
}
\usage{
SPARE(
  obj.null,
  Geno.mtx,
  missing.cutoff = 0.05,
  min.maf = 0.05,
  p.cutoff = 0.001
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.null}{
    an R object returned from 'Null_model()'
}
  \item{Geno.mtx}{
    Genotype matrix, with subjects IDs in row names and SNP names in column names
}
  \item{missing.cutoff}{
  Cut-off rate for number of missing genotypes
  }
  \item{min.maf}{
  Minimum minor allele frequency
  }
  \item{p.cutoff}{
  Cut-off P value for implementing a saddle-point approximation.
  }
}
\details{
To run SPARE for a matrix of genotype data, the user should first fit a null model based on the phenotype data. This null model is used to compute martingale residuals, and their distribution is empirically estimated using a the saddle-point approximation.

Subsequently, this function can be used to compute association statistics between SNPs in the genotype matrix and martingale residuals.
}
\value{
The output file contains a matrix with the following columns
\item{SNP}{
  Name of SNP
}
\item{Beta}{
  Test statistic
}
\item{pMG}{
  P value obtained from testing the test statistic without SPA correction
}
\item{pSPA}{
  P value obtained from testing the test statistic with SPA (recommended)
}
\item{SE}{
  Standard error of test statistic
}
\item{Z}{
  Z-statistic corresponding to the test statistic
}
\item{MAF}{
  Minor allele frequency of SNP
}
\item{Missing}{
  Fraction of missing SNPs
}
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

# Read phenotype data
datafile = system.file('extdata', 'Recurrencedata', package = 'SPARE')
data <- read.csv2(datafile, stringsAsFactors = FALSE)
data = as.data.frame(data)

# Generate genotype data in matrix
N = length(unique(data$subject))
m_loci = 10000
MAF = 0.3
Geno.mtx = matrix(rbinom(N * m_loci, 2, prob = MAF), nrow = N)
rownames(Geno.mtx) = unique(data$subject)
colnames(Geno.mtx) = paste0('rs',1:m_loci)

# Fit a null model
fitme = coxme(Surv(tstart, tstop, Status) ~ X1 + X2 + X3 + (1|subject), data = data, ties = 'breslow')

# Compute the null model object
obj.null = Null_model(fitme, data, IDs = unique(data$subject))

# Carry out the saddlepoint approximation for recurrent events
Outcome = SPARE(obj.null, Geno.mtx)

}
