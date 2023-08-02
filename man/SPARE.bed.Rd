\name{SPARE.bed}
\alias{SPARE.bed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Saddlepoint approximation for recurrent events for .bed files
}
\description{
A fast and accurate method for a recurrent event GWAS
}
\usage{
SPARE.bed(
  bedfile,
  gIDs,
  obj.null,
  output.file = NULL,
  chr = NULL,
  missing.cutoff = 0.05,
  min.maf = 0.05,
  p.cutoff = 0.001,
  memory = 512,
  maxchunksize = 5e4
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bedfile}{
    character, the name of the .bed file, without the '.bed' extension
}
  \item{gIDs}{
    character vector of genotype IDs, in the same order as the .bgen file
}
  \item{obj.null}{
    an R object returned from 'Null_model()'
}
  \item{output.file}{
    character, the name of the output file for GWAS summary statistics
  }
  \item{chr}{
  chromosome number
  }
  \item{missing.cutoff}{
  cut-off rate for number of missing genotypes
  }
  \item{min.maf}{
  minimum minor allele frequency
  }
  \item{p.cutoff}{
  cut-off P value for implementing a saddle-point approximation
  }
  \item{memory}{
  memory of the device in MB, used to compute chunk size
  }
  \item{maxchunksize}{
  maximum size of chunks of SNPs to be simultaneously analyzed in SPARE
  }
}
\details{
To run SPARE for .bed files, the user should first fit a null model based on the phenotype data. This null model computes martingale residuals, and subsequently the distribution of martingale residuals is empirically estimated using a the saddle-point approximation.

After computing the null model, SPARE.bed() computes association statistics between SNPs in the .bed file and martingale residuals.
}
\value{
The output file contains a matrix with the following columns
\item{Chr}{
  Chromosome number
}
\item{SNP}{
  Name of SNP
}
\item{A1}{
  A1 allele
}
\item{A2}{
  A2 allele
}
\item{BP}{
  Position on chromosome in base pair
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
\item{log_HR_approx}{
  Approximate log hazard ratio of SNP
}
\item{SE2}{
  Standard error of test statistic, corrected using the SPA-based P value
}
\item{log_HR_approx_SE}{
  Standard error of log hazard ratio
}
\item{log_HR_approx_SE2}{
  SPA-corrected standard error of log hazard ratio
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
# Load required packages
library(coxme)
library(seqminer)

# Load genotype and phenotype data
famfile = system.file('extdata', 'bedfile.fam', package = 'SPARE')
fam = read.table(famfile)

bedfile = gsub('.fam', '', famfile)

datafile = system.file('extdata', 'Recurrencedata', package = 'SPARE')
data <- read.csv2(datafile, stringsAsFactors = FALSE)
data = as.data.frame(data)

# Fit a non-stratified calendar model as null model
fitme = coxme(Surv(tstart, tstop, Status) ~ X1 + X2 + X3 + (1|subject), data = data, ties = 'breslow')

# Compute the martingale residuals and saddle-point approximation
obj.null = Null_model(fitme, data, IDs = unique(data$subject))

# Perform association tests between .bed file and martingale residuals. Output is written to 'bedoutput'
SPARE.bed(bedfile, fam$v2 , obj.null, 'bedoutput')
}
