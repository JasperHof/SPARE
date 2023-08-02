\name{SPARE.bgen}
\alias{SPARE.bgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Saddlepoint approximation for recurrent events for .bgen files.
}
\description{
A fast and accurate method for a recurrent event GWAS
}
\usage{
SPARE.bgen(
  bgenfile,
  gIDs,
  obj.null,
  output.file = NULL,
  chr = NULL,
  missing.cutoff = 0.05,
  min.maf = 0.05,
  p.cutoff = 0.001,
  memory = 512,
  maxchunksize = 5e4,
  backingdir = 'Connections',
  backingfile = 'backfile'
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bgenfile}{
    character, the name of the .bgen file, without the '.bgen' extension
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
  \item{backingdir}{
  character, name of the directory used for the backingfiles
  }
  \item{backingfile}{
  name of the backingfiles that are used to load the .bgen data
  }
}
\details{
To run SPARE for .bgen files, the user should first fit a null model based on the phenotype data. This null model computes martingale residuals, and subsequently the distribution of martingale residuals is empirically estimated using a the saddle-point approximation.

After computing the null model, SPARE.bgen() computes association statistics between SNPs in the .bed file and martingale residuals.
}

\value{
The output file contains a matrix with the following columns:
\item{rsid}{
  The rs-id of the SNP
}
\item{SNP}{
  Identifier of SNP in 'CHR_BP_A1_A2' format
}
\item{Chr}{
  Chromosome
}
\item{BP}{
  Position on the chromosome, in base pair
}
\item{A1}{
  A1 allele
}
\item{A2}{
  A2 allele
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
Application is similar to the SPARE.bed() function for .bed files, however the sample IDs of the .bgen file should be exactly the same as the 'gIDs' argument.
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

