# dependencies:
library(coxme)
library(bigsnpr)
library(RSQLite)
library(seqminer)

#
#
#           Functions to be included in the SPARE package
#
#

MGres_ph <- function (fitph = NULL, data = NULL)
{
  if (is.null(fitph))  stop("no coxph object included")
  if (!"coxph" %in% class(fitph)) stop("object not of class coxph")
  if (is.null(data))  stop("no data object included")
  if (!"subject" %in% colnames(data))  stop("please include individuals as \"subject\" in dataframe")
  if (length(grep("subject", attr(fitph$terms, "term.labels"))) < 1)  warning("subject not included as frailty")
  mgres = rowsum(fitph$residuals, data$subject)
  mgres = mgres[match(unique(data$subject), rownames(mgres)),]

  count = rep(0,length(unique(data$subject)))
  for(i in 1:length(count)){
    count[i] = sum(as.matrix(fitph$y[which(data$subject == unique(data$subject)[i]),3]))
  }

  cumhaz = count - mgres

  names(mgres) = unique(data$subject)
  names(cumhaz) = unique(data$subject)

  object = list(resids = mgres, cumhaz = cumhaz, frail = as.numeric(unlist(fitph$history)[1]))

  return(object)
}

MGres <- function (fitme = NULL, data = NULL)
{
  MGres_check(fitme, data)

  ### Find the strata included in the coxme fit
  check_strat = strsplit(as.character(fitme$formulaList$fixed)[3],
                         "strata")[[1]]
  if (length(check_strat) > 1) {
    name = substring(strsplit(check_strat[2], ")")[[1]][1],
                     2)
    strats = data[, which(colnames(data) == name)]
  } else strats = rep(0, dim(data)[1])
  strat_list = unique(strats)

  cumhaz = rep(0, length(unique(data$subject)))

  ### Compute baseline hazards for all strata
  events = as.data.frame(as.matrix(fitme$y))
  for (strat in strat_list) {
    in_strat = which(strats == strat)

    ### Distinguish models where 'start' time of risk interval is specified
    if (dim(events)[2] > 2) {
      basehaz = as.data.frame(cbind(rep(0, length = length(unique(events$stop))), unique(events$stop)))
      colnames(basehaz) = c('Haz', 'Time')
      eventtimes = events$stop

      for (i in 1:length(unique(events$stop))) {
        time = basehaz$Time[i]
        risk.set = which(events$start < time & events$stop >=
                           time & strats == strat)
        denom = sum(exp(fitme$linear.predictor[risk.set]))
        nom = sum(events$stop == time & events$status ==
                    1 & strats == strat)
        basehaz$Haz[i] = nom/denom
      }
    } else {
      basehaz = as.data.frame(cbind(rep(0, length = length(unique(events$time))), unique(events$time)))
      colnames(basehaz) = c('Haz', 'Time')
      eventtimes = events$time
      for (i in 1:length(unique(events$time))) {
        time = basehaz$Time[i]
        risk.set = which(events$time >= time & strats ==
                           strat)
        denom = sum(exp(fitme$linear.predictor[risk.set]))
        nom = sum(events$time == time & events$status ==
                    1 & strats == strat)
        basehaz$Haz[i] = nom/denom
      }
    }

    ### Compute individual cumulative hazards from baseline hazard and linear predictor
    for (i in 1:length(unique(data$subject))) {
      rows = which(data$subject == unique(data$subject)[i] &
                     strats == strat)
      if (length(rows) > 0) {
        for (r in 1:length(rows)) {
          if (dim(events)[2] > 2) {
            haz = sum(basehaz$Haz[which(events$start[rows[r]] <
                                          basehaz$Time & basehaz$Time <= events$stop[rows[r]])])
            cumhaz[i] = cumhaz[i] + haz * exp(fitme$linear.predictor[rows[r]])
          } else {
            haz = sum(basehaz$Haz[which(basehaz$Time <= events$time[rows[r]])])
            cumhaz[i] = cumhaz[i] + haz * exp(fitme$linear.predictor[rows[r]])
          }
        }
      }
    }
  }
  individual_prop_haz = rep(0, length(unique(data$subject)))
  count = rep(0, length(unique(data$subject)))
  for (i in 1:length(unique(data$subject))) {
    set = which(data$subject == unique(data$subject)[i])
    count[i] = sum(events$status[set] == 1)
  }
  resids = count - cumhaz

  names(resids) = unique(data$subject)
  names(cumhaz) = unique(data$subject)

  object = list(resids = resids, cumhaz = cumhaz, frail = as.numeric(fitme$vcoef))

  return(object)
}

#
# FUNCTION FOR ASSESSING NULL MODEL + SPA
#

Null_model <- function(fitme,
                       data,
                       IDs=NULL,
                       range=c(-20,20),
                       length.out = 50000)
{
  Call = match.call()

  # Compute martingale residuals and cumulative hazards
  if ("coxme" %in% class(fitme)){

    MG = MGres(fitme, data)

    mresid = MG$resids
    cumhaz = MG$cumhaz
    frail = MG$frail
  }
  if ("coxph" %in% class(fitme)){

    MG = MGres_ph(fitme, data)

    mresid = MG$resids
    cumhaz = MG$cumhaz
    frail = MG$frail
  }

  ### Check input arguments
  obj.check = check_input(data, IDs, mresid, range)

  ### Calculate empirical CGF for martingale residuals
  idx0 = qcauchy(1:length.out/(length.out+1))
  idx1 = idx0 * max(range) / max(idx0)

  cumul = NULL
  print("Compute empirical CGF for martingale residuals...")
  c = 0
  for(i in idx1){
    c = c+1
    t = i
    e_resid = exp(mresid*t)
    M0 = mean(e_resid)
    M1 = mean(mresid*e_resid)
    M2 = mean(mresid^2*e_resid)
    K0 = log(M0)
    K1 = M1/M0
    K2 = (M0*M2-M1^2)/M0^2
    cumul = rbind(cumul, c(t, K0, K1, K2))
    if(c %% 5000 == 0) print(paste0("Complete ",c,"/",length.out,"."))
  }

  K_org_emp = approxfun(cumul[,1], cumul[,2], rule=2)
  K_1_emp = approxfun(cumul[,1], cumul[,3], rule=2)
  K_2_emp = approxfun(cumul[,1], cumul[,4], rule=2)

  re = list(resid = mresid, cumhaz = cumhaz, frail = frail, K_org_emp = K_org_emp, K_1_emp = K_1_emp,
            K_2_emp = K_2_emp, Call = Call, IDs = IDs)

  class(re)<-"NULL_Model"
  return(re)
}

#
# Apply a SaddlePoint Approximation for Recurrent Event analysis (SPARE)
#

SPARE = function(obj.null,
               Geno.mtx,
               missing.cutoff = 0.05,
               min.maf = 0.05,
               p.cutoff = 0.001)
{
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf)

  ### Check input
  check_input_SPARE(obj.null, Geno.mtx, par.list)

  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  IDs = rownames(Geno.mtx)
  SNPs = colnames(Geno.mtx)

  #
  ### -------------- Quality Control ---------------
  #

  ### Only select genotypes that have a phenotype
  mresid = obj.null$resid
  Complete = intersect(names(mresid), IDs)
  Geno.mtx = Geno.mtx[which(rownames(Geno.mtx) %in% Complete),]

  ### Filter on call rate, minor allele frequency, and match entries of genotype and phenotype
  Geno.mtx[Geno.mtx == -9] = NA                  # for plink input
  SNP.callrate = colMeans(is.na(Geno.mtx))
  if(any(SNP.callrate > missing.cutoff)){ Geno.mtx = Geno.mtx[,-which(SNP.callrate > missing.cutoff)] }

  No_Geno <- is.na(Geno.mtx)
  Geno.mtx = na_mean(Geno.mtx)                # Impute missing values to mean for estimating MAF

  G = Geno.mtx[,which(colMeans(Geno.mtx) > 2*min.maf & colMeans(Geno.mtx) < 2*(1 - min.maf))]
  No_Geno = No_Geno[,which(colMeans(Geno.mtx) > 2*min.maf & colMeans(Geno.mtx) < 2*(1 - min.maf))]  # Also update the No_geno matrix

  if(dim(G)[2] == 0) return(NULL)               # No SNPs left after QC

  mresid = mresid[which(names(mresid) %in% Complete)]
  MG = mresid[match(rownames(G), names(mresid))]

  ### Perform Linear regression for all SNPs simultaneously using matrix multiplication

  print(paste0("SPARE started at ", Sys.time()))

  ### Define outcome Y and intercept X
  n = dim(G)[1]
  X = as.matrix(rep(1,dim(G)[1]), nrow = dim(G)[1],ncol=1)
  Y = as.matrix(MG)

  ### Scale genotypes to have mean 0 and compute slopes
  Str = G - as.matrix(X %*% colMeans(G))
  if(sum(No_Geno) > 0)  Str[No_Geno] <- 0
  b = as.numeric(crossprod(Y, Str))/as.numeric(colSums(Str ^ 2))

  ### Compute P values in SPARE
  S_sq = colSums(Str ^ 2)
  RSS = crossprod(Y^2, !No_Geno) - b ^ 2 * S_sq
  sigma_hat = RSS/(n - 2)
  error = sqrt(sigma_hat/ S_sq)
  pval.mg = as.numeric(2 * pnorm(-abs (b / error)))

  ### Prepare the output file
  b = as.numeric(b); error = as.numeric(error)

  # Compute approximate Hazard Ratio
  HR = rep(NA, length(b))
  for(i in 1:length(HR)){
    HR[i] = Beta_to_HR(b[i], G[,i], mresid, obj.null$cumhaz, obj.null$frail)
  }

  outcome = cbind(SNP = colnames(G), MAF = colMeans(G)/2, Missing = colSums(is.na(G)),
                  pSPA = pval.mg, pMG = pval.mg, log_HR_approx = HR ,Beta = b, SE = error,
                  Z = b/error)
  outcome = as.data.frame(outcome, stringsAsFactors = F)
  outcome = transform(outcome, Beta = as.numeric(Beta), log_HR_approx = as.numeric(log_HR_approx), pMG = as.numeric(pMG),
                      pSPA = as.numeric(pSPA), SE = as.numeric(SE), Z = as.numeric(Z),
                      MAF = as.numeric(MAF), Missing = as.integer(Missing))
  rownames(outcome) = NULL
  low.p = which(pval.mg < p.cutoff)

  # Perform SPA for SNPs with low P values
  for (i in low.p) {
    g = G[, i]
    g = g - mean(g)
    s = sum(g * MG)/sum(g^2)
    k0 = function(x) sum(obj.null$K_org_emp(g/sum(g^2) *
                                              x))
    k1 = function(x) sum(g/sum(g^2) * obj.null$K_1_emp(g/sum(g^2) *
                                                         x))
    k2 = function(x) sum((g/sum(g^2))^2 * obj.null$K_2_emp(g/sum(g^2) *
                                                             x))
    get_p = function(s, tail) {
      k1_root = function(x) k1(x) - s
      zeta = uniroot(k1_root, c(-2000, 2000))$root
      w = sign(zeta) * sqrt(2 * (zeta * s - k0(zeta)))
      v = zeta * sqrt(k2(zeta))
      pval = pnorm(w + 1/w * log(v/w), lower.tail = tail)
      return(pval)
    }
    p_SPA = get_p(abs(s), tail = F) + get_p(-abs(s), tail = T)
    outcome$pSPA[i] = p_SPA
  }

  # Compute SPA - corrected standard errors. Effects have to be slightly way from the null to avoid dividing by 0
  SNP_set = which(outcome$pMG < 0.75)
  SE2 = outcome$SE
  SE2[SNP_set] = sqrt(outcome$Beta[SNP_set]^2 / qchisq(outcome$pSPA[SNP_set],df=1, lower.tail = F))

  # Compute (SPA - corrected) standard errors for the approximate HRs
  log_HR_approx_SE = outcome$SE * (outcome$log_HR_approx / outcome$Beta)
  log_HR_approx_SE2 = SE2 * (outcome$log_HR_approx / outcome$Beta)

  outcome = cbind(outcome, SE2, log_HR_approx_SE, log_HR_approx_SE2)

  print(paste0("SPARE completed at ", Sys.time()))
  return(outcome)
}

#
# SPARE for .bgen files
#

SPARE.bgen <- function(bgenfile, gIDs,
                     obj.null,
                     output.file = NULL,
                     chr = NULL,
                     missing.cutoff = 0.05,
                     min.maf = 0.05,
                     p.cutoff = 0.001,
                     memory = 512,
                     maxchunksize = 5e4,
                     backingdir = 'Connections',
                     backingfile = 'backfile')
{

  bgenfile = paste0(bgenfile, '.bgen')
  bgifile = paste0(bgenfile,'.bgi')

  if(!file.exists(bgenfile)) stop("Could not find .bgen file")
  if(!file.exists(bgifile)) stop("Could not find .bgen.bgi file")
  if(!dir.exists(backingdir)){
    print(paste0('Creating directory ',paste0(getwd(),'/',backingdir),' for the backingfiles of the .bgen file'))
    dir.create(backingdir)
  }
  if(is.null(output.file)) stop('please provide name of output file')

  ### Create 'myid' variable for reading .bgen SNPs
  db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgifile)
  on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
  infos <- dplyr::collect(dplyr::tbl(db_con, "Variant"))
  infos$myid <- with(infos, paste(chromosome, position, allele1,
                                  allele2, sep = "_"))
  info = list()
  info[[1]] = infos$myid

  ### Set up chunk sizes
  size = dim(infos)[1]

  bytes_per_genotype <- (3 * 8) + 4       # 3 values of 8 bytes, + 4 bytes for ploidy
  bytes_per_variant = bytes_per_genotype * length(gIDs)
  memory_available = memory * 1024^2
  chunksize = floor(min(maxchunksize, memory_available/bytes_per_variant))

  reps = ceiling(size/chunksize)
  outcome = NULL

  print(paste0("Split all markers into ", reps, " chunks."))
  print(paste0("Each chunk includes less than ", chunksize, " markers."))

  ### Analyze .bgen files per chunk
  for(r in 1:reps){

    print(paste0('Reading chunk ', r, ' of ',reps))
    indices = c(((r-1)*chunksize + 1):min(r*chunksize, size))

    if(is.null(chr)) chr = as.numeric(infos$chromosome[indices[1]])
    backfile = paste0(getwd(),'/',backingdir,'/',backingfile,'_',chr,'_',r)

    snps = list(); snps[[1]] = infos$myid[indices]         # The names of the SNPs in this chunk
    bgen = snp_readBGEN(bgenfiles = bgenfile, list_snp_id = snps, backingfile = backfile, ncores = 1)

    file = paste0(backfile,'.rds')
    genotype <- snp_attach(file)
    G <- genotype$genotypes

    Geno.mtx = G[,]
    colnames(Geno.mtx) = infos$myid[indices]

    if(length(gIDs) != dim(Geno.mtx)[1]) stop("Length of gIDs not equal to number of samples in .bgen file.")
    rownames(Geno.mtx) = gIDs

    # Executing the linear regression
    outcome = SPARE(obj.null = obj.null,
                  Geno.mtx = Geno.mtx,
                  missing.cutoff = missing.cutoff,
                  min.maf = min.maf,
                  p.cutoff = p.cutoff)


    new_outcome = cbind(SNP = outcome$SNP, t(matrix(unlist(strsplit(outcome$SNP, split = '_')), nrow = 4)), outcome[,-1])
    colnames(new_outcome)[2:5] = c('Chr','BP','A1','A2')


    if(r == 1){
      data.table::fwrite(new_outcome, paste0(output.file,'.txt'), sep = "\t", append = F, row.names = F, col.names = T)
    } else {
      data.table::fwrite(new_outcome, paste0(output.file,'.txt'), sep = "\t", append = T, row.names = F, col.names = F)
    }

    ### Clean up directory by removing previous connections
    for(k in 1:r){
      prev.file = paste0(getwd(),'/Connections/tmpfile',chr,'_',k)
      unlink(prev.file)
      prev.file = paste0(getwd(),'/Connections/tmpfile',chr,'_',k)
      unlink(prev.file)
    }
  }
}

#
# SPARE for .bed files
#

SPARE.bed <- function(bedfile, gIDs,
                    obj.null,
                    output.file = NULL,
                    chr = NULL,
                    missing.cutoff = 0.05,
                    min.maf = 0.05,
                    p.cutoff = 0.001,
                    memory = 512,
                    maxchunksize = 5e4)
{


  bim.file = paste0(bedfile, ".bim")
  fam.file = paste0(bedfile, ".fam")
  bed.file = paste0(bedfile, ".bed")

  if(!file.exists(bim.file)) stop("Could not find paste0(bedfile,'.bim')")
  if(!file.exists(bed.file)) stop("Could not find paste0(bedfile,'.bed')")
  if(!file.exists(fam.file)) stop("Could not find paste0(bedfile,'.fam')")
  if(is.null(output.file)) stop('please provide name of output file')

  fam.data = read.table(fam.file, stringsAsFactors = F)
  bim.data = read.table(bim.file, stringsAsFactors = F)

  N = nrow(fam.data)
  M = nrow(bim.data)

  print(paste0("Totally ", M, " markers in plink files."))
  if(!any(obj.null$IDs %in% fam.data$V2)){
    stop("None of the subject IDs from null model were found in the .fam file")
  } else { print(paste0('In total, ', length(intersect(obj.null$IDs, fam.data$V2)),' samples found in genotype and phenotype'))}

  total.samples = intersect(obj.null$IDs, fam.data$V2)

  size = dim(bim.data)[1]

  bytes_per_genotype <- (3 * 8) + 4       # 3 values of 8 bytes, + 4 bytes for ploidy
  bytes_per_variant = bytes_per_genotype * N
  memory_available = memory * 1024^2
  chunksize = floor(min(maxchunksize, memory_available/bytes_per_variant))

  reps = ceiling(size/chunksize)
  outcome = NULL

  print(paste0("Split all markers into ", reps, " chunks."))
  print(paste0("Each chunk includes less than ", chunksize, " markers."))

  outcome = NULL

  for(r in 1:reps){

    print(paste0('Reading chunk ', r, ' of ',reps))
    indices = c(((r-1)*chunksize + 1):min(r*chunksize, size))

    Geno.mtx = seqminer::readPlinkToMatrixByIndex(bedfile, 1:N, indices)
    colnames(Geno.mtx) = bim.data$V2[indices]

    outcome = SPARE(obj.null = obj.null,
                  Geno.mtx = Geno.mtx,
                  missing.cutoff = missing.cutoff,
                  min.maf = min.maf,
                  p.cutoff = p.cutoff)
    new_outcome = cbind(Chr = bim.data[match(outcome$SNP,bim.data[,2]),1],
                        SNP = outcome$SNP,
                        A1 = bim.data[match(outcome$SNP,bim.data[,2]),5],
                        A2 = bim.data[match(outcome$SNP,bim.data[,2]),6],
                        BP = bim.data[match(outcome$SNP,bim.data[,2]),4],
                        outcome[,-1])

    if(r == 1){
      data.table::fwrite(new_outcome, paste0(output.file,'.txt'), sep = "\t", append = F, row.names = F, col.names = T)
    } else {
      data.table::fwrite(new_outcome, paste0(output.file,'.txt'), sep = "\t", append = T, row.names = F, col.names = F)
    }
  }
}

#
# CHECKING SPARE input
#

check_input_SPARE = function(obj.null, Geno.mtx, par.list)
{
  if(class(obj.null)!="NULL_Model")
    stop("obj.null should be a returned outcome from 'Null_model'")

  if(is.null(rownames(Geno.mtx))) stop("Row names of 'Geno.mtx' should be given.")
  if(is.null(colnames(Geno.mtx))) stop("Column names of 'Geno.mtx' should be given.")
  if(!is.numeric(Geno.mtx)|!is.matrix(Geno.mtx)) stop("Input 'Geno.mtx' should be a numeric matrix.")

  if(length(intersect(obj.null$IDs,rownames(Geno.mtx))) == 0) stop("None of 'IDs' are included in rownames(Geno.mtx).")
  print(paste0('In total, ', length(intersect(obj.null$IDs,rownames(Geno.mtx))),' samples with phenotype and genotype information'))

  if(!is.numeric(par.list$min.maf)|par.list$min.maf<0|par.list$min.maf>0.5) stop("Argument 'min.maf' should be a numeric value >= 0 and <= 0.5.")
  if(!is.numeric(par.list$missing.cutoff)|par.list$missing.cutoff<0|par.list$missing.cutoff>1) stop("Argument 'missing.cutoff' should be a numeric value between 0 and 1.")
}

#
# CHECKING INPUT NULL MODEL
#

check_input = function(data, IDs, mresid, range)
{
  if(is.null(IDs))
    stop("Argument 'IDs' is required in case of potential errors. For more information, please refer to 'Details'.")
  IDs = as.character(IDs)

  if(any(!is.element(IDs, data$subject)))
    stop("All elements in IDs should be also in data as 'subject'")

  if(anyDuplicated(IDs)!=0)
    stop("Argument 'IDs' should not have a duplicated element.")

  if(range[2]!=-1*range[1])
    stop("range[2] should be -1*range[1]")

  if(length(mresid)!=length(IDs))
    stop("length(mresid)!=length(IDs) where mresid are the martingale residuals.")

}

#
# FUNCTION FOR CHECKING INPUT COXME OBJECT AND DATA OBJECT
#


MGres_check = function(fitme, data){
  if(is.null(fitme)) stop('no coxme object included')
  if(class(fitme) != 'coxme') stop('object not of class coxme')
  if(is.null(data)) stop('no data object included')
  if(class(data) != 'data.frame') stop('data is not a data frame object')
  if(!'subject' %in% colnames(data)) stop('please include individuals as "subject" in dataframe')

  check_strat = strsplit(as.character(fitme$formulaList$fixed)[3],'strata')[[1]]
  if(length(check_strat) > 2){
    stop('do not include stratum/covariates with name strata')
  }

  if(length(check_strat) > 1){        # In this case, strata have been specified..
    name = substring(strsplit(check_strat[2], ')')[[1]][1],2)

    strats = data[,which(colnames(data) == name)]
    try(if(length(strats) != dim(data)[1]) stop('please include the strata once in the data frame'))
  } else {
    strats = rep(0, dim(data)[1])
  }
}

#
# Function for converting score statistics to Hazard Ratios
#

Beta_to_HR = function(beta, g, resids, cumhaz, frail){
  constant = sum(g^2 * cumhaz/(frail * cumhaz + 1))/sum(g^2)
  HR = beta/constant
  return(HR)
}

#
# Function for imputing missing values to mean
#

na_mean = function (x, option = "mean", maxgap = Inf) {
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
