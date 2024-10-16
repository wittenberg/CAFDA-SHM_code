## Variance components## load libraries
suppressPackageStartupMessages({
  library(spc)
  library(mgcv)
  library(dplyr)
  library(parallel, warn.conflicts = FALSE)
  RNGkind("L'Ecuyer-CMRG")
})


## -------------------------------------------------------------------------- ##
## Initialize setup ##

## -------------------------------------------------------------------------- ##

## functional data generation process
generate_profile_data <- function(hseq, delta, efi, meanU, fb1) {
  ################################
  # Step 4.1 new confounder data #
  ################################  
  ## 1. simulate structural component with shift "delta" in eigenfunction "efi"
  plain <- funData::simFunData(hseq, M=3, eFunType="Poly", eValType="exponential", N=1)
  ME    <- array(stats::rnorm(prod(dim(plain$simData@X)), mean=0, sd=.2), dim(plain$simData@X)) # White noise
  noisy <- plain$simData + funData::funData(plain$simData@argvals, ME) + plain$trueFuns@X[efi,]*delta
  # 2. simulate_covariate_mean_profile
  argvals     <- 2*pi*hseq/24
  covmeanprof <- sin(.3+argvals)
  # 3. simulate global intercept
  iglobal <- outer(-.05, sin(argvals*1/2)+cos(argvals*2))
  aldf    <- iglobal-mean(iglobal)
  ## 4. Simulate complete dataset
  z <- -1*runif(1, 0, 4) * covmeanprof + runif(1, 2, 10)
  w <- sapply(1:24, function(j) noisy@X[, j]) 
  data.frame(intercept=t(aldf), td01=hseq, meanU=meanU, z, fz=fb1(z), w) |>
    dplyr::mutate(u=intercept+meanU+fz+w)
}

## generation of scores
yaoetal <- function(data, phi, s2, sig, nu){
  idx    <- is.na(data)
  newphi <- phi[(idx==FALSE),]
  SIGMA  <- newphi%*%s2%*%t(newphi)+sig*diag(nrow(newphi))
  nu*c(t(newphi)%*%chol2inv(chol(SIGMA))%*%data[(idx==FALSE)])
}

## MEWMA run-length computation
MEWMA_arl <- function(mu0, h4, r, phi, hseq, delta, efi, meanU, fb1, gam2d, invSig, s2, sig, nu) {
  rl <- 0
  Z  <- NULL
  T2 <- 0
  while( T2 < h4 ){
    rl <- rl + 1
    simdf <- generate_profile_data(hseq, delta, efi, meanU, fb1)
    pred  <- predict(gam2d$gam, newdata=simdf)
    error <- simdf$u - pred
    X     <- yaoetal(error, phi, s2, sig, nu)
    if ( rl==1 ) {
      Z <- r*X + (1-r)*mu0
    } else {
      Z <- (r*X) + (1-r)*Z
    }
    T2 <- t(Z-mu0) %*% invSig %*% (Z-mu0)
  }
  rl
}

n  <- 1e4    ## Number of replications
r  <- 1     ## MEMWA smooting parameter change to .3 or 1
L0 <- 100    ## In-control ARL
p  <- 3      ## Number of eigenfunctions
h4 <- spc::mewma.crit(r, L0, p, hs=0, r=20) ## MEWMA control limit
RUNS <- 1:100
nc <- 120
## adapt this path to your local setup to load the estimated models and eigenfunctions
DIR <- "data/simulations/"
sapply(RUNS, function(ii) {

  # ## ----------------------------------------------------------------------- ###
  ## Start Cluster
  cl <- parallel::makeCluster(nc);
  
  # ## ----------------------------------------------------------------------- ###
  
  # ## ----------------------------------------------------------------------- ###
  
  fpcaResSimModel <- readRDS(paste0(DIR, "fpcaResSim_simulation_runlength_run_", ii, ".RDS"))
  gam2d   <- readRDS(paste0(DIR, "gam2d_simulation_runlength_run_", ii, ".RDS"))
  nu      <- nlme::VarCorr(gam2d$lme)[c("ef01", "ef02", "ef03"),][,"Variance"] |> ## Variance components
             t() |>
             as.data.frame() |>
             as.numeric()
  s2      <- diag(nu)                   ## Covariance matrix of the random effects
  invSig  <- diag(1/nu) / (r/(2-r))      ## Inverse of covariance matrix
  sig     <- gam2d$gam$sig2             ## Variance of the error term
  phi     <- fpcaResSimModel$efunctions ## Eigenfunctions from basic model
  hseq    <- seq(.5, 23.5, 1)
  delta_vec  <-  unique(c(seq(0, .5, .1), seq(.5, 4, .2), 4))
  meanU   <- 5                          ## Sensor mean level
  fb1     <- function(z) exp(-z/2.2)-.5 ## Transformation function f(z)
  mu0     <- c(rep(0, p))               ## Expected value in-control
  resEFI1 <- list()                     ## Results for different shifts in first eigenfunction
  resEFI2 <- list()                     ## Results for different shifts in second eigenfunction
  resEFI3 <- list()                     ## Results for different shifts in third eigenfunction


  ## Step 2:run-length computation

  # open cluster
  parallel::clusterSetRNGStream(cl, iseed=42)
  parallel::clusterEvalQ(cl, {library(mgcv)})
  parallel::clusterExport(cl, varlist=c("mu0", "h4", "r", "phi", "hseq", "delta_vec", 
                              "meanU", "fb1", "gam2d", "MEWMA_arl", "invSig", "resEFI1",
                              "resEFI2", "resEFI3", "generate_profile_data", "s2", "sig", 
                              "nu", "yaoetal"), envir=environment())
  # Shift in first eigenfunction
  resEFI1 <- lapply(seq_along(delta_vec), function(j) {
    list(
      parSapplyLB(cl, 1:n, function(i) {
        MEWMA_arl(mu0=mu0, h4=h4, r=r, phi=phi, hseq=hseq, delta=delta_vec[j], 
                  efi=1, meanU=meanU, fb1=fb1, gam2d=gam2d, invSig=invSig, s2=s2, 
                  sig=sig, nu=nu)
      })
    )
  })

  # Shift in second eigenfunction
  resEFI2 <- lapply(seq_along(delta_vec), function(j) {
    list(
      parSapplyLB(cl, 1:n, function(i) {
        MEWMA_arl(mu0=mu0, h4=h4, r=r, phi=phi, hseq=hseq, delta=delta_vec[j],
                  efi=2, meanU=meanU, fb1=fb1, gam2d=gam2d, invSig=invSig, s2=s2, 
                  sig=sig, nu=nu)
      })
    )
  })

  # Shift in third eigenfunction
  resEFI3 <- lapply(seq_along(delta_vec), function(j) {
    list(
      parSapplyLB(cl, 1:n, function(i) {
        MEWMA_arl(mu0=mu0, h4=h4, r=r, phi=phi, hseq=hseq, delta=delta_vec[j],
                  efi=3, meanU=meanU, fb1=fb1, gam2d=gam2d, invSig=invSig, s2=s2, 
                  sig=sig, nu=nu)
      })
    )
  })

  # First eigenfunction: Compute ARL and standard error
  res2  <- sapply(seq_along(delta_vec), function(i) resEFI1[[i]][[1]])
  ARL   <- apply(res2, 2, mean)
  ARLSE <- apply(res2, 2, sd)/sqrt(n)
  print(cbind(delta_vec, ARL, ARLSE))
  resA  <- data.frame(r, delta_vec, ARL, ARLSE)
  resA

  # Second eigenfunction: Compute ARL and standard error
  res2  <- sapply(seq_along(delta_vec), function(i) resEFI2[[i]][[1]])
  ARL   <- apply(res2, 2, mean)
  ARLSE <- apply(res2, 2, sd)/sqrt(n)
  print(cbind(delta_vec, ARL, ARLSE))
  resB  <- data.frame(r, delta_vec, ARL, ARLSE)
  resB

  # Third eigenfunction: Compute ARL and standard error
  res2  <- sapply(seq_along(delta_vec), function(i) resEFI3[[i]][[1]])
  ARL   <- apply(res2, 2, mean)
  ARLSE <- apply(res2, 2, sd)/sqrt(n)
  print(cbind(delta_vec, ARL, ARLSE))
  resC  <- data.frame(r, delta_vec, ARL, ARLSE)
  resC

  # Save results
  saveRDS(bind_rows(resA, resB, resC, .id="EFI"), paste0(DIR, "MEWMA_ARL_r_", r, "_N_300_RUN_", ii, ".rds"))

  ii
  # ## ----------------------------------------------------------------------- ###
  ## Stop Cluster
  parallel::stopCluster(cl)
  # ## ----------------------------------------------------------------------- ###
  
})

