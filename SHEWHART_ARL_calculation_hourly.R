## Variance components## load libraries
suppressPackageStartupMessages({
  library(spc)
  library(mgcv)
  library(dplyr)
  library(segmented)
  library(funData)
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

# Shewhart ARL computation based on simulation
SHEWHART_arl_hourly <- function(h, hseq, delta, efi, meanU, fb1, sig, modlm) {
  rl <- 0
  T2 <- 0
  while(all(T2 < h)){
    simdf <- generate_profile_data(hseq, delta, efi, meanU, fb1)
    pred  <- predict(modlm, newdata=simdf)
    error <- simdf$u - pred
    T2 <- (error/sig)^2
    if(any(T2>h)) {rl <- rl + min(which(T2>h)); break} else {rl <- rl + 24}
  }
  rl
}

n  <- 1e2    ## Number of replications
L0 <- 100*24 ## In-control ARL (rescaled to hours)
p  <- 1      ## Number of parameters
h <- spc::mewma.crit(l=1, L0, p, hs=0, r=20) ## Shewhart control limit
RUNS <- 1:100
nc <- 28#120
## adapt this path to your local setup to load the estimated models and eigenfunctions
DIR <- "data/simulations/"
sapply(RUNS, function(ii) {
  # ## ----------------------------------------------------------------------- ###
  ## Start Cluster
  cl <- parallel::makeCluster(nc);
  
  # ## ----------------------------------------------------------------------- ###
  
  # ## ----------------------------------------------------------------------- ###
  
  gam2d   <- readRDS(paste0(DIR, "gam2d_simulation_runlength_run_", ii, ".RDS"))
  dtai    <- gam2d$gam$model[c("u", "z", "ind_dayf", "td01")]
  modlm   <- segmented::segmented(lm(u~z, data=dtai)) # fit piecewise linear model
  # compute residuals, group by day, average hourly data to 24h and compute mean and sd
  fres    <- dtai |>
              mutate(res=u-predict(modlm, newdata=dtai), i=row_number()) |>
              summarise(meanres=mean(res), sdres=sd(res))
  sig     <- fres$sdres
  hseq    <- seq(.5, 23.5, 1)
  delta_vec  <- unique(c(seq(0, .5, .1), seq(.5, 4, .2), 4))
  meanU   <- 5                          ## Sensor mean level
  fb1     <- function(z) exp(-z/2.2)-.5 ## Transformation function f(z)
  resEFI1 <- list()                     ## Results for different shifts in first eigenfunction
  resEFI2 <- list()                     ## Results for different shifts in second eigenfunction
  resEFI3 <- list()                     ## Results for different shifts in third eigenfunction

  ## Step 2:run-length computation

  # open cluster
  parallel::clusterSetRNGStream(cl, iseed=42)
  parallel::clusterEvalQ(cl, {library(mgcv); library(dplyr); library(segmented); library(funData)})
  parallel::clusterExport(cl, varlist=c("h", "hseq", "delta_vec", "modlm",
                              "meanU", "fb1", "SHEWHART_arl_hourly", "resEFI1", "resEFI2", 
                              "resEFI3", "generate_profile_data", "sig"), envir=environment())
  # Shift in first eigenfunction
  resEFI1 <- lapply(seq_along(delta_vec), function(j) {
    list(
      parSapplyLB(cl, 1:n, function(i) {
        SHEWHART_arl_hourly(h=h, hseq=hseq, delta=delta_vec[j], efi=1, meanU=meanU, 
                     fb1=fb1, sig=sig, modlm=modlm)
      })
    )
  })

  # Shift in second eigenfunction
  resEFI2 <- lapply(seq_along(delta_vec), function(j) {
    list(
      parSapplyLB(cl, 1:n, function(i) {
        SHEWHART_arl_hourly(h=h,hseq=hseq, delta=delta_vec[j], efi=2, meanU=meanU, 
                     fb1=fb1, sig=sig, modlm=modlm)
      })
    )
  })

  # Shift in third eigenfunction
  resEFI3 <- lapply(seq_along(delta_vec), function(j) {
    list(
      parSapplyLB(cl, 1:n, function(i) {
        SHEWHART_arl_hourly(h=h,hseq=hseq, delta=delta_vec[j], efi=3, meanU=meanU, 
                     fb1=fb1, sig=sig, modlm=modlm)
      })
    )
  })

  # First eigenfunction: Compute ARL and standard error
  res2  <- sapply(seq_along(delta_vec), function(i) resEFI1[[i]][[1]])
  ARL   <- apply(res2, 2, mean)
  ARLSE <- apply(res2, 2, sd)/sqrt(n)
  print(cbind(delta_vec, ARL, ARLSE))
  resA  <- data.frame(delta_vec, ARL, ARLSE)
  resA

  # Second eigenfunction: Compute ARL and standard error
  res2  <- sapply(seq_along(delta_vec), function(i) resEFI2[[i]][[1]])
  ARL   <- apply(res2, 2, mean)
  ARLSE <- apply(res2, 2, sd)/sqrt(n)
  print(cbind(delta_vec, ARL, ARLSE))
  resB  <- data.frame(delta_vec, ARL, ARLSE)
  resB

  # Third eigenfunction: Compute ARL and standard error
  res2  <- sapply(seq_along(delta_vec), function(i) resEFI3[[i]][[1]])
  ARL   <- apply(res2, 2, mean)
  ARLSE <- apply(res2, 2, sd)/sqrt(n)
  print(cbind(delta_vec, ARL, ARLSE))
  resC  <- data.frame(delta_vec, ARL, ARLSE)
  resC

  # Save results
  saveRDS(bind_rows(resA, resB, resC, .id="EFI"), paste0(DIR, "SHEWHART_ARL_hourly_r_", 1, "_N_300_RUN_", ii, ".rds"))

  ii
  # ## ----------------------------------------------------------------------- ###
  ## Stop Cluster
  parallel::stopCluster(cl)
  # ## ----------------------------------------------------------------------- ###
  
})
