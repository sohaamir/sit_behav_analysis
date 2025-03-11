sit_run <- function(modelStr, test = TRUE, fitObj = NA, saveFit = FALSE, suffix=NULL,
                       nSubj = NULL, site = NULL, seed = NULL,
                       nIter = 2000, nwarmup = NULL, nCore = NULL,
                       adapt = 0.8, treedepth = 10, nThinn = 1 ) {
    
    # estimating only the choice, but not the bet / confidence
    
    library(rstan); library(loo); library(hBayesDM)
    L <- list()
    
    #### prepare data #### ===========================================================================
    load("_data/subset_data.rdata")
        
    #### preparation for running stan #### ============================================================
    # model string in a separate .stan file
    modelFile <- paste0("_scripts/stan_models/",modelStr,".stan")
    
    # setup up Stan configuration
    if (test == TRUE) {
        options(mc.cores = 1)
        nSamples <- 4
        nChains  <- 1 
        nBurnin  <- 0
        nThin    <- 1
        est_seed = sample.int(.Machine$integer.max, 1)
        
    } else {
        if (is.null(nCore)) {
            options(mc.cores = 4) 
        } else {
            options(mc.cores = nCore)  
        }
        
        if (is.null(seed)) {
            est_seed = sample.int(.Machine$integer.max, 1) 
        } else {
            est_seed = seed  
        }
        
        nSamples <- nIter
        nChains  <- 4
        if (is.null(nwarmup)) {
            nBurnin <- floor(nSamples/2)
        } else {
            nBurnin = nwarmup
        }
        nThin    <- nThinn
    }
    
    # parameter of interest (this could save both memory and space)
    poi <- create_pois(modelStr)
    
    #### run stan ####  ==============================================================================
    cat("\nDetails:\n")
    cat("Estimating", modelStr, "model for site", site, "... \n")
    cat(" # of MCMC samples (per chain) = ", nSamples, "\n")
    cat(" # of burn-in samples          = ", nBurnin, "\n")
    cat("Variant: ", modelStr, suffix, " ... \n", sep = "")
    startTime = Sys.time(); print(startTime)
    cat("Calling", nChains, "simulations in Stan... \n")
    cat("Running with", dataList$nSubjects, "participants... \n")
    
    rstan_options(auto_write = TRUE)
    
    stanfit <- stan(modelFile,
                    fit     = fitObj,
                    data    = dataList,
                    pars    = poi,
                    chains  = nChains,
                    iter    = nSamples,
                    warmup  = nBurnin,
                    thin    = nThin,
                    init    = "random",
                    seed    = est_seed,
                    control = list(adapt_delta = adapt, max_treedepth = treedepth),
                    verbose = FALSE)
    
    cat("Finishing", modelStr, "model simulation ... \n")
    endTime = Sys.time(); print(endTime)  
    cat("It took",as.character.Date(endTime - startTime), "\n")
    
    L$data <- dataList
    L$fit  <- stanfit
    
    class(L) <- "hBayesDM"
    
    if (saveFit == TRUE) {
        saveRDS(L, file = paste0('_stanfits/', '_', modelStr, suffix, '_PPC.RData'))
    }
    
    cat(' # --- rhat range:', range(rhat(L), na.rm=T), '\n')
    cat(' # --- LOOIC: \n')
    
    print(extract_looic_C(L,4,nSamples, nBurnin))
    # suppressWarnings(round(extract_looic_C(f)$looicC1C2))
    
    return(L)
}  # function run_model()


#### nested functions #### ===========================================================================

# ------------------------------------------------------------------------------------------------------
create_pois <- function(model){
    pois <- list()
    
    if (model == "RevLearn_m1a_C" || model == "RevLearn_m1b_C" ||
        model == "RevLearn_m2a_C" || model == "RevLearn_m2b_C") {
        pois <- c("lr_mu", "beta_mu", 
                  "lr_sd", "beta_sd", 
                  "lr", "beta", 
                  "c1_rep_forw", "c2_rep_forw",
                  "log_likc1", "log_likc2",
                  "lp__")
        
    } else if ( model == "RevLearn_m3a_C" || model == "RevLearn_m3b_C") {
        pois <- c("lr_my_mu", "lr_oth_mu", "beta_mu", "tau_oth_mu", 
                  "lr_my_kappa0", "lr_oth_kappa0", "beta_sd", "tau_oth_sd",
                  "lr_my", "lr_oth", "beta", "tau_oth",
                  "c1_rep_forw", "c2_rep_forw", 
                  "log_likc1", "log_likc2", 
                  "lp__")
        
    } else if ( model == "RevLearn_m4a_C" || model == "RevLearn_m4b_C" || 
                model == "RevLearn_m5a_C" || model == "RevLearn_m5b_C") {
        pois <- c("lr_mu", "beta_mu", 
                  "lr_kappa0", "beta_sd", 
                  "lr", "beta", 
                  "c1_rep_forw", "c2_rep_forw", 
                  "log_likc1", "log_likc2", 
                  "lp__")
        
    } else if ( model == "RevLearn_m6a_C" || model == "RevLearn_m6b_C" ) {
        pois <- c("lr_mu", "beta_mu", "disc_mu",
                  "lr_kappa0", "beta_sd", "disc_kappa0",
                  "lr", "beta", "disc",
                  "c1_rep_forw", "c2_rep_forw", 
                  "log_likc1", "log_likc2", 
                  "lp__")
        
    }
    return(pois)
} # function

#### end of function ####
