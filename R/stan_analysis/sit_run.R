library(rstan)
library(here)
library(conflicted)

# Stan-specific functions
conflicts_prefer(rstan::extract)
conflicts_prefer(rstan::loo)

# Matrix operations (important for Stan)
conflicts_prefer(Matrix::expand)
conflicts_prefer(Matrix::pack)
conflicts_prefer(Matrix::unpack)

# Stats and modeling functions
conflicts_prefer(stats::lag)
conflicts_prefer(stats::step)

# Data manipulation
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::combine)
conflicts_prefer(dplyr::recode)

# Plotting
conflicts_prefer(ggplot2::`%+%`)
conflicts_prefer(ggplot2::alpha)

# Mixed models
conflicts_prefer(lmerTest::lmer)

# Other functions
conflicts_prefer(car::logit)
conflicts_prefer(purrr::some)
conflicts_prefer(tidyr::smiths)
conflicts_prefer(psych::phi)
conflicts_prefer(psych::lookup)
conflicts_prefer(effectsize::cohens_d)
conflicts_prefer(effectsize::eta_squared)

# =============================================================================
#### Running Stan #### 
# =============================================================================
rstan_options(auto_write = TRUE)
options(mc.cores = 2)

rm(list=ls())

# Load data
load(here("data", "rdata", "testdata.RData"))

modelFile <- here("stan", "sit_stan_models", "new", "model1a.stan")


# Add Tsubj - assuming all subjects completed all trials
testdataList$Tsubj <- rep(testdataList$nTrials, testdataList$nSubjects)

# Stan sampling parameters - reduced for testing
nIter     <- 100  # Reduced iterations for testing
nChains   <- 2    # Reduced chains for testing
nWarmup   <- floor(nIter/2)
nThin     <- 1

cat("Estimating", modelFile, "model... \n")
startTime = Sys.time(); print(startTime)
cat("Calling", nChains, "simulations in Stan... \n")

fit_rl <- stan(modelFile, 
               data    = testdataList, 
               chains  = nChains,
               iter    = nIter,
               warmup  = nWarmup,
               thin    = nThin,
               init    = "random",
               seed    = 1450154626)

cat("Finishing", modelFile, "model simulation ... \n")
endTime = Sys.time(); print(endTime)  
cat("It took", as.character.Date(endTime - startTime), "\n")

# Basic parameter check
print(fit_rl, pars = c('lr_mu', 'beta_mu[1]', 'thres_diff_mu'))


# =============================================================================
#### Plots #### 
# =============================================================================

# Density plots for individual and group parameters
# Group level parameters
plot_dens_group <- stan_plot(fit_rl, 
                             pars=c('lr_mu', 'beta_mu', 'thres_diff_mu'), 
                             show_density=TRUE, 
                             fill_color='skyblue')

# Individual level parameters
plot_dens_indiv <- stan_plot(fit_rl, 
                             pars=c('lr', 'beta', 'thres_diff'), 
                             show_density=TRUE, 
                             fill_color='skyblue')

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# Extract posterior means using the 'mean-all chains' column (column 3)
lr_values <- get_posterior_mean(fit_rl, pars='lr')[,3]
beta_values <- get_posterior_mean(fit_rl, pars='beta')[,3]
thres_values <- get_posterior_mean(fit_rl, pars='thres_diff')[,3]

# Combine all values
pars_value <- c(lr_values, beta_values, thres_values)

# Create parameter names
n_subjects <- length(lr_values)
pars_name <- c(
  rep('lr', n_subjects),
  rep(paste0('beta', 1:7), each=n_subjects),
  rep('thres', n_subjects)
)

# Create dataframe for plotting
df <- data.frame(pars_value=pars_value, pars_name=factor(pars_name))

# Define plotting theme
myconfig <- theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Create violin plot
g1 <- ggplot(df, aes(x=pars_name, y=pars_value, color=pars_name, fill=pars_name)) +
  geom_violin(trim=TRUE, size=1) +
  stat_summary(fun.data=function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }, geom="pointrange", color="black", size=0.5) +
  scale_fill_brewer(palette="Set2") +
  scale_color_brewer(palette="Set2") +
  myconfig +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = 'Parameter Value')

print(g1)
print(plot_dens_group)

# Create directory if it doesn't exist
dir.create(here("output", "stan", "plots"), recursive = TRUE, showWarnings = FALSE)

# Save plots to the specified directory
ggsave(here("output", "stan", "plots", "parameter_distributions_violin.png"), 
       g1, width = 12, height = 6)

ggsave(here("output", "stan", "plots", "group_parameters_density.png"), 
       plot_dens_group, width = 8, height = 6)

ggsave(here("output", "stan", "plots", "individual_parameters_density.png"), 
       plot_dens_indiv, width = 8, height = 6)

# =============================================================================


