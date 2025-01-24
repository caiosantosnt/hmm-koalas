# Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
    library(tidyverse)
    library(momentuHMM)
    library(terra)
    library(lubridate)
})


# Load data ---------------------------------------------------------------
data <- readRDS("Rfdat_dists.rds")


# Remove unrealistic locations --------------------------------------------
library(atlastools)

# create a list
list.ind <- split(data, by = "ID")

# iterate over individuals
list.ind <- lapply(list.ind, function(dt) {
  # first process time to seconds
  dt <- dt %>%
    mutate(time = as.integer(
      as.numeric(timestamp)
    ))
  
  # calculate speed_in and speed_out
  dt <- dt %>%
    mutate(speed_in = atl_get_speed(dt,
                                    x = "x", y = "y",
                                    time = "time",
                                    type = "in"),
           speed_out = atl_get_speed(dt,
                                     x = "x", y = "y",
                                     time = "time",
                                     type = "out"))
  
  #Get 95th percentile of speed and angle
  # use sapply
  sp.thresholds <-
    sapply(dt[, list(speed_in, speed_out)],
           quantile,
           probs = 0.95, na.rm = T)
  
  # Filter by quantile
  q <- paste("speed_in <", sp.thresholds[1], "&",
             "speed_out <", sp.thresholds[2])
  
  dt <- atl_filter_covariates(
    data = dt,
    filters = c(q))
})

data.clean <- do.call(rbind, list.ind)
rm(list.ind)


# Split in shorter tracks -------------------------------------------------
dat <- data.clean |> 
  utilTM::split_ts(max_interval = 18,
                   min_ts_length = 96,
                   units = "hour") |> 
  mutate(ID_old = ID,
         ID = ID_split) |> 
  dplyr::select(-"ID_split")

# Prepare for regularisation ---------------------------------------------
# crawl initial states and parameters
theta <- fixPar <- ln.prior <- list()
for(i in unique(dat$ID)){
  theta[[i]]<-c(7, 1)
  fixPar[[i]]<-c(1,1,NA,NA)
  ln.prior[[i]] <- function(theta){-abs(theta[2]+3)/0.5} # laplace prior with location = -3 and scale = 0.5 
}

# Crawl - Regularisation -------------------------------------------------
crwOut <- crawlWrap(dat,
                    retryFits = 10,
                    attempts = 100,
                    ncores = 14,
                    retryParallel = T,
                    timeStep = "6 hours",
                    method = "Nelder-Mead",
                    theta = theta,
                    prior = ln.prior,
                    initialSANN = list(
                      maxit = 10000,
                      trace = 1))


# Include metrics ---------------------------------------------------------
## Calculate step-lengths and turning angles ----
fdat <- prepData(crwOut)


# Prepare to fit ----------------------------------------------------------
# determine knownStates based on map visualization
knownStates<-rep(NA,nrow(fdat))

# Specify state 3
knownStates[c(3016:3017,
              501:506, 570, 577, 607,
              7753, 7761, 7937,
              6877,
              5755:5757, 5764:5766, 5684:5686, 5723, 5724, 5678, 5679, 5684,
              10278:10279, 10286:10288, 10293:10295, 10321, 10322)] <- 3

# Specify state 2
knownStates[c(3013:3015,
              498, 507, 508, 610,
              7762, 7862:7864,
              6876, 6854,
              5773, 5777, 5778, 5781, 5782, 5719:5722, 5711, 5712,
              10293:10296, 10315, 10316)] <- 2


# Specify state 1
knownStates[c(3054, 3055, 3125,
              493:496, 499, 500, 601:605,
              7929:7932, 7838:7844,
              6860, 6861,
              5775, 5776, 5659:5661, 5672, 5673,
              10297:10303, 10308:10314)] <- 1

# define number of states
nbStates <- 3

# define the distribution of step lengths and turning angles
dist <- list(step = "gamma", angle = "wrpcauchy")

# define initial parameters
stepPar0.3s <- c(10, 30, 90, 5, 15, 45)
anglePar0.3s <- c(pi, pi, 0, 0.1, 0.01, 0.7)
beta0.3s <- matrix(c(-2, -100, -2, -2, -2, -2), nrow = 1)
beta.fixpar <- matrix(c(NA, -100, NA, NA, NA, NA), nrow = 1)

# fit model without covariates
mod.3s <- fitHMM(data = fdat, nbStates = 3, dist = dist, 
                 estAngleMean = list(angle = TRUE),
                 Par0 = list(step = stepPar0.3s,
                             angle = anglePar0.3s), 
                 knownStates = knownStates,
                 fixPar = list(beta = beta.fixpar),
                 nlmPar = list(print.level = 2))

#############
## Model 1 ##
#############
f1 <- ~ cosinor(hr, period = 24)
Par1 <- getPar0(mod.3s, formula = f1)
Par1$beta[1,2] <- -100
Par1$beta[-1,2] <- 0
beta.fixpar <- Par1$beta
beta.fixpar[,-2] <- NA

mod1 <- fitHMM(data = fdat, nbStates = 3, dist = dist, 
               estAngleMean = list(angle = TRUE),
               formula = f1,
               Par0 = list(step = Par1$Par$step,
                           angle = Par1$Par$angle),
               beta0 = Par1$beta,
               knownStates = knownStates, 
               fixPar = list(beta = beta.fixpar),
               nlmPar = list(print.level = 2))

#############
## Model 2 ##
#############
f2 <- ~ cosinor(hr, period = 24) + cosinor(yday, period = 365.25)
Par2 <- getPar0(mod1, formula = f2)
Par2$beta[1,2] <- -100
Par2$beta[-1,2] <- 0
beta.fixpar <- Par2$beta
beta.fixpar[,-2] <- NA

mod2 <- fitHMM(data = fdat, nbStates = 3, dist = dist, 
               estAngleMean = list(angle = TRUE),
               formula = f2,
               Par0 = list(step = Par2$Par$step,
                           angle = Par2$Par$angle),
               beta0 = Par2$beta,
               knownStates = knownStates, 
               fixPar = list(beta = beta.fixpar),
               nlmPar = list(print.level = 2))

#############
## Model 3 ##
#############
f3 <- ~ cosinor(hr, period = 24) + cosinor(yday, period = 365.25) + age
Par3 <- getPar0(mod2, formula = f3)
Par3$beta[1,2] <- -100
Par3$beta[-1,2] <- 0
beta.fixpar <- Par3$beta
beta.fixpar[,-2] <- NA

mod3 <- fitHMM(data = fdat, nbStates = 3, dist = dist, 
               estAngleMean = list(angle = TRUE),
               formula = f3,
               Par0 = list(step = Par3$Par$step,
                           angle = Par3$Par$angle),
               beta0 = Par3$beta,
               knownStates = knownStates, 
               fixPar = list(beta = beta.fixpar),
               nlmPar = list(print.level = 2))

#############
## Model 4 ##
#############
f4 <- ~ cosinor(hr, period = 24) + cosinor(yday, period = 365.25) + age + sex
Par4 <- getPar0(mod3, formula = f4)
Par4$beta[1,2] <- -100
Par4$beta[-1,2] <- 0
beta.fixpar <- Par4$beta
beta.fixpar[,-2] <- NA

mod4 <- fitHMM(data = fdat, nbStates = 3, dist = dist, 
               estAngleMean = list(angle = TRUE),
               formula = f4,
               Par0 = list(step = Par4$Par$step,
                           angle = Par4$Par$angle),
               beta0 = Par4$beta,
               knownStates = knownStates, 
               fixPar = list(beta = beta.fixpar),
               nlmPar = list(print.level = 2))

#############
## Model 5 ##
#############
f5 <- ~ cosinor(hr, period = 24) + cosinor(yday, period = 365.25) + age + sex + origin
Par5 <- getPar0(mod4, formula = f5)
Par5$beta[1,2] <- -100
Par5$beta[-1,2] <- 0
beta.fixpar <- Par5$beta
beta.fixpar[,-2] <- NA

mod5 <- fitHMM(data = fdat, nbStates = 3, dist = dist, 
               estAngleMean = list(angle = TRUE),
               formula = f5,
               Par0 = list(step = Par5$Par$step,
                           angle = Par5$Par$angle),
               beta0 = Par5$beta,
               knownStates = knownStates, 
               fixPar = list(beta = beta.fixpar),
               nlmPar = list(print.level = 2))

#############
## Model 6 ##
#############
f6 <- ~ cosinor(hr, period = 24) + cosinor(yday, period = 365.25) + age + sex + origin + urb
Par6 <- getPar0(mod5, formula = f6)
Par6$beta[1,2] <- -100
Par6$beta[-1,2] <- 0
beta.fixpar <- Par6$beta
beta.fixpar[,-2] <- NA

mod6 <- fitHMM(data = fdat, nbStates = 3, dist = dist, 
               estAngleMean = list(angle = TRUE),
               formula = f6,
               Par0 = list(step = Par6$Par$step,
                           angle = Par6$Par$angle),
               beta0 = Par6$beta,
               knownStates = knownStates, 
               fixPar = list(beta = beta.fixpar),
               nlmPar = list(print.level = 2))


#############
## Model 7 ##
#############
f7 <- ~ cosinor(hr, period = 24) + cosinor(yday, period = 365.25) + age + sex + origin + urb + rur
Par7 <- getPar0(mod5, formula = f6)
Par7$beta[1,2] <- -100
Par7$beta[-1,2] <- 0
beta.fixpar <- Par7$beta
beta.fixpar[,-2] <- NA

mod7 <- fitHMM(data = fdat, nbStates = 3, dist = dist, 
               estAngleMean = list(angle = TRUE),
               formula = f6,
               Par0 = list(step = Par7$Par$step,
                           angle = Par7$Par$angle),
               beta0 = Par7$beta,
               knownStates = knownStates, 
               fixPar = list(beta = beta.fixpar),
               nlmPar = list(print.level = 2))


################################################################################
########################   END   ###############################################
################################################################################