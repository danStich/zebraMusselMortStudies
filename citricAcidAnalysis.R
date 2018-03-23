# Front-end needs ---------------------------------------------------------
# Package install and load
  #install.packages("R2jags")
  library(R2jags)
  library(plyr)
  library(lubridate)

# Make a fuunction ot invert the logit transformation
  inv.logit = function(x){
    exp(x)/(1+exp(x))
  }

# Data --------------------------------------------------------------------
# Read it in
  vd1 = read.csv('caRound1.csv', stringsAsFactors = FALSE) # Round 1 of CA experiment  
  vd2 =read.csv('caRound2.csv', stringsAsFactors = FALSE)   # Round 2 of CA experiment

# Row bind data from the two rounds
  vd = rbind(vd1, vd2)

# Need to make a copy of the data frame with the zero time point of
# 100% survival to tighten up the intercepts.
  vd2 = vd
  vd2$time.pt=0
  #vd2$mort=0
  vd2$dead=0
  
# Add the zero time point to the data
  vd = rbind(vd, vd2)
  vd$start = 1
  
# Summarize the data by concentration and dose
  vd = ddply(vd, c('tank', 'time.pt', 'dose'), summarize, 
             dead=sum(dead), N=sum(start))
  
# Adult Model specification -----------------------------------------------
# Write the binomial model to a file in the working directory
# Base model
  modelString = "
    model{

    # Binomial likelihood
      for(i in 1:nruns){
        C[i] ~ dbinom(p[i], N[i])
        logit(p[i]) <- alpha + beta.time[dose[i]]*time.pt[i]
      }

    # Intercept term
      alpha ~ dunif(-10, -1)

    # Random intercepts and slopes for effects of dose and time on mortality
      for(t in 1:ndose){
        beta.time[t] ~ dnorm(mu.betaTime, tau.betaTime)
      }

    # Hyperparameters for random slopes for effect of time by dose
      mu.betaTime ~ dnorm(0, 0.001)
      tau.betaTime <- 1/(sigma.betaTime * sigma.betaTime)
      sigma.betaTime ~ dunif(0, 10)

    }"

  # Write the model to a file
    writeLines(modelString, con='vd_baseModel.txt')

# Adult model calibration --------------------------------------------------
# Package the data so JAGS can use it
  # Data for the base model
    vd.data = list(
      C = vd$dead,
      nruns = nrow(vd),
      ndose = length(unique(vd$dose)),
      dose = as.numeric(as.factor(vd$dose)),
      N = vd$N,
      time.pt = vd$time.pt#as.vector(scale(vd$time.pt))
    )

# Supply initial values for Markov Chains in the Gibbs sampler
  # Initial values for the base model
  inits = function(){
    list(
      alpha = runif(1, -10, 0),
      mu.betaTime = rnorm(1, 0, 1),
      sigma.betaTime = runif(1,0,10)
    )
  }

# Tell JAGS which parameters you want to trace during posterior sampling
  # Base model
    params = c(
      "alpha",
      "beta.time",
      "mu.betaTime",
      "sigma.betaTime"
    )

# MCMC settings
  ni = 33000
  nb = 3000
  nt = 3
  nc = 3

# Run the model
  # Base model
    vd.model = jags(
      data = vd.data,
      inits = inits,
      parameters.to.save = params,
      model.file = "vd_baseModel.txt",
      n.chains = nc,
      n.iter = ni,
      n.burnin = nb,
      n.thin = nt,
      working.directory = getwd(),
      DIC = TRUE
    )
    
# Print the model results
  print(vd.model, digits=3)

# Results -----------------------------------------------------------------
# Store the posterior estimates to an R object and write them out to a file for
# safe keeping
  # Extract results from model object and make into a dataframe
    vd.res = data.frame(vd.model$BUGSoutput$sims.list)

# Individual dosage predictions with uncertainty
  new.time=seq(0,72,1)
  # Predict mortality based on each of the coefficients from the posterior
  # Plot the posterior predictions
    par(mar=c(5,5,1,1))

    # First, make a blank plot to which we will add the others
      plot(0:0, xlim=c(0,72), ylim=c(0,1), xlab='Exposure time (minutes)',
           ylab = "", cex.lab=1.25, cex.axis=1.15, yaxt = 'n', type='l')

    # Now add the predictions for each draw from the posterior. Change the
    # number that follows 'beta.time.' to access the various posteriors for
    # doses 1-6
      for(i in 1:nrow(vd.res)){
        lines(new.time,
              inv.logit(vd.res$alpha[i]+vd.res$beta.time.7[i]*new.time),
              lwd=2,
              lty=1,
              col='gray87'
        )
      }
    # Add the mean of the posterior predictions, making sure to use the same
    # beta for the dose above
      lines(new.time,
            inv.logit(mean(vd.res$alpha)+mean(vd.res$beta.time.7)*new.time),
            col = 'red', lty = 1, lwd = 2
            )
    # Add lines for the 95% credible interval
      lines(new.time,
            inv.logit(
              quantile(vd.res$alpha,0.025)+quantile(vd.res$beta.time.7,0.025)*new.time),
              col = 'red', lty = 2, lwd = 1
            )
      lines(new.time,
            inv.logit(
              quantile(vd.res$alpha,0.975)+quantile(vd.res$beta.time.7,0.975)*new.time),
              col = 'red', lty = 2, lwd = 1
            )
    # Add the y-axis
      axis(2, las=2)
      mtext("Probability of mortality", 2, line=3, cex=1.25)




