# Front-end needs ---------------------------------------------------------
# Package install and load
  #install.packages("R2jags")
  library(R2jags)

# Set working directory
  setwd("P:/Research/2016_Interns/JoePerry/veliger_Bayes")

# Read in the data
  vd = read.csv('veligerData.csv')
  #vd = read.csv(file.choose())

# Change the names to match earlier data simulation used in code development
  names(vd)[c(4,5,16,17)] = c('dose', 'time.pt','mort','dead')

# Round off the number of dead individuals
  vd$dead = round(vd$dead)

# Remove lines with NA
  vd = na.omit(vd)

# Need to make a copy of the data frame with the zero time point of
# 100% survival to tighten up the intercepts.
  vd2 = vd
  vd2$time.pt=0
  vd2$mort=0
  vd2$dead=0

# Add the zero time point to the data
  vd = rbind(vd, vd2)

# Data simulation ---------------------------------------------------------
# NOTE: THIS CODE BLOCK WAS USED TO SIMULATE DATA DURING AIRPLANE TRAVEL WHEN
#       THE ORIGINAL DATA WERE INACCESSIBLE. IT IS FOR TESTING PURPOSES ONLY.

# # Simulate virkon concentrations
#   dose = rep(sort( rep( c( 0, .25, .5, 1, 2, 3 ), 10 ) ), 5)
#
# # Simulate time points corresponding to virkon concentrations
#   time.pt = rep( rep(seq(1, 10), length(unique(dose))), 5)
#
# # Simulate a relationship between dose, time, and mortality
#   mortality = c(
#     runif(10, 0.00, .15),
#     c(0, .10, .25, .50, .75, .80, .85, .90, .90, .90),
#     c(0, .25, .50, .75, .80, .85, .90, .90, .95, 1.00),
#     c(0, .50, .75, .85, .90, .95, 1.00, 1.00, 1.00, 1.00),
#     c(0, .50, .85, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00),
#     c(0, .75, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00)
#   )
#
# # Add some noise to the mortality rates
#   stoch.mort = round(rbeta(length(mortality), 500, 5), 2)
#   stoch.mort1 = round(rbeta(length(mortality), 500, 10), 2)
#
# # Now apply the stochastic parameters to mortality rates and get mortality for
# # all doses and time points.
#   mort = round( c(
#     mortality,
#     mortality * stoch.mort,
#     mortality * stoch.mort1,
#     mortality,
#     mortality * stoch.mort
#   ), 2)
#
# # Randomly draw starting number of veligers for each dose and time-point trial
#   set.seed(626)
#   N = round( runif( length( mort ), 50, 200 ) )
#
# # Calculate number dead at the end of each trial by multiplying starting number
# # by one minus the mortality rate for each combination
#   dead = round( (1-mort) * N )
#
# # Put all of the data together in a dataframe
#   vd = data.frame(dose, time.pt, mort, N, dead)
#
# # Make dose into a factor
#   vd$dose = as.numeric(as.factor(vd$dose))

# Model specification -----------------------------------------------------
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
      alpha ~ dunif(-10, 0)

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

# Covariate model
  modelString = "
    model{

    # Binomial likelihood
      for(i in 1:nruns){
        C[i] ~ dbinom(p[i], N[i])
        logit(p[i]) <- alpha + beta.time[dose[i]]*time.pt[i] +
                        beta.temp*temp[i] + beta.cond*cond[i] +
                        beta.pH*pH[i] + beta.O2*O2[i]
      }

    # Intercept term
      alpha ~ dunif(-10, 0)

    # Fixed effect of water quality parameters on mortality
      beta.temp ~ dnorm(0, 0.001)
      beta.cond ~ dnorm(0, 0.001)
      beta.pH ~ dnorm(0, 0.001)
      beta.O2 ~ dnorm(0, 0.001)

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
    writeLines(modelString, con='vd_covModel.txt')

# Model calibration --------------------------------------------------------
# Package the data so JAGS can use it
  # Data for the base model
    vd.data = list(
      C = vd$dead,
      nruns = nrow(vd),
      ndose = length(unique(vd$dose)),
      dose = as.numeric(as.factor(vd$dose)),
      N = vd$N,
      time.pt = vd$time.pt
    )

  # Data for the covariate model
    vd.Covdata = list(
      C = vd$dead,
      nruns = nrow(vd),
      ndose = length(unique(vd$dose)),
      dose = as.numeric(as.factor(vd$dose)),
      N = vd$N,
      time.pt = vd$time.pt,
      temp = vd$lakeTemp,
      pH = vd$lakePH,
      O2 = vd$lake02,
      cond = vd$lakeCond
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

  # Inits for the covariate model
    initsCov = function(){
      list(
        alpha = runif(1, -10, 0),
        mu.betaTime = rnorm(1, 0, 1),
        sigma.betaTime = runif(1,0,10),
        beta.temp = rnorm(1, 0, 1),
        beta.cond = rnorm(1, 0, 1),
        beta.pH = rnorm(1, 0, 1),
        beta.O2 = rnorm(1, 0, 1)
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

  # Covariate model
    paramsCov = c(
      "alpha",
      "beta.time",
      "mu.betaTime",
      "sigma.betaTime",
      "beta.temp",
      "beta.cond",
      "beta.pH",
      "beta.O2"
    )

# MCMC settings
  ni = 75000
  nb = 2500
  nt = 30
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

  # # Covariate model
  #   vd.model = jags(
  #     data = vd.Covdata,
  #     inits = initsCov,
  #     parameters.to.save = paramsCov,
  #     model.file = "vd_covModel.txt",
  #     n.chains = nc,
  #     n.iter = ni,
  #     n.burnin = nb,
  #     n.thin = nt,
  #     working.directory = getwd(),
  #     DIC = TRUE
  #   )

# Results -----------------------------------------------------------------
# Print the summary for coefficient estimates
  print(vd.model)

# Store the posterior estimates to an R object and write them out to a file for
# safe keeping
  # Extract results from model object and make into a dataframe
    vd.res = data.frame(vd.model$BUGSoutput$sims.list)
  # Write the data to a comma-separated file
    write.table(vd.res, 'vd_posteriors.txt', sep=',', row.names=FALSE,
                quote = FALSE)

# Make some histograms to examine coefficient estimates
# None of the coefficients overlaps zero, which means they are significant
# Because the values are all positive, this means that p(mortality) increases
# over time for all doses (excluding control, which decreases)

  # Histograms for each of the six contrations in order from 0 to 3.0 %
    hist(vd.res$beta.time.1, main='', col='gray87')
    hist(vd.res$beta.time.2, main='', col='gray87')
    hist(vd.res$beta.time.3, main='', col='gray87')
    hist(vd.res$beta.time.4, main='', col='gray87')
    hist(vd.res$beta.time.5, main='', col='gray87')
    hist(vd.res$beta.time.6, main='', col='gray87')

  # Histogram of the average response. There is considerably more overlap with
  # zero here. This is a good thing, since we are hoping we've sampled a wide
  # parameter space, and we know that there is more going on than just 'the
  # average'. But, you can still see that the general trend is still positive.
  # The reason that the HDI overlaps zero is because of the estimate for the
  # control, which is negative! That's why we want separate slopes.
    hist(vd.res$mu.betaTime, main='', col='gray87', xlab = "Grand mean of slope")
    abline( v=quantile(vd.res$mu.betaTime, probs=c(0.025, 0.975)),  col='red',
            lwd=2)

# Make predictions for each of the dose-response curves
  # Have another look at all of the doses
    sort( unique(vd$dose) )

  # Make a fuunction ot invert the logit transformation
    inv.logit = function(x){
      exp(x)/(1+exp(x))
    }

  # Make predictions from the posterior estimates to visualize mortality by dose
  # and time. First, we will look at the mean values for all of the combos
    # Create an empty list to hold the model predictions
      mean.preds = vector(mode = "list", length=length(unique(vd$dose)))

    # Create a vector to hold new values for time
      new.time = seq(0, 60, 1)

    # Make the predictions for dose response using a loop
      for(i in 2:7){
        mean.preds[[i-1]] = mean(vd.res[,1]) + mean(vd.res[ , i])*new.time
      }

    # Plot the posterior predictions
      par(mar=c(5,5,1,1))
      # First, make a blank plot to which we will add the others
        plot(0:0, xlim=c(0,60), ylim=c(0,1), xlab='Exposure time (minutes)',
             ylab = "", cex.lab=1.25, cex.axis=1.15, yaxt = 'n', type='l')

      # Add raw data points
        cols = c("red", "orange", "yellow", "green", "blue", "purple")
        points(vd$time.pt, vd$mort, pch = c(1, 2, 3, 4, 5, 6)[c(as.factor(vd$dose))],
               col=rev(cols)[c(as.factor(vd$dose))])

      # Now add the actual mean predictions
        for(i in 1:length(unique(vd$dose))){
          lines(new.time, inv.logit(mean.preds[[i]]),
                lwd=2,
                lty=seq(1,6,1)[i],
                col=rev(cols)[i]
          )
        }
      # Add the y-axis
        axis(2, las=2)
        mtext("Probability of mortality", 2, line=3, cex=1.25)
        legend(x=47, y=.8,
               legend = c('0.00', '0.25', '0.50', '1.00', '2.00', '3.00'),
               lty = seq(1,6),
               lwd=2,
               col=rev(cols)[1:6],
               bty='n',
               title = 'Concentration (%)'#,
               #title.adj = -2
               )


  # Individual dosage predictions with uncertainty
    # Predict mortality based on each of the coefficients from the posterior
    # Plot the posterior predictions
      par(mar=c(5,5,1,1))

      # First, make a blank plot to which we will add the others
        plot(0:0, xlim=c(0,60), ylim=c(0,1), xlab='Exposure time (minutes)',
             ylab = "", cex.lab=1.25, cex.axis=1.15, yaxt = 'n', type='l')

      # Now add the predictions for each draw from the posterior. Change the
      # number that follows 'beta.time.' to access the various posteriors for
      # doses 1-6
        for(i in 1:nrow(vd.res)){
          lines(new.time,
                inv.logit(vd.res$alpha[i]+vd.res$beta.time.5[i]*new.time),
                lwd=2,
                lty=1,
                col='gray87'
          )
        }
      # Add the mean of the posterior predictions, making sure to use the same
      # beta for the dose above
        lines(new.time,
              inv.logit(mean(vd.res$alpha)+mean(vd.res$beta.time.5)*new.time),
              col = 'red', lty = 1, lwd = 2
              )
      # Add lines for the 95% credible interval
        lines(new.time,
              inv.logit(
                quantile(vd.res$alpha,0.025)+quantile(vd.res$beta.time.5,0.025)*new.time),
                col = 'red', lty = 2, lwd = 1
              )
        lines(new.time,
              inv.logit(
                quantile(vd.res$alpha,0.975)+quantile(vd.res$beta.time.5,0.975)*new.time),
                col = 'red', lty = 2, lwd = 1
              )
      # Add the y-axis
        axis(2, las=2)
        mtext("Probability of mortality", 2, line=3, cex=1.25)

# Adult data --------------------------------------------------------------
# Read it in
  vd =read.csv('adults.csv')

# Need to make a copy of the data frame with the zero time point of
# 100% survival to tighten up the intercepts.
  vd2 = vd
  vd2$time.pt=0
  vd2$mort=0
  vd2$dead=0

# Add the zero time point to the data
  vd = rbind(vd, vd2)

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
      time.pt = vd$time.pt
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
  nt = 30
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

# Results -----------------------------------------------------------------
# Print the summary for coefficient estimates
  print(vd.model)

# Store the posterior estimates to an R object and write them out to a file for
# safe keeping
  # Extract results from model object and make into a dataframe
    vd.res = data.frame(vd.model$BUGSoutput$sims.list)
  # Write the data to a comma-separated file
    write.table(vd.res, 'adult_posteriors.txt', sep=',', row.names=FALSE,
                quote = FALSE)

# Make some histograms to examine coefficient estimates
# None of the coefficients overlaps zero, which means they are significant
# Because the values are all positive, this means that p(mortality) increases
# over time for all doses (excluding control, which decreases)

  # Histograms for each of the six contrations in order from 0 to 3.0 %
    hist(vd.res$beta.time.1, main='', col='gray87')
    hist(vd.res$beta.time.2, main='', col='gray87')
    hist(vd.res$beta.time.3, main='', col='gray87')
    hist(vd.res$beta.time.4, main='', col='gray87')
    hist(vd.res$beta.time.5, main='', col='gray87')
    hist(vd.res$beta.time.6, main='', col='gray87')

  # Histogram of the average response. There is considerably more overlap with
  # zero here. This is a good thing, since we are hoping we've sampled a wide
  # parameter space, and we know that there is more going on than just 'the
  # average'. But, you can still see that the general trend is still positive.
  # The reason that the HDI overlaps zero is because of the estimate for the
  # control, which is negative! That's why we want separate slopes.
    hist(vd.res$mu.betaTime, main='', col='gray87', xlab = "Grand mean of slope")
    abline( v=quantile(vd.res$mu.betaTime, probs=c(0.025, 0.975)),  col='red',
            lwd=2)

# Make predictions for each of the dose-response curves
  # Have another look at all of the doses
    sort( unique(vd$dose) )

  # Make a fuunction ot invert the logit transformation
    inv.logit = function(x){
      exp(x)/(1+exp(x))
    }

  # Make predictions from the posterior estimates to visualize mortality by dose
  # and time. First, we will look at the mean values for all of the combos
    # Create an empty list to hold the model predictions
      mean.preds = vector(mode = "list", length=length(unique(vd$dose)))

    # Create a vector to hold new values for time
      new.time = seq(0, 72, 1)

    # Make the predictions for dose response using a loop
      for(i in 2:6){
        mean.preds[[i-1]] = mean(vd.res[ , 1]) + mean(vd.res[ , i])*new.time
      }

    # Plot the posterior predictions
      par(mar=c(5,5,1,1))
      # First, make a blank plot to which we will add the others
        plot(0:0, xlim=c(0,72), ylim=c(0,1), xlab='Exposure time (h)',
             ylab = "", cex.lab=1.25, cex.axis=1.15, yaxt = 'n', type='l')

      # Add raw data points
        cols = c("red", "orange", "green", "blue", "purple")
        points(vd$time.pt, vd$mort, pch = c(1, 2, 3, 4, 5)[c(as.factor(vd$dose))],
               col=rev(cols)[c(as.factor(vd$dose))])

      # Now add the actual mean predictions
        for(i in 1:length(unique(vd$dose))){
          lines(new.time, inv.logit(mean.preds[[i]]),
                lwd=2,
                lty=seq(1,5,1)[i],
                col=rev(cols)[i]
          )
        }
      # Add the y-axis
        axis(2, las=2)
        mtext("Probability of mortality", 2, line=3, cex=1.25)
        legend(x=47, y=.8,
               legend = c('0.00', '0.10', '0.25', '0.50', '1.00'),
               lty = seq(1,5),
               lwd=2,
               col=rev(cols)[1:5],
               bty='n',
               title = 'Concentration (%)'#,
               #title.adj = -2
               )


  # Individual dosage predictions with uncertainty
    # Predict mortality based on each of the coefficients from the posterior
    # Plot the posterior predictions
      par(mar=c(5,5,1,1))

      # First, make a blank plot to which we will add the others
        plot(0:0, xlim=c(0,60), ylim=c(0,1), xlab='Exposure time (minutes)',
             ylab = "", cex.lab=1.25, cex.axis=1.15, yaxt = 'n', type='l')

      # Now add the predictions for each draw from the posterior. Change the
      # number that follows 'beta.time.' to access the various posteriors for
      # doses 1-6
        for(i in 1:nrow(vd.res)){
          lines(new.time,
                inv.logit(vd.res$alpha[i]+vd.res$beta.time.2[i]*new.time),
                lwd=2,
                lty=1,
                col='gray87'
          )
        }
      # Add the mean of the posterior predictions, making sure to use the same
      # beta for the dose above
        lines(new.time,
              inv.logit(mean(vd.res$alpha)+mean(vd.res$beta.time.2)*new.time),
              col = 'red', lty = 1, lwd = 2
              )
      # Add lines for the 95% credible interval
        lines(new.time,
              inv.logit(
                quantile(vd.res$alpha,0.025)+quantile(vd.res$beta.time.2,0.025)*new.time),
                col = 'red', lty = 2, lwd = 1
              )
        lines(new.time,
              inv.logit(
                quantile(vd.res$alpha,0.975)+quantile(vd.res$beta.time.2,0.975)*new.time),
                col = 'red', lty = 2, lwd = 1
              )
      # Add the y-axis
        axis(2, las=2)
        mtext("Probability of mortality", 2, line=3, cex=1.25)




