# Front-end needs ---------------------------------------------------------
# Package install and load
  #install.packages("R2jags")
  library(R2jags)
  library(akima)

# Read in the data
  vd = read.csv("eliseVeligers.csv")
  #vd = read.csv("eliseAdults.csv")
  
# Need to make a copy of the data frame with the zero time point of
# 100% survival to tighten up the intercepts.
  vd2 = vd
  vd2$time.pt=0
  vd2$mort=0
  vd2$dead=0

# Add the zero time point to the data
  vd = rbind(vd, vd2)

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
  ni = 7500
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
    hist(vd.res$beta.time.7, main='', col='gray87')    

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

  # # Make predictions from the posterior estimates to visualize mortality by dose
  # # and time. First, we will look at the mean values for all of the combos
  #   # Create an empty list to hold the model predictions
  #     mean.preds = vector(mode = "list", length=length(unique(vd$dose)))
  # 
  #   # Create a vector to hold new values for time
  #     new.time = seq(0, 60, 1)
  # 
  #   # Make the predictions for dose response using a loop
  #     for(i in 2:7){
  #       mean.preds[[i-1]] = mean(vd.res[,1]) + mean(vd.res[ , i])*new.time
  #     }
  # 
  #   # Plot the posterior predictions
  #     par(mar=c(5,5,1,1))
  #     # First, make a blank plot to which we will add the others
  #       plot(0:0, xlim=c(0,60), ylim=c(0,1), xlab='Exposure time (minutes)',
  #            ylab = "", cex.lab=1.25, cex.axis=1.15, yaxt = 'n', type='l')
  # 
  #     # Add raw data points
  #       cols = c("red", "orange", "yellow", "green", "blue", "purple")
  #       points(vd$time.pt, vd$mort, pch = c(1, 2, 3, 4, 5, 6)[c(as.factor(vd$dose))],
  #              col=rev(cols)[c(as.factor(vd$dose))])
  # 
  #     # Now add the actual mean predictions
  #       for(i in 1:length(unique(vd$dose))){
  #         lines(new.time, inv.logit(mean.preds[[i]]),
  #               lwd=2,
  #               lty=seq(1,6,1)[i],
  #               col=rev(cols)[i]
  #         )
  #       }
  #     # Add the y-axis
  #       axis(2, las=2)
  #       mtext("Probability of mortality", 2, line=3, cex=1.25)
  #       legend(x=47, y=.8,
  #              legend = c('0.00', '0.25', '0.50', '1.00', '2.00', '3.00'),
  #              lty = seq(1,6),
  #              lwd=2,
  #              col=rev(cols)[1:6],
  #              bty='n',
  #              title = 'Concentration (%)'#,
  #              #title.adj = -2
  #              )

  # Individual dosage predictions with uncertainty
    # Predict mortality based on each of the coefficients from the posterior
    # Plot the posterior predictions
      par(mar=c(5,5,1,1))

      # First, make a blank plot to which we will add the others
        plot(0:0, xlim=c(0, max(vd$time.pt)), ylim=c(0,1),
             xlab='Exposure time (hours)',
             ylab = "", cex.lab=1.25, cex.axis=1.15, yaxt = 'n', type='l')

      # Now add the predictions for each draw from the posterior. Change the
      # number that follows 'beta.time.' to access the various posteriors for
      # doses 1-7
        # Create a vector to hold new values for time
          new.time = seq(0, max(vd$time.pt), 1)
        # Make predictions and plot
          for(i in 1:nrow(vd.res)){
            lines(new.time,
                  inv.logit(vd.res$alpha[i]+vd.res$beta.time.1[i]*new.time),
                  lwd=2,
                  lty=1,
                  col='gray87'
            )
          }
        
      # Add the mean of the posterior predictions, making sure to use the same
      # beta for the dose above
        lines(new.time,
              inv.logit(mean(vd.res$alpha)+mean(vd.res$beta.time.1)*new.time),
              col = 'red', lty = 1, lwd = 2
              )
      # Add lines for the 95% credible interval
        lines(new.time,
              inv.logit(
                quantile(vd.res$alpha,0.025)+quantile(vd.res$beta.time.1,0.025)*new.time),
                col = 'red', lty = 2, lwd = 1
              )
        lines(new.time,
              inv.logit(
                quantile(vd.res$alpha,0.975)+quantile(vd.res$beta.time.1,0.975)*new.time),
                col = 'red', lty = 2, lwd = 1
              )
      # Add the y-axis
        axis(2, las=2)
        mtext("Probability of mortality", 2, line=3, cex=1.25)
        
        
  
# Model specification- continuous -----------------------------------------------------
# Write the binomial model to a file in the working directory
# Base model
  modelString = "
    model{

    # Binomial likelihood
      for(i in 1:nruns){
        C[i] ~ dbinom(p[i], N[i])
        logit(p[i]) <- alpha + beta.time*time.pt[i] + beta.dose*dose[i]
      }

    # Intercept term
      alpha ~ dunif(-10, 0)

    # Effect of time as continuous
      beta.time ~ dnorm(0, 0.001)

    # Effect of dose as continuous
      beta.dose ~ dnorm(0, 0.001)

    }"

  # Write the model to a file
    writeLines(modelString, con='vd_doseContinuous.txt')

# Model calibration- continuous --------------------------------------------------------
# Package the data so JAGS can use it
  # Data for the base model
    vd.data = list(
      C = vd$dead,
      nruns = nrow(vd),
      dose = vd$dose,
      N = vd$N,
      time.pt = vd$time.pt
    )

# Supply initial values for Markov Chains in the Gibbs sampler
  # Initial values for the base model
  inits = function(){
    list(
      alpha = runif(1, -10, 0),
      beta.time = rnorm(1, 0, 1),
      beta.dose = rnorm(1, 0, 1)      
    )
  }

# Tell JAGS which parameters you want to trace during posterior sampling
  # Base model
    params = c(
      "alpha",
      "beta.time",
      "beta.dose"
    )

# MCMC settings
  ni = 7500
  nb = 2500
  nt = 30
  nc = 3

# Run the model
  # Base model
    vd.modelC = jags(
      data = vd.data,
      inits = inits,
      parameters.to.save = params,
      model.file = "vd_doseContinuous.txt",
      n.chains = nc,
      n.iter = ni,
      n.burnin = nb,
      n.thin = nt,
      working.directory = getwd(),
      DIC = TRUE
    )
    
# Print the model summary
  print(vd.modelC, digits = 3)     
  
# Contour plots -----------------------------------------------------------
# Make new values of dose and time to predict with
  newdose = sample(unique(vd$dose), 25000, replace = T) #sample(seq(from=0, to=4, by=0.25), 25000, replace = T)  
  newtime = sample(unique(vd$time.pt), 25000, replace = T)#sample(seq(from=min(vd$time.pt), to=max(vd$time.pt), 1),
                   #25000, replace = T)
  
# Collect posterior means
  alpha = vd.modelC$BUGSoutput$mean$alpha
  bdose = vd.modelC$BUGSoutput$mean$beta.dose  
  btime = vd.modelC$BUGSoutput$mean$beta.time 
  
# Predict mortality based on these things
  newmort = inv.logit( alpha + bdose*newdose + btime*newtime )
  
# Make a new dataframe with columns of interest
	persp.test = data.frame(x=newtime,	y=newdose,	z=newmort)
  
# Order the data frame and remove NA values
	persp.test = na.omit(persp.test[with(persp.test, order(x, y)), ])

# Interpolate values of 'z' at evenly spaced values of x and y.
	im = with(persp.test, interp(x, y, z, nx=15, ny=7, duplicate = 'mean'))

# Make the contour plot
    par(mar=c(5, 5.2, 1, 10))
  # filled.contour is the function that actually makes the contour plot
	  filled.contour(
	    im$x,                                 # The variable to be displayed on the x-axis
	    im$y,                                 # The variable to be displayed on the y-axis
	    im$z,                                 # The response variable you wish to plot
	    col=rev(gray.colors(20)),             # Could also choose 'grey.colors' or 'topo.colors'. If you want the ramp to go the other way, just delete the 'rev'. Note that you will need to change the 20 in parentheses to match the number of levels that you actually have or want to display.
	    main = '',                            # I don't like in-figure titles. You can add one, though. You will, however, need to change the 'mar' argument in the call to par above.
	    ylim=c(0,4),                          # Set max and min of y-axis to your data range
	    xlim=c(0,max(vd$time.pt)),            # Set max and min of x-axis to your data range
	    xlab="Time (h)",                      # Change the words in the quotes to change the x-axis label
	    cex.lab=1.5,                          # This makes the labels 1.5x larger than default
	    plot.axes = {                         # This argument tells R to print the axes, but increas the size
	      contour(                            # This is the line that adds the contour lines
	        im$x,                             # The variable to be displayed on the x-axis
	        im$y,                             # The variable to be displayed on the y-axis
	        im$z,                             # The response variable you wish to plot
	        nlevels = 20,                     # This number needs to match the one in 'col' on line 102
	        levels =.5,                       # Draw line for LD50
	        drawlabels = FALSE,               # The labels are realy ugly
	        col = 'black',
	        lwd = 2,
	        lty = 2,
	        add = TRUE                        # Add the lines to the current plot
	      );                                  # Close the call to the contour line function
	      axis(1, cex.axis=1.25);             # X axis tick marks & tick labels
	      axis(2, cex.axis=1.25)              # Y axis tick marks & tick labels
	    }                                     # Close the argument plot.axes
	  )                                       # Close the call to filled.contour
	# Finally, add a label for the y-axis
    mtext(side = 2, "EarthTec (ppm)", line=4, cex.lab=1.5, cex=1.5)
    mtext(side = 4, "Mortality", line=8, cex.lab=1.5, cex=1.5)    
