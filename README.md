# zebraMusselMortStudies
R code for estimating zebra mussel mortality from Bayesian
logistic regression models.

<br>
 
## File descriptions

`virkon_Bayes_noIntercept.R`
  - Bayesian hierarchical models for estimating mortality
    of zebra mussel veligers following exposure to Virkon
    at different time points and concentrations.
  - Response is binomial (number dead out of number initial)
  - Concentration as a categorical effect
 
 `citricAcidAnalysis.R`
  - Bayesian hierarchical models for estimating mortality
    of zebra mussel veligers following exposure to C.A.
    at different time points and concentrations.
  - Essentially identical to `virkon_Bayes_noIntercept.R`
    with differences in data manipulation and screening
 
 `earthTech_Bayes.R`
  - Bayesian logistic regression models for estimating 
    probability of zebra mussel mortality following long-term
    exposure to EarthTec
  - Response is binomial (number dead out of number initial)
  - Concentration and time points both as continuous covariates
  
  
