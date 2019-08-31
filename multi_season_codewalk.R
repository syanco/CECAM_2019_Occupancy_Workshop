#####################################################
##  SINGLE-SPECIES, MULTI-SEASON OCCUPANCY MODELS  ##
#####################################################

#make sure you load the library unmarked
library(unmarked)

#############################
##  Create simulated data  ##
#############################

#create a site index
site <- 1:100

#create 'structure' to the data by randomizing detection as functuion of site
multi_dat <- data.frame("visit1" = rbinom(100, 1, site/100), 
                        "visit2" = rbinom(100, 1, site/100),
                        "visit3" = rbinom(100, 1, site/100),
                        "visit4" = rbinom(100, 1, site/100),
                        "visit5" = rbinom(100, 1, site/100),
                        "visit6" = rbinom(100, 1, site/100),
                        "visit7" = rbinom(100, 1, site/100),
                        "visit8" = rbinom(100, 1, site/100),
                        "visit9" = rbinom(100, 1, site/100))

#simulate site-level covars not that cov1 should be correlated with occupancy
#cov2 should not...
site_covars <- data.frame("cov1" = site/100, "cov2" = runif(100, 0, 1))

#simulate year-site-level covariates - we'll only do dummy var for year itself
yearly_site_covars <- data.frame("year" = rep(c("01", "02", "03"), 100))

#################################
## Load into `unmarked` format ##
#################################

#let's look at the documentation for `unmarkedMultFrame`
?unmarkedMultFrame

#create the frame
umf.mult <- unmarkedMultFrame(y = multi_dat, siteCovs = site_covars, 
                              yearlySiteCovs = yearly_site_covars, 
                              numPrimary = 3)

## Time to fit models ##

#have a peek at the colext function help file
?colext

#intercept only model
fm1 <- colext(~1, ~1, ~1, ~1, umf.mult)

#year-dependent col and ext
fm2 <- colext(~1, ~year, ~year, ~1, umf.mult)

#fully time dependent model
fm3 <- colext(~1, ~year, ~year, ~year, umf.mult)

#occu cov1, time invariant
fm4 <- colext(~cov1, ~1, ~1, ~1, umf.mult)

#occu cov2, time invariant
fm5 <- colext(~cov2, ~1, ~1, ~1, umf.mult)

#occu cov1, col anmd ext by year
fm6 <- colext(~cov1, ~year, ~year, ~1, umf.mult)

#occu cov2, col anmd ext by year
fm7 <- colext(~cov2, ~year, ~year, ~1, umf.mult)

#let's compare models
models <- fitList('fm1: psi(.)gam(.)eps(.)p(.)' = fm1,
                  'fm2: psi(.)gam(Y)eps(Y)p(.)' = fm2,
                  'fm3: psi(.)gam(Y)eps(Y)p(Y)' = fm3,
                  'fm4: psi(Cov1)gam(.)eps(.)p(.)' = fm4,
                  'fm5: psi(Cov2)gam(.)eps(.)p(.)' = fm5,
                  'fm6: psi(Cov1)gam(Y)eps(Y)p(Y)' = fm6,
                  'fm7: psi(Cov2)gam(Y)eps(Y)p(Y)' = fm7)

#summarize model selction
(ms_table <- modSel(models))

#pull out components of models
coef(ms_table) # Estimates only
SE(ms_table) # Standard errors only

#also can store everything into a df
toExport <- as(ms_table, "data.frame") # Everything
str(toExport)

#let's examine the top model for performance
summary(fm4)
#in fact there is even more stored in our model objects
names(fm4)

#remember our estimates are logit scaled, so let's back transform
backTransform(fm4, type="det") #gets coeff off the logit scale
confint(backTransform(fm4, type="det"))

#the above only wored becausr there was no covariate, for psi we need to pick
#a value of cov1 in this case
backTransform(linearComb(fm4, c(1,0), type="psi"))
#this requires a pretty good understanding of linear models to make sense...

# plot effect of mixed forest from mf8
newdata <- data.frame("cov1" = seq(0, 1, by=.01))
model_predicted <- predict(fm4, type="psi", newdata=newdata, appendData=TRUE)
with(model_predicted, {
  plot(cov1, Predicted, 
       ylim=c(0,1.3), 
       type="l",
       xlab="cov1",
       ylab=expression(hat(psi)), cex.lab=0.8, cex.axis=0.8)
  lines(cov1, Predicted+1.96*SE, col=gray(0.7)) #get model CIs
  lines(cov1, Predicted-1.96*SE, col=gray(0.7))
})
