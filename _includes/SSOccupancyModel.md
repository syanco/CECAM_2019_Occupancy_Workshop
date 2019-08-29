Single season occupancy model tutorial
================
Allison K. Pierce
8/19/2019

Before we get started we need to load in the R packages needed for this
coding tutorial. Make sure you install these packages first before you
try to load
them.

``` r
### Code in this chunk loads some custom function and the R libraries we need for the rest of the tutorial
library(unmarked)
```

    ## Loading required package: lattice

    ## Loading required package: parallel

    ## Loading required package: Rcpp

    ## Loading required package: reshape2

# Exploring study design with simulated data

When designing an occupancy survey you need to decide how many sites or
plots to survey and how many times you will visit each site in order to
get a good estimate of occupancy that addresses your goals for the
study. This might depend on practical matters such as number of
available surveyors, time available to do surveys, money for pay and
travel, etc. However, there are other considerations as well. Is the
area of interest very large or small? Is the organism you are studying
easy to detect or very cryptic? Is the organism common or rare? At what
spatial or temporal scale can you reasonably assume the organism is not
moving in and out of the site (closure)?

These are just a few of the important questions you should think about
when designing an occupancy study. One way to help inform these
decisions is to do a pilot survey which can help give you an idea on the
detectability of the organism and how feasible it is to visit N number
of plots across the study area. This can be especially useful if you are
unfamiliar with the organism and/or study area but does require
additional time.

Alternatively, or in addition to a pilot survey, we can use R to
simulate data using pre-set occupancy probabilities and detection
probabilities. This way we can evaluate how well we can expect to
estimate occupancy and detection under different scenarios. Not only
will this help to design more effcient surveys it will also help us
understand the modelling process.

## Simulating occupancy surveys

First we will start with a simple presence/absence only example where we
are just estimating occupancy and detection probabilities. The nice
thing about a simulation is we can set and know the TRUE occupancy and
detection probabilities and compare them to the estimates from the
model. We can then do experiments ‘in silico’ (in the computer) to
change the number of sites, visits, or any other aspect of the
simulation to see how they affect the model’s ability to accurately and
precisely estimate the true preset occupancy and detection
probabilities.

We will be simulating data many times in this tutorial. To save a lot of
typing I wrote a custom function. Writing functions is beyond the scope
of this tutorial so all you really need to know is how to use it and
what it outputs. The ‘simoccudata’ function will generate stochastic
detection histories based on preset values for number of sites, number
of visits, true occupancy (psi), and true detection (p). Even if you
keep all the presets the same you will get a different set of detection
histories each time you run the function but this is consistent with the
variation seen from sampling error which occurs anytime we only measure
a sample of the whole population (or in this case only visit a few sites
instead of the entire study
area).

``` r
# custom function that simulates detection histories for a specfied number of sites and visits from known occupancy and detection probabilities
simoccudata <- function(nsites, nvisits, psi, p){
  studyarea <- matrix(NA, 100, 100)
  
  #y is the observed detections
  y <- matrix(NA, nsites, nvisits)
  
  #z is the true occupancy at each site
  z <- rbinom(nsites, 1, psi)
  
  for(i in 1:nsites){
    y[i,] <- rbinom(nvisits, 1, z[i] * p)
  }
  y
}
```

Let’s simulate some data\! Let’s assume we have an organism that’s
pretty easy to detect like a unicorn. It has a big horn, white, sparkly,
and is pretty large.

![Unicorns\!](unicorn-2001367_1920.png)

Maybe if it’s laying down we might not see it so let’s say 80% of the
time we will detect it if its there. Let’s also assume we are surveying
a unicorn preserve, eventhough the area is protected unicorns are pretty
rare and only occupy 20% of the sites in the preserve. The unicorns
breed in meadows and avoid the forested areas because thier long horns
get tangled in branches. There are over 1000 meadows in the preserve and
it would be nearly impossible to survey them all before the 5 week
breeding season ends and they leave the preserve. So let’s start out by
simulating data for visiting 100 meadows once a week.

``` r
#constant p and psi
p <- 0.8
psi <- 0.2

set.seed(42)

y <- simoccudata(nsites = 100, nvisits = 5, psi = psi, p = p)
head(y)
```

    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]    1    1    1    1    0
    ## [2,]    0    1    1    1    1
    ## [3,]    0    0    0    0    0
    ## [4,]    1    0    1    1    1
    ## [5,]    0    0    0    0    0
    ## [6,]    0    0    0    0    0

In the simulated data, 0 means we saw no unicorns in the meadow in that
visit and 1 means we saw at least 1 unicorn in the meadow during the
visit. Each row is the full detection history for a meadow and each
column is one visit/survey.

## Fitting a single season occupancy model to simulated data

The ‘unmarked’ package in R is a popular tool for modeling occupancy in
R that has a lot of documentation and support resources but there are
other tools availble both in and outside of the R program. This package
also has functions for other kinds of models to estimate presence and
abundance that is beyond the scope of the workshop. To use the occupancy
modeling functions in ‘unmarked’ we need to format our data in a special
way using the ‘unmarkedFrameOccu’ function and save it as a new object.
We can get a summary of the new object to check our data was reformatted
correctly. This will also tell us how many detections we had.

``` r
occdata <- unmarkedFrameOccu(y=y)
summary(occdata)
```

    ## unmarkedFrame Object
    ## 
    ## 100 sites
    ## Maximum number of observations per site: 5 
    ## Mean number of observations per site: 5 
    ## Sites with at least one detection: 23 
    ## 
    ## Tabulation of y observations:
    ##   0   1 
    ## 405  95

Looks good so let’s fit some models\! For a single season model we will
use the ‘occu’ function in ‘unmarked’. With a single season occupancy
model we are estimating 2 state variables, occupancy probability and
detection probability. We need to supply a formula for each state
variable that describes how mean detection and occupancy varies. First
we will start with the simplest model that mean detection and occupancy
is constant and doesn’t vary as a function of site, visit, etc. This is
denoted with the ‘~1’ formula repeated twice in the function call. In
the ‘occu’ function the first formula describes detection and the second
occupancy. We can get a brief summary of the model results using the
‘summary’ function.

``` r
#Detection and occupancy probability is constant
dotmodel <- occu(~1 ~1, data = occdata)
summary(dotmodel)
```

    ## 
    ## Call:
    ## occu(formula = ~1 ~ 1, data = occdata)
    ## 
    ## Occupancy (logit-scale):
    ##  Estimate    SE     z P(>|z|)
    ##     -1.21 0.238 -5.08 3.7e-07
    ## 
    ## Detection (logit-scale):
    ##  Estimate    SE    z  P(>|z|)
    ##      1.56 0.246 6.32 2.62e-10
    ## 
    ## AIC: 218.1164 
    ## Number of sites: 100
    ## optim convergence code: 0
    ## optim iterations: 13 
    ## Bootstrap iterations: 0

How is detection probability over 1? Notice that these estimates are on
the logit scale and represent log odds not probabilities. To get
probabilities we need to back transform the estimates. The
‘backtransform’ function is one easy way to do this. We can also use
the ‘confint’ function to get 95% confidence intervals for the back
transformed estimates.

``` r
#back transforms the occupancy estimate
psi_estimate <- backTransform(dotmodel, type = "state")
psi_estimate
```

    ## Backtransformed linear combination(s) of Occupancy estimate(s)
    ## 
    ##  Estimate     SE LinComb (Intercept)
    ##      0.23 0.0421   -1.21           1
    ## 
    ## Transformation: logistic

``` r
confint(psi_estimate)
```

    ##      0.025     0.975
    ##  0.1579111 0.3224925

``` r
#back transforms the detection estimate
p_estimate <- backTransform(dotmodel, type = "det")
p_estimate
```

    ## Backtransformed linear combination(s) of Detection estimate(s)
    ## 
    ##  Estimate     SE LinComb (Intercept)
    ##     0.826 0.0354    1.56           1
    ## 
    ## Transformation: logistic

``` r
confint(p_estimate)
```

    ##      0.025     0.975
    ##  0.7454128 0.8849479

With 5 visits over 50 sites we’ve done a fairly decent job estimating
the true occupancy of 0.2 and detection of 0.8. Maybe we can lighten the
workload a little and reduce our number of visits. Let’s run the same
model but simulate data across different numbers of visits ranging from
2 to 10. For now we will keep the TRUE detection and occupancy
probabilities the same. To do this we will repeat the same thing we did
above 9 times changing the number of visits. To make it less tedious we
will use a for loop in our code.

``` r
nvisits <- 2:10
simresults <- data.frame(psi = rep(NA, length(nvisits)), 
                         psi_SE = rep(NA, length(nvisits)),
                         p = rep(NA, length(nvisits)), 
                         p_SE = rep(NA, length(nvisits)))

for(i in 1:length(nvisits)){
  y <- simoccudata(nsites = 100, nvisits = nvisits[i], psi = psi, p = p)
  occdata <- unmarkedFrameOccu(y=y)
  model <- occu(~1 ~1, data = occdata)
  psi_estimate <- backTransform(model, type = "state")
  p_estimate <- backTransform(model, type = "det")
  simresults$psi[i] <- coef(psi_estimate)
  simresults$psi_SE[i] <- SE(psi_estimate)
  simresults$p[i] <- coef(p_estimate)
  simresults$p_SE[i] <- SE(p_estimate)
}

head(simresults)
```

    ##         psi     psi_SE         p       p_SE
    ## 1 0.2112506 0.04991438 0.6153836 0.11227046
    ## 2 0.1829841 0.03911909 0.7469344 0.06349727
    ## 3 0.1801140 0.03844132 0.8467637 0.04270123
    ## 4 0.2300920 0.04210014 0.7909887 0.03805225
    ## 5 0.1600226 0.03666588 0.7707108 0.04296430
    ## 6 0.1899967 0.03922998 0.8195403 0.03334964

Let’s plot the results of our simulation using the estimates and the
standard errors to get an idea of the accuracy and precision of our
results compared to the TRUE occupancy and detection probabilities.

``` r
plot(nvisits, simresults$p, 
    ylim=range(c(0,1)),
    pch=19)
# horizontal error bars
arrows(nvisits, simresults$p - simresults$p_SE,  nvisits, simresults$p + simresults$p_SE, length=0.05, angle=90, code=3)

abline(h = p, col = "red")
title("Detection probability")
```

![](SSOccupancyModel_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
plot(nvisits, simresults$psi, 
    ylim=range(c(0,1)),
    pch=19)
# horizontal error bars
arrows(nvisits, simresults$psi - simresults$psi_SE,  nvisits, simresults$psi + simresults$psi_SE, length=0.05, angle=90, code=3)

abline(h = psi, col = "red")
title("Occupancy probability")
```

![](SSOccupancyModel_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Looks like we could reduce the number of our visits and still do a
decent job of estimating both parameters\! Visiting 100 meadows is still
a lot of work so can we reduce that as well? We will use the same code
as above but modify it to step through number of sites.

``` r
nvisits <- 5
nsites <- seq(from = 5, to = 100, by = 5)
simresults_sites <- data.frame(psi = rep(NA, length(nsites)), 
                         psi_SE = rep(NA, length(nsites)),
                         p = rep(NA, length(nsites)), 
                         p_SE = rep(NA, length(nsites)))

for(i in 1:length(nsites)){
  y <- simoccudata(nsites = nsites[i], nvisits = nvisits, psi = psi, p = p)
  occdata <- unmarkedFrameOccu(y=y)
  model <- occu(~1 ~1, data = occdata)
  psi_estimate <- backTransform(model, type = "state")
  p_estimate <- backTransform(model, type = "det")
  simresults_sites$psi[i] <- coef(psi_estimate)
  simresults_sites$psi_SE[i] <- SE(psi_estimate)
  simresults_sites$p[i] <- coef(p_estimate)
  simresults_sites$p_SE[i] <- SE(p_estimate)
}

head(simresults_sites)
```

    ##         psi     psi_SE         p       p_SE
    ## 1 0.4000038 0.21909122 0.8999911 0.09489298
    ## 2 0.2000017 0.12649229 0.8999904 0.09489327
    ## 3 0.3333258 0.12171772 0.8799714 0.06502788
    ## 4 0.2500258 0.09683472 0.8399152 0.07343443
    ## 5 0.2417017 0.08604773 0.6288743 0.09065688
    ## 6 0.2334085 0.07724520 0.7997413 0.06785283

Let’s plot\!

``` r
plot(nsites, simresults_sites$p, 
    ylim=range(c(0,1)),
    pch=19)
# horizontal error bars
arrows(nsites, simresults_sites$p - simresults_sites$p_SE,  nsites, simresults_sites$p + simresults_sites$p_SE, length=0.05, angle=90, code=3)

abline(h = p, col = "red")
title("Detection probability")
```

![](SSOccupancyModel_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
plot(nsites, simresults_sites$psi, 
    ylim=range(c(0,1)),
    pch=19)
# horizontal error bars
arrows(nsites, simresults_sites$psi - simresults_sites$psi_SE,  nsites, simresults_sites$psi + simresults_sites
       $psi_SE, length=0.05, angle=90, code=3)

abline(h = psi, col = "red")
title("Occupancy probability")
```

![](SSOccupancyModel_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

We could probably reduce the number of sites we visit but probably not
by more than
half.

# Fitting models with covariates

### interactive scripting session - will have short writeup here for reference

``` r
unicorndata <- read.csv("unicorns.csv")
sitecovs <- read.csv("unicorn_sitecovs.csv")
obscovs <- read.csv("unicorn_obscovs.csv")

occdata <- unmarkedFrameOccu(y=unicorndata, siteCovs = sitecovs, obsCovs = obscovs)
summary(occdata)
```

    ## unmarkedFrame Object
    ## 
    ## 100 sites
    ## Maximum number of observations per site: 3 
    ## Mean number of observations per site: 3 
    ## Sites with at least one detection: 16 
    ## 
    ## Tabulation of y observations:
    ##   0   1 
    ## 264  36 
    ## 
    ## Site-level covariates:
    ##  grassheight   shrubcover          distance_w  
    ##  Mid  :53    Min.   :0.0000361   Min.   :0.00  
    ##  Short:23    1st Qu.:0.0843713   1st Qu.:1.00  
    ##  Tall :24    Median :0.3128667   Median :2.00  
    ##              Mean   :0.3628649   Mean   :2.29  
    ##              3rd Qu.:0.5405953   3rd Qu.:3.00  
    ##              Max.   :0.9850743   Max.   :6.00  
    ## 
    ## Observation-level covariates:
    ##  observer
    ##  A:154   
    ##  B:146

``` r
m1 <- occu(~1 ~1, data = occdata)
m2 <- occu(~1 ~grassheight, data = occdata)
m3 <- occu(~1 ~shrubcover, data = occdata)
m4 <- occu(~1 ~distance_w, data = occdata)
m5 <- occu(~observer ~1, data = occdata)
m6 <- occu(~observer ~grassheight, data = occdata)
m7 <- occu(~observer ~shrubcover, data = occdata)
m8 <- occu(~observer ~distance_w, data = occdata)
m9 <- occu(~grassheight ~1, data = occdata)
m10 <- occu(~grassheight ~shrubcover, data = occdata)
m11 <- occu(~grassheight ~distance_w, data = occdata)
```

``` r
fitmodels <- fitList("p(.)psi(.)" = m1,
               "p(.)psi(grassheight)" = m2,
               "p(.)psi(shrubcover)" = m3,
               "p(.)psi(distance_w)" = m4,
               "p(observer)psi(.)" = m5,
               "p(observer)psi(grassheight)" = m6,
               "p(observer)psi(shrubcover)" = m7,
               "p(observer)psi(distance_w)" = m8,
               "p(grassheight)psi(.)" = m9,
               "p(grassheight)psi(shrubcover)" = m10,
               "p(grassheight)psi(distance_w)" = m11)
```

``` r
modSel(fitmodels)
```

    ##                               nPars    AIC delta  AICwt cumltvWt
    ## p(grassheight)psi(distance_w)     5 138.81  0.00 0.6154     0.62
    ## p(.)psi(distance_w)               3 141.66  2.85 0.1483     0.76
    ## p(grassheight)psi(.)              4 142.65  3.83 0.0904     0.85
    ## p(observer)psi(distance_w)        4 143.66  4.84 0.0546     0.91
    ## p(grassheight)psi(shrubcover)     5 144.26  5.44 0.0405     0.95
    ## p(.)psi(.)                        2 145.37  6.56 0.0232     0.97
    ## p(.)psi(shrubcover)               3 146.92  8.10 0.0107     0.98
    ## p(observer)psi(.)                 3 147.37  8.56 0.0085     0.99
    ## p(observer)psi(shrubcover)        4 148.92 10.10 0.0039     1.00
    ## p(.)psi(grassheight)              4 149.30 10.49 0.0033     1.00
    ## p(observer)psi(grassheight)       5 151.30 12.48 0.0012     1.00
