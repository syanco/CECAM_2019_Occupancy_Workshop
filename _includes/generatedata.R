nvisits <- 3
nsites <- 100

#site level covariates
grassheight <- sample(c("Short","Mid", "Tall"), nsites, replace = T, prob = c(0.25, 0.5, 0.25))
shrubcover <- rbeta(nsites, 0.6, 1)
distance_w <- rpois(nsites, lambda = 2)

#observation level covariates
observer <- sample(c("A","B"), nvisits * nsites, replace = T)

psi <- plogis(-0.5 + (-0.5 * distance_w))

grasscov <- NA
for(i in 1:nsites){
  grasscov[i] <- switch(grassheight[i], "Tall" = -3, "Mid" = -1, "Short" = 0.5)
}

#grasscov <- rep(grasscov, each = 3)

# obscov <- NA
# for(i in 1:(nsites*nvisits)){
#   obscov[i] <- switch(observer[i], "A" = 0.5, "B" = 0)
# }

p <- plogis(2 + grasscov)
#p <- matrix(p, ncol = nvisits)


simoccudata <- function(nsites, nvisits, psi, p){
  #y is the observed detections
  y <- matrix(NA, nsites, nvisits)
  
  #z is the true occupancy at each site
  z <- NA
  for(i in 1:nsites){
    z[i] <- rbinom(1, 1, psi[i])
    y[i,] <- rbinom(nvisits, 1, z[i] * p[i])
  }
  y
}

data <- simoccudata(nsites, nvisits, psi, p)
sitecovs <- data.frame(grassheight, shrubcover, distance_w)
obscovs <- data.frame(observer)

write.csv(data, file = "unicorns.csv", row.names = F)
write.csv(sitecovs, file = "unicorn_sitecovs.csv", row.names = F)
write.csv(obscovs, file = "unicorn_obscovs.csv", row.names = F)
