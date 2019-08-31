nvisits <- 4
nsites <- 150

#site level covariates
grassheight <- sample(c("Short","Mid", "Tall"), nsites, replace = T, prob = c(0.35, 0.4, 0.25)) #covaries with p
shrubcover <- runif(nsites,0,0.6) #does not covary with p or psi
distance_w <- rpois(nsites, lambda = 2.5) + runif(nsites,0,1) #covaries with psi

#observation level covariates
observer <- sample(c("A","B"), nvisits * nsites, replace = T) #does not covary with p or psi

#simulate true psi with negative relationship with distance to water
psi <- plogis(0.3 + (-0.8 * distance_w))

#create vector of effects for categorical data for generating p
grasscov <- NA
for(i in 1:nsites){
  grasscov[i] <- switch(grassheight[i], "Tall" = -1, "Mid" = 1, "Short" = 1.5)
}

#simulate true detection at different levels of the site level categorical variable
p <- plogis(grasscov)


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

#make sure relationships are as we intended
plot(psi, distance_w)
boxplot(p~grassheight)

#output data for analysis
data <- simoccudata(nsites, nvisits, psi, p)
data <- as.data.frame(data)
names(data) <- c("week1", "week2", "week3", "week4")
sitecovs <- data.frame(grassheight, shrubcover, distance_w)
obscovs <- data.frame(observer)

write.csv(data, file = "unicorns.csv", row.names = F)
write.csv(sitecovs, file = "unicorn_sitecovs.csv", row.names = F)
write.csv(obscovs, file = "unicorn_obscovs.csv", row.names = F)


