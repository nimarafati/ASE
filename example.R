source("broken.R")

#create some fake data
indat <- data.frame(latitude=0:12)
indat$allele1 <- 0.5+rnorm(dim(indat)[1],sd=0.05) #no relation between latitude and freq
c=6
w=10
indat$allele2 <- rnorm(dim(indat)[1],sd=0.05) + (1+tanh(2*(indat$latitude-c)/w))/2

# fit the broken stick to allele1
bs.a1 <- bs(indat$latitude,indat$allele1)
#fit a broken stick to allele2
bs.a2 <- bs(indat$latitude,indat$allele2)

#The estimate of slope and midpoint can be extracted by:
extract.bs(bs.a1$par)
extract.bs(bs.a2$par)

