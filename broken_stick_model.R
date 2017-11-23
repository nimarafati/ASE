#
#
# broken stick regression model (modified from script provided by Allan E. Strand)
#
#
bs <- function(X,Y) 
  {

#    initial values
    p <- rep(0,4)
    p[1:2] <- c(quantile(X,0.25),mean(Y)) #point 1
    p[3:4] <- c(quantile(X,0.75),mean(Y)) #point 2
    optim(p,llike.bs,X=X,Y=Y)
  }

flt <- function(X,Y)
  {
    optim(p=c(0.5),llike.flt,X=X,Y=Y)
  }


llike.bs <- function(p,X,Y,weight=0.1)# current objective function but sets sd to
                              #the sd of allele freqs across cline
  {
    p1 <- p[1:2]
    p2 <- p[3:4]
    sdY <- sd(Y)
#    sd=0.2
    if ((p1[1]<p2[1])&(sdY>0)) {ok <- T} else {ok <- F}# points not swapped, sd ok
    if (ok & ((p1[1]>min(X))&(p2[1]<max(X)))) {ok <- T} else {ok <- F}# points not at end of line
    if (ok & (min(c(p1[2],p2[2])>=0))) {ok <- T} else {ok <- F} #range finding
    if (ok & (max(c(p1[2],p2[2])<=1))) {ok <- T} else {ok <- F} #range finding
    
    if (ok)
      {
        middleX <- X[which((X>p1[1])&(X<p2[1]))]
        middleY <- Y[which((X>p1[1])&(X<p2[1]))]
#        print(middleX)
        if (length(middleX)<=1) #applies penalties for slopes that either lack intermediate points or
                                #have x-param values that differ from actual populations
                                # weight determines how important the left-right penalty will be
          {
            leftpenalty <- weight*(p1[1]-max(X[X<=p1[1]]))
            rightpenalty <- weight*(min(X[X>=p2[1]])-p2[1])
            penalty <- -1*(leftpenalty+rightpenalty)
#            cat(paste("[",leftpenalty,", ",rightpenalty,"]; "))
          }
        else
          {
            penalty <- 0
          }
        
        leftX <- X[which(X<=p1[1])]
        leftY <- Y[which(X<=p1[1])]
        rightX <- X[which(X>=p2[1])]
        rightY <- Y[which(X>=p2[1])]

        resleft <- p1[2]-leftY
        resright <- p2[2]-rightY
        
        slp <- (p2-p1)[2]/(p2-p1)[1] #rise over run!
        mnp <- (p2+p1)/2
        int = mnp[2]-mnp[1]*slp
        resmid = middleY-(middleX*slp + int)

        -sum(c(penalty,dnorm(c(resleft),mean=0,sd=sdY,log=T),dnorm(c(resmid),mean=0,sd=sdY,log=T),dnorm(c(resright),mean=0,sd=sdY,log=T)))
      } else {Inf}
  }


llike.flt <- function(p,X,Y)
  {
    sdY <- sd(Y)
    -sum(dnorm(Y-p[1],mean=0,sd=sdY,log=T))
  }




extract.bs <- function(bs.par)
  {
    slope <- (bs.par[2]-bs.par[4])/(bs.par[1]-bs.par[3])
    mid <- mean(bs.par[c(1,3)])
    c(slope=slope,mid=mid)
  }

