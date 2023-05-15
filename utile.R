
VonNeumann <- function(n, p=1, graine)
{
  x <-  rep(graine,n*p+1)
  for(i in 2:(n*p+1))
  {
    numbers <- strsplit(format(x[i-1]^2,scientific=FALSE),'')[[1]]
    while(length(numbers)>4){ 
      numbers <- numbers[2:(length(numbers)-1)] 
    }
    x[i] <- as.numeric(numbers)%*%(10^seq(length(numbers)-1,0,-1))
  }
  x <- matrix(x[2:(n*p+1)],nrow=n,ncol=p)
  return(x)
}


MersenneTwister <- function(n, p=1, graine)
{
  set.seed(graine,kind='Mersenne-Twister')
  x <- sample.int(2^32-1,n*p)
  x <- matrix(x,nrow=n,ncol=p)
  return(x)
}


RANDU <- function(n, p=1, graine)
{
  a<-65539
  b<-0
  m<-2^31
  x <-  rep(graine,n*p+1)
  for(i in 2:(n*p+1))
  {
    x[i] <- (a*x[i-1]+b)%%m
  }
  x <- matrix(x[2:(n*p+1)],nrow=n,ncol=p)
  return(x)
}


StandardMinimal <- function(n, p=1, graine)
{
  a<-16807
  b<-0
  m<-2^31-1
  x <-  rep(graine,n*p+1)
  for(i in 2:(n*p+1))
  {
    x[i] <- (a*x[i-1]+b)%%m
  }
  x <- matrix(x[2:(n*p+1)],nrow=n,ncol=p)
  return(x)
}

binary <- function(x)
{
  if((x<2^31)&(x>=0))
    return( as.integer(intToBits(as.integer(x))) )
  else{
    if((x<2^32)&(x>0))
      return( c(binary(x-2^31)[1:31], 1) )
    else{
      cat('Erreur dans binary : le nombre etudie n est pas un entier positif en 32 bits.\n')
      return(c())
    }
  }
}

Frequency <- function(x, nb)
{
  sum <- 0
  for(xElem in x)
  {
    eps <- binary(xElem)
    for(i in 1:nb)
    {
      sum <- sum + (2*eps[i]-1)
    }
  }
  
  sObs = abs(sum)/sqrt(nb*length(x))
  
  Pvaleur = 2*(1-pnorm(sObs))
  
  return(Pvaleur)
}


Runs <-  function(x,nb)
{
  cpt <- 0
  n <- nb*length(x)
  t <- 2/sqrt(n)
  for(xElem in x)
  {
    eps <- binary(xElem)
    
    for(i in 1:nb)
    {
      if (eps[i]==1) {
        cpt<-cpt+1;
      }
    }
    
  }
  
  pi <- cpt/n
  
  if (abs(pi-1/2) < t) {
    vObs <- 1
    lastBit <- 0
    
    for (k in 1:length(x)) {
      eps <- binary(x[k])
      
      for (i in 1:nb)
      {
        if (k!=1) {
          if(i==1) {
            if (eps[i] != lastBit) {
              vObs <- vObs + 1
            }
          } else {
            if (eps[i] != eps[i-1]) {
              vObs <- vObs + 1
            }
          }
        } else {
          if(i!=1) {
            if (eps[i] != eps[i-1]) {
              vObs <- vObs + 1
            }
            if (i==nb) {
              lastBit <- eps[i]
            }
          }
        }
      }
    }
    
    Pvaleur <- 2*( 1-pnorm( abs(vObs-2*n*pi*(1-pi)) / (2*sqrt(n)*pi*(1-pi)) ) )
  } else {
    Pvaleur <- 0.0
  }
  return(Pvaleur)
}










