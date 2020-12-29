convolve <- function(dens1, dens2,
                     cdf1=Vectorize(function(x){integrate(dens1,-Inf,x)$value}),
                     cdf2=Vectorize(function(x){integrate(dens2,-Inf,x)$value}),
                     delta=0.01, epsilon=0.0001)
{
  stopifnot(is.function(dens1), is.function(dens2),
            is.function(cdf1), is.function(cdf2),
            delta>0, epsilon>0)
  
  symKL <- function(d)
  {
    stopifnot(d>=0)
    func1 <- function(x)
    {
      d2md <- dens2(x-d)
      if (is.element("log", names(formals(dens2)))) {
        logd2   <- dens2(x, log=TRUE)
        logd2md <- dens2(x-d, log=TRUE)
      } else {
        logd2   <- log(dens2(x))
        logd2md <- log(d2md)
      }
      return(ifelse((logd2>-Inf) & (logd2md>-Inf), (logd2md-logd2)*d2md, 0.0))
    }
    
    func2 <- function(x)
    {
      d2 <- dens2(x)
      if (is.element("log", names(formals(dens2)))) {
        logd2   <- dens2(x, log=TRUE)
        logd2md <- dens2(x-d, log=TRUE)
      } else {
        logd2   <- log(d2)
        logd2md <- log(dens2(x-d))
      }
      return(ifelse((logd2>-Inf) & (logd2md>-Inf), (logd2-logd2md)*d2, 0.0))
    }
    int1 <- integrate(func1, -Inf, Inf)
    if (int1$message != "OK")
      warning(paste0("Problem computing KL-divergence (1): \"", int1$message,"\""))
    int2 <- integrate(func2, -Inf, Inf)
    if (int2$message != "OK")
      warning(paste0("Problem computing KL-divergence (2): \"", int2$message,"\""))
    return(int1$value + int2$value)
  }
  
  step <- sqrt(delta)
  while (symKL(step) < delta) step <- 2*step
  ur <- uniroot(function(X){return(symKL(X)-delta)}, lower=0, upper=step)
  step <- ur$root
  mini <- -1
  while (cdf1(mini) > epsilon/2) mini <- 2*mini
  maxi <- 1
  while (cdf1(maxi) < 1-(epsilon/2)) maxi <- 2*maxi
  ur <- uniroot(function(X){return(cdf1(X)-epsilon/2)},
                lower=mini, upper=maxi)
  mini <- ur$root
  ur <- uniroot(function(X){return(cdf1(X)-(1-epsilon/2))},
                lower=mini, upper=maxi)
  maxi <- ur$root
  k <- ceiling((maxi-mini)/(2*step))+1
  support <- mini - ((k*2*step)-(maxi-mini))/2 + (0:(k-1))*2*step
  
  margins <- support[-1]-step
  
  weight <- rep(NA, length(support))
  for (i in 1:(k-1))
    weight[i] <- cdf1(margins[i])
  weight[k] <- 1
  for (i in k:2)
    weight[i] <- weight[i]-weight[i-1]
  grid <- cbind("lower"=c(-Inf, margins),
                "upper"=c(margins, Inf),
                "reference"=support,
                "prob"=weight)
  
  density <- function(x)
  {
    return(apply(matrix(x,ncol=1), 1,
                 function(x){sum(grid[,"prob"]*dens2(x-grid[,"reference"]))}))
  }
  
  cdf <- function(x)
  {
    return(apply(matrix(x,ncol=1), 1,
                 function(x){sum(grid[,"prob"]*cdf2(x-grid[,"reference"]))}))
  }
  
  quantile <- function(p)
  {
    quant <- function(pp)
    {
      mini <- -1
      while (cdf(mini) > pp) mini <- 2*mini
      maxi <- 1
      while (cdf(maxi) < pp) maxi <- 2*maxi
      ur <- uniroot(function(x){return(cdf(x)-pp)}, lower=mini, upper=maxi)
      return(ur$root)      
    }
    proper <- ((p>0) & (p<1))
    result <- rep(NA,length(p))
    if (any(proper)) result[proper] <- apply(matrix(p[proper],ncol=1), 1, quant)
    return(result)
  }
  
  return(list("delta"    = delta,     
              "epsilon"  = epsilon,   
              "binwidth" = 2*step,    
              "bins"     = k,         
              "support"  = grid,      
              "density"  = density,   
              "cdf"      = cdf,       
              "quantile" = quantile))
}