

Break = "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"

notice = "    \"metaling\", A suite of R functions for Modern Meta-Analysis in Second Language Research.
    Copyright (C) 2019-present  Reza Norouzian, rnorouzian@gmail.com\n"

message(Break, notice, Break)


#==================================================================================================================

#palette('R3')

#==================================================================================================================


autoreg <- function(steps, r){
  
  steps <- if(length(steps) == 1) steps+1 else length(steps)+1
  x <- diag(steps)
  r <- data.frame(r^abs(row(x)-col(x)))
  rownames(r) <- colnames(r) <- c("pre", paste0("post", 1:(steps-1)))
  return(r)
} 


#===============================================================================================================================

t2d <- function(t, n1, n2 = NA){
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  t/sqrt(N)
}

#===============================================================================================================================

d2t <- function(d, n1, n2 = NA){
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  d*sqrt(N)
}

#===============================================================================================================================

trim <- function(X){
  X <- setNames(X, trimws(names(X)))
  y <- sapply(names(X), function(x) is.character(as.vector(X[[x]])))
  X[y] <- lapply(X[y], trimws)
  return(X)
}
 
#===============================================================================================================================

isFALSE <- function(x) identical(FALSE, x)              
              
#===============================================================================================================================

sdif <- function(n = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, t.pair = NA, df = NA,
                 sdp = NA){
  
  n <- ifelse(!is.na(n), n, ifelse(is.na(n) & !is.na(df), df + 1, NA))
  
  ifelse(!is.na(r) & !is.na(sdpre) & !is.na(sdpos), sqrt(sdpre^2+sdpos^2-2*r*sdpre*sdpos),
         ifelse(!is.na(n) & is.na(r) & !is.na(t.pair) & !is.na(mpre) & !is.na(mpos), sqrt((n*(mpos - mpre)^2)/t.pair^2), 
                ifelse(!is.na(r) & !is.na(sdp), sqrt(2*sdp^2*(1-r)), NA)))
}

#===============================================================================================================================

rdif <- function(n = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, t.pair = NA, df = NA, sdif = NA, 
                 sdp = NA) {
  
  sdif <- ifelse(is.na(sdif), sdif(sdpre = sdpre, sdpos = sdpos, t.pair = t.pair, n = n, mpos = mpos, mpre = mpre, df = df, sdp = sdp), sdif)
  
  ifelse(!is.na(sdif) & is.na(sdp) & !is.na(sdpre) & !is.na(sdpos), (sdpre^2 + sdpos^2 - sdif^2)/(2*sdpre*sdpos), 
         ifelse(!is.na(sdp) & !is.na(sdif), 1 - (sdif^2/(2*sdp^2)), NA))
}

#===============================================================================================================================

cfactor <- function(df) exp(lgamma(df/2)-log(sqrt(df/2)) - lgamma((df-1)/2))

#===============================================================================================================================

reget <- function(List, what){
  
  s <- substitute(what)  
 
  if(class(List)[1] != "list") List <- list(List)  
    
  h <- lapply(List, function(x) do.call("subset", list(x, s)))
  
  res <- Filter(NROW, h)
  
  if(length(res) == 0) NULL else res
}

#===============================================================================================================================
              
odiag <- function(x) x[(n <- nrow(x))^2-(1:n)*(n-1)]
              
#===============================================================================================================================
              
get.uni <- function(data, what){
  
  data$study.name <- trimws(data$study.name)
  m <- split(data, data$study.name)
  m <- Filter(NROW, rm.allrowNA2(m)) 
  
  G <- substitute(what)
  E <- quote(x$x)
  E[[3]] <- G[[2]]
  G[[2]] <- E
  
  f <- sapply(m, function(x) sum(eval(G)) == nrow(x))
  
  h <- m[names(f)[f]]
  
  res <- Filter(NROW, h)
  
  if(length(res) == 0) NULL else res
}


#===============================================================================================================================


get.gen <- function(data, what){
  
  s <- substitute(what)  
  data$study.name <- trimws(data$study.name)
  m <- split(data, data$study.name)
  m <- Filter(NROW, rm.allrowNA2(m)) 
  
  h <- lapply(m, function(x) do.call("subset", list(x, s)))
  
  res <- Filter(NROW, h)
  
  if(length(res) == 0) NULL else res
}

#===============================================================================================================================

find.stud <- function(data, what, timevar = TRUE){
  
  s <- substitute(what)
  data$study.name <- trimws(data$study.name)
  if(!timevar) { unique(as.vector(subset(data, eval(s))$study.name))
    
  } else {
    chep <- sort(unique(na.omit(data$time)))
    G <- lapply(chep, function(x) bquote(.(s) & time == .(x)))
    setNames(lapply(seq_along(G), function(j) unique(as.vector(subset(data, eval(G[[j]]))$study.name))), as.character(chep))
  }
}       

#===============================================================================================================================
                    
find.miss <- function(data, space = FALSE, all = FALSE){    
 
res <- Filter(length, lapply(data, function(x) which(if(!space & !all) is.na(x) else if(space & !all) x == "" else is.na(x) | x == "")))
if(length(res) == 0) NA else res
}                    
                    
#===============================================================================================================================
                  
mod.level <- function(data, what){
  
sort(unique(na.omit(unlist(data[paste0(substitute(what))]))))
  
}

#===============================================================================================================================
                                                     
mods.level <- function(data){
  
f <- function(data, what) sort(unique(na.omit(unlist(data[what]))))

  ar <- c(formalArgs(d.prepos)[-c(21, 22)], c("dint", "SD", "id"))
  
  dot.names <- names(data)[!names(data) %in% ar]
  
setNames(lapply(seq_along(dot.names), function(i) f(data = data, what = dot.names[i])), dot.names)
}                  
                  
#===============================================================================================================================
              
              
cor.mat <- function(r, dim) { 
  
  m <- matrix(r, dim, dim)
  diag(m) <- 1
  m
}              
              
#===============================================================================================================================
              
option1 <- function(ds, sds, r = .5){ 
  
V <- sds^2  
 
m <- length(V)

r <- cor.mat(r, m)

SD <- sqrt((1/m)^2 * sum(sqrt(outer(V, V)) * r))

D <- mean(ds)

return(c(D, SD))
}

#===============================================================================================================================
              
decimal <- function(x, k = 3) format(round(x, k), nsmall = k)              

#===============================================================================================================================
                
roundi <- function(x, digits = 7){

if(!inherits(x, "data.frame")) stop("Only used for a 'data.frame'.", call. = FALSE)
  
num <- sapply(x, is.numeric)

x[num] <- lapply(x[num], round, digits)

return(x)
}                
                
#===============================================================================================================================
              
cov.dint <- function(sds, r = .5, no.names = TRUE){

m <- length(sds) 

r <- cor.mat(r, m)

D <- diag(sds)

m <- D%*%r%*%D

if(!no.names) rownames(m) <- colnames(m) <- paste0("d", 1:length(sds))
return(m)
}              

              
#===============================================================================================================================
              
option2 <- function(ds, sds, r = .5){
    
    d <- matrix(ds)
    
    if(length(d) == 1) { return(c(d, sds)) }
    
    r <- cor.mat(r, length(d))
    
    e <- matrix(rep(1, length(d)))
    
    A <- cov.dint(sds, r)
    
    w <- t((solve(A)%*%e)%*%solve((t(e)%*%solve(A)%*%e)))

    se <- as.vector(sqrt(solve(t(e)%*%solve(A)%*%e)))
  
    ave.d <- as.vector(w%*%d)
    
    return(c(ave.d, se))
}              
              
              
#===============================================================================================================================
              
fuse <- function(..., per.study){
  
  ll <- per.study
  
  L <- list(...)
  
  if(all(sapply(list(...), inherits, "list"))){
    
  g <- lapply(1:length(L), function(i) split(L[[i]], rep(seq_along(ll), ll)))
  
  h <- lapply(1:length(L), function(i) lapply(g[[i]], function(x) do.call(rbind, x)))
  
  lapply(1:length(h), function(i) Filter(NROW, h[[i]]))
  
  } else {
    
    g <- split(L, rep(seq_along(ll), ll))
    
    h <- lapply(g, function(x) do.call(rbind, x))
    
    Filter(NROW, h)
  }
}              
              
#===============================================================================================================================
               
pair1 <- function(j, k){
  lapply(seq_along(j), function(i) {x1 <- expand.grid(d1 = j[[i]]$d, d2 = k[[i]]$d); 
  row.names(x1) <- c(outer(row.names(j[[i]]), row.names(k[[i]]), FUN = paste)); 
  setNames(split(as.matrix(x1), row(x1)), paste(names(k[i]), row.names(x1), sep = ""))})
} 

#===============================================================================================================================                
                
pair2 <- function(j, k){ 
  f1 <- function(x, y) {
    dat <- expand.grid(x, y)
    split(as.matrix(dat), row(dat))
  }
  do.call(Map, c(f = f1, list(lapply(j, `[[`, "d"), lapply(k, `[[`, "d"))))
}
                
#===============================================================================================================================
                
pair <- function(j, k){
  lapply(seq_along(j), function(i) {
    x1 <- expand.grid(d1 = j[[i]]$d, d2 = k[[i]]$d);
    split(as.matrix(x1), row(x1))})
}                
                
#===============================================================================================================================

                
dit1 <- Vectorize(function(dppc, dppt, nc, nt, n.sim = 1e5){
  
  like1 <- function(x) dt(dppc*sqrt(nc), df = nc - 1, ncp = x*sqrt(nc))
  like2 <- function(x) dt(dppt*sqrt(nt), df = nt - 1, ncp = x*sqrt(nt))
  
  d1 <- AbscontDistribution(d = like1)
  d2 <- AbscontDistribution(d = like2)
  
  dif <- distr::r(d2 - d1)(n.sim)
  
  din <- dppt - dppc
   SD <- sd(dif)
  
  return(c(dint = din, SD = SD))
})                      
                
#===============================================================================================================================
                
                
dit2 <- Vectorize(function(dppc, dppt, nc, nt, rev.sign = FALSE){
  
  a <- dppc
  b <- dppt
  
  din <- b - a  
  
  test <- if(!rev.sign || rev.sign & b < 0 & a < 0 & abs(b) < abs(a)) FALSE else TRUE
  
  like1 <- function(x) dt(dppc*sqrt(nc), df = nc - 1, ncp = x*sqrt(nc))
  like2 <- function(x) dt(dppt*sqrt(nt), df = nt - 1, ncp = x*sqrt(nt))
  
  d1 <- AbscontDistribution(d = like1, low1 = -15, up1 = 15, withStand = TRUE)
  d2 <- AbscontDistribution(d = like2, low1 = -15, up1 = 15, withStand = TRUE)
  
  like.dif <- function(x) distr::d(if(test) -(d2 - d1) else d2 - d1)(x)
  
  din <- if(test) -din else din
  
  Mean <- integrate(function(x) x*like.dif(x), -Inf, Inf)[[1]]
  SD <- sqrt(integrate(function(x) x^2*like.dif(x), -Inf, Inf)[[1]] - Mean^2)
  
  return(c(dint = din, SD = SD))
})

             
#===============================================================================================================================
             
dit3 <- Vectorize(function(dppc, dppt, nc, nt, n.sim = 1e5, rev.sign = FALSE){
  
  like1 <- function(x) dt(dppc*sqrt(nc), df = nc - 1, ncp = x*sqrt(nc))
  like2 <- function(x) dt(dppt*sqrt(nt), df = nt - 1, ncp = x*sqrt(nt))
  
  d1 <- AbscontDistribution(d = like1, low1 = -15, up1 = 15, withStand = TRUE)
  d2 <- AbscontDistribution(d = like2, low1 = -15, up1 = 15, withStand = TRUE)
  
  dif <- distr::r(d2 - d1)(n.sim)
  
  a <- dppc
  b <- dppt
  
  din <- b - a 
  
  di <- ifelse(!rev.sign || rev.sign & b < 0 & a < 0 & abs(b) < abs(a), din, -din)
  
  SD <- sd(dif)
  
  return(c(dint = di, SD = SD))
})                     

#===============================================================================================================================
             
dit4 <- Vectorize(function(dppc, dppt, nc, nt, n.sim = 1e5, rev.sign = FALSE){
  
  a <- dppc
  b <- dppt
  
  din <- b - a 
  
  test <- if(!rev.sign || rev.sign & b < 0 & a < 0 & abs(b) < abs(a)) FALSE else TRUE
  
  like1 <- function(x) dt(dppc*sqrt(nc), df = nc - 1, ncp = x*sqrt(nc))
  like2 <- function(x) dt(dppt*sqrt(nt), df = nt - 1, ncp = x*sqrt(nt))
  
  d1 <- AbscontDistribution(d = like1, low1 = -15, up1 = 15, withStand = TRUE)
  d2 <- AbscontDistribution(d = like2, low1 = -15, up1 = 15, withStand = TRUE)
  
  dif <- distr::r(if(test) -(d2 - d1) else d2 - d1)(n.sim)
  
  din <- if(test) -din else din
  
  SD <- sd(dif)
  
  return(c(dint = din, SD = SD))
})
  
#===============================================================================================================================
    
dit5 <- Vectorize(function(dppc, dppt, nc, nt, n.sim = 1, rev.sign = FALSE){
  
  a <- dppc
  b <- dppt
  
  din <- b - a 
  
  test <- if(!rev.sign || rev.sign & b < 0 & a < 0 & abs(b) < abs(a)) FALSE else TRUE
  
  vc <- 1/nc + dppc^2/(2*nc)
  
  vt <- 1/nt + dppt^2/(2*nt)
  
  SD <- sqrt(vt + vc)
  
  din <- if(test) -din else din
  
  return(c(dint = din, SD = SD))
})       
    
#===============================================================================================================================    
    
dit <- Vectorize(function(dppc, dppt, nc, nt, n.sim = 1, rev.sign = FALSE){
  
  a <- dppc
  b <- dppt
  
  din <- b - a 
  
  test <- if(!rev.sign || rev.sign & b < 0 & a < 0 & abs(b) < abs(a)) FALSE else TRUE
  
  vc <- 1/nc + (1 - ( (nc-3) / ((nc-1)*cfactor(nc-1)^2) )) * dppc^2
  vt <- 1/nt + (1 - ( (nt-3) / ((nt-1)*cfactor(nt-1)^2) )) * dppt^2
  
  SD <- sqrt(vt + vc)
  
  din <- if(test) -din else din
  return(c(dint = din, SD = SD))
})            
    
#===============================================================================================================================
             
pairup <- function(rv){ rv[1:(length(rv)/2)] }             
             
#===============================================================================================================================

var.d <- function(d, n1, n2 = NA, g = FALSE, r = .37, cont.grp = FALSE){
  
  v <- if(is.na(n2) & !cont.grp) (1/n1) + ((d^2)/(2*n1)) 
  else if(is.na(n2) & cont.grp) ((2*(1-r))/n1) + ((d^2)/(2*n1)) 
  else ((n1+n2)/(n1*n2)) + ((d^2)/(2*(n1+n2)) )
  df <- ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)
  
  ifelse(g == TRUE, cfactor(df)^2 * v, v)
}
  
#===============================================================================================================================

se.d <- function(d, n1, n2 = NA, g = FALSE, r = .37, cont.grp = FALSE) sqrt( var.d(d, n1, n2 = n2, g = g, r = r, cont.grp = cont.grp) )
                
#===============================================================================================================================

t.testb <- function(m1, m2, s1, s2, n1, n2 = NA, m0 = 0, var.equal = FALSE, sdif = NA, r = NA, digits = 6){
  
  if(var.equal & !is.na(n2))
    {
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  } else if(!var.equal & !is.na(n2))
    {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    df <- ((s1^2/n1 + s2^2/n2)^2)/((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))
  }else
  {
    se <- if(!is.na(sdif)) sdif/sqrt(n1) else sdif(sdpre = s1, sdpos = s2, r = r)/sqrt(n1)
    df <- n1 - 1
  }
  
  t <- (m2-m1-m0)/se
  
  a <- round(data.frame(mean.dif = m2-m1, std.error = se, t.value = t, p.value = 2*pt(-abs(t),df)), digits)
  a$paired <- if(is.na(n2)) TRUE else FALSE
  a    
}        
  
#===============================================================================================================================
 
fill <- function(refdf, ...) 
  {
  L <- list(...)
  res <- lapply(L, function(x) {
    toadd <- setdiff(names(refdf), names(x))
    x[toadd] <- refdf[toadd]
    x
  })
  c(list(refdf), res)
}         
             
#===============================================================================================================================
             
opt1 <- function(sds, r = .5){ 
  
  V <- sds^2  
  
  m <- length(V)
  
  r <- cor.mat(r, m)
  
  sqrt((1/m)^2 * sum(sqrt(outer(V, V)) * r))
}

#===============================================================================================================================
 
hdir <- function(sample, level = .95, digits = 1e2){
  
  if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  sorted <- sort(sample)
  index <- ceiling(level*length(sorted))
  n <- length(sorted)- index
  width <- numeric(n)
  for(i in 1:n){
    width[i] <- sorted[i+ index]- sorted[i]
  }
  lower <- sorted[which.min(width)]
  upper <- sorted[which.min(width)+ index]
  return(round(c(lower, upper), digits))
}
             
#===============================================================================================================================
             
             
study.plot <- function(fit, adjust = 1, na.rm = TRUE, n = 1e4, hdi = TRUE, level = .95, xlab = "dint (short-long)", ylim = NA, xlim = NA, labels = NA, bottom = 1, top = 1, scale = 1){
  
  lab <- fit$labels
  
  L <- lapply(1:fit$k, function(x) fit$rposterior(1e4, tau.sample = FALSE, ind = x))
  
  loop <- length(L)
  soop <- seq_len(loop)
  
  a <- lapply(L, function(x) density(x, adjust = adjust, na.rm = na.rm, n = n))
  
  from <- numeric(loop)
  to <- numeric(loop)
  hi <- numeric(loop)
  if(hdi) CI <- matrix(NA, loop, 2)
  mode <- numeric(loop)
  
  for(i in soop){
    from[i] <- min(a[[i]]$x)
    to[i] <- max(a[[i]]$x)
    hi[i] <- max(a[[i]]$y)
    if(hdi) CI[i,] <- hdir(L[[i]], level = level)
    mode[i] <- a[[i]]$x[which.max(a[[i]]$y)]
  }
  
  f = hi + soop
  m = scale*hi + soop
  
  plot(rep(soop, 2), rep(soop, 2), type = "n", xlim = if(is.na(xlim)) c(min(from), max(to)) else xlim, ylim = if(is.na(ylim)) c(bottom*1, top*max(f)) else ylim, ylab = NA, yaxt = "n", xlab = xlab, mgp = c(2, .3, 0))
  axis(2, at = soop, labels = if(is.na(labels)) lab else labels, font = 2, las = 1, cex.axis = .8, tck = -.012, mgp = c(2, .3, 0), padj = rep(.35, loop))
  abline(h = soop, col = 8, lty = 3)
  
  for(i in soop){
    polygon(x = a[[i]]$x, y = scale*a[[i]]$y +i, col = adjustcolor(i, .4), border = NA, xpd = NA)
  }
  
  if(hdi){   
    segments(CI[, 1], soop, CI[, 2], soop, lend = 1, lwd = 4, col = soop, xpd = NA)                            
    segments(mode, soop, mode, m, lty = 3, xpd = NA, lend = 1)  
    points(mode, soop, pch = 22, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
    I = decimal(CI, 2); o = decimal(mode, 2)
    text(c(CI[,1], o, CI[,2]), soop, c(I[,1], o, I[,2]), pos = 3, font = 2, cex = .8, xpd = NA)
  }  
  return(invisible(a))
}             
             
#===============================================================================================================================
             
one.rm <- function(List){
  
  nms <- unlist(lapply(List, names))  
  keep <- nms[duplicated(nms)]
  res <- lapply(List, function(x) x[names(x) %in% keep])
  Filter(NROW, res)
}             

#===============================================================================================================================
       
comb.dif.mean <- function(List, na.rm = TRUE){

pt <- unlist(List[2:3])
List[[2]] <- tapply(pt, names(pt), FUN = mean, na.rm = na.rm)
List[[3]] <- NULL
List[[1]] - List[[2]]
}       

#===============================================================================================================================
       
comb.dif.sd <- function(List, r = .5){
  
  pt <- unlist(List[2:3])
  List[[2]] <- tapply(pt, names(pt), FUN = opt1, r = r)
  List[[3]] <- NULL
  sdif(sdpre = List[[1]], sdpos = List[[2]], r = r)
}       
       
#===============================================================================================================================              

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
  
 
#===============================================================================================================================

  
funnel.bayesmeta <- function(x,
                             main = deparse(substitute(x)),
                             xlab = "Effect Size (dint)", xlim = NULL,
                             ylab = "Standard Error (SE)", study.name = TRUE,
                             FE = FALSE, legend = FALSE, shrink = FALSE, show.mu = TRUE, pt.cex = 1, mu.cex = .8, mu.pos = 2, ...)
{
  
  
  stopifnot(is.element("bayesmeta", class(x)))
  
  yrange <- c(0.0, max(x$sigma))
  
  sevec <- seq(from=0, yrange[2]*1.04, le=27)
  
  intRE <- matrix(NA_real_, nrow=length(sevec), ncol=2,
                  dimnames=list(NULL, c("lower","upper")))
  intRE[1,] <- x$qposterior(theta.p=c(0.025, 0.975), predict=TRUE)
  for (i in 2:length(sevec)){
    conv <- try(convolve(dens1=function(a, log=FALSE){return(x$dposterior(theta=a, predict=TRUE, log=log))},
                         dens2=function(b, log=FALSE){return(dnorm(x=b, mean=0, sd=sevec[i], log=log))},
                         cdf1 =function(a){return(x$pposterior(theta=a, predict=TRUE))},
                         cdf2 =function(b){return(pnorm(q=b, mean=0, sd=sevec[i]))}))
    if (all(class(conv)!="try-error")) {
      intRE[i,] <- conv$quantile(p=c(0.025, 0.975))
    }
  }
  
  intFE <- matrix(NA_real_, nrow=length(sevec), ncol=2,
                  dimnames=list(NULL, c("lower","upper")))
  cm <- x$cond.moment(tau=0)
  for (i in 1:length(sevec)){
    intFE[i,] <- qnorm(c(0.025, 0.975), mean=cm[1,"mean"], sd=sqrt(cm[1,"sd"]^2+sevec[i]^2))
  }
  FEcol="red3"
  REcol="blue3"
  
  xl <- range(intRE)
  xlim <- if(is.null(xlim)) range(c(xl, range(x$y, finite = TRUE)), finite = TRUE) else xlim
  
  plot(xl, -yrange, type="n", xlim = xlim,
       ylab=ylab, xlab=xlab, main=main, axes = FALSE, ...)
  
  polygon(c(intRE[,1], rev(intRE[,2])), c(-sevec, rev(-sevec)), col="grey90", border=NA)
  if (FE) polygon(c(intFE[,1], rev(intFE[,2])), c(-sevec, rev(-sevec)), col="grey80", border=NA)
  
  lines(c(intRE[1,1], intRE[1,1], NA, intRE[1,2], intRE[1,2]),
        c(0,-max(sevec), NA, 0, -max(sevec)), col="grey75", lty="dashed")
  abline(h=0, col="darkgrey")
  yticks <- pretty(yrange)
  abline(h=-yticks[yticks>0], col="grey75", lty="15")
  
  matlines(intRE, cbind(-sevec, -sevec), col=REcol, lty="dashed")
  if (FE) matlines(intFE, cbind(-sevec, -sevec), col=FEcol, lty="dotted")
  lines(rep(x$summary["median","mu"], 2), range(-sevec), col= REcol, lty= 2)
  if (FE) lines(rep(cm[1,"mean"], 2), range(-sevec), col=FEcol, lty="dotted")
  
  if(show.mu) text(x$summary["median","mu"], mean(par('usr')[3:4])*.1, bquote(mu == .(round(x$summary["median","mu"], 3))), font = 2, col = REcol, srt = 90, pos = mu.pos, cex = mu.cex)
  
  lines(c(0, 0), c(-1,1)*max(sevec), col = "darkgrey")
  
  points(x$y, -x$sigma, pch=21, col="magenta", bg="cyan", cex= pt.cex)
  
  if(shrink) points(x$theta[5,], -x$sigma, pch=21, col=adjustcolor("gray40", .5), bg= adjustcolor("gray40", .5), cex=pt.cex)
  #if(shrink) {
  # segments(x$y, -x$sigma, x$theta[5,], -x$sigma, col="gray40", lty = 3)
  # points(x$theta[5,], -x$sigma, pch=21, col="gray60", bg= "gray60", cex=1.2)
  #}
  
  if(study.name)text(x$y, -x$sigma, x$labels, cex = .65, font = 2, pos = 3, xpd = NA)
  
  if (FE && legend)
    legend("topleft", c("RE model", "FE model"),
           col=c(REcol, FEcol), lty=c("dashed", "dotted"), bg="white")
  axis(1, ...) ; axis(2, at=-yticks, labels=yticks, ...); box()
  invisible()
}
  
  
  
#===============================================================================================================================
              
              
d.prepos1 <- function(d = NA, study.name = NA, group.name = NA, n = NA, mdif = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, rev.sign = FALSE, autoreg = FALSE, t.pair = NA, df = NA, sdif = NA, post, control, outcome, ...) 
{
  
  if(missing(control) || missing(post) || missing(outcome)) stop("'post', 'outcome' and/or 'control' missing in the EXCEL sheet.", call. = FALSE)  
  
  r <- ifelse(autoreg == TRUE & !is.na(r), autoreg(max(post, na.rm = TRUE), r)[,1][-1][post], r)
  
  n <- ifelse(!is.na(n), n, ifelse(is.na(n) & !is.na(df), df + 1, NA))
  d <- ifelse(!is.na(d), d, ifelse(!is.na(t.pair) & !is.na(n), t2d(t.pair, n), NA))
  mdif <- ifelse(!is.na(mdif), mdif, ifelse(!is.na(mpre) & !is.na(mpre) & is.na(mdif), mpos - mpre, NA))
  sdif <- ifelse(is.na(sdif), sdif(sdpre = sdpre, sdpos = sdpos, t = t.pair, r = r, n = n, mpos = mpos, mpre = mpre), sdif)
  cor. <- ifelse(is.na(r), rdif(n = n, mpre = mpre, mpos = mpos, sdpre = sdpre, sdpos = sdpos, sdif = sdif), r)
  d <- ifelse(!is.na(mdif) & is.na(d) & !is.na(sdif), mdif/sdif, d)
  d <- ifelse(rev.sign == TRUE, -d, d)*cfactor(n-1)
  #se <- se.d(d, n1 = n, g = TRUE)

  out <- data.frame(d = d, n = n, sdif = sdif, rpr.po = cor., post, control, outcome, ...)
  
  #if(!anyNA(group.name) & length(group.name) == nrow(out)) row.names(out) <- as.character(group.name) else if(!anyNA(group.name) & length(group.name) != nrow(out)) stop("'group.name' incorrectly specified.", call. = FALSE)
  
  if(all(is.na(out$d))) stop("\ninsufficient info. to calculate effect size(s).", call. = FALSE)
  
  return(out) 
}


#================================================================================================================================
             
             
d.prepos2 <- function(d = NA, study.name = NA, group.name = NA, n = NA, mdif = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, rev.sign = FALSE, autoreg = FALSE, t.pair = NA, df = NA, sdif = NA, post, control, outcome, ...) 
{
  
  if(missing(control) || missing(post) || missing(outcome)) stop("'post', 'outcome' and/or 'control' missing in the EXCEL sheet.", call. = FALSE)  
 
  r <- ifelse(autoreg == TRUE & !is.na(r), autoreg(max(post, na.rm = TRUE), r)[,1][-1][post], r)
  
  n <- ifelse(!is.na(n), n, ifelse(is.na(n) & !is.na(df), df + 1, NA))
  d <- ifelse(!is.na(d), d, ifelse(!is.na(t.pair) & !is.na(n), t2d(t.pair, n), NA))
  mdif <- ifelse(!is.na(mdif), mdif, ifelse(!is.na(mpre) & !is.na(mpre) & is.na(mdif), mpos - mpre, NA))
  sdif <- ifelse(is.na(sdif), sdif(sdpre = sdpre, sdpos = sdpos, t = t.pair, r = r, n = n, mpos = mpos, mpre = mpre), sdif)
  r <- ifelse(is.na(r), rdif(n = n, mpre = mpre, mpos = mpos, sdpre = sdpre, sdpos = sdpos, sdif = sdif), r)
  d <- ifelse(!is.na(mdif) & is.na(d) & !is.na(sdif), mdif/sdif, d)*cfactor(n-1)
  if(anyNA(d) & anyNA(r)) stop("'r' must be defined. If none available, we suggest '.6'.", call. = FALSE)
  
  out <- data.frame(d, n, sdif, r, rev.sign, post, control, outcome, ...)
  
  if(all(is.na(out$d))) stop("\ninsufficient info. to calculate effect size(s).", call. = FALSE)
  
  return(out) 
}             
             
#================================================================================================================================
             
 
d.prepos3 <- function(d = NA, study.name = NA, group.name = NA, n = NA, mdif = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, rev.sign = FALSE, rev.group = FALSE, autoreg = FALSE, t.pair = NA, df = NA, sdif = NA, post, control, outcome, ...) 
{
  
  if(missing(control) || missing(post) || missing(outcome)) stop("'post', 'outcome' and/or 'control' missing in the EXCEL sheet.", call. = FALSE)  
  
  r <- ifelse(autoreg == TRUE & !is.na(r), autoreg(max(post, na.rm = TRUE), r)[,1][-1][post], r)
  
  n <- ifelse(!is.na(n), n, ifelse(is.na(n) & !is.na(df), df + 1, NA))
  d <- ifelse(!is.na(d), d, ifelse(!is.na(t.pair) & !is.na(n), t2d(t.pair, n), NA))
  mdif <- ifelse(!is.na(mdif), mdif, ifelse(!is.na(mpre) & !is.na(mpre) & is.na(mdif), mpos - mpre, NA))
  sdif <- ifelse(is.na(sdif), sdif(sdpre = sdpre, sdpos = sdpos, t.pair = t.pair, r = r, n = n, mpos = mpos, mpre = mpre), sdif)
  r <- ifelse(is.na(r), rdif(n = n, mpre = mpre, mpos = mpos, t.pair = t.pair, sdpre = sdpre, sdpos = sdpos, sdif = sdif), r)
  d <- ifelse(!is.na(mdif) & is.na(d) & !is.na(sdif), mdif/sdif, d)*cfactor(n-1)
  if(anyNA(d) & anyNA(r)) stop("'r' must be defined. If none available, we suggest '.6'.", call. = FALSE)
  d <- ifelse(rev.group, -d, d)
  
  out <- data.frame(d, n, sdif, r, rev.sign, post, control, outcome, ...)
  
  if(all(is.na(out$d))) stop("\ninsufficient info. to calculate effect size(s).", call. = FALSE)
  
  return(out) 
}                          
             

#===============================================================================================================================
             
             
d.prepos4 <- function(d = NA, study.name = NA, group.name = NA, n = NA, mdif = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, rev.sign = FALSE, rev.group = FALSE, autoreg = FALSE, t.pair = NA, df = NA, sdif = NA, post, control, outcome, ...) 
{
  
  if(missing(control) || missing(post) || missing(outcome)) stop("'post', 'outcome' and/or 'control' missing in the EXCEL sheet.", call. = FALSE)
  
  r <- ifelse(autoreg == TRUE, autoreg(max(post, na.rm = TRUE), r)[,1][-1][post], r)
  
  n <- ifelse(!is.na(n), n, ifelse(is.na(n) & !is.na(df), df + 1, NA))
  d <- ifelse(!is.na(d), d, ifelse(!is.na(t.pair) & !is.na(n), t2d(t.pair, n), NA))
  mdif <- ifelse(!is.na(mdif), mdif, ifelse(!is.na(mpre) & !is.na(mpre) & is.na(mdif), mpos - mpre, NA))
  sdif <- ifelse(is.na(sdif), sdif(sdpre = sdpre, sdpos = sdpos, t.pair = t.pair, r = r, n = n, mpos = mpos, mpre = mpre), sdif)
  r <- ifelse(is.na(r), rdif(n = n, mpre = mpre, mpos = mpos, t.pair = t.pair, sdpre = sdpre, sdpos = sdpos, sdif = sdif), r)
  r <- ifelse(is.na(sdif) & is.na(r), .6, r)
  sdif <- sdif(sdpre = sdpre, sdpos = sdpos, t.pair = t.pair, r = r, n = n, mpos = mpos, mpre = mpre)
  d <- ifelse(!is.na(mdif) & is.na(d) & !is.na(sdif), mdif/sdif, d)*cfactor(n-1)
  d <- ifelse(rev.group, -d, d)
  
  out <- data.frame(d, n, sdif, r, rev.sign, post, control, outcome, ...)
  
  if(all(is.na(out$d))) stop("\ninsufficient info. to calculate effect size(s).", call. = FALSE)
  
  return(out) 
}                     
            

#================================================================================================================================       
       
       
d.prepos <- function(d = NA, study.name = NA, group.name = NA, n = NA, mdif = NA, stder = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, rev.sign = FALSE, rev.group = FALSE, autoreg = FALSE, t.pair = NA, df = NA, sdif = NA, post = NA, control = NA, outcome = NA, time = NA, ...) 
{
  
  if(anyNA(control) || anyNA(post) || anyNA(outcome) || anyNA(time)) stop("'post', 'outcome', 'time', or 'control' missing in the EXCEL sheet.", call. = FALSE)
  
  rev.sign <- ifelse(is.na(rev.sign), FALSE, rev.sign)
  rev.group <- ifelse(is.na(rev.group), FALSE, rev.group)
  autoreg <- ifelse(is.na(autoreg), FALSE, autoreg)
  control <- ifelse(is.na(control), FALSE, control)
  
  r <- ifelse(autoreg == TRUE, autoreg(max(post, na.rm = TRUE), r)[,1][-1][post], r)
  
  n <- ifelse(!is.na(n), n, ifelse(is.na(n) & !is.na(df), df + 1, NA))
  mdif <- ifelse(!is.na(mdif), mdif, ifelse(!is.na(mpre) & !is.na(mpre) & is.na(mdif), mpos - mpre, NA))
  t.pair <- ifelse(!is.na(t.pair), t.pair, ifelse(is.na(t.pair) & !is.na(mdif) & !is.na(stder), mdif/stder, NA))
  d <- ifelse(!is.na(d), d, ifelse(!is.na(t.pair) & !is.na(n), t2d(t.pair, n), ifelse(!is.na(mdif) & !is.na(sdif), mdif/sdif, NA)))
  sdif <- ifelse(is.na(sdif), sdif(sdpre = sdpre, sdpos = sdpos, t.pair = t.pair, r = r, n = n, mpos = mpos, mpre = mpre), sdif)
  r <- ifelse(is.na(r), rdif(n = n, mpre = mpre, mpos = mpos, t.pair = t.pair, sdpre = sdpre, sdpos = sdpos, sdif = sdif), r)
  r <- ifelse(is.na(r) & is.na(sdif), .6, r)
  sdif <- sdif(sdpre = sdpre, sdpos = sdpos, t.pair = t.pair, r = r, n = n, mpos = mpos, mpre = mpre)
  d <- ifelse(!is.na(mdif) & is.na(d) & !is.na(sdif), mdif/sdif, d)*cfactor(n-1)
  d <- ifelse(rev.group, -d, d)
  
  out <- data.frame(d, n, sdif, r, rev.sign, post, control, outcome, time, ...)
  
  if(all(is.na(out$d))) stop("insufficient info. to calculate effect size(s).", call. = FALSE)
  
  return(out) 
}

       
       
#================================================================================================================================


dint1 <- function(..., per.study = NULL, study.name = NA, group.name = NA, n.sim = 1e5, by, data = NULL)
{
  
  L <- if(!is.null(data)){
    
    m <- split(data, data$study.name)        
    
    m[[1]] <- NULL                          
    
    if(is.null(reget(m, control))) stop("Required 'control' group not found.", call. = FALSE)
    
    ar <- head(formalArgs(d.prepos), -1)
    
    dot.names <- names(m[[1]])[!names(m[[1]]) %in% ar]
    
    args <- lapply(m, function(x) unclass(x[c(ar, dot.names)]))
    
    argsT <- setNames(lapply(names(args[[1]]), 
             function(i) lapply(args, `[[`, i)), names(args[[1]]))
    
    do.call(Map, c(f = d.prepos, argsT))
    
  } else { fuse(... = ..., per.study = per.study) }
  
  
  if(!missing(by)){
    
    s <- substitute(by)    
    
    k <- as.list(s)
    
    if("control" %in% k || "!control" %in% k) stop("'control' can't be a moderating variable either alone or with other variables.", call. = FALSE)
    
    H <- lapply(L, function(x) do.call("subset", list(x, s)))
    
    res <- Filter(NROW, H)
    
    L <- if(length(res) == 0) stop("No study with the requested moderators found.", call. = FALSE) else res
  }
  
  G <- function(m, n.sim)
  {
    
    cdel1 <- reget(m, control & post == 2 & outcome == 1)
    cdel2 <- reget(m, control & post == 3 & outcome == 1)
    cs <- reget(m, control & post == 1 & outcome == 1)
    
    tdel1 <- reget(m, !control & post == 2 & outcome == 1)
    tdel2 <- reget(m, !control & post == 3 & outcome == 1)
    ts <- reget(m, !control & post == 1 & outcome == 1) 
    
    if(all(sapply(list(cdel1, cdel2, tdel1, tdel2, ts, cs), is.null))) stop("Either 'control' or 'post' incorrectly coded.", call. = FALSE)
    
    short <- all(sapply(list(cs, ts), function(x) !is.null(x)))
    
    del1 <- all(sapply(list(cdel1, tdel1), function(x) !is.null(x)))
    
    del2 <- all(sapply(list(cdel2, tdel2), function(x) !is.null(x)))
    
    
    cdel1..2 <- reget(m, control & post == 2 & outcome == 2)
    cdel2..2 <- reget(m, control & post == 3 & outcome == 2)
    cs..2 <- reget(m, control & post == 1 & outcome == 2)
    
    tdel1..2 <- reget(m, !control & post == 2 & outcome == 2)
    tdel2..2 <- reget(m, !control & post == 3 & outcome == 2)
    ts..2 <- reget(m, !control & post == 1 & outcome == 2)


    short..2 <- all(sapply(list(cs..2, ts..2), function(x) !is.null(x)))
    
    del1..2 <- all(sapply(list(cdel1..2, tdel1..2), function(x) !is.null(x)))
    
    del2..2 <- all(sapply(list(cdel2..2, tdel2..2), function(x) !is.null(x)))
    
    
    if(short){
      nc1 <- m$n[m$control & m$post == 1 & m$outcome == 1]
      nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 1]
      dps <- pair(cs, ts)  
      dppc1 <- dppcs <- sapply(1:length(dps), function(i) dps[[i]][[1]][1])
      dppt1 <- dppts <- sapply(1:length(dps), function(i) dps[[i]][[1]][2])
     # group.name1 <- unlist(lapply(1:length(dps), function(i) names(dps[[i]])))
      SHORT <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim)))
     # row.names(SHORT) <- group.name1
    }
    
    
    if(short..2){
      nc1 <- m$n[m$control & m$post == 1 & m$outcome == 2]
      nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 2]
      dps <- pair(cs..2, ts..2)  
      dppc1 <- sapply(1:length(dps), function(i) dps[[i]][[1]][1])
      dppt1 <- sapply(1:length(dps), function(i) dps[[i]][[1]][2])
      #group.name1 <- unlist(lapply(1:length(dps), function(i) names(dps[[i]])))
      SHORT..2 <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim)))
     # row.names(SHORT..2) <- group.name1
    }
    
    
    if(del1){
      nc2 <- m$n[m$control & m$post == 2 & m$outcome == 1]
      nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 1]
      dpdel1 <- pair(cdel1, tdel1)
      dppc2 <- dppcdel1 <- sapply(1:length(dpdel1), function(i) dpdel1[[i]][[1]][1])
      dppt2 <- dpptdel1 <- sapply(1:length(dpdel1), function(i) dpdel1[[i]][[1]][2])
      #group.name2 <- unlist(lapply(1:length(dpdel1), function(i) names(dpdel1[[i]])))
      DEL1 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim)))
      #row.names(DEL1) <- group.name2
    }
    
    
    if(del1..2){
      nc2 <- m$n[m$control & m$post == 2 & m$outcome == 2]
      nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 2]
      dpdel1 <- pair(cdel1..2, tdel1..2)
      dppc2 <- sapply(1:length(dpdel1), function(i) dpdel1[[i]][[1]][1])
      dppt2 <- sapply(1:length(dpdel1), function(i) dpdel1[[i]][[1]][2])
     # group.name2 <- unlist(lapply(1:length(dpdel1), function(i) names(dpdel1[[i]])))
      DEL1..2 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim)))
      #row.names(DEL1..2) <- group.name2
    }
    
    
    if(del2){
      nc3 <- m$n[m$control & m$post == 3 & m$outcome == 1]
      nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 1]
      dpdel2 <- pair(cdel2, tdel2)
      dppc3 <- dppcdel2 <- sapply(1:length(dpdel2), function(i) dpdel2[[i]][[1]][1])
      dppt3 <- dpptdel2 <- sapply(1:length(dpdel2), function(i) dpdel2[[i]][[1]][2])
    #  group.name3 <- unlist(lapply(1:length(dpdel2), function(i) names(dpdel2[[i]])))
      DEL2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim)))
      #row.names(DEL2) <- group.name3
    }
    
    if(del2..2){
      nc3 <- m$n[m$control & m$post == 3 & m$outcome == 2]
      nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 2]
      dpdel2 <- pair(cdel2..2, tdel2..2)
      dppc3 <- sapply(1:length(dpdel2), function(i) dpdel2[[i]][[1]][1])
      dppt3 <- sapply(1:length(dpdel2), function(i) dpdel2[[i]][[1]][2])
      #group.name3 <- unlist(lapply(1:length(dpdel2), function(i) names(dpdel2[[i]])))
      DEL2..2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim)))
      #row.names(DEL2..2) <- group.name3
    }
    
    
    list(SHORT = if(short) SHORT else NULL, SHORT..2 = if(short..2) SHORT..2 else NULL, DEL1 = if(del1) DEL1 else NULL, DEL1..2 = if(del1..2) DEL1..2 else NULL, DEL2 = if(del2) DEL2 else NULL, DEL2..2 = if(del2..2) DEL2..2 else NULL) 
  }
  
  h <- lapply(1:length(L), function(i) G(m = L[[i]], n.sim = n.sim))
  
  if(!is.null(data)) study.name <- names(L)
  
  names(h) <- if(anyNA(study.name)) paste0("Study", seq_along(h)) else if(!anyNA(study.name) & length(study.name) == length(h)) study.name else if(!anyNA(study.name) & length(study.name) != length(h)) stop("'study.name' incorrectly specified.", call. = FALSE)
  
  return(h)
  
}
              
#================================================================================================================================
              
              
dint2 <- function(..., per.study = NULL, study.name = NA, n.sim = 1e5, by, data = NULL)
{
  
  L <- if(!is.null(data)){
    
    m <- split(data, data$study.name)        
    
    m[[1]] <- NULL                          
    
    if(is.null(reget(m, control))) stop("Required 'control' group not found.", call. = FALSE)
    
    ar <- head(formalArgs(d.prepos), -1)
    
    dot.names <- names(m[[1]])[!names(m[[1]]) %in% ar]
    
    args <- lapply(m, function(x) unclass(x[c(ar, dot.names)]))
    
    argsT <- setNames(lapply(names(args[[1]]), function(i) lapply(args, `[[`, i)), names(args[[1]]))
    
    do.call(Map, c(f = d.prepos, argsT))
    
  } else { fuse(... = ..., per.study = per.study) }
  
  
  if(!missing(by)){ 
    
    s <- substitute(by)
    k <- as.list(s)
    
    if("control" %in% k || "!control" %in% k) stop("'control' can't be a moderating variable either alone or with other variables.", call. = FALSE)
    
    H <- lapply(L, function(x) do.call("subset", list(x, s)))
    
    H <- Filter(NROW, H)
    h <- if(length(H) == 0) stop("No study with the requested moderators found.", call. = FALSE) else H
    
    g <- lapply(L[names(h)], function(x) subset(x, control))
    L <- Map(rbind, h, g)   
  }
  
  G <- function(m, n.sim)
  {
    
    cdel1 <- reget(m, control & post == 2 & outcome == 1)
    cdel2 <- reget(m, control & post == 3 & outcome == 1)
    cs <- reget(m, control & post == 1 & outcome == 1)
    
    tdel1 <- reget(m, !control & post == 2 & outcome == 1)
    tdel2 <- reget(m, !control & post == 3 & outcome == 1)
    ts <- reget(m, !control & post == 1 & outcome == 1) 
    
    if(all(sapply(list(cdel1, cdel2, tdel1, tdel2, ts, cs), is.null))) stop("Either 'control' or 'post' incorrectly coded.", call. = FALSE)
    
    short <- all(sapply(list(cs, ts), function(x) !is.null(x)))
    
    del1 <- all(sapply(list(cdel1, tdel1), function(x) !is.null(x)))
    
    del2 <- all(sapply(list(cdel2, tdel2), function(x) !is.null(x)))
    
    
    cdel1..2 <- reget(m, control & post == 2 & outcome == 2)
    cdel2..2 <- reget(m, control & post == 3 & outcome == 2)
    cs..2 <- reget(m, control & post == 1 & outcome == 2)
    
    tdel1..2 <- reget(m, !control & post == 2 & outcome == 2)
    tdel2..2 <- reget(m, !control & post == 3 & outcome == 2)
    ts..2 <- reget(m, !control & post == 1 & outcome == 2)
    
    
    short..2 <- all(sapply(list(cs..2, ts..2), function(x) !is.null(x)))
    
    del1..2 <- all(sapply(list(cdel1..2, tdel1..2), function(x) !is.null(x)))
    
    del2..2 <- all(sapply(list(cdel2..2, tdel2..2), function(x) !is.null(x)))
    
    
    if(short){
      nc1 <- m$n[m$control & m$post == 1 & m$outcome == 1]
      nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 1]
      dps <- pair(cs, ts)  
      dppc1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][1])
      dppt1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][2])
      # group.name1 <- unlist(lapply(1:length(dps), function(i) names(dps[[i]])))
      SHORT <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim)))
      # row.names(SHORT) <- group.name1
    }
    
    
    if(short..2){
      nc1 <- m$n[m$control & m$post == 1 & m$outcome == 2]
      nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 2]
      dps <- pair(cs..2, ts..2)  
      dppc1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][1])
      dppt1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][2])
      #group.name1 <- unlist(lapply(1:length(dps), function(i) names(dps[[i]])))
      SHORT..2 <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim)))
      # row.names(SHORT..2) <- group.name1
    }
    
    
    if(del1){
      nc2 <- m$n[m$control & m$post == 2 & m$outcome == 1]
      nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 1]
      dpdel1 <- pair(cdel1, tdel1)
      dppc2 <- dppcdel1 <- sapply(1:length(dpdel1), function(i) dpdel1[[i]][[1]][1])
      dppt2 <- dpptdel1 <- sapply(1:length(dpdel1), function(i) dpdel1[[i]][[1]][2])
      #group.name2 <- unlist(lapply(1:length(dpdel1), function(i) names(dpdel1[[i]])))
      DEL1 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim)))
      #row.names(DEL1) <- group.name2
    }
    
    
    if(del1..2){
      nc2 <- m$n[m$control & m$post == 2 & m$outcome == 2]
      nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 2]
      dpdel1 <- pair(cdel1..2, tdel1..2)
      dppc2 <- sapply(1:length(dpdel1), function(i) dpdel1[[i]][[1]][1])
      dppt2 <- sapply(1:length(dpdel1), function(i) dpdel1[[i]][[1]][2])
      # group.name2 <- unlist(lapply(1:length(dpdel1), function(i) names(dpdel1[[i]])))
      DEL1..2 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim)))
      #row.names(DEL1..2) <- group.name2
    }
    
    
    if(del2){
      nc3 <- m$n[m$control & m$post == 3 & m$outcome == 1]
      nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 1]
      dpdel2 <- pair(cdel2, tdel2)
      dppc3 <- dppcdel2 <- sapply(1:length(dpdel2), function(i) dpdel2[[i]][[1]][1])
      dppt3 <- dpptdel2 <- sapply(1:length(dpdel2), function(i) dpdel2[[i]][[1]][2])
      #  group.name3 <- unlist(lapply(1:length(dpdel2), function(i) names(dpdel2[[i]])))
      DEL2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim)))
      #row.names(DEL2) <- group.name3
    }
    
    if(del2..2){
      nc3 <- m$n[m$control & m$post == 3 & m$outcome == 2]
      nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 2]
      dpdel2 <- pair(cdel2..2, tdel2..2)
      dppc3 <- sapply(1:length(dpdel2), function(i) dpdel2[[i]][[1]][1])
      dppt3 <- sapply(1:length(dpdel2), function(i) dpdel2[[i]][[1]][2])
      #group.name3 <- unlist(lapply(1:length(dpdel2), function(i) names(dpdel2[[i]])))
      DEL2..2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim)))
      #row.names(DEL2..2) <- group.name3
    }
    
    list(SHORT = if(short) SHORT else NULL, SHORT..2 = if(short..2) SHORT..2 else NULL, DEL1 = if(del1) DEL1 else NULL, DEL1..2 = if(del1..2) DEL1..2 else NULL, DEL2 = if(del2) DEL2 else NULL, DEL2..2 = if(del2..2) DEL2..2 else NULL) 
  }
  
  h <- lapply(1:length(L), function(i) G(m = L[[i]], n.sim = n.sim))
  
  if(!is.null(data)) study.name <- names(L)
  
  names(h) <- if(anyNA(study.name)) paste0("Study", seq_along(h)) else if(!anyNA(study.name) & length(study.name) == length(h)) study.name else if(!anyNA(study.name) & length(study.name) != length(h)) stop("'study.name' incorrectly specified.", call. = FALSE)
  
  return(h) 
}              
              
#================================================================================================================================
              
              
dinti <- function(data = NULL, by, impute = FALSE, n.sim = 1e5)
{   
    m <- split(data, data$study.name)        
    
    m[[1]] <- NULL                          
    
    if(is.null(reget(m, control))) stop("Required 'control' group not found.", call. = FALSE)
    
 if(impute) {  
    ar <- formalArgs(rdif)[c(-7, -9)]

   args <- lapply(m, function(x) unclass(x[ar]))

   argsT <- setNames(lapply(names(args[[1]]), function(i) lapply(args, `[[`, i)), names(args[[1]]))

   f <- do.call(Map, c(f = rdif, argsT))

   f <- lapply(f, na.locf0)
                             
   m <- Map(function(x, y) transform(x, r = na.locf0(y, fromLast = TRUE)), m, f) 
   }   
      
    ar <- head(formalArgs(d.prepos), -1)
    
    dot.names <- names(m[[1]])[!names(m[[1]]) %in% ar]
    
    args <- lapply(m, function(x) unclass(x[c(ar, dot.names)]))
    
    argsT <- setNames(lapply(names(args[[1]]), function(i) lapply(args, `[[`, i)), names(args[[1]]))
    
  L <- do.call(Map, c(f = d.prepos, argsT))
    
  
  
  if(!missing(by)){ 
    
    s <- substitute(by)
    k <- as.list(s)
    
    if("control" %in% k || "!control" %in% k) stop("'control' can't be a moderating variable either alone or with other variables.", call. = FALSE)
    
    H <- lapply(L, function(x) do.call("subset", list(x, s)))
    
    H <- Filter(NROW, H)
    h <- if(length(H) == 0) stop("No study with the requested moderators found.", call. = FALSE) else H
    
    g <- lapply(L[names(h)], function(x) subset(x, control))
    L <- Map(rbind, h, g)   
  }
  
  G <- function(m, n.sim)
  {
    
    cdel1 <- reget(m, control & post == 2 & outcome == 1)
    cdel2 <- reget(m, control & post == 3 & outcome == 1)
    cs <- reget(m, control & post == 1 & outcome == 1)
    
    tdel1 <- reget(m, !control & post == 2 & outcome == 1)
    tdel2 <- reget(m, !control & post == 3 & outcome == 1)
    ts <- reget(m, !control & post == 1 & outcome == 1) 
    
    if(all(sapply(list(cdel1, cdel2, tdel1, tdel2, ts, cs), is.null))) stop("Either 'control' or 'post' incorrectly coded.", call. = FALSE)
    
    short <- all(sapply(list(cs, ts), function(x) !is.null(x)))
    
    del1 <- all(sapply(list(cdel1, tdel1), function(x) !is.null(x)))
    
    del2 <- all(sapply(list(cdel2, tdel2), function(x) !is.null(x)))
    
    
    cdel1..2 <- reget(m, control & post == 2 & outcome == 2)
    cdel2..2 <- reget(m, control & post == 3 & outcome == 2)
    cs..2 <- reget(m, control & post == 1 & outcome == 2)
    
    tdel1..2 <- reget(m, !control & post == 2 & outcome == 2)
    tdel2..2 <- reget(m, !control & post == 3 & outcome == 2)
    ts..2 <- reget(m, !control & post == 1 & outcome == 2)
    
    
    short..2 <- all(sapply(list(cs..2, ts..2), function(x) !is.null(x)))
    
    del1..2 <- all(sapply(list(cdel1..2, tdel1..2), function(x) !is.null(x)))
    
    del2..2 <- all(sapply(list(cdel2..2, tdel2..2), function(x) !is.null(x)))
    
    
    if(short){
      nc1 <- m$n[m$control & m$post == 1 & m$outcome == 1]
      nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 1]
      dps <- pair(cs, ts)  
      dppc1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][1])
      dppt1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][2])
      rv <- m$rev.sign[m$post == 1 & m$outcome == 1]
      
      SHORT <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim, rev.sign = pairup(rv))))
      
    }
    
    
    if(short..2){
      nc1 <- m$n[m$control & m$post == 1 & m$outcome == 2]
      nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 2]
      dps <- pair(cs..2, ts..2)  
      dppc1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][1])
      dppt1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][2])
      rv <- m$rev.sign[m$post == 1 & m$outcome == 2]
      
      SHORT..2 <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim, rev.sign = pairup(rv))))
    }
    
    
    if(del1){
      nc2 <- m$n[m$control & m$post == 2 & m$outcome == 1]
      nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 1]
      dpdel1 <- pair(cdel1, tdel1)
      dppc2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][1])
      dppt2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][2])
      rv <- m$rev.sign[m$post == 2 & m$outcome == 1]
      
      DEL1 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim, rev.sign = pairup(rv))))
    }
    
    
    if(del1..2){
      nc2 <- m$n[m$control & m$post == 2 & m$outcome == 2]
      nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 2]
      dpdel1 <- pair(cdel1..2, tdel1..2)
      dppc2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][1])
      dppt2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][2])
      rv <- m$rev.sign[m$post == 2 & m$outcome == 2]
      
      DEL1..2 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim, rev.sign = pairup(rv))))
    }
    
    
    if(del2){
      nc3 <- m$n[m$control & m$post == 3 & m$outcome == 1]
      nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 1]
      dpdel2 <- pair(cdel2, tdel2)
      dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
      dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
      rv <- m$rev.sign[m$post == 3 & m$outcome == 1]
      
      DEL2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
    }
    
    if(del2..2){
      nc3 <- m$n[m$control & m$post == 3 & m$outcome == 2]
      nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 2]
      dpdel2 <- pair(cdel2..2, tdel2..2)
      dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
      dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
      rv <- m$rev.sign[m$post == 3 & m$outcome == 2]
      
      DEL2..2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
    }
    
    list(SHORT = if(short) SHORT else NULL, SHORT..2 = if(short..2) SHORT..2 else NULL, DEL1 = if(del1) DEL1 else NULL, DEL1..2 = if(del1..2) DEL1..2 else NULL, DEL2 = if(del2) DEL2 else NULL, DEL2..2 = if(del2..2) DEL2..2 else NULL) 
  }
  
  h <- lapply(1:length(L), function(i) G(m = L[[i]], n.sim = n.sim))
  
  if(!is.null(data)) study.name <- names(L)
  
  names(h) <- if(anyNA(study.name)) paste0("Study", seq_along(h)) else if(!anyNA(study.name) & length(study.name) == length(h)) study.name else if(!anyNA(study.name) & length(study.name) != length(h)) stop("'study.name' incorrectly specified.", call. = FALSE)
  
  return(h) 
}               
              
              
#================================================================================================================================
              
              
              
meta.withini <- function(data = NULL, by, tau.prior = function(x){dhalfnormal(x)}, impute = FALSE, n.sim = 1e5){
  
  L <- eval(substitute(dinti(data = data, by = by, impute = impute, n.sim = n.sim)))
  
  study.name <- names(L)
  
  G <- function(object, tau.prior)
  {
    
    d1 <- object$SHORT$dint
    d1..2 <- object$SHORT..2$dint
    sd1 <- object$SHORT$SD
    sd1..2 <- object$SHORT..2$SD
    
    
    d2 <- object$DEL1$dint
    d2..2 <- object$DEL1..2$dint
    sd2 <- object$DEL1$SD
    sd2..2 <- object$DEL1..2$SD
    
    
    d3 <- object$DEL2$dint
    d3..2 <- object$DEL2..2$dint
    sd3 <- object$DEL2$SD
    sd3..2 <- object$DEL2..2$SD
    
    Short <- all(sapply(list(d1), function(x) !is.null(x)))
    
    Del1 <- all(sapply(list(d2), function(x) !is.null(x)))
    
    Del2 <- all(sapply(list(d3), function(x) !is.null(x)))
    
    Short..2 <- all(sapply(list(d1..2), function(x) !is.null(x)))
    
    Del1..2 <- all(sapply(list(d2..2), function(x) !is.null(x)))
    
    Del2..2 <- all(sapply(list(d3..2), function(x) !is.null(x)))
    
    
    if(Short & length(d1) == 1 & !Short..2) { 
      
      short <- c(Mean.dint.short = d1, SD.dint.short = sd1)
    }
    
    
    if(Short & length(d1) == 1 & Short..2 & length(d1..2) == 1) { 
      
      
     res1 <- bayesmeta(              y = c(d1, d1..2),
                                 sigma = c(sd1, sd1..2),
                                labels = NULL, tau.prior = tau.prior)
   res1$call <- match.call(expand.dots = FALSE)
   short <- c(Mean.dint.short = res1$summary["mean","mu"], SD.dint.short = res1$summary["sd","mu"])
      
    }
    
    
    if(Short & length(d1) > 1){
      result1 <- bayesmeta(                y = d1,
                                       sigma = sd1,
                                      labels = NULL, tau.prior = tau.prior)
      result1$call <- match.call(expand.dots = FALSE)
      
      short <- c(result1$summary["mean","mu"], result1$summary["sd","mu"])
      
    }
    
    
    if(Short & length(d1) > 1 & Short..2 & length(d1..2) > 1) { 
      
      
      res1 <- bayesmeta(                y = d1..2,
                                    sigma = sd1..2,
                                   labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      
      ds <- c(short[1], res1$summary["mean","mu"])
      sds <- c(short[2], res1$summary["sd","mu"])
      
      
      resi1 <- bayesmeta(                y = ds,
                                     sigma = sds,
                                    labels = NULL, tau.prior = tau.prior)
      resi1$call <- match.call(expand.dots = FALSE)
      
      short <- c(Mean.dint.short = resi1$summary["mean","mu"], SD.dint.short = resi1$summary["sd","mu"])
      
    }
    
    #####
    
    
    if(Del1 & length(d2) == 1 & !Del1..2) { 
      
      del1 <- c(Mean.dint.del1 = d2, SD.dint.del1 = sd2)
      
    }
    
    
    if(Del1 & length(d2) == 1 & Del1..2 & length(d2..2) == 1) { 
      
      res2 <- bayesmeta(                   y = c(d2, d2..2),
                                       sigma = c(sd2, sd2..2),
                                      labels = NULL, tau.prior = tau.prior)
         res2$call <- match.call(expand.dots = FALSE)
         del1 <- c(Mean.dint.del1 = res2$summary["mean","mu"], SD.dint.del1 = res2$summary["sd","mu"])
    }
    
    
    if(Del1 & length(d2) > 1){
      
      result2 <- bayesmeta(                y = d2,
                                       sigma = sd2,
                                      labels = NULL, tau.prior = tau.prior)
      result2$call <- match.call(expand.dots = FALSE)
      
      del1 <- c(result2$summary["mean","mu"], result2$summary["sd","mu"])
      
    }
    
    if(Del1 & length(d2) > 1 & Del1..2 & length(d2..2) > 1) { 
      
      
      res2 <- bayesmeta(                y = d2..2,
                                    sigma = sd2..2,
                                   labels = NULL, tau.prior = tau.prior)
      res2$call <- match.call(expand.dots = FALSE)
      
      
      ds <- c(del1[1], res2$summary["mean","mu"])
      sds <- c(del1[2], res2$summary["sd","mu"])
      
      
      resi2 <- bayesmeta(                y = ds,
                                     sigma = sds,
                                    labels = NULL, tau.prior = tau.prior)
      resi2$call <- match.call(expand.dots = FALSE)
      
      del1 <- c(Mean.dint.del1 = resi2$summary["mean","mu"], SD.dint.del1 = resi2$summary["sd","mu"])
      
    }
    
    
  ####
    
    if(Del2 & length(d3) == 1 & !Del2..2) { 
      
      del2 <- c(Mean.dint.del2 = d3, SD.dint.del2 = sd3)
    }
    
    
    if(Del2 & length(d3) == 1 & Del2..2 & length(d3..2) == 1) { 
      
      res3 <- bayesmeta(                        y = c(d3, d3..2),
                                            sigma = c(sd3, sd3..2),
                                           labels = NULL, tau.prior = tau.prior)
              res3$call <- match.call(expand.dots = FALSE)
     del2 <- c(Mean.dint.del2 = res3$summary["mean","mu"], SD.dint.del2 = re3$summary["sd","mu"])
 
    }
    
    
    if(Del2 & length(d3) > 1){
      
      result3 <- bayesmeta(                     y = d3,
                                            sigma = sd3,
                                           labels = NULL, tau.prior = tau.prior)
           result3$call <- match.call(expand.dots = FALSE)
      
      del2 <- c(result3$summary["mean","mu"], result3$summary["sd","mu"])
      
    }
    
    if(Del2 & length(d3) > 1 & Del2..2 & length(d3..2) > 1) { 
      
      
      res3 <- bayesmeta(                y = d3..2,
                                    sigma = sd3..2,
                                   labels = NULL, tau.prior = tau.prior)
      res3$call <- match.call(expand.dots = FALSE)
      
      
      ds <- c(del2[1], res3$summary["mean","mu"])
      sds <- c(del2[2], res3$summary["sd","mu"])
      
      
      resi3 <- bayesmeta(                y = ds,
                                     sigma = sds,
                                    labels = NULL, tau.prior = tau.prior)
      resi3$call <- match.call(expand.dots = FALSE)
      
      del2 <- c(Mean.dint.del2 = resi3$summary["mean","mu"], SD.dint.del2 = resi3$summary["sd","mu"])
      
    }
    
    
    out <- data.frame(Mean.dint.short = if(Short)short[1] else NA, SD.dint.short = if(Short) short[2]else NA, Mean.dint.del1 = if(Del1)del1[1]else NA, SD.dint.del1 = if(Del1)del1[2]else NA, Mean.dint.del2 = if(Del2)del2[1]else NA, SD.dint.del2 = if(Del2)del2[2]else NA, row.names = NULL) 
    return(out)
  }             
  
  h <- lapply(1:length(L), function(i) G(object = L[[i]], tau.prior = tau.prior))
  names(h) <- study.name
  
  return(h)
}              
              
              
              
#================================================================================================================================
              
              

meta.bayesi <- function(data = NULL, by, tau.prior = function(x){dhalfnormal(x)}, long = FALSE, impute = FALSE, n.sim = 1e5)
{
  

j <- eval(substitute(meta.withini(data = data, by = by, tau.prior = tau.prior, impute = impute, n.sim = n.sim)))
  
study.name <- names(j)

L <- lapply(c('Mean.dint.short', 'SD.dint.short', 'Mean.dint.del1', 'SD.dint.del1', 'Mean.dint.del2',
              'SD.dint.del2'), function(i) {V <- unlist(sapply(j, `[[`, i)); V[!is.na(V)]})

d <- list(L[[1]], L[[3]], L[[5]])

ds <- Filter(NROW, d)
sds <- Filter(NROW, list(L[[2]], L[[4]], L[[6]]))

test <- sapply(d, function(x) length(x) >= 2)

if(all(!test)) stop("Insufficient studies to meta-analyze either 'short-' or 'long-term' effects.", call. = FALSE)


if(test[1]) { result1 <- bayesmeta(     y = ds[[1]],
                                    sigma = sds[[1]],
                                   labels = names(ds[[1]]), tau.prior = tau.prior)
   result1$call <- match.call(expand.dots = FALSE)
} 
  

if(test[2]) { result2 <- bayesmeta(     y = ds[[2]],
                                    sigma = sds[[2]],
                                   labels = names(ds[[2]]), tau.prior = tau.prior)
   result2$call <- match.call(expand.dots = FALSE)
}  
  

if(test[3]) { result3 <- bayesmeta(     y = ds[[3]],
                                    sigma = sds[[3]],
                                   labels = names(ds[[3]]), tau.prior = tau.prior)
   result3$call <- match.call(expand.dots = FALSE)
}  


if(test[2] & test[3] & long){
  
ddelys <- c(result2$summary["mean","mu"], result3$summary["mean","mu"])
sdelys <- c(result2$summary["sd","mu"], result3$summary["sd","mu"])

             result4 <- bayesmeta(      y = ddelys,
                                    sigma = sdelys,
                                   labels = c("Delay1", "Delay2"), tau.prior = tau.prior)
   result4$call <- match.call(expand.dots = FALSE)

   if(test[1])return(list(SHORT = result1, LONG = result4))
   if(!test[1])return(list(LONG = result4))
}

if(!test[1]) message("NOTE: No or insufficient studies to meta-analyze 'short-term' effects.")
if(!test[2]) message("NOTE: No or insufficient studies to meta-analyze 'delayed 1' effects.")
if(!test[3]) message("NOTE: No or insufficient studies to meta-analyze 'delayed 2' effects.")


if(!long || long & !test[2] || long & !test[3]){ 

if(all(test)) return(list(SHORT = result1, DEL1 = result2, DEL2 = result3))
if(test[1] & test[2] & !test[3]) return(list(SHORT = result1, DEL1 = result2))
if(test[1] & !test[2] & !test[3]) return(list(SHORT = result1))
if(!test[1] & test[2] & !test[3]) return(list(DEL1 = result2))
if(!test[1] & !test[2] & test[3]) return(list(DEL2 = result3))
if(!test[1] & test[2] & test[3]) return(list(DEL1 = result2, DEL2 = result3))
   }             
}
               
#===============================================================================================================================                  
                 
dint.plot2 <- function(..., main = NULL, ylab = "Effect Size (dint)", labels = NULL, percent = FALSE, lwd = 1){
  
  m <- Filter(NROW, lapply(list(...), function(x) x[!is.na(x)]))
  L <- length(m)
  n <- substitute(...())
  graphics.off()
  org.par <- par(no.readonly = TRUE)
  on.exit(par(org.par))
  
  
  if(L > 1L) { par(mfrow = n2mfrow(L)) ; set.margin() ; if(percent) par(mar = c(1.5, 2.6, 1.8, 1.6)) }
  
  G <- function(fit, main, labels){  
    
    L <- length(fit)  
    
    bs <- all(sapply(fit, inherits, "bayesmeta"))
    
    if(bs){
      
      mu <- sapply(1:L, function(i) fit[[i]]$summary["mean","mu"])
      lo <- sapply(1:L, function(i) fit[[i]]$summary["95% lower","mu"])
      hi <- sapply(1:L, function(i) fit[[i]]$summary["95% upper","mu"])
      k <- sapply(1:L, function(i) fit[[i]]$k)
      
    } else {
      
      mu <- sapply(1:L, function(i) fit[[i]]$reg_table$b.r[[1]])
      lo <- sapply(1:L, function(i) fit[[i]]$reg_table$CI.L[[1]])            
      hi <- sapply(1:L, function(i) fit[[i]]$reg_table$CI.U[[1]])
      k <- sapply(1:L, function(i) fit[[i]]$N)
    }
    
    x <- 0:(L-1)
    
    plot(x, mu, type = "l", xlim = range(x)+c(-.05, .05), ylim = range(lo, hi), ylab = ylab, lwd = lwd, lty = 2, lend = 1, font.lab = 2,
         xaxt = "n", xlab = NA, panel.last = axis(1, at = x, labels = labels), main = main, las = 1, cex.axis = .9, padj = .3)
    
    invisible(lapply(seq_len(L), function(i) if(!is.na(mu[i])) lines(c(i-1, i-1), c(lo[i], hi[i]), lwd = 4, lend = 1, col = if(bs) 2 else 4)))
    
    text(x, .98*hi, paste0("(k = ", k,")"), cex = .65, font = 2, xpd = NA, srt = 90, pos = 2)
    
    points(x, mu, pch = 22, cex = 6.3, bg = "cyan", col = "magenta", xpd = NA)
    
    text(x, mu,
         round(mu, 3), cex = .9, font = 2, xpd = NA)
    
    text(x, lo,
         round(lo, 3), cex = .9, font = 2, xpd = NA, pos = 1)
    
    text(x, hi,
         round(hi, 3), cex = .9, font = 2, xpd = NA, pos = 3)
    
    if(percent){
      
      text(x*1.02, c(lo, .95*mu, hi),
           paste0("[", dint.norm(c(lo, mu, hi)),"]"), cex = .7, font = 2, xpd = NA, pos = 4, col = "magenta")
      }
    
    mu
  }
  
  res <- invisible(lapply(seq_len(L), function(i) G(m[[i]], main = if(is.null(main)) n[[i]] else if(is.na(main)) NA else main[i], labels = if(is.null(labels)) names(m[[i]]) else labels[[i]])))
  
  z <- if(is.null(main)) as.character(n) else main
  
  mu <- unlist(res)
  data.frame(plot.name = rep(z, each = lengths(res)), mu = round(mu, 4), percent.mu = dint.norm(mu))
}
   
#===============================================================================================================================                                                    
                                                    
dint.plot <- function(..., main = NULL, ylab = "Effect Size (dint)", labels = NULL, percent = FALSE, lwd = 1
                     , reset = TRUE){
  
  m <- Filter(NROW, lapply(list(...), function(x) x[!is.na(x)]))
  L <- length(m)
  n <- substitute(...())

if(reset){
  graphics.off()
  org.par <- par(no.readonly = TRUE)
  on.exit(par(org.par))
}                          
  
  
  if(L > 1L) { par(mfrow = n2mfrow(L)) ; set.margin() ; if(percent) par(mar = c(1.5, 2.6, 1.8, 1.6)) }
  
  G <- function(fit, main, labels){  
    
    L <- length(fit)  
    
    bs <- all(sapply(fit, inherits, "bayesmeta"))
    
    if(bs){
      
      mu <- sapply(1:L, function(i) fit[[i]]$summary["mean","mu"])
      lo <- sapply(1:L, function(i) fit[[i]]$summary["95% lower","mu"])
      hi <- sapply(1:L, function(i) fit[[i]]$summary["95% upper","mu"])
      k <- sapply(1:L, function(i) fit[[i]]$k)
      
    } else {
      
      mu <- sapply(1:L, function(i) fit[[i]]$reg_table$b.r[[1]])
      lo <- sapply(1:L, function(i) fit[[i]]$reg_table$CI.L[[1]])            
      hi <- sapply(1:L, function(i) fit[[i]]$reg_table$CI.U[[1]])
      k <- sapply(1:L, function(i) fit[[i]]$N)
      dfs <- sapply(1:L, function(i) fit[[i]]$reg_table$dfs[[1]] < 4)
    }
    
    x <- 0:(L-1)
    
    plot(x, mu, type = "l", xlim = range(x)+c(-.05, .05), ylim = range(lo, hi), ylab = ylab, lwd = lwd, lty = 2, lend = 1, font.lab = 2,
         xaxt = "n", xlab = NA, panel.last = axis(1, at = x, labels = labels), main = main, las = 1, cex.axis = .9, padj = .3)
    
    invisible(lapply(seq_len(L), function(i) if(!is.na(mu[i])) lines(c(i-1, i-1), c(lo[i], hi[i]), lwd = 4, lend = 1, col = if(bs) 2 else if(!bs & dfs[i]) 8 else 4)))
    
    text(x, .98*hi, paste0("(k = ", k,")"), cex = .65, font = 2, xpd = NA, srt = 90, pos = 2)
    
   # rec <- matrix(rep(c(0.3, 0.22), each = length(x)), ncol = 2)
   # symbols(x, mu, rectangles = rec, inches = FALSE, add = TRUE, bg = "cyan", fg = "magenta")
    
    points(x, mu, pch = 22, cex = 6.3, bg = "cyan", col = "magenta", xpd = NA)
    
    text(x, mu, round(mu, 3), cex = .9, font = 2, xpd = NA)
    
    text(x, lo, round(lo, 3), cex = .9, font = 2, xpd = NA, pos = 1)
    
    text(x, hi, round(hi, 3), cex = .9, font = 2, xpd = NA, pos = 3)
    
    if(percent){
      
      text(x*1.04, c(lo, mu, hi),
           paste0("[", dint.norm(c(lo, mu, hi)),"]"), cex = .7, font = 2, xpd = NA, pos = 4, col = "magenta")
    }
    mu
  }
  
  res <- invisible(lapply(seq_len(L), function(i) G(m[[i]], main = if(is.null(main)) n[[i]] else if(is.na(main)) NA else main[i], labels = if(is.null(labels)) names(m[[i]]) else labels[[i]])))
  
  z <- if(is.null(main)) as.character(n) else main
  
  mu <- unlist(res)
  data.frame(plot.name = rep(z, lengths(res)), mu = round(mu, 4), percent.mu = dint.norm(mu))
}                            
                                                    
#===============================================================================================================================
                  
                  
 dint <- function(data = NULL, by, impute = FALSE, n.sim = 1e4)
 {
   
   data <- rm.allrowNA(trim(data))
   check <- "study.name" %in% names(data)
   if(!check) stop("Add a new column named 'study.name'.", call. = FALSE) 

   
   m <- split(data, data$study.name)
   
   if(is.null(reget(m, control))) stop("Required 'control' group not found.", call. = FALSE)

   if(!missing(by)){ 
     
     s <- if(is.call(by)) by else substitute(by)
     k <- as.list(s)
     
     if("control" %in% k || "!control" %in% k) stop("'control' can't be a moderating variable either alone or with other variables.", call. = FALSE)
     
     H <- lapply(m, function(x) do.call("subset", list(x, s)))
     
     H <- Filter(NROW, H)
     h <- if(length(H) == 0) stop("No study with the requested moderators found.", call. = FALSE) else H
     
     nms <- names(which(!sapply(h, function(x) any(x$control))))
     
     if(length(nms) > 0) {
       h[nms]  <- Map(function(x, y) rbind(y, x[x$control, ]), m[nms], h[nms])
     }
     m <- h
   }
   
   
   if(impute) { 
     
     ar <- formalArgs(rdif)[c(-7, -9)]
     
     args <- lapply(m, function(x) unclass(x[ar]))
     
     argsT <- setNames(lapply(names(args[[1]]), function(i) lapply(args, `[[`, i)), names(args[[1]]))
     
     f <- do.call(Map, c(f = rdif, argsT))
     
     f <- lapply(f, na.locf0)
     
     m <- Map(function(x, y) transform(x, r = na.locf0(y, fromLast = TRUE)), m, f) 
   }     
   
   
   ar <- formalArgs(d.prepos)[-c(21, 22)]
   
   dot.names <- names(data)[!names(data) %in% ar]
   
   args <- lapply(m, function(x) unclass(x[c(ar, dot.names)]))
   
   argsT <- setNames(lapply(names(args[[1]]), function(i) lapply(args, `[[`, i)), names(args[[1]]))
   
   L <- do.call(Map, c(f = d.prepos, argsT))
   
   
   G <- function(m, n.sim)
   {
     
     cdel3 <- reget(m, control & post == 4 & outcome == 1)
     cdel1 <- reget(m, control & post == 2 & outcome == 1)
     cdel2 <- reget(m, control & post == 3 & outcome == 1)
     cs <- reget(m, control & post == 1 & outcome == 1)
     
     tdel3 <- reget(m, !control & post == 4 & outcome == 1)
     tdel1 <- reget(m, !control & post == 2 & outcome == 1)
     tdel2 <- reget(m, !control & post == 3 & outcome == 1)
     ts <- reget(m, !control & post == 1 & outcome == 1) 
     
     if(all(sapply(list(cdel1, cdel2, cdel3, tdel1, tdel2, tdel3, ts, cs), is.null))) stop("Either 'control' or 'post' incorrectly coded.", call. = FALSE)
     
     short <- all(sapply(list(cs, ts), function(x) !is.null(x)))
     
     del1 <- all(sapply(list(cdel1, tdel1), function(x) !is.null(x)))
     
     del2 <- all(sapply(list(cdel2, tdel2), function(x) !is.null(x)))
     
     del3 <- all(sapply(list(cdel3, tdel3), function(x) !is.null(x)))
     
     cdel3..2 <- reget(m, control & post == 4 & outcome == 2)
     cdel1..2 <- reget(m, control & post == 2 & outcome == 2)
     cdel2..2 <- reget(m, control & post == 3 & outcome == 2)
     cs..2 <- reget(m, control & post == 1 & outcome == 2)
     
     tdel3..2 <- reget(m, !control & post == 4 & outcome == 2)
     tdel1..2 <- reget(m, !control & post == 2 & outcome == 2)
     tdel2..2 <- reget(m, !control & post == 3 & outcome == 2)
     ts..2 <- reget(m, !control & post == 1 & outcome == 2)
     
     
     short..2 <- all(sapply(list(cs..2, ts..2), function(x) !is.null(x)))
     
     del1..2 <- all(sapply(list(cdel1..2, tdel1..2), function(x) !is.null(x)))
     
     del2..2 <- all(sapply(list(cdel2..2, tdel2..2), function(x) !is.null(x)))
     
     del3..2 <- all(sapply(list(cdel3..2, tdel3..2), function(x) !is.null(x)))
     
     
     cdel3..3 <- reget(m, control & post == 4 & outcome == 3)
     cdel1..3 <- reget(m, control & post == 2 & outcome == 3)
     cdel2..3 <- reget(m, control & post == 3 & outcome == 3)
     cs..3 <- reget(m, control & post == 1 & outcome == 3)
     
     tdel3..3 <- reget(m, !control & post == 4 & outcome == 3)
     tdel1..3 <- reget(m, !control & post == 2 & outcome == 3)
     tdel2..3 <- reget(m, !control & post == 3 & outcome == 3)
     ts..3 <- reget(m, !control & post == 1 & outcome == 3)
     
     short..3 <- all(sapply(list(cs..3, ts..3), function(x) !is.null(x)))
     
     del1..3 <- all(sapply(list(cdel1..3, tdel1..3), function(x) !is.null(x)))
     
     del2..3 <- all(sapply(list(cdel2..3, tdel2..3), function(x) !is.null(x))) 
     
     del3..3 <- all(sapply(list(cdel3..3, tdel3..3), function(x) !is.null(x)))
     
     
     cdel3..4 <- reget(m, control & post == 4 & outcome == 4)
     cdel1..4 <- reget(m, control & post == 2 & outcome == 4)
     cdel2..4 <- reget(m, control & post == 3 & outcome == 4)
     cs..4 <- reget(m, control & post == 1 & outcome == 4)
     
     tdel3..4 <- reget(m, !control & post == 4 & outcome == 4)
     tdel1..4 <- reget(m, !control & post == 2 & outcome == 4)
     tdel2..4 <- reget(m, !control & post == 3 & outcome == 4)
     ts..4 <- reget(m, !control & post == 1 & outcome == 4)
     
     short..4 <- all(sapply(list(cs..4, ts..4), function(x) !is.null(x)))
     
     del1..4 <- all(sapply(list(cdel1..4, tdel1..4), function(x) !is.null(x)))
     
     del2..4 <- all(sapply(list(cdel2..4, tdel2..4), function(x) !is.null(x)))         
     
     del3..4 <- all(sapply(list(cdel3..4, tdel3..4), function(x) !is.null(x)))
     
     
     if(short){
       nc1 <- m$n[m$control & m$post == 1 & m$outcome == 1]
       nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 1]
       dps <- pair(cs, ts)  
       dppc1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][1])
       dppt1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 1 & m$outcome == 1]
       
       SHORT <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(short..2){
       nc1 <- m$n[m$control & m$post == 1 & m$outcome == 2]
       nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 2]
       dps <- pair(cs..2, ts..2)  
       dppc1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][1])
       dppt1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 1 & m$outcome == 2]
       
       SHORT..2 <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(short..3){
       nc1 <- m$n[m$control & m$post == 1 & m$outcome == 3]
       nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 3]
       dps <- pair(cs..3, ts..3)  
       dppc1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][1])
       dppt1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 1 & m$outcome == 3]
       
       SHORT..3 <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(short..4){
       nc1 <- m$n[m$control & m$post == 1 & m$outcome == 4]
       nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 4]
       dps <- pair(cs..4, ts..4)  
       dppc1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][1])
       dppt1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 1 & m$outcome == 4]
       
       SHORT..4 <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del1){
       nc2 <- m$n[m$control & m$post == 2 & m$outcome == 1]
       nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 1]
       dpdel1 <- pair(cdel1, tdel1)
       dppc2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][1])
       dppt2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 2 & m$outcome == 1]
       
       DEL1 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del1..2){
       nc2 <- m$n[m$control & m$post == 2 & m$outcome == 2]
       nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 2]
       dpdel1 <- pair(cdel1..2, tdel1..2)
       dppc2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][1])
       dppt2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 2 & m$outcome == 2]
       
       DEL1..2 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del1..3){
       nc2 <- m$n[m$control & m$post == 2 & m$outcome == 3]
       nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 3]
       dpdel1 <- pair(cdel1..3, tdel1..3)
       dppc2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][1])
       dppt2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 2 & m$outcome == 3]
       
       DEL1..3 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del1..4){
       nc2 <- m$n[m$control & m$post == 2 & m$outcome == 4]
       nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 4]
       dpdel1 <- pair(cdel1..4, tdel1..4)
       dppc2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][1])
       dppt2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 2 & m$outcome == 4]
       
       DEL1..4 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del2){
       nc3 <- m$n[m$control & m$post == 3 & m$outcome == 1]
       nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 1]
       dpdel2 <- pair(cdel2, tdel2)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 3 & m$outcome == 1]
       
       DEL2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     if(del2..2){
       nc3 <- m$n[m$control & m$post == 3 & m$outcome == 2]
       nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 2]
       dpdel2 <- pair(cdel2..2, tdel2..2)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 3 & m$outcome == 2]
       
       DEL2..2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del2..3){
       nc3 <- m$n[m$control & m$post == 3 & m$outcome == 3]
       nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 3]
       dpdel2 <- pair(cdel2..3, tdel2..3)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 3 & m$outcome == 3]
       
       DEL2..3 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del2..4){
       nc3 <- m$n[m$control & m$post == 3 & m$outcome == 4]
       nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 4]
       dpdel2 <- pair(cdel2..4, tdel2..4)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 3 & m$outcome == 4]
       
       DEL2..4 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     if(del3){
       nc3 <- m$n[m$control & m$post == 4 & m$outcome == 1]
       nt3 <- m$n[m$control == FALSE & m$post == 4 & m$outcome == 1]
       dpdel2 <- pair(cdel3, tdel3)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 4 & m$outcome == 1]
       
       DEL3 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     if(del3..2){
       nc3 <- m$n[m$control & m$post == 4 & m$outcome == 2]
       nt3 <- m$n[m$control == FALSE & m$post == 4 & m$outcome == 2]
       dpdel2 <- pair(cdel3..2, tdel3..2)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 4 & m$outcome == 2]
       
       DEL3..2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del3..3){
       nc3 <- m$n[m$control & m$post == 4 & m$outcome == 3]
       nt3 <- m$n[m$control == FALSE & m$post == 4 & m$outcome == 3]
       dpdel2 <- pair(cdel3..3, tdel3..3)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 4 & m$outcome == 3]
       
       DEL3..3 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del3..4){
       nc3 <- m$n[m$control & m$post == 4 & m$outcome == 4]
       nt3 <- m$n[m$control == FALSE & m$post == 4 & m$outcome == 4]
       dpdel2 <- pair(cdel3..4, tdel3..4)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 4 & m$outcome == 4]
       
       DEL3..4 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     list(SHORT = if(short) SHORT else NULL, SHORT..2 = if(short..2) SHORT..2 else NULL, SHORT..3 = if(short..3) SHORT..3 else NULL, SHORT..4 = if(short..4) SHORT..4 else NULL,
          DEL1 = if(del1) DEL1 else NULL, DEL1..2 = if(del1..2) DEL1..2 else NULL, DEL1..3 = if(del1..3) DEL1..3 else NULL, DEL1..4 = if(del1..4) DEL1..4 else NULL, 
          DEL2 = if(del2) DEL2 else NULL, DEL2..2 = if(del2..2) DEL2..2 else NULL, DEL2..3 = if(del2..3) DEL2..3 else NULL, DEL2..4 = if(del2..4) DEL2..4 else NULL,
          DEL3 = if(del3) DEL3 else NULL, DEL3..2 = if(del3..2) DEL3..2 else NULL, DEL3..3 = if(del3..3) DEL3..3 else NULL, DEL3..4 = if(del3..4) DEL3..4 else NULL) 
   }
   
   setNames(lapply(L, G, n.sim = n.sim), names(L))
 }

#=======================================================================================================================================
             
             
 dintA <- function(data = NULL, by, impute = FALSE, n.sim = 1e4)
 {
     
   names(data) <- trimws(names(data))
   check <- "study.name" %in% names(data)
   if(!check) stop("Add a new column named 'study.name'.", call. = FALSE) 
   
   data$study.name <- trimws(data$study.name)
   data <- rm.allrowNA(data) 
   
   m <- split(data, data$study.name)         
   m <- Filter(NROW, rm.allrowNA2(m)) 
   if(!(length(unique(data$study.name)) == length(m))) stop("Each 'study.name' must be distinct.", call. = FALSE)
   
   if(is.null(reget(m, control))) stop("Required 'control' group not found.", call. = FALSE)
   
   
   if(!missing(by)){ 
     
     s <- substitute(by)
     k <- as.list(s)
     
     if("control" %in% k || "!control" %in% k) stop("'control' can't be a moderating variable either alone or with other variables.", call. = FALSE)
     
     H <- lapply(m, function(x) do.call("subset", list(x, s)))
     
     H <- Filter(NROW, H)
     h <- if(length(H) == 0) stop("No study with the requested moderators found.", call. = FALSE) else H
     
     nms <- names(which(!sapply(h, function(x) any(x$control))))
     
     if(length(nms) > 0) {
       h[nms]  <- Map(function(x, y) rbind(y, x[x$control, ]), m[nms], h[nms])
     }
     m <- h
   }
   
   
   if(impute) { 
     
     ar <- formalArgs(rdif)[c(-7, -9)]
     
     args <- lapply(m, function(x) unclass(x[ar]))
     
     argsT <- setNames(lapply(names(args[[1]]), function(i) lapply(args, `[[`, i)), names(args[[1]]))
     
     f <- do.call(Map, c(f = rdif, argsT))
     
     f <- lapply(f, na.locf0)
     
     m <- Map(function(x, y) transform(x, r = na.locf0(y, fromLast = TRUE)), m, f) 
   }     
   
   
   ar <- formalArgs(d.prepos)[-c(21, 22)]
   
   dot.names <- names(data)[!names(data) %in% ar]
   
   args <- lapply(m, function(x) unclass(x[c(ar, dot.names)]))
   
   argsT <- setNames(lapply(names(args[[1]]), function(i) lapply(args, `[[`, i)), names(args[[1]]))
   
   L <- do.call(Map, c(f = d.prepos, argsT))
   
   
   G <- function(m, n.sim)
   {
     
     cdel3 <- reget(m, control & post == 4 & outcome == 1)
     cdel1 <- reget(m, control & post == 2 & outcome == 1)
     cdel2 <- reget(m, control & post == 3 & outcome == 1)
     cs <- reget(m, control & post == 1 & outcome == 1)
     
     tdel3 <- reget(m, !control & post == 4 & outcome == 1)
     tdel1 <- reget(m, !control & post == 2 & outcome == 1)
     tdel2 <- reget(m, !control & post == 3 & outcome == 1)
     ts <- reget(m, !control & post == 1 & outcome == 1) 
     
     if(all(sapply(list(cdel1, cdel2, cdel3, tdel1, tdel2, tdel3, ts, cs), is.null))) stop("Either 'control' or 'post' incorrectly coded.", call. = FALSE)
     
     short <- all(sapply(list(cs, ts), function(x) !is.null(x)))
     
     del1 <- all(sapply(list(cdel1, tdel1), function(x) !is.null(x)))
     
     del2 <- all(sapply(list(cdel2, tdel2), function(x) !is.null(x)))
     
     del3 <- all(sapply(list(cdel3, tdel3), function(x) !is.null(x)))
     
     cdel3..2 <- reget(m, control & post == 4 & outcome == 2)
     cdel1..2 <- reget(m, control & post == 2 & outcome == 2)
     cdel2..2 <- reget(m, control & post == 3 & outcome == 2)
     cs..2 <- reget(m, control & post == 1 & outcome == 2)
     
     tdel3..2 <- reget(m, !control & post == 4 & outcome == 2)
     tdel1..2 <- reget(m, !control & post == 2 & outcome == 2)
     tdel2..2 <- reget(m, !control & post == 3 & outcome == 2)
     ts..2 <- reget(m, !control & post == 1 & outcome == 2)
     
     
     short..2 <- all(sapply(list(cs..2, ts..2), function(x) !is.null(x)))
     
     del1..2 <- all(sapply(list(cdel1..2, tdel1..2), function(x) !is.null(x)))
     
     del2..2 <- all(sapply(list(cdel2..2, tdel2..2), function(x) !is.null(x)))
     
     del3..2 <- all(sapply(list(cdel3..2, tdel3..2), function(x) !is.null(x)))
     
     
     cdel3..3 <- reget(m, control & post == 4 & outcome == 3)
     cdel1..3 <- reget(m, control & post == 2 & outcome == 3)
     cdel2..3 <- reget(m, control & post == 3 & outcome == 3)
     cs..3 <- reget(m, control & post == 1 & outcome == 3)
     
     tdel3..3 <- reget(m, !control & post == 4 & outcome == 3)
     tdel1..3 <- reget(m, !control & post == 2 & outcome == 3)
     tdel2..3 <- reget(m, !control & post == 3 & outcome == 3)
     ts..3 <- reget(m, !control & post == 1 & outcome == 3)
     
     short..3 <- all(sapply(list(cs..3, ts..3), function(x) !is.null(x)))
     
     del1..3 <- all(sapply(list(cdel1..3, tdel1..3), function(x) !is.null(x)))
     
     del2..3 <- all(sapply(list(cdel2..3, tdel2..3), function(x) !is.null(x))) 
     
     del3..3 <- all(sapply(list(cdel3..3, tdel3..3), function(x) !is.null(x)))
     
     
     cdel3..4 <- reget(m, control & post == 4 & outcome == 4)
     cdel1..4 <- reget(m, control & post == 2 & outcome == 4)
     cdel2..4 <- reget(m, control & post == 3 & outcome == 4)
     cs..4 <- reget(m, control & post == 1 & outcome == 4)
     
     tdel3..4 <- reget(m, !control & post == 4 & outcome == 4)
     tdel1..4 <- reget(m, !control & post == 2 & outcome == 4)
     tdel2..4 <- reget(m, !control & post == 3 & outcome == 4)
     ts..4 <- reget(m, !control & post == 1 & outcome == 4)
     
     short..4 <- all(sapply(list(cs..4, ts..4), function(x) !is.null(x)))
     
     del1..4 <- all(sapply(list(cdel1..4, tdel1..4), function(x) !is.null(x)))
     
     del2..4 <- all(sapply(list(cdel2..4, tdel2..4), function(x) !is.null(x)))         
     
     del3..4 <- all(sapply(list(cdel3..4, tdel3..4), function(x) !is.null(x)))
     
     
     if(short){
       nc1 <- m$n[m$control & m$post == 1 & m$outcome == 1]
       nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 1]
       dps <- pair(cs, ts)  
       dppc1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][1])
       dppt1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 1 & m$outcome == 1]
       
       SHORT <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(short..2){
       nc1 <- m$n[m$control & m$post == 1 & m$outcome == 2]
       nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 2]
       dps <- pair(cs..2, ts..2)  
       dppc1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][1])
       dppt1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 1 & m$outcome == 2]
       
       SHORT..2 <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(short..3){
       nc1 <- m$n[m$control & m$post == 1 & m$outcome == 3]
       nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 3]
       dps <- pair(cs..3, ts..3)  
       dppc1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][1])
       dppt1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 1 & m$outcome == 3]
       
       SHORT..3 <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(short..4){
       nc1 <- m$n[m$control & m$post == 1 & m$outcome == 4]
       nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 4]
       dps <- pair(cs..4, ts..4)  
       dppc1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][1])
       dppt1 <- sapply(1:lengths(dps), function(i) dps[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 1 & m$outcome == 4]
       
       SHORT..4 <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del1){
       nc2 <- m$n[m$control & m$post == 2 & m$outcome == 1]
       nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 1]
       dpdel1 <- pair(cdel1, tdel1)
       dppc2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][1])
       dppt2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 2 & m$outcome == 1]
       
       DEL1 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del1..2){
       nc2 <- m$n[m$control & m$post == 2 & m$outcome == 2]
       nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 2]
       dpdel1 <- pair(cdel1..2, tdel1..2)
       dppc2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][1])
       dppt2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 2 & m$outcome == 2]
       
       DEL1..2 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del1..3){
       nc2 <- m$n[m$control & m$post == 2 & m$outcome == 3]
       nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 3]
       dpdel1 <- pair(cdel1..3, tdel1..3)
       dppc2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][1])
       dppt2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 2 & m$outcome == 3]
       
       DEL1..3 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del1..4){
       nc2 <- m$n[m$control & m$post == 2 & m$outcome == 4]
       nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 4]
       dpdel1 <- pair(cdel1..4, tdel1..4)
       dppc2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][1])
       dppt2 <- sapply(1:lengths(dpdel1), function(i) dpdel1[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 2 & m$outcome == 4]
       
       DEL1..4 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del2){
       nc3 <- m$n[m$control & m$post == 3 & m$outcome == 1]
       nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 1]
       dpdel2 <- pair(cdel2, tdel2)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 3 & m$outcome == 1]
       
       DEL2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     if(del2..2){
       nc3 <- m$n[m$control & m$post == 3 & m$outcome == 2]
       nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 2]
       dpdel2 <- pair(cdel2..2, tdel2..2)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 3 & m$outcome == 2]
       
       DEL2..2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del2..3){
       nc3 <- m$n[m$control & m$post == 3 & m$outcome == 3]
       nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 3]
       dpdel2 <- pair(cdel2..3, tdel2..3)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 3 & m$outcome == 3]
       
       DEL2..3 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del2..4){
       nc3 <- m$n[m$control & m$post == 3 & m$outcome == 4]
       nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 4]
       dpdel2 <- pair(cdel2..4, tdel2..4)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 3 & m$outcome == 4]
       
       DEL2..4 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     if(del3){
       nc3 <- m$n[m$control & m$post == 4 & m$outcome == 1]
       nt3 <- m$n[m$control == FALSE & m$post == 4 & m$outcome == 1]
       dpdel2 <- pair(cdel3, tdel3)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 4 & m$outcome == 1]
       
       DEL3 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     if(del3..2){
       nc3 <- m$n[m$control & m$post == 4 & m$outcome == 2]
       nt3 <- m$n[m$control == FALSE & m$post == 4 & m$outcome == 2]
       dpdel2 <- pair(cdel3..2, tdel3..2)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 4 & m$outcome == 2]
       
       DEL3..2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del3..3){
       nc3 <- m$n[m$control & m$post == 4 & m$outcome == 3]
       nt3 <- m$n[m$control == FALSE & m$post == 4 & m$outcome == 3]
       dpdel2 <- pair(cdel3..3, tdel3..3)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 4 & m$outcome == 3]
       
       DEL3..3 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     
     if(del3..4){
       nc3 <- m$n[m$control & m$post == 4 & m$outcome == 4]
       nt3 <- m$n[m$control == FALSE & m$post == 4 & m$outcome == 4]
       dpdel2 <- pair(cdel3..4, tdel3..4)
       dppc3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][1])
       dppt3 <- sapply(1:lengths(dpdel2), function(i) dpdel2[[1]][[i]][2])
       rv <- m$rev.sign[m$post == 4 & m$outcome == 4]
       
       DEL3..4 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, rev.sign = pairup(rv))))
     }
     
     list(SHORT = if(short) SHORT else NULL, SHORT..2 = if(short..2) SHORT..2 else NULL, SHORT..3 = if(short..3) SHORT..3 else NULL, SHORT..4 = if(short..4) SHORT..4 else NULL,
          DEL1 = if(del1) DEL1 else NULL, DEL1..2 = if(del1..2) DEL1..2 else NULL, DEL1..3 = if(del1..3) DEL1..3 else NULL, DEL1..4 = if(del1..4) DEL1..4 else NULL, 
          DEL2 = if(del2) DEL2 else NULL, DEL2..2 = if(del2..2) DEL2..2 else NULL, DEL2..3 = if(del2..3) DEL2..3 else NULL, DEL2..4 = if(del2..4) DEL2..4 else NULL,
          DEL3 = if(del3) DEL3 else NULL, DEL3..2 = if(del3..2) DEL3..2 else NULL, DEL3..3 = if(del3..3) DEL3..3 else NULL, DEL3..4 = if(del3..4) DEL3..4 else NULL) 
   }
   
   setNames(lapply(L, G, n.sim = n.sim), names(L))
 }
             
             
#=======================================================================================================================================

              
meta.within5 <- function(data = NULL, by, tau.prior = function(x){dhalfnormal(x)}, impute = FALSE, n.sim = 1e5){
  
  L <- eval(substitute(dint(data = data, by = by, impute = impute, n.sim = n.sim)))
  
  study.name <- names(L)
  
  G <- function(m, tau.prior)
  {
    
    d1 <- m$SHORT$dint
    d1..2 <- m$SHORT..2$dint
    d1..3 <- m$SHORT..3$dint
    d1..4 <- m$SHORT..4$dint
    
    sd1 <- m$SHORT$SD
    sd1..2 <- m$SHORT..2$SD
    sd1..3 <- m$SHORT..3$SD
    sd1..4 <- m$SHORT..4$SD
    
    d2 <- m$DEL1$dint
    d2..2 <- m$DEL1..2$dint
    d2..3 <- m$DEL1..3$dint
    d2..4 <- m$DEL1..4$dint
    
    sd2 <- m$DEL1$SD
    sd2..2 <- m$DEL1..2$SD
    sd2..3 <- m$DEL1..3$SD
    sd2..4 <- m$DEL1..4$SD
    
    d3 <- m$DEL2$dint
    d3..2 <- m$DEL2..2$dint
    d3..3 <- m$DEL2..3$dint
    d3..4 <- m$DEL2..4$dint
    
    sd3 <- m$DEL2$SD
    sd3..2 <- m$DEL2..2$SD
    sd3..3 <- m$DEL2..3$SD
    sd3..4 <- m$DEL2..4$SD
    
Short <- !is.null(d1)    ; Del1 <- !is.null(d2)       ; Del2 <- !is.null(d3)
Short..2 <- !is.null(d1..2) ; Del1..2 <- !is.null(d2..2) ; Del2..2 <- !is.null(d3..2)
Short..3 <- !is.null(d1..3) ; Del1..3 <- !is.null(d2..3) ; Del2..3 <- !is.null(d3..3)
Short..4 <- !is.null(d1..4) ; Del1..4 <- !is.null(d2..4) ; Del2..4 <- !is.null(d3..4)
    
    
    if(Short & length(d1) == 1 & !Short..2 & !Short..3 & !Short..4) { 
      
      short <- c(Mean.dint.short = d1, SD.dint.short = sd1)
    }
    
    
    if(Short & length(d1) == 1 & Short..2 & length(d1..2) == 1 & !Short..3 & !Short..4) {
      
      res1 <- bayesmeta(                  y = c(d1, d1..2),
                                      sigma = c(sd1, sd1..2),
                                     labels = NULL, tau.prior = tau.prior)
        res1$call <- match.call(expand.dots = FALSE)
        
      short <- c(Mean.dint.short = res1$summary["mean","mu"], SD.dint.short = res1$summary["sd","mu"])
    }
    
    
    if(Short & length(d1) == 1 & Short..2 & length(d1..2) == 1 & Short..3 & length(d1..3) == 1 & !Short..4) { 
      
      res1 <- bayesmeta(                   y = c(d1, d1..2, d1..3),
                                       sigma = c(sd1, sd1..2, sd1..3),
                                      labels = NULL, tau.prior = tau.prior)
         res1$call <- match.call(expand.dots = FALSE)
         
      short <- c(Mean.dint.short = res1$summary["mean","mu"], SD.dint.short = res1$summary["sd","mu"])
      
    }
    
  
    if(Short & length(d1) == 1 & Short..2 & length(d1..2) == 1 & Short..3 & length(d1..3) == 1 & Short..4 & length(d1..4) == 1) { 
      
      
      res1 <- bayesmeta(                     y = c(d1, d1..2, d1..3, d1..4),
                                         sigma = c(sd1, sd1..2, sd1..3, sd1..4),
                                        labels = NULL, tau.prior = tau.prior)
           res1$call <- match.call(expand.dots = FALSE)
           
      short <- c(Mean.dint.short = res1$summary["mean","mu"], SD.dint.short = res1$summary["sd","mu"])
    }
      
    
    
    if(Short & length(d1) > 1){
      res1 <- bayesmeta(                        y = d1,
                                            sigma = sd1,
                                           labels = NULL, tau.prior = tau.prior)
              res1$call <- match.call(expand.dots = FALSE)
      
      short1 <- c(res1$summary["mean","mu"], res1$summary["sd","mu"])
      
    }
    
    
   if(Short..2 & length(d1..2) > 1){
      
            res2 <- bayesmeta(                y = d1..2,
                                          sigma = sd1..2,
                                         labels = NULL, tau.prior = tau.prior)
            res2$call <- match.call(expand.dots = FALSE)
    
     short2 <- c(res2$summary["mean","mu"], res2$summary["sd","mu"])
    
  }
    
    
    if(Short..3 & length(d1..3) > 1){
      
      res3 <- bayesmeta(                     y = d1..3,
                                         sigma = sd1..3,
                                        labels = NULL, tau.prior = tau.prior)
           res3$call <- match.call(expand.dots = FALSE)
      
      short3 <- c(res3$summary["mean","mu"], res3$summary["sd","mu"])
      
    }   
    
    
    
    if(Short..4 & length(d1..4) > 1){
      
      res4 <- bayesmeta(                    y = d1..4,
                                        sigma = sd1..4,
                                       labels = NULL, tau.prior = tau.prior)
          res4$call <- match.call(expand.dots = FALSE)
      
      short4 <- c(res4$summary["mean","mu"], res4$summary["sd","mu"])
      
    }   
    
    
    if(Short & length(d1) > 1 & !Short..2 & !Short..3 & !Short..4) {
      
      short <- c(Mean.dint.short = short1[1], SD.dint.short = short1[2])
      
    }
    
    
    if(Short & length(d1) > 1 & Short..2 & length(d1..2) > 1 & !Short..3 & !Short..4) {
      
      
      ds <- c(short1[1], short2[1])
      sds <- c(short1[2], short2[2])
      
      resi1 <- bayesmeta(                     y = ds,
                                          sigma = sds,
                                         labels = NULL, tau.prior = tau.prior)
           resi1$call <- match.call(expand.dots = FALSE)
      
      short <- c(Mean.dint.short = resi1$summary["mean","mu"], SD.dint.short = resi1$summary["sd","mu"])
      
    }
    
    
    if(Short & length(d1) > 1 & Short..2 & length(d1..2) > 1 & Short..3 & length(d1..3) > 1 & !Short..4) {
      
      
      ds <- c(short1[1], short2[1], short3[1])
      sds <- c(short1[2], short2[2], short3[2])
      
      resi1 <- bayesmeta(                     y = ds,
                                          sigma = sds,
                                         labels = NULL, tau.prior = tau.prior)
           resi1$call <- match.call(expand.dots = FALSE)
      
      short <- c(Mean.dint.short = resi1$summary["mean","mu"], SD.dint.short = resi1$summary["sd","mu"])
      
    }

    
    
    if(Short & length(d1) > 1 & Short..2 & length(d1..2) > 1 & Short..3 & length(d1..3) > 1 & Short..4 & length(d1..4) > 1 ) {
      
      
      ds <- c(short1[1], short2[1], short3[1], short4[1])
      sds <- c(short1[2], short2[2], short3[2], short4[2])
      
      resi1 <- bayesmeta(                     y = ds,
                                          sigma = sds,
                                         labels = NULL, tau.prior = tau.prior)
           resi1$call <- match.call(expand.dots = FALSE)
      
      short <- c(Mean.dint.short = resi1$summary["mean","mu"], SD.dint.short = resi1$summary["sd","mu"])
      
    }
    
    #####
    
    
    if(Del1 & length(d2) == 1 & !Del1..2 & !Del1..3 & !Del1..4) { 
      
      del1 <- c(Mean.dint.del1 = d2, SD.dint.del1 = sd2)
    }
    
    
    if(Del1 & length(d2) == 1 & Del1..2 & length(d2..2) == 1 & !Del1..3 & !Del1..4) {
      
      res1 <- bayesmeta(                       y = c(d2, d2..2),
                                           sigma = c(sd2, sd2..2),
                                          labels = NULL, tau.prior = tau.prior)
             res1$call <- match.call(expand.dots = FALSE)
      
      del1  <- c(Mean.dint.del1 = res1$summary["mean","mu"], SD.dint.del1 = res1$summary["sd","mu"])
    }
    
    
    if(Del1 & length(d2) == 1 & Del1..2 & length(d2..2) == 1 & Del1..3 & length(d2..3) == 1 & !Del1..4) { 
      
      res1 <- bayesmeta(                        y = c(d2, d2..2, d2..3),
                                            sigma = c(sd2, sd2..2, sd2..3),
                                           labels = NULL, tau.prior = tau.prior)
              res1$call <- match.call(expand.dots = FALSE)
      
      del1 <- c(Mean.dint.del1 = res1$summary["mean","mu"], SD.dint.del1 = res1$summary["sd","mu"])
      
    }
    
    
    if(Del1 & length(d2) == 1 & Del1..2 & length(d2..2) == 1 & Del1..3 & length(d2..3) == 1 & Del1..4 & length(d2..4) == 1) { 
      
      
      res1 <- bayesmeta(                          y = c(d2, d2..2, d2..3, d2..4),
                                              sigma = c(sd2, sd2..2, sd2..3, sd2..4),
                                             labels = NULL, tau.prior = tau.prior)
                res1$call <- match.call(expand.dots = FALSE)
      
      del1 <- c(Mean.dint.del1 = res1$summary["mean","mu"], SD.dint.del1 = res1$summary["sd","mu"])
    }
    
    
    
    if(Del1 & length(d2) > 1){
      
      res1 <- bayesmeta(                        y = d2,
                                            sigma = sd2,
                                           labels = NULL, tau.prior = tau.prior)
              res1$call <- match.call(expand.dots = FALSE)
      
       del11 <- c(res1$summary["mean","mu"], res1$summary["sd","mu"])
      
    }
    
    
    if(Del1..2 & length(d2..2) > 1){
      
      res2 <- bayesmeta(                     y = d2..2,
                                         sigma = sd2..2,
                                        labels = NULL, tau.prior = tau.prior)
           res2$call <- match.call(expand.dots = FALSE)
      
      del12 <- c(res2$summary["mean","mu"], res2$summary["sd","mu"])
      
    }
    
    
    if(Del1..3 & length(d2..3) > 1){
      
      res3 <- bayesmeta(                     y = d2..3,
                                         sigma = sd2..3,
                                        labels = NULL, tau.prior = tau.prior)
           res3$call <- match.call(expand.dots = FALSE)
      
     del13 <- c(res3$summary["mean","mu"], res3$summary["sd","mu"])
      
    }   
    
    
    
    if(Del1..4 & length(d2..4) > 1){
      
      res4 <- bayesmeta(                    y = d2..4,
                                        sigma = sd2..4,
                                       labels = NULL, tau.prior = tau.prior)
          res4$call <- match.call(expand.dots = FALSE)
      
     del14 <- c(res4$summary["mean","mu"], res4$summary["sd","mu"])
      
    }   
    
    
    if(Del1 & length(d2) > 1 & !Del1..2 & !Del1..3 & !Del1..4) {
      
      del1 <- c(Mean.dint.del1 = del11[1], SD.dint.del1 = del11[2])
      
    }
    
    
    if(Del1 & length(d2) > 1 & Del1..2 & length(d2..2) > 1 & !Del1..3 & !Del1..4) {
      
      
      ds <- c(del11[1], del12[1])
      sds <- c(del11[2], del12[2])
      
      resi1 <- bayesmeta(                         y = ds,
                                              sigma = sds,
                                             labels = NULL, tau.prior = tau.prior)
               resi1$call <- match.call(expand.dots = FALSE)
      
      del1 <- c(Mean.dint.del1 = resi1$summary["mean","mu"], SD.dint.del1 = resi1$summary["sd","mu"])
      
    }
    
    
    if(Del1 & length(d2) > 1 & Del1..2 & length(d2..2) > 1 & Del1..3 & length(d2..3) > 1 & !Del1..4) {
      
      
      ds <- c(del11[1], del12[1], del13[1])
      sds <- c(del11[2], del12[2], del13[2])
      
      resi1 <- bayesmeta(                          y = ds,
                                               sigma = sds,
                                              labels = NULL, tau.prior = tau.prior)
                resi1$call <- match.call(expand.dots = FALSE)
      
      del1 <- c(Mean.dint.del1 = resi1$summary["mean","mu"], SD.dint.del1 = resi1$summary["sd","mu"])
      
    }
    
    
    
    if(Del1 & length(d2) > 1 & Del1..2 & length(d2..2) > 1 & Del1..3 & length(d2..3) > 1 & Del1..4 & length(d2..4) > 1 ) {
      
      
      ds <- c(del11[1], del12[1], del13[1], del14[1])
      sds <- c(del11[2], del12[2], del13[2], del14[2])
      
      resi1 <- bayesmeta(                          y = ds,
                                               sigma = sds,
                                              labels = NULL, tau.prior = tau.prior)
                resi1$call <- match.call(expand.dots = FALSE)
      
      del1 <- c(Mean.dint.del1 = resi1$summary["mean","mu"], SD.dint.del1 = resi1$summary["sd","mu"])
      
    }
    
    
    ####
    
    if(Del2 & length(d3) == 1 & !Del2..2 & !Del2..3 & !Del2..4) { 
      
      del2 <- c(Mean.dint.del2 = d3, SD.dint.del2 = sd3)
    }
    
    
    if(Del2 & length(d3) == 1 & Del2..2 & length(d3..2) == 1 & !Del2..3 & !Del2..4) {
      
      res1 <- bayesmeta(                            y = c(d3, d3..2),
                                                sigma = c(sd3, sd3..2),
                                               labels = NULL, tau.prior = tau.prior)
                  res1$call <- match.call(expand.dots = FALSE)
      
      del2  <- c(Mean.dint.del2 = res1$summary["mean","mu"], SD.dint.del2 = res1$summary["sd","mu"])
    }
    
    
    if(Del2 & length(d3) == 1 & Del2..2 & length(d3..2) == 1 & Del2..3 & length(d3..3) == 1 & !Del2..4) { 
      
      res1 <- bayesmeta(                             y = c(d3, d3..2, d3..3),
                                                 sigma = c(sd3, sd3..2, sd3..3),
                                                labels = NULL, tau.prior = tau.prior)
                   res1$call <- match.call(expand.dots = FALSE)
      
      del2 <- c(Mean.dint.del2 = res1$summary["mean","mu"], SD.dint.del2 = res1$summary["sd","mu"])
      
    }
    
    
    if(Del2 & length(d3) == 1 & Del2..2 & length(d3..2) == 1 & Del2..3 & length(d3..3) == 1 & Del2..4 & length(d3..4) == 1) { 
      
      
      res1 <- bayesmeta(                               y = c(d3, d3..2, d3..3, d3..4),
                                                   sigma = c(sd3, sd3..2, sd3..3, sd3..4),
                                                  labels = NULL, tau.prior = tau.prior)
                     res1$call <- match.call(expand.dots = FALSE)
      
      del2 <- c(Mean.dint.del2 = res1$summary["mean","mu"], SD.dint.del2 = res1$summary["sd","mu"])
    }
    
    
    
    if(Del2 & length(d3) > 1){
      res1 <- bayesmeta(                             y = d3,
                                                 sigma = sd3,
                                                labels = NULL, tau.prior = tau.prior)
                   res1$call <- match.call(expand.dots = FALSE)
      
      del21 <- c(res1$summary["mean","mu"], res1$summary["sd","mu"])
      
    }
    
    
    if(Del2..2 & length(d3..2) > 1){
      
      res2 <- bayesmeta(                          y = d3..2,
                                              sigma = sd3..2,
                                             labels = NULL, tau.prior = tau.prior)
                res2$call <- match.call(expand.dots = FALSE)
      
      del22 <- c(res2$summary["mean","mu"], res2$summary["sd","mu"])
      
    }
    
    
    if(Del2..3 & length(d3..3) > 1){
      
      res3 <- bayesmeta(                          y = d3..3,
                                              sigma = sd3..3,
                                             labels = NULL, tau.prior = tau.prior)
                res3$call <- match.call(expand.dots = FALSE)
      
      del23 <- c(res3$summary["mean","mu"], res3$summary["sd","mu"])
      
    }   
    
    
    
    if(Del2..4 & length(d3..4) > 1){
      
      res4 <- bayesmeta(                        y = d3..4,
                                            sigma = sd3..4,
                                           labels = NULL, tau.prior = tau.prior)
              res4$call <- match.call(expand.dots = FALSE)
      
      del24 <- c(res4$summary["mean","mu"], res4$summary["sd","mu"])
      
    }   
    
    
    if(Del2 & length(d3) > 1 & !Del2..2 & !Del2..3 & !Del2..4) {
      
      del2 <- c(Mean.dint.del2 = del21[1], SD.dint.del2 = del21[2])
      
    }
    
    
    if(Del2 & length(d3) > 1 & Del2..2 & length(d3..2) > 1 & !Del2..3 & !Del2..4) {  ##### START HERE change "del1" or "Del1" to "del2" or "Del2".
      
      
      ds <- c(del21[1], del22[1])
      sds <- c(del21[2], del22[2])
      
      resi1 <- bayesmeta(                              y = ds,
                                                   sigma = sds,
                                                  labels = NULL, tau.prior = tau.prior)
                    resi1$call <- match.call(expand.dots = FALSE)
      
      del2 <- c(Mean.dint.del2 = resi1$summary["mean","mu"], SD.dint.del2 = resi1$summary["sd","mu"])
      
    }
    
    
    if(Del2 & length(d3) > 1 & Del2..2 & length(d3..2) > 1 & Del2..3 & length(d3..3) > 1 & !Del2..4) {
      
      
      ds <- c(del21[1], del22[1], del23[1])
      sds <- c(del21[2], del22[2], del23[2])
      
      resi1 <- bayesmeta(                          y = ds,
                                                   sigma = sds,
                                                   labels = NULL, tau.prior = tau.prior)
      resi1$call <- match.call(expand.dots = FALSE)
      
      del2 <- c(Mean.dint.del2 = resi1$summary["mean","mu"], SD.dint.del2 = resi1$summary["sd","mu"])
      
    }
    
    
    if(Del2 & length(d3) > 1 & Del2..2 & length(d3..2) > 1 & Del2..3 & length(d3..3) > 1 & Del2..4 & length(d3..4) > 1 ) {
      
      
      ds <- c(del21[1], del22[1], del23[1], del24[1])
      sds <- c(del21[2], del22[2], del23[2], del24[2])
      
      resi1 <- bayesmeta(                          y = ds,
                                                   sigma = sds,
                                                   labels = NULL, tau.prior = tau.prior)
      resi1$call <- match.call(expand.dots = FALSE)
      
      del2 <- c(Mean.dint.del2 = resi1$summary["mean","mu"], SD.dint.del2 = resi1$summary["sd","mu"])
      
    }
###
    
    out <- data.frame(Mean.dint.short = if(Short)short[1] else NA, SD.dint.short = if(Short) short[2]else NA, Mean.dint.del1 = if(Del1)del1[1]else NA, SD.dint.del1 = if(Del1)del1[2]else NA, Mean.dint.del2 = if(Del2)del2[1]else NA, SD.dint.del2 = if(Del2)del2[2]else NA, row.names = NULL) 
    return(out)
  }             
  
  h <- lapply(1:length(L), function(i) G(m = L[[i]], tau.prior = tau.prior))
  names(h) <- study.name
  
  return(h)
}              

#=======================================================================================================================================
              
meta.withiniii <- function(data = NULL, by, tau.prior = function(x){dhalfnormal(x)}, impute = FALSE, n.sim = 1e5, option = 2, r = .5){
  
  L <- eval(substitute(dint(data = data, by = by, impute = impute, n.sim = n.sim)))
  
  study.name <- names(L)

  G <- function(m, tau.prior)
  {
    
    d1 <- m$SHORT$dint
    d1..2 <- m$SHORT..2$dint
    d1..3 <- m$SHORT..3$dint
    d1..4 <- m$SHORT..4$dint
    
    sd1 <- m$SHORT$SD
    sd1..2 <- m$SHORT..2$SD
    sd1..3 <- m$SHORT..3$SD
    sd1..4 <- m$SHORT..4$SD
    
    d2 <- m$DEL1$dint
    d2..2 <- m$DEL1..2$dint
    d2..3 <- m$DEL1..3$dint
    d2..4 <- m$DEL1..4$dint
    
    sd2 <- m$DEL1$SD
    sd2..2 <- m$DEL1..2$SD
    sd2..3 <- m$DEL1..3$SD
    sd2..4 <- m$DEL1..4$SD
    
    d3 <- m$DEL2$dint
    d3..2 <- m$DEL2..2$dint
    d3..3 <- m$DEL2..3$dint
    d3..4 <- m$DEL2..4$dint
    
    sd3 <- m$DEL2$SD
    sd3..2 <- m$DEL2..2$SD
    sd3..3 <- m$DEL2..3$SD
    sd3..4 <- m$DEL2..4$SD
    
    Short <- !is.null(d1)    ; Del1 <- !is.null(d2)       ; Del2 <- !is.null(d3)
    Short..2 <- !is.null(d1..2) ; Del1..2 <- !is.null(d2..2) ; Del2..2 <- !is.null(d3..2)
    Short..3 <- !is.null(d1..3) ; Del1..3 <- !is.null(d2..3) ; Del2..3 <- !is.null(d3..3)
    Short..4 <- !is.null(d1..4) ; Del1..4 <- !is.null(d2..4) ; Del2..4 <- !is.null(d3..4)
    
    
    if(Short & length(d1) == 1 & !Short..2 & !Short..3 & !Short..4) { 
      
      short <- c(Mean.dint.short = d1, SD.dint.short = sd1)
    }
    
    
    if(Short & length(d1) == 1 & Short..2 & length(d1..2) == 1 & !Short..3 & !Short..4) {
      
      
      if(option == 2){
        
        res <- option2(c(d1, d1..2), c(sd1, sd1..2), r = r)
        
        short <- c(Mean.dint.short = res[1], SD.dint.short = res[2])
      }
      
      if(option == 1){
        
        res <- option1(c(d1, d1..2), c(sd1, sd1..2), r = r)
        
        short <- c(Mean.dint.short = res[1], SD.dint.short = res[2])
      }
      
      if(option == 3){
      
      res1 <- bayesmeta(                  y = c(d1, d1..2),
                                          sigma = c(sd1, sd1..2),
                                          labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      short <- c(Mean.dint.short = res1$summary["mean","mu"], SD.dint.short = res1$summary["sd","mu"])
      }
    }
    
    
    if(Short & length(d1) == 1 & Short..2 & length(d1..2) == 1 & Short..3 & length(d1..3) == 1 & !Short..4) { 
      
      
      if(option == 2){
        
        res1 <- option2(c(d1, d1..2, d1..3), c(sd1, sd1..2, sd1..3), r = r)
        
        short <- c(Mean.dint.short = res1[1], SD.dint.short = res1[2])
        
      }
      
      if(option == 1){
        
        res1 <- option1(c(d1, d1..2, d1..3), c(sd1, sd1..2, sd1..3), r = r)
        
        short <- c(Mean.dint.short = res1[1], SD.dint.short = res1[2])
      }
      
      if(option == 3){
        
      res1 <- bayesmeta(                   y = c(d1, d1..2, d1..3),
                                           sigma = c(sd1, sd1..2, sd1..3),
                                           labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      short <- c(Mean.dint.short = res1$summary["mean","mu"], SD.dint.short = res1$summary["sd","mu"])
      }
    }
    
    
    if(Short & length(d1) == 1 & Short..2 & length(d1..2) == 1 & Short..3 & length(d1..3) == 1 & Short..4 & length(d1..4) == 1) { 
      
      ds <- c(d1, d1..2, d1..3, d1..4)
      sds <- c(sd1, sd1..2, sd1..3, sd1..4)
  
      if(option == 2){
        
        res1 <- option2(ds, sds, r = r)
        
        short <- c(Mean.dint.short = res1[1], SD.dint.short = res1[2])
      }
      
      if(option == 1){
        
        res1 <- option1(ds, sds, r = r)
        
        short <- c(Mean.dint.short = res1[1], SD.dint.short = res1[2])
      }
      
      if(option == 3){
        
      res1 <- bayesmeta(                     y = ds,
                                             sigma = sds,
                                             labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      short <- c(Mean.dint.short = res1$summary["mean","mu"], SD.dint.short = res1$summary["sd","mu"])
      }
    }
    
    
    
    if(Short & length(d1) > 1){
      
      
      if(option == 2){
        
        res1 <- option2(d1, sd1, r = r)
        
        short1 <- c(res1[1], res1[2])
        
      }
      
      if(option == 1){
        
        res1 <- option1(d1, sd1, r = r)
        
        short1 <- c(res1[1], res1[2])
        
      }
      
      if(option == 3){
      res1 <- bayesmeta(                        y = d1,
                                                sigma = sd1,
                                                labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      short1 <- c(res1$summary["mean","mu"], res1$summary["sd","mu"])
      }
    }
    
    
    if(Short..2 & length(d1..2) > 1){
      
      if(option == 2){
        
        res2 <- option2(d1..2, sd1..2, r = r)
        
        short2 <- c(res2[1], res2[2])
        
      }
      
      if(option == 1){
        
        res2 <- option1(d1..2, sd1..2, r = r)
        
        short2 <- c(res2[1], res2[2])
        
      }
      
      if(option == 3){
        
      res2 <- bayesmeta(                y = d1..2,
                                        sigma = sd1..2,
                                        labels = NULL, tau.prior = tau.prior)
      res2$call <- match.call(expand.dots = FALSE)
      
      short2 <- c(res2$summary["mean","mu"], res2$summary["sd","mu"])
      }
    }
    
    
    if(Short..3 & length(d1..3) > 1){
      
      
      if(option == 2){
        
        res3 <- option2(d1..3, sd1..3, r = r)
        
        short3 <- c(res3[1], res3[2])
        
      }
      
      if(option == 1){
        
        res3 <- option1(d1..3, sd1..3, r = r)
        
        short3 <- c(res3[1], res3[2])
        
      }
      
      if(option == 3){
        
      res3 <- bayesmeta(                     y = d1..3,
                                             sigma = sd1..3,
                                             labels = NULL, tau.prior = tau.prior)
      res3$call <- match.call(expand.dots = FALSE)
      
      short3 <- c(res3$summary["mean","mu"], res3$summary["sd","mu"])
      }
    }   
    
    
    
    if(Short..4 & length(d1..4) > 1){
      
      
      if(option == 2){
        
        res4 <- option2(d1..4, sd1..4, r = r)
        
        short4 <- c(res4[1], res4[2])
        
      }
      
      if(option == 1){
        
        res4 <- option1(d1..4, sd1..4, r = r)
        
        short4 <- c(res4[1], res4[2])
      }
      
      if(option == 3){
      res4 <- bayesmeta(                    y = d1..4,
                                            sigma = sd1..4,
                                            labels = NULL, tau.prior = tau.prior)
      res4$call <- match.call(expand.dots = FALSE)
      
      short4 <- c(res4$summary["mean","mu"], res4$summary["sd","mu"])
      }
    }   
    
    
    if(Short & length(d1) > 1 & !Short..2 & !Short..3 & !Short..4) {
      
      short <- c(Mean.dint.short = short1[1], SD.dint.short = short1[2])
      
    }
    
    
    if(Short & length(d1) > 1 & Short..2 & length(d1..2) > 1 & !Short..3 & !Short..4) {
      
      
      ds <- c(short1[1], short2[1])
      sds <- c(short1[2], short2[2])
      
      
      if(option == 2){
        
        resi1 <- option2(ds, sds, r = r)
        
        short <- c(Mean.dint.short = resi1[1], SD.dint.short = resi1[2])
        
      }
      
      if(option == 1){
        
        resi1 <- option1(ds, sds, r = r)
        
        short <- c(Mean.dint.short = resi1[1], SD.dint.short = resi1[2])
        
      }
      
      
      if(option == 3){
      resi1 <- bayesmeta(                     y = ds,
                                              sigma = sds,
                                              labels = NULL, tau.prior = tau.prior)
      resi1$call <- match.call(expand.dots = FALSE)
      
      short <- c(Mean.dint.short = resi1$summary["mean","mu"], SD.dint.short = resi1$summary["sd","mu"])
      }
    }
    
    
    if(Short & length(d1) > 1 & Short..2 & length(d1..2) > 1 & Short..3 & length(d1..3) > 1 & !Short..4) {
      
      
      ds <- c(short1[1], short2[1], short3[1])
      sds <- c(short1[2], short2[2], short3[2])
      
      
      if(option == 2){
        
        resi1 <- option2(ds, sds, r = r)
        
        short <- c(Mean.dint.short = resi1[1], SD.dint.short = resi1[2])
        
      }
      
      if(option == 1){
        
        resi1 <- option1(ds, sds, r = r)
        
        short <- c(Mean.dint.short = resi1[1], SD.dint.short = resi1[2])
      }
      
      if(option == 3){
        
      resi1 <- bayesmeta(                     y = ds,
                                              sigma = sds,
                                              labels = NULL, tau.prior = tau.prior)
      resi1$call <- match.call(expand.dots = FALSE)
      
      short <- c(Mean.dint.short = resi1$summary["mean","mu"], SD.dint.short = resi1$summary["sd","mu"])
      }
    }
    
    
    
    if(Short & length(d1) > 1 & Short..2 & length(d1..2) > 1 & Short..3 & length(d1..3) > 1 & Short..4 & length(d1..4) > 1 ) {
      
      
      ds <- c(short1[1], short2[1], short3[1], short4[1])
      sds <- c(short1[2], short2[2], short3[2], short4[2])
      
      if(option == 2){
        
        resi1 <- option2(ds, sds, r = r)
        
        short <- c(Mean.dint.short = resi1[1], SD.dint.short = resi1[2])
        
      }
      
      if(option == 1){
        
        resi1 <- option1(ds, sds, r = r)
        
        short <- c(Mean.dint.short = resi1[1], SD.dint.short = resi1[2])
        
      }
      
      if(option == 3){
        
      resi1 <- bayesmeta(                     y = ds,
                                              sigma = sds,
                                              labels = NULL, tau.prior = tau.prior)
      resi1$call <- match.call(expand.dots = FALSE)
      
      short <- c(Mean.dint.short = resi1$summary["mean","mu"], SD.dint.short = resi1$summary["sd","mu"])
      }
    }
    
    #####
    
    
    if(Del1 & length(d2) == 1 & !Del1..2 & !Del1..3 & !Del1..4) { 
      
      del1 <- c(Mean.dint.del1 = d2, SD.dint.del1 = sd2)
    }
    
    
    if(Del1 & length(d2) == 1 & Del1..2 & length(d2..2) == 1 & !Del1..3 & !Del1..4) {
      
      
      if(option == 2){
        
        res1 <- option2(c(d2, d2..2), c(sd2, sd2..2), r = r)
        
        del1  <- c(Mean.dint.del1 = res1[1], SD.dint.del1 = res1[2])
        
      }
      
      if(option == 1){
        
        res1 <- option1(c(d2, d2..2), c(sd2, sd2..2), r = r)
        
        del1  <- c(Mean.dint.del1 = res1[1], SD.dint.del1 = res1[2])
      }
      
      if(option == 3){
        
      res1 <- bayesmeta(                       y = c(d2, d2..2),
                                               sigma = c(sd2, sd2..2),
                                               labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      del1  <- c(Mean.dint.del1 = res1$summary["mean","mu"], SD.dint.del1 = res1$summary["sd","mu"])
      }

    }
    
    
    if(Del1 & length(d2) == 1 & Del1..2 & length(d2..2) == 1 & Del1..3 & length(d2..3) == 1 & !Del1..4) { 
      
      ds <- c(d2, d2..2, d2..3)
      sds <- c(sd2, sd2..2, sd2..3)
      
      
      if(option == 2){
        
        res1 <- option2(ds, sds, r = r)
        
        del1 <- c(Mean.dint.del1 = res1[1], SD.dint.del1 = res1[2])
        
      }
      
      if(option == 1){
        
        res1 <- option1(ds, sds, r = r)
        
        del1 <- c(Mean.dint.del1 = res1[1], SD.dint.del1 = res1[2])
      }
      
      if(option == 3){
        
      res1 <- bayesmeta(                        y = ds,
                                                sigma = sds,
                                                labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      del1 <- c(Mean.dint.del1 = res1$summary["mean","mu"], SD.dint.del1 = res1$summary["sd","mu"])
      }

    }
    
    
    if(Del1 & length(d2) == 1 & Del1..2 & length(d2..2) == 1 & Del1..3 & length(d2..3) == 1 & Del1..4 & length(d2..4) == 1) { 
      
      ds <- c(d2, d2..2, d2..3, d2..4)
      sds <- c(sd2, sd2..2, sd2..3, sd2..4)
      
      
      if(option == 2){
        
        res1 <- option2(ds, sds, r = r)
        
        del1 <- c(Mean.dint.del1 = res1[1], SD.dint.del1 = res1[2])
      }
      
      if(option == 1){
        
        res1 <- option1(ds, sds, r = r)
        
        del1 <- c(Mean.dint.del1 = res1[1], SD.dint.del1 = res1[2])
      }
      
      if(option == 3){
      res1 <- bayesmeta(                          y = ds,
                                                  sigma = sds,
                                                  labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      del1 <- c(Mean.dint.del1 = res1$summary["mean","mu"], SD.dint.del1 = res1$summary["sd","mu"])
      }
    }
    
    if(Del1 & length(d2) > 1){
      
      if(option == 2){
        
        res1 <- option2(d2, sd2, r = r)
        
        del11 <- c(res1[1], res1[2])
        
      }
      
      if(option == 1){
        
        res1 <- option1(d2, sd2, r = r)
        
        del11 <- c(res1[1], res1[2])
      }
      
      if(option == 3){
      res1 <- bayesmeta(                        y = d2,
                                                sigma = sd2,
                                                labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      del11 <- c(res1$summary["mean","mu"], res1$summary["sd","mu"])
      }
    }
    
    
    if(Del1..2 & length(d2..2) > 1){
      
      
      if(option == 2){
        
        res2 <- option2(d2..2, sd2..2, r = r)
        
        del12 <- c(res2[1], res2[2])
        
      }
      
      if(option == 1){
        
        res2 <- option1(d2..2, sd2..2, r = r)
        
        del12 <- c(res2[1], res2[2])
        
      }
      
      if(option == 3){
        
      res2 <- bayesmeta(                     y = d2..2,
                                             sigma = sd2..2,
                                             labels = NULL, tau.prior = tau.prior)
      res2$call <- match.call(expand.dots = FALSE)
      
      del12 <- c(res2$summary["mean","mu"], res2$summary["sd","mu"])
      }

    }
    
    
    if(Del1..3 & length(d2..3) > 1){
      
      
      if(option == 2){
        
        res3 <- option2(d2..3, sd2..3, r = r)
        
        del13 <- c(res3[1], res3[2])
        
      }
      
      if(option == 1){
        
        res3 <- option1(d2..3, sd2..3, r = r)
        
        del13 <- c(res3[1], res3[2])
      }
      
      if(option == 3){
        
      res3 <- bayesmeta(                     y = d2..3,
                                             sigma = sd2..3,
                                             labels = NULL, tau.prior = tau.prior)
      res3$call <- match.call(expand.dots = FALSE)
      
      del13 <- c(res3$summary["mean","mu"], res3$summary["sd","mu"])
      }
    }   
 
    
    if(Del1..4 & length(d2..4) > 1){
      
      
      if(option == 2){
        
        res4 <- option2(d2..4, sd2..4, r = r)
        
        del14 <- c(res4[1], res4[2])
        
      }
      
      if(option == 1){
        
        res4 <- option1(d2..4, sd2..4, r = r)
        
        del14 <- c(res4[1], res4[2])
      }
      
      if(option == 3){
      res4 <- bayesmeta(                    y = d2..4,
                                            sigma = sd2..4,
                                            labels = NULL, tau.prior = tau.prior)
      res4$call <- match.call(expand.dots = FALSE)
      
      del14 <- c(res4$summary["mean","mu"], res4$summary["sd","mu"])
      }
    }   
    
    
    if(Del1 & length(d2) > 1 & !Del1..2 & !Del1..3 & !Del1..4) {
      
      del1 <- c(Mean.dint.del1 = del11[1], SD.dint.del1 = del11[2])
      
    }
    
    
    if(Del1 & length(d2) > 1 & Del1..2 & length(d2..2) > 1 & !Del1..3 & !Del1..4) {
      
      
      ds <- c(del11[1], del12[1])
      sds <- c(del11[2], del12[2])

      if(option == 2){
        
        resi1 <- option2(ds, sds, r = r)
        
        del1 <- c(Mean.dint.del1 = resi1[1], SD.dint.del1 = resi1[2])
        
      }
      
      if(option == 1){
        
        resi1 <- option1(ds, sds, r = r)
        
        del1 <- c(Mean.dint.del1 = resi1[1], SD.dint.del1 = resi1[2])
      } 
      
      if(option == 3){
      resi1 <- bayesmeta(                         y = ds,
                                                  sigma = sds,
                                                  labels = NULL, tau.prior = tau.prior)
      resi1$call <- match.call(expand.dots = FALSE)
      
      del1 <- c(Mean.dint.del1 = resi1$summary["mean","mu"], SD.dint.del1 = resi1$summary["sd","mu"])
      }
    }
    
    
    if(Del1 & length(d2) > 1 & Del1..2 & length(d2..2) > 1 & Del1..3 & length(d2..3) > 1 & !Del1..4) {
      
      
      ds <- c(del11[1], del12[1], del13[1])
      sds <- c(del11[2], del12[2], del13[2])
      
 
      if(option == 2){
        
        resi1 <- option2(ds, sds, r = r)
        
        del1 <- c(Mean.dint.del1 = resi1[1], SD.dint.del1 = resi1[2])
        
      }
      
      if(option == 1){
        
        resi1 <- option1(ds, sds, r = r)
        
        del1 <- c(Mean.dint.del1 = resi1[1], SD.dint.del1 = resi1[2])
      } 
      
      if(option == 3){
      resi1 <- bayesmeta(                          y = ds,
                                                   sigma = sds,
                                                   labels = NULL, tau.prior = tau.prior)
      resi1$call <- match.call(expand.dots = FALSE)
      
      del1 <- c(Mean.dint.del1 = resi1$summary["mean","mu"], SD.dint.del1 = resi1$summary["sd","mu"])
      }
    }
    
    
    
    if(Del1 & length(d2) > 1 & Del1..2 & length(d2..2) > 1 & Del1..3 & length(d2..3) > 1 & Del1..4 & length(d2..4) > 1) {
      
      
      ds <- c(del11[1], del12[1], del13[1], del14[1])
      sds <- c(del11[2], del12[2], del13[2], del14[2])
      
      
      if(option == 2){
        
        resi1 <- option2(ds, sds, r = r)
        
        del1 <- c(Mean.dint.del1 = resi1[1], SD.dint.del1 = resi1[2])
        
      }
      
      if(option == 1){
        
        resi1 <- option1(ds, sds, r = r)
        
        del1 <- c(Mean.dint.del1 = resi1[1], SD.dint.del1 = resi1[2])
      }
      
      if(option == 3){
      resi1 <- bayesmeta(                          y = ds,
                                                   sigma = sds,
                                                   labels = NULL, tau.prior = tau.prior)
      resi1$call <- match.call(expand.dots = FALSE)
      
      del1 <- c(Mean.dint.del1 = resi1$summary["mean","mu"], SD.dint.del1 = resi1$summary["sd","mu"])
      }
 
    }
    
    
    
    if(Del2 & length(d3) == 1 & !Del2..2 & !Del2..3 & !Del2..4) { 
      
      del2 <- c(Mean.dint.del2 = d3, SD.dint.del2 = sd3)
    }
    
    
    if(Del2 & length(d3) == 1 & Del2..2 & length(d3..2) == 1 & !Del2..3 & !Del2..4) {
      
      if(option == 2){
        
        res1 <- option2(c(d3, d3..2), c(sd3, sd3..2), r = r)
        
        del2  <- c(Mean.dint.del2 = res1[1], SD.dint.del2 = res1[2])
        
      }
      
      if(option == 1){
        
        res1 <- option1(c(d3, d3..2), c(sd3, sd3..2), r = r)
        
        del2  <- c(Mean.dint.del2 = res1[1], SD.dint.del2 = res1[2])
      }
      
      if(option == 3){
        
      res1 <- bayesmeta(                            y = c(d3, d3..2),
                                                    sigma = c(sd3, sd3..2),
                                                    labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      del2  <- c(Mean.dint.del2 = res1$summary["mean","mu"], SD.dint.del2 = res1$summary["sd","mu"])
      }
    }
    
    
    if(Del2 & length(d3) == 1 & Del2..2 & length(d3..2) == 1 & Del2..3 & length(d3..3) == 1 & !Del2..4) { 
      
      
      ds <- c(d3, d3..2, d3..3)
      sds <- c(sd3, sd3..2, sd3..3)
      
      if(option == 2){
        
        res1 <- option2(ds, sds, r = r)
        
        del2  <- c(Mean.dint.del2 = res1[1], SD.dint.del2 = res1[2])
        
      }
      
      if(option == 1){
        
        res1 <- option1(ds, sds, r = r)
        
        del2  <- c(Mean.dint.del2 = res1[1], SD.dint.del2 = res1[2])
      }
      
      if(option == 3){
        
      res1 <- bayesmeta(                             y = ds,
                                                     sigma = sds,
                                                     labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      del2 <- c(Mean.dint.del2 = res1$summary["mean","mu"], SD.dint.del2 = res1$summary["sd","mu"])
      }
    }
    
    
    if(Del2 & length(d3) == 1 & Del2..2 & length(d3..2) == 1 & Del2..3 & length(d3..3) == 1 & Del2..4 & length(d3..4) == 1) { 
      
      
      ds <- c(d3, d3..2, d3..3, d3..4)
      sds <- c(sd3, sd3..2, sd3..3, sd3..4)
      
      if(option == 2){
        
        res1 <- option2(ds, sds, r = r)
        
        del2  <- c(Mean.dint.del2 = res1[1], SD.dint.del2 = res1[2])
      }
      
      if(option == 1){
        
        res1 <- option1(ds, sds, r = r)
        
        del2  <- c(Mean.dint.del2 = res1[1], SD.dint.del2 = res1[2])
      }
      
      if(option == 3){
        
      res1 <- bayesmeta(                               y = ds,
                                                       sigma = sds,
                                                       labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      del2 <- c(Mean.dint.del2 = res1$summary["mean","mu"], SD.dint.del2 = res1$summary["sd","mu"])
      }
  }  
    
    
    if(Del2 & length(d3) > 1){
      
      if(option == 1){
        
        res1 <- option1(d3, sd3, r = r)
        
        del21 <- c(res1[1], res1[2])
      }
      
      if(option == 2){
        
        res1 <- option2(d3, sd3, r = r)
        
        del21 <- c(res1[1], res1[2])
      }
      
      if(option == 3){
      res1 <- bayesmeta(                             y = d3,
                                                     sigma = sd3,
                                                     labels = NULL, tau.prior = tau.prior)
      res1$call <- match.call(expand.dots = FALSE)
      
      del21 <- c(res1$summary["mean","mu"], res1$summary["sd","mu"])
      }
    }
    
    
    if(Del2..2 & length(d3..2) > 1){
      
      if(option == 1){
        
        res2 <- option1(d3..2, sd3..2, r = r)
        
        del22 <- c(res2[1], res2[2])
      }
      
      if(option == 2){
        
        res2 <- option2(d3..2, sd3..2, r = r)
        
        del22 <- c(res2[1], res2[2])
      }
      
      
      if(option == 3){
      res2 <- bayesmeta(                          y = d3..2,
                                                  sigma = sd3..2,
                                                  labels = NULL, tau.prior = tau.prior)
      res2$call <- match.call(expand.dots = FALSE)
      
      del22 <- c(res2$summary["mean","mu"], res2$summary["sd","mu"])
      }
    }
    
    
    if(Del2..3 & length(d3..3) > 1){
      
      if(option == 1){
        
        res3 <- option1(d3..3, sd3..3, r = r)
        
        del23 <- c(res3[1], res3[2])
      }
      
      if(option == 2){
        
        res3 <- option2(d3..3, sd3..3, r = r)
        
        del23 <- c(res3[1], res3[2])
      }
      
      if(option == 3){
        
      res3 <- bayesmeta(                          y = d3..3,
                                                  sigma = sd3..3,
                                                  labels = NULL, tau.prior = tau.prior)
      res3$call <- match.call(expand.dots = FALSE)
      
      del23 <- c(res3$summary["mean","mu"], res3$summary["sd","mu"])
      }
    }   
    
    
    
    if(Del2..4 & length(d3..4) > 1){
      
      
      if(option == 1){
        
        res4 <- option1(d3..4, sd3..4, r = r)
        
        del24 <- c(res4[1], res4[2])
      }
      
      if(option == 2){
        
        res4 <- option2(d3..4, sd3..4, r = r)
        
        del24 <- c(res4[1], res4[2])
        
      }
      
      
      if(option == 3){
      res4 <- bayesmeta(                        y = d3..4,
                                                sigma = sd3..4,
                                                labels = NULL, tau.prior = tau.prior)
      res4$call <- match.call(expand.dots = FALSE)
      
      del24 <- c(res4$summary["mean","mu"], res4$summary["sd","mu"])
      }
    }   
    
    
    if(Del2 & length(d3) > 1 & !Del2..2 & !Del2..3 & !Del2..4) {
      
      del2 <- c(Mean.dint.del2 = del21[1], SD.dint.del2 = del21[2])
      
    }
    
    
    if(Del2 & length(d3) > 1 & Del2..2 & length(d3..2) > 1 & !Del2..3 & !Del2..4) {  ##### START HERE change "del1" or "Del1" to "del2" or "Del2".
      
      
      ds <- c(del21[1], del22[1])
      sds <- c(del21[2], del22[2])
      
      if(option == 1){
        
        resi1 <- option1(ds, sds, r = r)
        
        del2 <- c(Mean.dint.del2 = resi1[1], SD.dint.del2 = resi1[2])
      }
      
      if(option == 2){
        
        resi1 <- option2(ds, sds, r = r)
        
        del2 <- c(Mean.dint.del2 = resi1[1], SD.dint.del2 = resi1[2])
        
      }
      
      if(option == 3){
      resi1 <- bayesmeta(                              y = ds,
                                                       sigma = sds,
                                                       labels = NULL, tau.prior = tau.prior)
      resi1$call <- match.call(expand.dots = FALSE)
      
      del2 <- c(Mean.dint.del2 = resi1$summary["mean","mu"], SD.dint.del2 = resi1$summary["sd","mu"])
      }
    }
    
    
    if(Del2 & length(d3) > 1 & Del2..2 & length(d3..2) > 1 & Del2..3 & length(d3..3) > 1 & !Del2..4) {
      
      
      ds <- c(del21[1], del22[1], del23[1])
      sds <- c(del21[2], del22[2], del23[2])
      
      if(option == 1){
        
        resi1 <- option1(ds, sds, r = r)
        
        del2 <- c(Mean.dint.del2 = resi1[1], SD.dint.del2 = resi1[2])
      }
      
      if(option == 2){
        
        resi1 <- option2(ds, sds, r = r)
        
        del2 <- c(Mean.dint.del2 = resi1[1], SD.dint.del2 = resi1[2])
        
      }
      
      if(option == 3){
      resi1 <- bayesmeta(                          y = ds,
                                                   sigma = sds,
                                                   labels = NULL, tau.prior = tau.prior)
      resi1$call <- match.call(expand.dots = FALSE)
      
      del2 <- c(Mean.dint.del2 = resi1$summary["mean","mu"], SD.dint.del2 = resi1$summary["sd","mu"])
      }
    }
    
    
    
    if(Del2 & length(d3) > 1 & Del2..2 & length(d3..2) > 1 & Del2..3 & length(d3..3) > 1 & Del2..4 & length(d3..4) > 1 ) {
      
      
      ds <- c(del21[1], del22[1], del23[1], del24[1])
      sds <- c(del21[2], del22[2], del23[2], del24[2])
      
      if(option == 1){
        
        resi1 <- option1(ds, sds, r = r)
        
        del2 <- c(Mean.dint.del2 = resi1[1], SD.dint.del2 = resi1[2])
      }
      
      if(option == 2){
        
        resi1 <- option2(ds, sds, r = r)
        
        del2 <- c(Mean.dint.del2 = resi1[1], SD.dint.del2 = resi1[2])
        
      }
      
      if(option == 3){
      resi1 <- bayesmeta(                          y = ds,
                                                   sigma = sds,
                                                   labels = NULL, tau.prior = tau.prior)
      resi1$call <- match.call(expand.dots = FALSE)
      
      del2 <- c(Mean.dint.del2 = resi1$summary["mean","mu"], SD.dint.del2 = resi1$summary["sd","mu"])
      }
    }
    ###
    
    out <- data.frame(Mean.dint.short = if(Short)short[1] else NA, SD.dint.short = if(Short) short[2]else NA, Mean.dint.del1 = if(Del1)del1[1]else NA, SD.dint.del1 = if(Del1)del1[2]else NA, Mean.dint.del2 = if(Del2)del2[1]else NA, SD.dint.del2 = if(Del2)del2[2]else NA, row.names = NULL) 
    return(out)
  }             
  
  h <- lapply(1:length(L), function(i) G(m = L[[i]], tau.prior = tau.prior))
  names(h) <- study.name
  
  return(h)
}             
              
              
#=======================================================================================================================================

              
meta.within <- function(data = NULL, by, impute = FALSE, n.sim = 1e5, option = 2, r = .5){
  
  L <- eval(substitute(dint(data = data, by = by, impute = impute, n.sim = n.sim)))
  
  study.name <- names(L)
  
  G <- function(m)
  {
    
    d1 <- m$SHORT$dint
    d1..2 <- m$SHORT..2$dint
    d1..3 <- m$SHORT..3$dint
    d1..4 <- m$SHORT..4$dint
    
    sd1 <- m$SHORT$SD
    sd1..2 <- m$SHORT..2$SD
    sd1..3 <- m$SHORT..3$SD
    sd1..4 <- m$SHORT..4$SD
    
    d2 <- m$DEL1$dint
    d2..2 <- m$DEL1..2$dint
    d2..3 <- m$DEL1..3$dint
    d2..4 <- m$DEL1..4$dint
    
    sd2 <- m$DEL1$SD
    sd2..2 <- m$DEL1..2$SD
    sd2..3 <- m$DEL1..3$SD
    sd2..4 <- m$DEL1..4$SD
    
    d3 <- m$DEL2$dint
    d3..2 <- m$DEL2..2$dint
    d3..3 <- m$DEL2..3$dint
    d3..4 <- m$DEL2..4$dint
    
    sd3 <- m$DEL2$SD
    sd3..2 <- m$DEL2..2$SD
    sd3..3 <- m$DEL2..3$SD
    sd3..4 <- m$DEL2..4$SD
    
    d1s <- c(d1, d1..2, d1..3, d1..4)
    sd1s <- c(sd1, sd1..2, sd1..3, sd1..4)
    
    d2s <- c(d2, d2..2, d2..3, d2..4)
    sd2s <- c(sd2, sd2..2, sd2..3, sd2..4)
    
    d3s <- c(d3, d3..2, d3..3, d3..4)
    sd3s <- c(sd3, sd3..2, sd3..3, sd3..4)
    
    Short <- !is.null(d1s)    ; Del1 <- !is.null(d2s)       ; Del2 <- !is.null(d3s)
    
    res <- if(option == 1 & Short)  option1(d1s, sd1s, r = r) else if(option == 2 & Short) option2(d1s, sd1s, r = r) else NA
    
    short <- c(Mean.dint.short = res[1], SD.dint.short = res[2])
    
    res1 <- if(option == 1 & Short)  option1(d2s, sd2s, r = r) else if(option == 2 & Del1) option2(d2s, sd2s, r = r) else NA
    
    del1  <- c(Mean.dint.del1 = res1[1], SD.dint.del1 = res1[2])
    
    res2 <- if(option == 1 & Short) option1(d3s, sd3s, r = r) else if(option == 2 & Del2) option2(d3s, sd3s, r = r) else NA
    
    del2  <- c(Mean.dint.del2 = res2[1], SD.dint.del2 = res2[2])
    
    out <- data.frame(Mean.dint.short = if(Short)short[1] else NA, SD.dint.short = if(Short) short[2]else NA,
                      Mean.dint.del1 = if(Del1)del1[1]else NA, SD.dint.del1 = if(Del1)del1[2]else NA,
                      Mean.dint.del2 = if(Del2)del2[1]else NA,
                      SD.dint.del2 = if(Del2)del2[2]else NA, row.names = NULL)
    return(out)
  }
setNames(lapply(L, G), study.name)
}              
              
              
#=======================================================================================================================================
              
              
meta.bayesii <- function(data = NULL, by, tau.prior = function(x){dhalfnormal(x)}, impute = FALSE, long = FALSE, option = 2, r = .5, dif = TRUE, n.sim = 1e5)
{
   
  j <- eval(substitute(meta.within(data = data, by = by, tau.prior = tau.prior, impute = impute, n.sim = n.sim, option = option, r = r)))
  
  study.name <- names(j)
  
  L <- lapply(c('Mean.dint.short', 'SD.dint.short', 'Mean.dint.del1', 'SD.dint.del1', 'Mean.dint.del2',
                'SD.dint.del2'), function(i) {V <- unlist(sapply(j, `[[`, i)); V[!is.na(V)]})
  
  d <- list(L[[1]], L[[3]], L[[5]])
  
  ds <- Filter(NROW, d)
  sds <- Filter(NROW, list(L[[2]], L[[4]], L[[6]]))
  
  test <- sapply(d, function(x) length(x) >= 2)
  
  if(all(!test)) stop("Insufficient studies to meta-analyze either 'short-' or 'long-term' effects.", call. = FALSE)
  
  
 if(dif){
   
   ds <- one.rm(ds)
   sds <- one.rm(sds)
   ds <- comb.dif.mean(ds)
   sds <- comb.dif.sd(sds, r = r)
   
   res <- bayesmeta(                y = ds,
                                sigma = sds,
                               labels = names(ds), tau.prior = tau.prior)
   res$call <- match.call(expand.dots = FALSE)
   
   return(res)
   
 } else {
  
  if(test[1]) { result1 <- bayesmeta(     y = ds[[1]],
                                          sigma = sds[[1]],
                                          labels = names(ds[[1]]), tau.prior = tau.prior)
  result1$call <- match.call(expand.dots = FALSE)
  } 
  
  
  if(test[2]) { result2 <- bayesmeta(     y = ds[[2]],
                                          sigma = sds[[2]],
                                          labels = names(ds[[2]]), tau.prior = tau.prior)
  result2$call <- match.call(expand.dots = FALSE)
  }  
  
  
  if(test[3]) { result3 <- bayesmeta(     y = ds[[3]],
                                          sigma = sds[[3]],
                                          labels = names(ds[[3]]), tau.prior = tau.prior)
  result3$call <- match.call(expand.dots = FALSE)
  }  
  
  
  if(test[2] & test[3] & long){
    
    ddelys <- c(result2$summary["mean","mu"], result3$summary["mean","mu"])
    sdelys <- c(result2$summary["sd","mu"], result3$summary["sd","mu"])
    
    result4 <- bayesmeta(      y = ddelys,
                               sigma = sdelys,
                               labels = c("Delay1", "Delay2"), tau.prior = tau.prior)
    result4$call <- match.call(expand.dots = FALSE)
    
    if(test[1])return(list(SHORT = result1, LONG = result4))
    if(!test[1])return(list(LONG = result4))
  }
  
  if(!test[1]) message("NOTE: No or insufficient studies to meta-analyze 'short-term' effects.")
  if(!test[2]) message("NOTE: No or insufficient studies to meta-analyze 'delayed 1' effects.")
  if(!test[3]) message("NOTE: No or insufficient studies to meta-analyze 'delayed 2' effects.")
  
  
  if(!long || long & !test[2] || long & !test[3]){ 
    
    if(all(test)) return(list(SHORT = result1, DEL1 = result2, DEL2 = result3))
    if(test[1] & test[2] & !test[3]) return(list(SHORT = result1, DEL1 = result2))
    if(test[1] & !test[2] & !test[3]) return(list(SHORT = result1))
    if(!test[1] & test[2] & !test[3]) return(list(DEL1 = result2))
    if(!test[1] & !test[2] & test[3]) return(list(DEL2 = result3))
    if(!test[1] & test[2] & test[3]) return(list(DEL1 = result2, DEL2 = result3))
     }             
   }                  
}                  

                 
#===============================================================================================================================
                 
                 
meta.bayes <- function(data = NULL, by, tau.prior = function(x){dhalfnormal(x)}, impute = FALSE, long = FALSE, option = 2, r = .5, dif = FALSE, n.sim = 1e5, combine = FALSE)
{
  
  j <- eval(substitute(meta.within(data = data, by = by, impute = impute, n.sim = n.sim, option = option, r = r)))
  
  L <- lapply(c('Mean.dint.short', 'SD.dint.short', 'Mean.dint.del1', 'SD.dint.del1', 'Mean.dint.del2',
                'SD.dint.del2'), function(i) {V <- unlist(sapply(j, `[[`, i)); V[!is.na(V)]})
  
  d <- list(L[[1]], L[[3]], L[[5]])
  
  ds <- Filter(NROW, d)
  sds <- Filter(NROW, list(L[[2]], L[[4]], L[[6]]))
  
  test <- sapply(d, function(x) length(x) >= 2)
  
  if(all(!test)) stop("Insufficient studies to meta-analyze either 'short-' or 'long-term' effects.", call. = FALSE)
  
  
  if(combine){
    
    pt <- unlist(ds)
    ds <- tapply(pt, names(pt), FUN = mean)
    
    pt <- unlist(sds)
    sds <- tapply(pt, names(pt), FUN = opt1, r = r)
    
    
    res <- bayesmeta(                y = ds,
                                     sigma = sds,
                                     labels = names(ds), tau.prior = tau.prior)
    res$call <- match.call(expand.dots = FALSE)
    
    return(res)
  }
  
  
  if(!combine){  
    
    if(dif){
      
      ds <- one.rm(ds)
      sds <- one.rm(sds)
      ds <- comb.dif.mean(ds)
      sds <- comb.dif.sd(sds, r = r)
      
      res <- bayesmeta(                y = ds,
                                       sigma = sds,
                                       labels = names(ds), tau.prior = tau.prior)
      res$call <- match.call(expand.dots = FALSE)
      
      return(res)
      
    } else {
      
      if(test[1]) { result1 <- bayesmeta(          y = ds[[1]],
                                                   sigma = sds[[1]],
                                                   labels = names(ds[[1]]), tau.prior = tau.prior)
      result1$call <- match.call(expand.dots = FALSE)
      } 
      
      
      if(test[2]) { result2 <- bayesmeta(          y = ds[[2]],
                                                   sigma = sds[[2]],
                                                   labels = names(ds[[2]]), tau.prior = tau.prior)
      result2$call <- match.call(expand.dots = FALSE)
      }  
      
      
      if(test[3]) { result3 <- bayesmeta(     y = ds[[3]],
                                              sigma = sds[[3]],
                                              labels = names(ds[[3]]), tau.prior = tau.prior)
      result3$call <- match.call(expand.dots = FALSE)
      }  
      
      
      if(test[2] & test[3] & long){
        
        ddelys <- c(result2$summary["mean","mu"], result3$summary["mean","mu"])
        sdelys <- c(result2$summary["sd","mu"], result3$summary["sd","mu"])
        
        result4 <- bayesmeta(      y = ddelys,
                                   sigma = sdelys,
                                   labels = c("Delay1", "Delay2"), tau.prior = tau.prior)
        result4$call <- match.call(expand.dots = FALSE)
        
        if(test[1])return(list(SHORT = result1, LONG = result4))
        if(!test[1])return(list(LONG = result4))
      }
      
      if(!test[1]) message("NOTE: No or insufficient studies to meta-analyze 'short-term' effects.")
      if(!test[2]) message("NOTE: No or insufficient studies to meta-analyze 'delayed 1' effects.")
      if(!test[3]) message("NOTE: No or insufficient studies to meta-analyze 'delayed 2' effects.")
      
      
      if(!long || long & !test[2] || long & !test[3]){ 
        
        if(all(test)) return(list(SHORT = result1, DEL1 = result2, DEL2 = result3))
        if(test[1] & test[2] & !test[3]) return(list(SHORT = result1, DEL1 = result2))
        if(test[1] & !test[2] & !test[3]) return(list(SHORT = result1))
        if(!test[1] & test[2] & !test[3]) return(list(DEL1 = result2))
        if(!test[1] & !test[2] & test[3]) return(list(DEL2 = result3))
        if(!test[1] & test[2] & test[3]) return(list(DEL1 = result2, DEL2 = result3))
      }             
    }
  }
}                 
                 
                 
#===============================================================================================================================
          
meta.robust <- function(f, data, group, w.model = "CORR", small = TRUE, rho = .8){ 
  
  data <- roundi(data)
  
  f <- if(missing(f)) formula(dint~1) else formula(f)
  
  if(!missing(group)) { 
    
    s <- substitute(group) 
    data <- subset(data, eval(s)) 
    f <- formula(bquote(.(f[[2]]) ~ 1))
  }
  
  m <- robu(f, data = data, studynum = study.name, var = SD^2, small = small, model = w.model, rho = rho)
  m$ml <- f
  m
}                                
    
#===============================================================================================================================
                                      
group.center <- function (var, grp) 
{
  grp <- as.factor(grp)
  grp <- as.numeric(grp)
  var <- as.numeric(var)
  return(var - tapply(var, grp, mean, na.rm = TRUE)[grp])
}

#===============================================================================================================================
                                      
group.mean <- function (var, grp) 
{
  grp <- as.factor(grp)
  grp <- as.numeric(grp)
  var <- as.numeric(var)
  return(tapply(var, grp, mean, na.rm = TRUE)[grp])
}                                      

#===================================Modern Inter-rater Reliability in Meta-Analysis=====================================================

rm.allrowNA <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) i[rowSums(is.na(i) | i == "") != ncol(i), , drop = FALSE])
    
  } else { X[rowSums(is.na(X) | X == "") != ncol(X), , drop = FALSE] }
}

#===============================================================================================================================
       
rm.allcolNA <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) i[, colSums(is.na(i) | i == "") != nrow(i), drop = FALSE])
    
  } else { X[, colSums(is.na(X) | X == "") != nrow(X), drop = FALSE] }
}

#===============================================================================================================================
           
rm.colrowNA <- function(X){

r <- rm.allrowNA(X)
rm.allcolNA(r)  

}                                      

           
#================================================================================================================================
           
           
full.clean2 <- function(X, which, keep.org.name = TRUE)
  {
  
  X <- rm.colrowNA(X)
  
  X <- lapply(X, function(x) setNames(x, sub("\\.\\d+$", "", names(x))))
  
  if(inherits(X, "list") & length(X) >= 2){
    
    X <- lapply(seq_along(X), function(i) X[[i]][ , !names(X[[i]]) %in% which])
    if(keep.org.name) lapply(X, function(x) setNames(X, sub("\\.\\d+$", "", names(X))))
    X
    
  } else {
    
    X <- X[[1]][ , !names(X[[1]]) %in% which]
    if(keep.org.name) names(X) <- sub("\\.\\d+$", "", names(X))
    list(X)
  }
}           
           

#================================================================================================================================
                             

full.clean1 <- function(X, which)
{
  
  f <- function(dat, vec) {
    i1 <- !names(dat) %in% vec
    setNames(dat[i1], names(dat)[i1])
  }
  
  X <- rm.colrowNA(X)
  
  X <- lapply(X, function(x) setNames(x, sub("\\.\\d+$", "", names(x))))
    
  lapply(X, f, vec = which)
}                             
                             
#================================================================================================================================
  
drop.col <- function(dat, vec){
  
  vec <- trimws(vec)
  names(dat) <- trimws(names(dat))
    
  f <- function(dat, vec) {
    i1 <- !names(dat) %in% vec
    setNames(dat[i1], names(dat)[i1])
  }
  
  if(inherits(dat, "list")) { lapply(dat, f, vec = vec)
  } else { f(dat, vec) }
}               
              
#================================================================================================================================              
              
full.clean <- function(X, omit, all = TRUE, omit.auto.suffix = TRUE)
{
  
  X <- rm.colrowNA(X)
  
  X <- if(inherits(X, "list") & omit.auto.suffix){ lapply(X, function(x) setNames(x, sub("\\.\\d+$", "", names(x)))) 
    
    } else if(inherits(X, "data.frame") & omit.auto.suffix) { setNames(X, sub("\\.\\d+$", "", names(X))) } else { X }
  
  if(all){ X } else { 
    
    drop.col(X, vec = omit)
  }
}              
              
#================================================================================================================================                                      
                                      
kap <- function (x, level = .95)
{
    
  x <- matrix(table(x[[1]], x[[2]]), 2)
  d  <- diag(x)
  n  <- sum(x)
  nc <- ncol(x)
  colFreqs <- colSums(x)/n
  rowFreqs <- rowSums(x)/n
  
  kappa <- function (po, pc) (po - pc) / (1 - pc)
  
  std  <- function (p, pc, kw, W = diag(1, ncol = nc, nrow = nc)) {
    sqrt((sum(p * sweep(sweep(W, 1, W %*% colSums(p) * (1 - kw)), 2, W %*% rowSums(p) * (1 - kw)) ^ 2) - (kw - pc * (1 - kw)) ^ 2) / crossprod(1 - pc) / n)
  }
  
  po <- sum(d) / n
  pc <- crossprod(colFreqs, rowFreqs)[1]
  k <- kappa(po, pc)
  s <- std(x / n, pc, k)
  
  p <- (1 - level) / 2
  q <- qnorm(c(p, 1-p))
  ci <- k + q*s
  
  return(c(
    KAPPA = k,
    lower = ci[1],
    upper = ci[2],
    conf.level = level))
}                                                                             
  
#===============================================================================================================================                                      
                                      
kappa <- function(X, Y, level = .95, raw.sheet = FALSE){
  
  if(raw.sheet){
    
    ar <- head(formalArgs(d.prepos), -1)
    dot.names <- names(X)[!names(X) %in% ar]
    X <- X[dot.names]
  }
  
  L <- Map(table, X, Y[names(X)])
  
  lapply(L, kap, level = level)
}    
    
#===============================================================================================================================
                                      

efa <- function(x, factors, data = NULL, covmat = NULL, n.obs = NA,
                subset, na.action = "na.omit", start = NULL, center = FALSE,
                scores = c("none", "regression", "Bartlett"),
                rotation = "varimax", control = NULL, ...)
{
  cc <- match.call(expand.dots = FALSE)
  cc[[1]] <- quote(factanal)
  fit <- eval.parent(cc)
  fit$call <- match.call(expand.dots = FALSE)
 
  noncent <- if(is.null(data)) scale(as.data.frame(x), center = center) else scale(as.data.frame(data), center = center)
  
  Rvv_1 <- solve(fit$correlation)
  Pvf <- fit$loadings
  Wvf <- Rvv_1%*%Pvf
  
  scores <- data.frame(noncent%*%Wvf)
  
  fit$scores <- scores
  
  return(fit)
}                                                         
                                                                          

#===============================================================================================================================
                                                          
detail2 <- function(X, useNA = "ifany"){
  
  nr <- nrow(X)
  nc <- ncol(X)
  tab <- table(row(X), unlist(X), useNA = useNA)
  pj <- apply(tab, 2, sum)/(nr * nc)
  pjk <- (apply(tab^2, 2, sum) - nr * nc * pj)/(nr * nc * (nc - 1) * pj)
  K <- (pjk - pj)/(1 - pj)
  h <- names(K)
  h[is.na(h)] <- "NA"
  setNames(K, h)
}         

#===============================================================================================================================
                                                          
detail <- function(X, useNA = "ifany") {
  X <- as.matrix(X)
  tab <- table(row(X), unlist(X), useNA = useNA)
  w <- diag(ncol(tab))
  rosum <- rowSums(tab)
  obs_oc <- tab * (t(w %*% t(tab)) - 1)
  obs_c <- colSums(obs_oc)
  max_oc <- tab * (rosum - 1)
  max_c <- colSums(max_oc)
  SA <- obs_c / max_c
  h <- names(SA)
  h[is.na(h)] <- "NA"
  setNames(SA, h)
}                                                       
                                                          
#===============================================================================================================================                                                          

set.margin <- function() 
{
  par(mgp = c(1.5, 0.14, 0), mar = c(2.5, 2.6, 1.8, .5), 
      tck = -0.02)
}                                                         

#===============================================================================================================================                                                          
 
splot <- function(y, main, lwd = 5, lend = 2, show.sa = FALSE, digits = 3, cex.sa = .9){
  
  ll <- length(y)
  
  x <- seq_len(ll)
  
  plot(x, y, type = "h", main = main, xlim = c(.95, 1.02*max(x)), ylim = 0:1,
       ylab = "SA%", xaxt = "n", xlab = "Category", lend = lend, lwd = lwd,
       col = colorRampPalette(c(4, 2))(ll), font.lab = 2, 
       panel.first = abline(h = 0, col = 8), las = 1, cex.axis = .9)
  
  if(show.sa) text(x[y != 0]-.015, .4, round(y[y != 0], digits), pos = 2, xpd = NA, srt = 90, font = 2, cex = cex.sa)
  
  axis(1, at = x, labels = names(y))
}
                                                       
                                                          
#===============================================================================================================================
      
irr <- int <- function (X, nsim = 1e3, useNA = "ifany", level = .95, digits = 6, raw = TRUE) 
{
  
  if(!inherits(X, c("data.frame", "matrix", "table"))) stop("Ratings must be 'data.frame', 'matrix', and if not raw, a 'table'.", call. = FALSE)
  
  if(raw) X <- table(row(X), unlist(X), useNA = useNA)
  
  X2 <- X * (X - 1)
  sumcol <- colSums(X)
  sumrow <- rowSums(X)
  nc <- ncol(X)
  nr <- nrow(X)
  tot <- sum(X)
  pij <- X2/(sumrow * (sumrow - 1))
  pi <- rowSums(pij)
  p <- mean(pi)
  pj <- sumcol/tot
  pj2 <- pj^2
  pe <- sum(pj2)
  KAPPA <- (p - pe)/(1 - pe)
  s <- (nc * p - 1)/(nc - 1)
  pi.v.boot <- replicate(nsim, pi.boot <- sample(pi, size = nr, replace = TRUE))
  p.boot <- colMeans(pi.v.boot)
  s.boot <- sapply(seq_len(nsim), function(i) (nc * p.boot[i] - 1)/(nc - 1))
  
  p <- (1 - level) / 2
  s.boot.ci <- quantile(s.boot, probs = c(p, 1-p), na.rm = TRUE)
  
  return(round(c(Fleiss_KAPPA = KAPPA, 
                 Sindex = s, 
                 lower = s.boot.ci[[1]], 
                 upper = s.boot.ci[[2]], 
                 conf.level = level), digits))
}                                      

                   
#===============================================================================================================================                   
                   
                                   
int2 <- function(X, level = .95, useNA = "ifany", nsim = 1e3, digits = 4, raw = TRUE){ 
  
  X <- table(row(X), unlist(X), useNA = useNA)
  
  agree.mat <- as.matrix(X) 
  n <- nrow(agree.mat) # number of studies or groups within studies
  q <- ncol(agree.mat) # number of categories
  f <- 0               # population correction 
  
    weights.mat <- diag(q)

  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more
  
  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat) * sum(pi.vec*(1-pi.vec)) / (q*(q-1))
  ac1 <- (pa-pe)/(1-pe)
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec == 0) # this replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec
  
  pe.r2 <- pe*(ri.vec>=2)
  ac1.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pe.ivec <- (sum(weights.mat)/(q*(q-1))) * (agree.mat%*%(1-pi.vec))/ri.vec
  ac1.ivec.x <- ac1.ivec - 2*(1-ac1) * (pe.ivec-pe)/(1-pe)
  
  var.ac1 <- ((1-f)/(n*(n-1))) * sum((ac1.ivec.x - ac1)^2)
  stderr <- sqrt(var.ac1)
  p.value <- 2*(1-pt(ac1/stderr,n-1))
  
  lower <- ac1 - stderr*qt(1-(1-level)/2,n-1)
  upper <- min(1,ac1 + stderr*qt(1-(1-level)/2,n-1))
  
  return(round(c(AC = ac1, lower = lower, upper = upper, conf.level = level), digits))
}                   
                   
#===============================================================================================================================
                       
is.constant <- function(x) length(unique(x)) == 1L 

#===============================================================================================================================                   
                   
drop.inner.list <- function(L, what, omit.auto.suffix = TRUE) {
  
  if(omit.auto.suffix) L <- lapply(L, function(x) setNames(x, sub("\\.\\d+$", "", names(x))))
  
  L[!names(L) %in% what]
}
                                   
#===============================================================================================================================
   
is.unique <- function(X, which){
  
  f <- function(X, which) { nrow(unique(X[which])) == nrow(X[which]) }
  
  test <- if(inherits(X, "list")) sapply(X, f, which) else if(inherits(X, "data.frame")) f(X, which) else {
    
    length(unique(X)) == length(X)
  }
  base::all(test)  
  }                                   
                                                                          
#===============================================================================================================================
                                 
meta_rate <- function(..., sub.name = "group.name", nsim = 1e3, level = .95, useNA = "ifany", type = c("s", "ac"), na.rm = FALSE, digits = 3, common = FALSE, all = TRUE, drop = NULL, plot = TRUE, lwd = 5, lend = 1, show.sa = TRUE, sub.level = NULL, study.level = NULL, file.name = NULL, reset = TRUE, rev.page = FALSE, cex.sa = .9)
{
  
  r <- list(...) 
  
  if(!all(sapply(r, inherits, c("data.frame", "matrix")))) stop("Coding-sheet(s) must be 'Excel CSV' files, 'data.frame' or 'matrix'.", call. = FALSE)
  
  n.df <- length(r)
  
  r <- lapply(r, as.data.frame)
  
  ar <- formalArgs(d.prepos)[-c(2, 22)]
  
  r <- full.clean(r, ar, all)
  
  r <- lapply(r, trim)
  
  check <- all(sapply(r, function(i) "study.name" %in% names(i)))
  
  if(!check) stop("Add a new column named 'study.name' to the coding sheet(s).", call. = FALSE)
  
  r <- lapply(r, function(x) do.call(rbind, c(split(x, x$study.name), make.row.names = FALSE)))
  
  drop <- trimws(drop)              
  drop <- drop[!drop %in% "study.name"]
  
  if(length(drop) != 0) r <- drop.col(r, drop)   
  
  r <- unname(r)
  
  sub.name <- trimws(sub.name)
  if(n.df == 1) tbl <- table(names(r[[1]])[!names(r[[1]]) %in% c("study.name", sub.name)])
  
  com.names <- if(n.df >= 2) { 
    
    ok <- is.constant(sapply(r, nrow))
    
    if(!ok) stop("The coding sheets don't have the same number of rows.", call. = FALSE)
    
    if(common) { Reduce(intersect, lapply(r, names)) 
      
    } else {
      
      vec <- names(unlist(r, recursive = FALSE))
      unique(vec[duplicated(vec)])
    }
    
  } else { 
    
    if(common) { 
      
      names(which(tbl == max(tbl)))
      
    } else {
      
      names(which(tbl >= 2))
    }
  }
  
  dot.names <- if(all) com.names else com.names[!com.names %in% ar]
  
  if(length(dot.names) == 0) stop("No 2 raters detected OR no two moderators names match.", call. = FALSE)
  
  if(n.df >= 2) { 
    
    r <- do.call(cbind, r)
    
    tbl <- table(names(r)[!names(r) %in% c("study.name", sub.name)]) 
    
  } else { r <- r[[1]]
  
  }
  
  n.coder <- if(common) { 
    
    tbl[tbl == max(tbl)] 
    
  } else {
    
    tbl[tbl >= 2]
  }
  
  i1 <- colnames(r) != 'study.name'
  st.level <- names(which(sapply(split.default(r[i1], names(r)[i1]), function(x) all(!colSums(!aggregate(.~ study.name, transform(x, study.name = r$study.name), FUN = is.constant)[-1])))))
  
  st.level <- st.level[st.level %in% dot.names]
  
  exclude <- trimws(sub.level)
  
  st.level <- st.level[!st.level %in% c(exclude,"study.name", sub.name)]
  
  L <- split.default(r[names(r) %in% dot.names], names(r)[names(r) %in% dot.names])
  
  if(length(st.level) != 0) L[st.level] <- lapply(L[st.level], function(x) x[ave(seq_along(x[[1]]), r$study.name, FUN = seq_along) == 1, ]) 
  
  L <- drop.inner.list(L, c("study.name", sub.name))
  
  if(na.rm) L <- lapply(L, na.omit)
  
  type <- trimws(type)
  f <- if(type == "s") int else int2
  out <- lapply(L, f, nsim = nsim, level = level, digits = digits, useNA = useNA, raw = TRUE)
  
  A <- lapply(L, detail, useNA = useNA)
  
  study.level <- sapply(seq_along(out), function(i) names(out)[[i]] %in% st.level)
  
  d <- data.frame(out)
  
  d[] <- lapply(d, as.list)
  
  if(plot){
    
    n <- length(L)
    
    if(reset){
      graphics.off()
      org.par <- par(no.readonly = TRUE)
      on.exit(par(org.par))
    }
    dev <- if(!rev.page) n2mfrow(n) else rev(n2mfrow(n))
    if(n > 1L) { par(mfrow = dev) ; set.margin() }
    
    invisible(mapply(splot, y = A, main = names(A), lwd = lwd, lend = lend, show.sa = show.sa, digits = digits, cex.sa = cex.sa))
  }
  
  res <- data.frame(t(rbind(d, row.comprd = sapply(L, nrow), min.cat = sapply(A, function(i) if(any(i < 1)) names(i)[which.min(i)] else "--"), 
                            n.coder = n.coder, study.level = ifelse(study.level, "Yes", "No"))))
  
  file.name <- trimws(file.name)
  
  if(length(file.name) != 0){
    output <- data.frame(lapply(res, unlist))
    nm <- paste0(file.name, ".csv")
    ur <- try(write.csv(output, nm), silent = TRUE)
    if(inherits(ur, "try-error")) stop(paste0("\nClose the Excel file '", nm, "' and try again OR pick another file name."), call. = FALSE)
    message(paste0("\nNote: Check folder '", basename(getwd()),"' for the Excel file '", nm, "'.\n"))
  }
  
  return(res)
}
                                               
#===============================================================================================================================
      
metal <- function(data = NULL, mod, mu.prior = c("mean" = NA, "sd" = NA), tau.prior = function(x){dhalfnormal(x)}, impute = FALSE, n.sim = 1e4, option = 1, r = .5){
  
  f1 <- function(data, zy, impute, n.sim, option, r){ 
    
    L <- eval(substitute(dint(data = data, by = zy, impute = impute, n.sim = n.sim)))
    
    ds <- Filter(Negate(is.null), lapply(L, function(x) do.call(rbind, x)$dint))
    sds <- Filter(Negate(is.null), lapply(L, function(x) do.call(rbind, x)$SD)) 
    
    f <- if(option == 1) option1 else option2
    
    mapply(f, ds = ds, sds = sds, r = r, SIMPLIFY = FALSE)
  }
  
  f2 <- function(j, tau.prior){  
    
    ds <- sapply(seq_along(j), function(i) j[[i]][1])
    sds <- sapply(seq_along(j), function(i) j[[i]][2])
    
    test <- length(ds) >= 2
    
    if(!test) return(NA)
    
    res <- bayesmeta(        y = ds,
                             sigma = sds,
                             labels = names(j), 
                             tau.prior = tau.prior,
                             mu.prior = mu.prior)
    res$call <- match.call(expand.dots = FALSE)
    
    return(res)
  }  
  
  chep <- sort(unique(na.omit(data$time)))
  
  G <- if(missing(mod)) { lapply(chep, function(y) bquote(time == .(y))) 
    
  } else {
    
    s <- substitute(mod)
    lapply(chep, function(x) bquote(.(s) & time == .(x)))
  }
  
  go <- length(G)
  
  k <- vector("list", go)
  
  for(w in seq_len(go)) k[[w]] <- try(f1(data = data, zy = G[[w]], impute = impute, n.sim = n.sim, option = option, r = r), silent = TRUE)
  
  so <- length(k)
  
  z <- vector("list", so)
  
  for(a in seq_len(so)) z[[a]] <- f2(j = k[[a]], tau.prior = tau.prior)
  
  setNames(z, chep)
}       
            

#===============================================================================================================================
         
meta.set <- function(data, file.name = NULL, na = ""){
  
  data <- rm.allrowNA(trim(data))
  
  # L <- dint(data)
  L <- suppressWarnings(dint(data))
    
  d <- do.call(rbind, 
               Map(cbind, pp <- Filter(Negate(is.null), lapply(L, function(x) 
                 do.call(rbind, x))), 
                 id = seq_along(pp)))
  
  h <- cbind(study.name = sub("(.*)\\.(SHORT|DEL(1|2|3))(\\.+\\d.*)?", "\\1", rownames(d)), d)
  rownames(h) <- NULL
  
  n <- split(h, h$study.name)  
  hh <- setNames(lapply(n, NROW), names(n))
  
  D <- data
  
  m <- split(D, D$study.name)
  
  ff <- setNames(lapply(seq_along(m), function(i) NROW(subset(m[[i]], !control))), names(m))
  
  eq <- identical(ff, hh)
  
  if(eq){
    
    message(paste0("'dints' successfully computed and reshaped into 'long form'."))
    
    ar <- formalArgs(d.prepos)[-(20:22)]
    
    mod.names <- names(D)[!names(D) %in% ar]
    
    mods <- subset(D[order(D$study.name), ], !control, select = mod.names)
    
    H <- roundi(cbind(h, mods))
    
    if(!is.null(file.name)) { 
      file.name <- trimws(file.name)
      nm <- paste0(file.name, ".csv")
      ur <- try(write.csv(H, nm, row.names = FALSE, na = na), silent = TRUE)
      if(inherits(ur, "try-error")) stop(paste0("\nClose the Excel file '", nm, "' and try again OR pick another file name."), call. = FALSE)
      message(paste0("\nNote: Check folder '", basename(getwd()),"' for the Excel file '", nm, "'.\n"))
      }
    
    return(H)
    
  } else {
    
    message(paste0("\nProblem in coding sheet detected. See error analysis below:\n"))
    ur <- try(test.sheet(D, metaset = TRUE), silent = TRUE)
    
    if(inherits(ur, "try-error")) { stop("\nIncomplete data: check descrptive columns ('n', 'mpre' etc.)", call. = FALSE) 
      
      } else {
    
    message("\nSo, pre-post effect sizes can be calculated but not 'dints'... why?")    
        
    ur <- try(dint(D), silent = TRUE)
    
    if(inherits(ur, "try-error")) { stop("\nDue to incorrect/incomplete data: Check columns 'post', 'control', & 'outcome'.", call. = FALSE) 
      
      } else if(is.null(unlist(ur))) {
      
        message("Due probably to column 'post' having inconsistent values.")
        
    } else { message("\nDue to coding inconsistencies/errors (ex. part of a study coded or entered into this program but not all of it).")}
       }
    }
}               
         
#===============================================================================================================================
                                      
                                      
fold <- function(x, at, labels = seq_len(length(at)+1), xlab = "x", ylab = "Frequency", lend = 1, na.rm = FALSE, plot = FALSE, ...){
 
  if(na.rm) x <- na.omit(x)
  
  cats <- cut(x, breaks = c(-Inf, at, Inf), include.lowest = TRUE, labels = labels)
  tab <- table(x, dnn = NULL)
  cattab <- table(cats, dnn = NULL)
  
if(plot){
  
  graphics.off()
  org.par <- par(no.readonly = TRUE)
  on.exit(par(org.par))   
  par(mfrow = c(2, 1))
  set.margin() 
  
  cols <- colorRampPalette(c(4, 2))(length(unique(cats)))
  
  grp <- cut(as.numeric(names(tab)), 
             breaks = c(-Inf, at, Inf), 
             include.lowest = TRUE)
  
  plot(tab, xlab = xlab, ylab = ylab, main = "Original", panel.f = abline(h = 0, col = 8), col = cols[grp], lend = lend, ...)
  
  plot(cattab, xlab = xlab, ylab = ylab, main = "Categorized", panel.f = abline(h = 0, col = 8), col = cols, lend = lend, ...)
  
  box()
} 
  list(Original = tab, Categorized = cattab, cats = as.numeric(cats))
} 
 
                                                      
                                      
#===============================================================================================================================
                                      
rm.allrowNA2 <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) if(NROW(i) != 0) Filter(NROW, i[rowSums(is.na(i) | i == "") != ncol(i), , drop = FALSE]) else Filter(NROW, i))
    
  } else { X[rowSums(is.na(X) | X == "") != ncol(X), , drop = FALSE] }
}                                      
  
#===============================================================================================================================
           
d_prepo <- function(data = NULL)
{
  
  if(is.null(reget(data, control))) stop("Required 'control' group not found.", call. = FALSE)
  
  ar <- formalArgs(d.prepos)[-c(21, 22)]
  
  dot.names <- names(data)[!names(data) %in% ar]
  
  args <- unclass(data[c(ar, dot.names)])
  
  argsT <- setNames(args, names(args))
  
  do.call(d.prepos, argsT)
}

#===============================================================================================================================

test.sheet <- function(data, metaset = FALSE){
  
  data <- rm.allrowNA(trim(data))
  
  check <- "study.name" %in% names(data)
  if(!check) stop("Add a new column named 'study.name'.", call. = FALSE)
  
  L <- split(data, data$study.name)         
  
  f <- function(number){
    
    ns <- names(L)[number]
    
    z <- try(d_prepo(L[[number]]), silent = TRUE)
    
    if(inherits(z, "try-error")) message("Error: pre-post coding problem in: *", toString(dQuote(ns)), "*") else message(if(metaset)"" else "Ok: ","No pre-post coding problem in ", toString(dQuote(ns)))
  }
  invisible(lapply(seq_along(L), function(i) f(i)))
}     
 
#===============================================================================================================================
                   
dint.norm <- function(dint) noquote(paste0(round(pnorm(dint) - pnorm(0), 4)*1e2, "%"))   
                   
#===============================================================================================================================                   
                   
do.factor <- function(data, exclude = NULL, char = TRUE, drop = NULL){
  
  colnames(data) <- trimws(colnames(data))
  
  if(!is.null(drop)) data <- drop.col(data, drop)
  
  data <- rm.allrowNA(data) 
  
  ar <- c(formalArgs(d.prepos)[-(20:22)], c("SD", "dint", "id"), exclude)
  
  dot.names <- names(data)[!names(data) %in% ar]
  
  data[dot.names] <- lapply(data[dot.names], if(char) as.character else as.factor)
  
  return(data)
}                  
                                  
#===============================================================================================================================

meta.in <- function(data = NULL, by, impute = FALSE, n.sim = 1e5, option = 1, r = .5){ 
  
  L <- eval(substitute(dintA(data = data, by = by, impute = impute, n.sim = n.sim)))
  
  ds <- Filter(Negate(is.null), lapply(L, function(x) do.call(rbind, x)$dint))
  sds <- Filter(Negate(is.null), lapply(L, function(x) do.call(rbind, x)$SD))   
  
  f <- if(option == 1) option1 else option2
  
  mapply(f, ds = ds, sds = sds, r = r, SIMPLIFY = FALSE)
}

#========================================================================================

meta.out <- function(data = NULL, by, impute = FALSE, n.sim = 1e5, option = 1, r = .5, mu.prior = mu.norm(-6, 6), tau.prior = function(x){dhalfnormal(x)}){  
  
  j <- eval(substitute(meta.in(data = data, by = by, impute = impute, n.sim = n.sim, option = option, r = r)))
  
  ds <- sapply(seq_along(j), function(i) j[[i]][1])
  sds <- sapply(seq_along(j), function(i) j[[i]][2])
  
  test <- length(ds) >= 2
  
  if(!test) return(NA)
  
  res <- bayesmeta(                y = ds,
                                   sigma = sds,
                                   labels = names(j), 
                                   tau.prior = tau.prior,
                                   mu.prior = mu.prior)
  res$call <- match.call(expand.dots = FALSE)
  
  return(res)
}                                         

#===============================================================================================================================================================================================================================
                
                
forest.rob <- function(x, zoom, xlab = "effect size (dint)", refline = NULL, cex, level = .95, col = NULL, main = NA, order.by = FALSE, col.by.cluster = TRUE, refit = FALSE, wsize = 1, slab = TRUE, summary = TRUE, ...)
{
    
  order.f <- if(order.by == "weight" || order.by == "effect") TRUE else FALSE
  
  check <- x$ml[[3]] != 1
  
  if(is.null(refline)) refline <- x$reg_table$b.r[[1]]
  
  if(check) message("Note: Overall effect displayed for intercept-only models.")
  
  mis <- missing(zoom)
  d <- cbind(x$data.full, x$data, orig.nm = as.vector(x$study_orig_id))
  s <- substitute(zoom)
  if(!mis) d <- subset(d, eval(s))
  d <- if(order.by == "effect") d[order(d$effect.size, decreasing = TRUE), ] else if(order.by == "weight") d[order(d$r.weights, decreasing = TRUE), ] else d
  grp <- d$study
  
  if(slab){
  s <- if(order.f & mis) unique(d$orig.nm)[grp] else d$orig.nm
  elab <- if(order.f & mis) s[seq_len(length(s))] else s
  }
  
  cols <- rep(c(1:2, "green4", 4, "magenta3", "gold4", "darkred", "blue3"), nrow(d))
  
  col <- if(!is.null(col)) col else
    if(order.f || is.null(col) & !col.by.cluster || order.f & col.by.cluster) 1 
  else cols[grp]
  
  y <- d$effect.size
  vi <- d$var.eff.size
  w <- d$r.weights
  
  m <- if(!refit) x else { fo <- formula(x$ml) 
  h <- robu(fo, var = d$var.eff.size, study = d$orig.nm, small = x$small, model = x$modelweights, data = d, rho = x$mod_info$rho) 
  h$ml <- fo
  h }
  
  if(refit & refline == x$reg_table$b.r[[1]]) refline <- m$reg_table$b.r[[1]]
  
  f <- forest.default(x = y, vi = vi, psize = wsize*(w/max(w)),
                      level = level,           
                      refline = refline,
                      xlab = xlab,
                      slab = if(order.f & slab) elab else NA,
                      cex = cex, efac = 0, col = col, mgp = c(1, .3, 0), ...)
  
  ES <- m$reg_table$b.r[[1]]
  ES.CI.L <- m$reg_table$CI.L[[1]]
  ES.CI.U <- m$reg_table$CI.U[[1]]
  
  rows <- f$rows
  
  if(!order.f & slab){
    
    grp <- rev(which(!duplicated(grp)))
    
    text(f$xlim[1], rev(rows)[-grp], elab[-grp], pos = 4, cex = f$cex)
    
    text(f$xlim[1], rev(rows)[grp], elab[grp], pos = 4, cex = f$cex, col = "red4", font = 4)
  }
  
  mtext(text = main, font = 2, line = -2)
  
  if(!check & summary) addpoly.default(ES, ci.lb = ES.CI.L, ci.ub = ES.CI.U, mlab = expression(bold("mean effect ("*mu*")")), 
                             level = level, cex = f$cex, col = "cyan", rows = par('usr')[3], font = 2, xpd = NA, border = "magenta")
  
  abline(h = max(rows)+1, lwd = 1, col = 0, xpd = NA)
  
  if(refit) return(m)
}

                      
#========================================================================================                      
                      
forest.bayes1 <- function(x, xlab = "effect size", refline = x$summary["median","mu"], cex = NULL, ...){

  f <- forest.default(x = x$y, sei = x$sigma,
                          showweight = FALSE,  
                          ylim = c(-x$k-4, 1),
                          level = 95,         
                          refline = refline,
                          xlab = xlab,
                          slab = x$labels,
                          rows = seq(-2, -x$k - 1, by = -1),
                          cex = cex, ...)
  
  if(is.null(cex)) cex <- f$cex
  
  abline(h = max(f$rows)+1, lwd = 1, col = 0, xpd = NA)
  
  addpoly(x$summary["median","mu"], ci.lb=x$summary["95% lower","mu"], ci.ub=x$summary["95% upper","mu"],
                   rows = -x$k-3, mlab= "mean effect", level=95, cex=cex, xpd = NA, col = "cyan", border = "magenta", font = 2, ...)
  
  addpoly(x$summary["median","theta"], ci.lb=x$summary["95% lower","theta"], ci.ub=x$summary["95% upper","theta"],
                   rows = -x$k-4.5, mlab= "prediction", level=95, cex=cex, xpd = NA, col = 8, font = 2, ...)
}  
    
                      
#========================================================================================

                      
forest.bayes <- function(x, zoom, xlab = "effect size", refline = x$summary["median","mu"], cex = NULL, order.by = FALSE, refit = FALSE, ...){
  
  
  mis <- missing(zoom)
  d <- data.frame(dint = x$y, sei = x$sigma, labels = x$labels)
  s <- substitute(zoom)
  if(!mis) d <- subset(d, eval(s))
  
  if(!mis & refit) { 
    
    x <- bayesmeta::bayesmeta(y = d$dint, sigma = d$sei, labels = d$labels)
    if(is.null(refline)) refline <- x$summary["median","mu"]
  }
    
  d <- if(order.by) d[order(d$dint, decreasing = TRUE), ] else d
  k <- nrow(d)
  
  f <- forest.default(x = d$dint, sei = d$sei,
                      showweight = FALSE,  
                      ylim = c(-k-4, 1),
                      level = 95,         
                      refline = refline,
                      xlab = xlab,
                      slab = d$labels,
                      rows = seq(-2, -k - 1, by = -1),
                      cex = cex, ...)
  
  if(is.null(cex)) cex <- f$cex
  
  abline(h = max(f$rows)+1, lwd = 1, col = 0, xpd = NA)
  
  addpoly.default(x$summary["median","mu"], ci.lb=x$summary["95% lower","mu"], ci.ub=x$summary["95% upper","mu"],
          rows = -k-3, mlab= expression(bold("mean effect ("*mu*")")), level=95, cex=cex, xpd = NA, col = "cyan", border = "magenta", font = 2, ...)
  
  addpoly.default(x$summary["median","theta"], ci.lb=x$summary["95% lower","theta"], ci.ub=x$summary["95% upper","theta"],
          rows = -k-4.5, mlab= expression(bold("prediction ("*mu*")")), level=95, cex=cex, xpd = NA, col = 8, font = 2, ...)
}  
                      
                      
#========================================================================================

forest.dint <- function(x, zoom, xlab = "Effect Size (dint)", refline = NULL, cex = NULL, level = .95, col = NULL, col.by.cluster = FALSE,  refit = FALSE, order.by = FALSE, wsize = 1, space = TRUE, slab = TRUE, summary = TRUE, reset = TRUE, ...){
  
  
  par.mgp <- par("mgp")
  par(mgp = c(1.8, .3, 0))  
  on.exit(par(mgp = par.mgp))
  
  if(space){
    par.mar <- par("mar")
    par(mar = c(2.7, 3, 0, 1))
    if(reset){
      on.exit(par(mar = par.mar))
    }
  }
  
  if(inherits(x, "robu")) {  
    
    eval(substitute(forest.rob(x = x, zoom = zoom, xlab = xlab, refline = refline, cex = cex, level = level, order.by = order.by, col.by.cluster = col.by.cluster, col = col, refit = refit, wsize = wsize, slab = slab, summary = summary, ...)))
    
  }else{
    
    if(is.null(refline)) refline <- x$summary["median","mu"]
    
    eval(substitute(forest.bayes(x = x, zoom = zoom, xlab = xlab, refline = refline, order.by = order.by, cex = cex, refit = refit, ...)))
  }
}   
                      
#========================================================================================
                
                
funnel.rob <- function(x, zoom, xlab = "effect size (dint)", ylab = "SD", refline = x$reg_table$b.r[[1]],
                       cex = 1, level = .95, col = "magenta", main = deparse(substitute(x)),
                       back = 8, shade = 0, hlines = NA,
                       pch = 21, bg = "cyan", refit = FALSE, ...){
  
  check <- x$ml[[3]] != 1
  if(check) message("Note: Overall effect line only used for intercept-only models.")
  
  d <- cbind(x$data.full, x$data, orig.nm = as.vector(x$study_orig_id))
  s <- substitute(zoom)
  if(!missing(zoom)) d <- subset(d, eval(s))
  
  y <- d$effect.size
  vi <- d$var.eff.size
  
  m <- if(!refit) x else { fo <- formula(x$ml)
  h <- robu(fo, var = d$var.eff.size, study = d$orig.nm, small = x$small, model = x$modelweights, data = d, rho = x$mod_info$rho)
  h$ml <- fo
  h }
  
  if(refit & refline == x$reg_table$b.r[[1]]) refline <- m$reg_table$b.r[[1]]
  
  if(refline == 0 & level == .95 & shade == 0){level <- c(.95, .99) ; shade = c("white", "orange")}
  
  funnel.default(x = y,
                 vi = vi,
                 level = level,          
                 refline = refline,
                 xlab = xlab,
                 cex = cex, col = col, ylab = ylab,
                 back = back, shade = shade, hlines = hlines,
                 pch = pch, bg = bg, mgp = c(1.7, .5, 0), ...)
  box()
  
  ref <- if(refline == 0) bquote(H[0]*":"*~mu == .(0)) else bquote(mu == .(round(refline, 3)))
  
  mtext(main, font = 2, ...)
  
  text(refline, 0, ref, col = "magenta", font = 2, pos = 4)
  
  if(refit) return(m)
}               
            
#========================================================================================

funnel.dint <- function(x, zoom, xlab = "Effect Size (dint)", ylab = "SD", refline = x$reg_table$b.r[[1]], 
                        cex = 1, level = .95, col = "magenta", main = deparse(substitute(x)),
                        back = 8, shade = 0, hlines = NA, 
                        pch = 21, bg = "cyan", study.name = TRUE,
                        FE = FALSE, legend = FALSE, shrink = FALSE, show.mu = TRUE,
                        refit = FALSE, ...){
  
  if(inherits(x, "robu")) { eval(substitute(funnel.rob(x = x, zoom = zoom, level = level,           
                                                       refline = refline,
                                                       xlab = xlab, main = main,
                                                       cex = cex, col = col, ylab = ylab, 
                                                       back = back, shade = shade, hlines = hlines, 
                                                       pch = pch, bg = bg, refit = refit, ...)))
    
  } else { funnel.bayesmeta(x = x, main = main, study.name = study.name, FE = FE, 
                  legend = legend, shrink = shrink, show.mu = show.mu, ...) 
  }
}
                      
#==============================================================================================
                      
                      
meta.bayes <- function(data = NULL, by, option = 1, r = .5, mu.prior = mu.norm(-6, 6), tau.prior = function(x){dhalfnormal(x)}){  
  
  data$study.name <- trimws(data$study.name)
  
  metain <- function(data = NULL, by, option = 1, r = .5){ 
    
    m <- split(data, data$study.name)
    m <- Filter(NROW, rm.allrowNA2(m)) 
    
    L <- if(missing(by)) m else { s <- substitute(by) ; h <- lapply(m, function(x) do.call("subset", list(x, s))) ;
    res <- Filter(NROW, h) ; if(length(res) == 0) NULL else res}
    
    ds <- Filter(Negate(is.null), lapply(seq_along(L), function(i) L[[i]]$dint))
    sds <- Filter(Negate(is.null), lapply(seq_along(L), function(i) L[[i]]$SD))
    
    f <- if(option == 1) option1 else option2
    
    setNames(mapply(f, ds = ds, sds = sds, r = r, SIMPLIFY = FALSE), names(L))
  }
  
  j <- eval(substitute(metain(data = data, by = by, option = option, r = r)))
  
  ds <-  sapply(seq_along(j), function(i) j[[i]][1])
  sds <- sapply(seq_along(j), function(i) j[[i]][2])
  
  test <- length(ds) >= 2
  
  if(!test) return(NA)
  
  res <- bayesmeta(                y = ds,
                                   sigma = sds,
                                   labels = names(j), 
                                   tau.prior = tau.prior,
                                   mu.prior = mu.prior)
  res$call <- match.call(expand.dots = FALSE)
  
  return(res)
}                
                
                
#==============================================================================================================================================   
                
 find.norms <- function(low, high, cover = .99, digits = 6){
   
   f <- Vectorize(mu.norm)
   data.frame(t(f(low = low, high = high, cover = cover, digits = digits)))
 }              
                
#====================================================================================================
                
mu.norm <- find.norm <- function(low, high, cover = .99, digits = 6){
  
  options(warn = -1)
  
  cover[cover >= 1] <- .999999999999999
  cover[cover <= 0] <- .000000000000001
  
  p1 <- (1 - cover) / 2 
  p2 <- 1 - p1
  
  q <- c(low, high)  
  alpha <- c(p1, p2)
  
  is.df <- function(a, b, sig = 4) (round(a, sig) != round(b, sig))
  
  if (p1 <= 0 || p2 >= 1 || q[1] >= q[2] || p1 >= p2) {
    
    stop("Incorrect 'low' and/or 'high' or 'cover' values.", call. = FALSE)
    
  } else {
    
    beta <- qnorm(alpha)
    
    parm <- solve(cbind(1, beta), q)
    
    q <- qnorm(c(p1, p2), parm[[1]], parm[[2]])
  }
  
  if(is.df(low, q[[1]]) || is.df(high, q[[2]])) {
    
    stop("Change 'low' and/or 'high' or 'cover' values.", call. = FALSE)
    
  } else {
    
    return(round(c(mean = parm[[1]], sd = parm[[2]]), digits = digits))
  }
}                               
  
                
#========================================================================================
                
meta.stats <- function(..., stat = "median", digits = 4){
  
  m <- list(...)
  cl <- class(...)
  len <- length(m)
  m <- if(cl == "list" & len == 1) as.list(...) else m
  n <- substitute(...())
  stat <- trimws(stat)
  
  f <- function(fit, stat = "median"){
    
    if(inherits(fit, "bayesmeta")){
      
      oo <- round(c(K = fit$k, mu = fit$summary["mean","mu"],low = fit$summary["95% lower","mu"], up = fit$summary["95% upper","mu"], I2 = fit$I2(fit$summary[stat,1]), tau = fit$summary[stat,1], BF01.mu = fit$bayesfactor[1,2], BF01.tau = fit$bayesfactor[1,1]),digits)
      
      c(oo, perc.mu = dint.norm(oo[2]))
      
    } else if(inherits(fit, "robu")){
      
      oo <- round(c(K = fit$N, mu = fit$reg_table$b.r[[1]], low = fit$reg_table$CI.L[[1]], up = fit$reg_table$CI.U[[1]], I2 = as.vector(fit$mod_info$I.2), tau = sqrt(as.vector(fit$mod_info$tau.sq))), digits)
      
      c(oo, perc.mu = dint.norm(oo[2]))
      
    } else {
      
      c(I2 = fit$I2, tau = sqrt(fit$tau2))
    }
  }
  m <- m[!is.na(m)]
  out <- setNames(lapply(m, f, stat = stat), if(cl == "list" & len == 1) names(m) else n)
  d <- data.frame(out)
  d[] <- lapply(d, as.list)
  data.frame(t(d))
}                              
                
#================================================================================================================================================================
              
metal.dint <- function(data = NULL, by, over = time, mu.prior = mu.norm(-6, 6), tau.prior = function(x){dhalfnormal(x)}, 
                       option = 1, r = .5, method = c("robust", "bayes")){
  
  over <- deparse(substitute(over))
  
  method <- match.arg(method)
  
  chep <- sort(unique(na.omit(data[[over]])))
  
  G <- if(missing(by)) { lapply(chep, function(y) bquote(.(as.name(noquote(over))) == .(y))) 
    
  } else {
    
    s <- substitute(by)
    lapply(chep, function(x) bquote(.(s) & .(as.name(noquote(over))) == .(x)))
  }
  
  if(method == "robust"){
    
    dd <- Filter(length, lapply(seq_along(G), function(j) subset(data, eval(G[[j]]))))
    
    if(length(dd) == 0) return(NA)
    
    f <- try(setNames(lapply(seq_along(dd), function(x) robu(dint~1, data = dd[[x]], studynum = as.vector(study.name), var = SD^2)), paste(over, chep)), silent = TRUE)
    if(inherits(f, "try-error")) NA else f[sapply(seq_along(f), function(i) is.finite(f[[i]]$reg_table$dfs[1]) & f[[i]]$reg_table$dfs[1] > 1)]
  } 
  
  else 
    
  {
    
    data$study.name <- trimws(data$study.name)
    
    f1 <- function(data = NULL, zy, option = 1, r = .5){ 
      
      data <- if(missing(zy)) data else { s <- substitute(zy) ; subset(data, eval(s)) }
      
      m <- split(data, data$study.name)
      L <- Filter(NROW, rm.allrowNA2(m)) 
      
      ds <- Filter(Negate(is.null), lapply(seq_along(L), function(i) L[[i]]$dint))
      sds <- Filter(Negate(is.null), lapply(seq_along(L), function(i) L[[i]]$SD))
      
      f <- if(option == 1) option1 else option2
      
      setNames(mapply(f, ds = ds, sds = sds, r = r, SIMPLIFY = FALSE), names(L))
    }
    
    f2 <- function(j, tau.prior, mu.prior){  
      
      ds <- sapply(seq_along(j), function(i) j[[i]][1])
      sds <- sapply(seq_along(j), function(i) j[[i]][2])
      
      test <- length(ds) >= 2
      
      if(!test) return(NA)
      
      res <- bayesmeta(        y = ds,
                               sigma = sds,
                               labels = names(j), 
                               tau.prior = tau.prior,
                               mu.prior = mu.prior)
      res$call <- match.call(expand.dots = FALSE)
      
      return(res)
    }  
    
    go <- length(G)
    
    k <- vector("list", go)
    
    for(w in seq_len(go)) k[[w]] <- f1(data = data, zy = eval(G[[w]]), option = option, r = r)
    
    so <- length(k)
    
    z <- vector("list", so)
    
    for(a in seq_len(so)) z[[a]] <- f2(j = k[[a]], tau.prior = tau.prior, mu.prior = mu.prior)
    
    setNames(z, paste(over, chep))
  }
}
           
#================================================================================================================================================================           
           
best.model <- function(mod.names, data, n.best = 50, small = FALSE, model = c("CORR", "HIER"), rho = .8){
   
   f <- unlist(lapply(seq_along(mod.names), function(n) combn(mod.names, n, FUN = function(i) paste0("dint~", paste0(i, collapse = "+")))))
   
   res <- sapply(f, function(j) as.double(robu2(formula(j), data = data, study = study.name, var = SD^2, small = small, model = model, rho = rho)$mod_info$tau.sq))
   
   len <- length(res)
   
   if(n.best > len) n.best <- len
   
   data.frame(Tau.sq = sort(res)[seq_len(n.best)])
 }           

#================================================================================================================================================================
                 
tplot <- function(y, main, lwd = 4, lend = 2, cat.level = 0, low = FALSE){
  
  if(!low) low <- NULL
  
  z <- length(y) 
  x <- seq_len(z)
  
  if(cat.level != 0 & z >= cat.level) { main <- bquote(bold(.(main)~symbol(("\326")))) ; col.main <- "magenta"} else { main <- main ; col.main <- 1} 
  plot(x, y, type = "h", main = main, xlim = c(.95, 1.02*max(x)),
       ylab = "Frequency", axes = FALSE, xlab = "Category", lwd = lwd,
       col = colorRampPalette(c(4, 2))(z), font.lab = 2, lend = lend, col.main = col.main)
  box()
  axis(1, at = which(!names(y) %in% names(low)), labels = names(y)[!names(y) %in% names(low)], cex.axis = .9)
  if(!is.null(low)) axis(1, at = which(names(y) %in% names(low)), labels = names(low), cex.axis = .9, col = "magenta", col.axis = "magenta", font = 2)
  if(!is.null(low)) text(which(names(y) %in% names(low)), low, low, font = 2, cex = .9, pos = 3, col = "magenta")
  axis(2, at = pretty(y), cex.axis = .85, las = 1, padj = .3)
}

#================================================================================================================================================================
                 
plot.mods <- function(data, exclude = NULL, lwd = 4, lend = 2, cat.level = 0, code = NULL, low = NULL){
  
  lo <- low
  low <- if(is.null(low)) FALSE else TRUE
  
  names(data) <- trimws(names(data))
  
  data <- rm.colrowNA(data) 
  
  ar <- c(formalArgs(d.prepos)[-(20:22)], c("SD", "dint", "id", "study.name"), exclude)
  
  mods <- names(data)[!names(data) %in% ar]
  
  A <- setNames(lapply(seq_along(mods), function(i) table(data[[mods[i]]], dnn = NULL)), mods)
  Ls <- lapply(A, length)
  
  bad <- Ls < 2
  bad.names <- names(A[bad])
  A <- A[!bad]
  
  A <- A[Ls >= cat.level]
  if(length(A) == 0) stop(paste("No variable with cat.level >=", if(cat.level != 0) cat.level else 2, "found."), call. = FALSE)  
  
  if(!is.null(code)){
    target <- sapply(seq_along(A), function(i) any(names(A[[i]]) == code))
    A <- A[target]
  }
  
  if(low){
    target <- sapply(seq_along(A), function(i) any(A[[i]] <= lo))
    A <- A[target]
    low <- if(length(A) == 0) FALSE else setNames(lapply(seq_along(A), function(i) A[[i]][which(A[[i]] <= lo)]), names(A))
  }
  
  n <- length(A)
  graphics.off()
  org.par <- par(no.readonly = TRUE)
  on.exit(par(org.par))
  if(n > 1L) { par(mfrow = n2mfrow(n)) ; set.margin() }
  
  if(length(bad.names) != 0) message("Note: ", toString(dQuote(bad.names)), "ignored due to insufficient category levels.")
  if(length(A) > 0) { invisible(mapply(tplot, y = A, main = names(A), lwd = lwd, lend = lend, cat.level = cat.level, low = low)) ; return(A) } else NA
}

#================================================================================================================================================================
                                             
study.level1 <- function(data, exclude = NULL){
  
  ar <- c(formalArgs(d.prepos)[-c(20:22)], c("SD", "dint", "id"), exclude)  
  
  mods <- names(data)[!names(data) %in% ar]
  
  tmp <- do.call(rbind, lapply(mods, function(x){
    d <- setNames(unique(data[c("study.name", x)]), c("study.name", "code"))
    transform(d, mod.name = x)
  }))
  
  res <- tmp[with(tmp, ave(code, code, mod.name, FUN = length) == 1),]
  
  h <- Filter(length, lapply(split(res, res$study.name, drop = TRUE), `row.names<-`, NULL))
  
  mix <- if(length(h) == 0) NA else h
  grp <- group.level(data = data, exclude = exclude)
  
  if(!is.na(grp) & !is.na(mix)) { 
    
    mix[names(grp)] <- Map(function(x, y) {y1 <- stack(y); 
    subset(x, !(mod.name %in% y1$ind & code %in% y1$values))}, mix[names(grp)], grp)
    h <- Filter(nrow, mix) 
    if(length(h) == 0) NA else h
    
  } else NA
} 
                       
#================================================================================================================================================================
                       
mix.level1 <- function(data, exclude = NULL){
  
  ar <- c(formalArgs(d.prepos)[-c(20:22)], c("SD", "dint", "id"), exclude)  
  
  mods <- names(data)[!names(data) %in% ar]
  
  tmp <- do.call(rbind, lapply(mods, function(x){
    d <- setNames(unique(data[c("study.name", x)]), c("study.name", "code"))
    transform(d, mod.name = x)
  }))
  
  res <- tmp[with(tmp, ave(code, code, mod.name, FUN = length) == 1),]
  
  h <- Filter(length, lapply(split(res, res$study.name, drop = TRUE), `row.names<-`, NULL))
  
  if(length(h) == 0) NA else h
}    

#================================================================================================================================================================
                       
mix.level <- function(data, exclude = NULL){
  
  ar <- c(formalArgs(d.prepos)[-c(20:22)], c("SD", "dint", "id"), exclude)  
  
  mods <- names(data)[!names(data) %in% ar]
  
  tmp <- do.call(rbind, lapply(mods, function(x){
    d <- setNames(unique(data[c("study.name", x)]), c("study.name", "code"))
    transform(d, mod.name = x)
  }))
  
  res <- tmp[with(tmp, ave(code, code, mod.name, FUN = length) == 1),]
  if(nrow(res) == 0) return(NULL)
  res <- res[order(res$study.name),]
  rownames(res) <- NULL
  res
}                           
                       
#================================================================================================================================================================
   
group.level1 <- function(data, exclude = NULL){

  ar <- c(formalArgs(d.prepos)[-c(2,20:22)], c("SD", "dint", "id"), exclude)
  
  d <- drop.col(data, ar)
  
  m <- split(d, d$study.name)
  m <- Filter(NROW, rm.allrowNA2(m))
  
  input.order <- unsplit(m, d$study.name)
  
  molten <- data.frame(input.order[, 1, drop = FALSE], stack(input.order[, -1]))
  
  res <- molten[as.logical(ave(molten[['values']], molten[['ind']], 
                               FUN = function(x) !duplicated(x) & !duplicated(x, fromLast = TRUE))), ]
  vec <- res$values
  names(vec) <- res$ind
  
  h <- Filter(NROW, split(vec, as.character(res$study.name)))
  
  if(length(h) == 0) NA else h
}                                   

                        
#================================================================================================================================================================
                               
group.level <- function(data, exclude = NULL){
  
  ar <- c(formalArgs(d.prepos)[-c(2,20:22)], c("SD", "dint", "id"), exclude)
  
  d <- drop.col(data, ar)
  
  molten <- data.frame(d[, 1, drop = FALSE], stack(d[, -1]))
  
  res <- molten[as.logical(ave(molten[['values']], molten[['ind']], 
                               FUN = function(x) !duplicated(x) & !duplicated(x, fromLast = TRUE))), ]
  
  if(nrow(res) == 0) return(NULL)
  res <- setNames(res[order(res$study.name),], c("study.name", "code", "mod.name"))
  rownames(res) <- NULL
  res
}                               

#================================================================================================================================================================
                               
study.level <- function(data, exclude = NULL){
  
  mix <- mix.level(data = data, exclude = exclude)
  grp <- group.level(data = data, exclude = exclude)
  
  if(!is.null(grp) & !is.null(mix)) { 
    
    res <- subset(mix, !((study.name %in% grp$study.name) & (code %in% grp$code) & (mod.name %in% grp$mod.name)))
    res <- res[order(res$study.name),]    
    rownames(res) <- NULL
    res
  } else if(is.null(grp) & !is.null(mix)){ mix } else { NULL }
}                                                                              

#================================================================================================================================================================
                               
print.labdf <- function(data) {
  dn <- dimnames(data)
  names(dn) <- attr(data, "rclab")
  data <- as.matrix(data)
  dimnames(data) <- dn
  print(data, quote = FALSE)
}                               
                               
#================================================================================================================================================================
                               
exam.code2 <- function(data, exclude = NULL, rule = 1, lwd = 4, lend = 2, cat.level = 6){
  
  names(data) <- trimws(names(data))
  check <- "study.name" %in% names(data)
  if(!check) stop("Add a new column named 'study.name'.", call. = FALSE)
  
  data$study.name <- trimws(data$study.name)
  data <- rm.allrowNA(data)
  z <- length(unique(data$study.name))
  if(z < 2) stop("At least two coded studies required.", call. = FALSE)
  
  m <- split(data, data$study.name)         
  m <- Filter(NROW, rm.allrowNA2(m)) 
  if(z != length(m)) stop("Each 'study.name' must be distinct.", call. = FALSE)
  
  exclude <- trimws(exclude)  
  excl <- setdiff(exclude, "study.name")
  
  exclude <- if(!is.null(excl) & length(excl) != 0) exclude else NULL
  
  h <- if(rule == 1) study.level(data = data, exclude = exclude) else
    if(rule == 2) group.level(data = data, exclude = exclude) else 
      plot.mods(data = data, exclude = exclude, lwd = lwd, lend = lend, cat.level = cat.level)
  
  if(rule == 1 & !is.null(h) || rule == 2 & !is.null(h)){    
    attr(h, "rclab") <- c("", paste0("Violations of Rule ", rule, ":"))
    class(h) <- c("labdf", class(h))
    h } else { h }
}   

#================================================================================================================================================================
    
suggest <- function(data, exclude = NULL){

  ar <- c(formalArgs(d.prepos)[-c(2,20:22)], c("SD", "dint", "id"), exclude)
  data <- drop.col(data, ar)
  
DF <- data[with(data, ave(as.numeric(study.name), study.name, FUN = length)) >= 3, ]
long.df <- cbind(DF[1], stack(DF[-1]))
res <- long.df[with(long.df, ave(values, study.name, ind, values, FUN = length) == 1 &
               ave(values, study.name, ind, FUN = function(x) length(unique(x))) == 2), ]

if(nrow(res) == 0) return(NULL)
res <- setNames(res[order(res$study.name),], c("study.name", "code", "mod.name"))   
rownames(res) <- NULL
res
}    
    
#================================================================================================================================================================
    
exam.code3 <- function(data, exclude = NULL, rule = 1, lwd = 4, lend = 2, cat.level = 0, code = NULL, low = 4, suggest = FALSE){
  
  names(data) <- trimws(names(data))
  check <- "study.name" %in% names(data)
  if(!check) stop("Add a new column named 'study.name'.", call. = FALSE)
  
  data$study.name <- trimws(data$study.name)
  data <- rm.colrowNA(data)
  
  if(length(unique(data$study.name)) < 2) stop("At least two coded studies required.", call. = FALSE)
  
  exclude <- trimws(exclude)  
  excl <- setdiff(exclude, "study.name")
  
  exclude <- if(!is.null(excl) & length(excl) != 0) exclude else NULL
  
h <- if(rule == 1 & !suggest) study.level(data = data, exclude = exclude) else
      if(rule == 2 & !suggest) group.level(data = data, exclude = exclude) else 
       if(rule == 3 & !suggest) plot.mods(data = data, exclude = exclude, lwd = lwd, lend = lend, cat.level = cat.level, code = code, low = low) else
        if(suggest) suggest(data = data, exclude = exclude)
  
  if(rule == 1 & !is.null(h) || rule == 2 & !is.null(h) || suggest & !is.null(h)){    
    
    attr(h, "rclab") <- c("", paste0(if(suggest) "Possible Mistakes" else paste("Violations of Rule", rule), ":"))
    class(h) <- c("labdf", class(h))
    h } else { h }
}      
 
#================================================================================================================================================================                                     
 
rule3 <- function(data, exclude = NULL, low = 4){
  
ar <- c(formalArgs(d.prepos)[-c(20:22)], c("SD", "dint", "id"), exclude)  

mods <- names(data)[!names(data) %in% ar]

A <- setNames(lapply(seq_along(mods), function(i) table(data[[mods[i]]], dnn = NULL)), mods)

target <- sapply(seq_along(A), function(i) any(A[[i]] <= low))
A <- A[target]

if(length(A) == 0) return(NULL)

lo <- setNames(lapply(seq_along(A), function(i) A[[i]][which(A[[i]] <= low)]), names(A))

lst <- Filter(length, lapply(split(data[names(lo)], data$study.name), 
                              function(dat) Filter(nrow, Map(function(x, y) 
                                merge(x, y[setdiff(names(y), "values")], by = "ind"), lapply(dat, 
                                function(x) stack(table(x))), lapply(lo, stack)))))

res <- do.call(rbind, c(Map(cbind, study.name = names(lst), lapply(lst, 
                 function(x) do.call(rbind, c(Map(cbind, x, mod.name = names(x)),
                 make.row.names = FALSE)))), make.row.names = FALSE))

r3 <- setNames(res, c("study.name","code","occurs","mod.name"))[, c(1,2,4,3)]
r3 <- r3[order(r3$mod.name),]    
rownames(r3) <- NULL

grp <- group.level(data = data, exclude = exclude)

if(is.null(grp)){ r3 } else {
res <- subset(r3, !((study.name %in% grp$study.name) & (code %in% grp$code) & (mod.name %in% grp$mod.name)))
res <- res[order(res$mod.name),]    
rownames(res) <- NULL
res
  } 
}                                     

#================================================================================================================================================================
                                                                   
exam.code <- function(data, exclude = NULL, rule = 1, lwd = 4, lend = 2, cat.level = 0, code = NULL, low = NULL, suggest = FALSE, plot = TRUE){
   
   data <- rm.colrowNA(trim(data))
   check <- "study.name" %in% names(data)
   if(!check) stop("Add a new column named 'study.name'.", call. = FALSE)
   
   if(length(unique(data$study.name)) < 2) stop("At least two coded studies required.", call. = FALSE)
   
   message("
Rule 1: An efficient coding sheet is one in which at least '2 studies' coded the same category of a moderator (ex. 2 studies with moderator 'X' coded 'y').

Rule 2: An efficient coding sheet is one in which at least '2 study groups' coded the same category of a moderator (ex. 2 groups with moderator 'X' coded 'y').

Rule 3: An efficient coding sheet is one in which there is no moderator whose categories occur less than ~4 times overall.
") 
    
   exclude <- trimws(exclude)  
   excl <- setdiff(exclude, "study.name")
   
   exclude <- if(!is.null(excl) & length(excl) != 0) exclude else NULL
   
   if(plot) invisible(plot.mods(data = data, exclude = exclude, lwd = lwd, lend = lend, cat.level = cat.level, code = code, low = low))
   
   h <- if(rule == 1 & !suggest) study.level(data = data, exclude = exclude) else
     if(rule == 2 & !suggest) group.level(data = data, exclude = exclude) else 
       if(rule == 3 & !suggest) rule3(data = data, exclude = exclude, low = if(is.null(low)) 4 else low) else
         if(suggest) suggest(data = data, exclude = exclude)
   
   if(!is.null(h)){    
     
     attr(h, "rclab") <- c("", paste0(if(suggest) "Possible Mistakes" else paste("Violations of Rule", rule), ":"))
     class(h) <- c("labdf", class(h)) 
      h } else { h }
 }                                                                           
  
                                      
#================================================================================================================================================================
                                      
robu.weight <- function(...){
  
  if(!all(sapply(list(...), inherits, "robu"))) stop("Non-robust variance model(s) detected.", call. = FALSE)
  
  m <- list(...)
  n <- substitute(...())
  
  f <- function(fit){
    fit$data.full$r.weights
  }
  setNames(lapply(m, f), n)
}                                     

#================================================================================================================================================================
                                      
maxs <- function(x, n = 2){
   
   len <- length(x)
   if(n >len){
     warning('n greater than length(x). Setting n = length(x)', call. = FALSE)
     n <- length(x)
   }
   sort(x, partial = len-n+1)[len-n+1]
 }
 
#================================================================================================================================================================
                                      
 mins <- function(x, n = 2){
   
   len <- length(x)
   if(n >len){
     warning('n greater than length(x). Setting n = length(x)', call. = FALSE)
     n <- length(x)
   }
   sort(x)[n]
 }                                      

                                      
#========================================================================================================================================================
                                      
outlier <- function(data, n = 5){
  
  data$study.name <- as.vector(data$study.name)  
  
  len <- seq_len(n)  
  
  up <- maxs(data$dint, len)
  lo <- mins(data$dint, len)
  
  one <- data.frame(t(sapply(len, function(i) {h <- subset(data, dint == eval(up[i])); data.frame(row = as.numeric(rownames(h)), h)})))[,c(2,1,3)]
  two <- data.frame(t(sapply(len, function(i) {h <- subset(data, dint == eval(lo[i])); data.frame(row = as.numeric(rownames(h)), h)})))[,c(2,1,3)]
  cbind(one, two)
  
}                                      

                           
#================================================================================================================================================================
                           
egger <- function(...){
  
  m <- list(...)
  if(!all(sapply(m, inherits, c("robu", "bayesmeta")))) stop("Unsupported meta-analytic model(s) detected.", call. = FALSE)
  L <- length(m)
  n <- substitute(...())  
  
  fe <- function(fit){
    
if(inherits(fit, "robu")){ 
  
    X <- cbind(1, fit$data.full$sd.eff.size)
    
    tmp <- lm((fit$data.full$effect.size) ~ X - 1)
    coef.na <- is.na(coef(tmp))
    if(any(coef.na)) stop("Model matrix not of full rank.", call. = FALSE) 
    
    f <- fit$ml[[2]]
    
    m <- robu(formula(bquote(.(f) ~ X - 1)), data = fit$data, var = fit$data.full$var.eff.size, study = fit$study_orig_id, rho = fit$mod_info$rho, small = TRUE, model = fit$modelweights)  
    
    h <- round(data.frame(b1 = m$reg_table$b.r[2], t.value = m$reg_table$t[2], p.value = m$reg_table$p[2], b1.lower = m$reg_table$CI.L[2], b1.upper = m$reg_table$CI.U[2]), 4)
    
  } else {
    
    sei <- fit$sigma
    X <- cbind(1, sei)
    y <- fit$y
    
    m <- rma.uni(y, sei = sei, mods = X, intercept = FALSE)
    
    h <- round(data.frame(b1 = m$b[[2]], z.value = m$zval[2], p.value = m$pval[2], b1.lower = m$ci.lb[2], b1.upper = m$ci.ub[2]), 4)
  }
  
    result <- symnum(h$p.value, cut = c(0, .001, .01, .05, .1, 1), na = FALSE, symbols = c("***", "**", "*", ":-)", ":-))"), corr = FALSE)
    h <- cbind(h, result)

    attr(h, "rclab") <- c("", "(H0: funnel is symmetric)\nEgger symmetry test:")
    class(h) <- c("labdf", class(h)) 
    return(h)
  }
  setNames(lapply(m, fe), n)
}                            

#================================================================================================================================================================
                                                      
egger.data <- function(data = NULL, by, r = .5, dependent = FALSE, dep.model = "CORR", ef.name = "dint", se.name = "SD"){  
  
  data <- trim(data)
  
  ef.name <- trimws(ef.name)
  se.name <- trimws(se.name)
  
  metain <- function(data = NULL, by, r = .5){ 
    
    m <- split(data, data$study.name)
    m <- Filter(NROW, rm.allrowNA2(m)) 
    
    L <- if(missing(by)) m else { s <- substitute(by) ; h <- lapply(m, function(x) do.call("subset", list(x, s))) ;
    res <- Filter(NROW, h) ; if(length(res) == 0) NULL else res}
    
    ds <- Filter(Negate(is.null), lapply(seq_along(L), function(i) L[[i]][[ef.name]]))
    sds <- Filter(Negate(is.null), lapply(seq_along(L), function(i) L[[i]][[se.name]]))
    
    setNames(mapply(option1, ds = ds, sds = sds, r = r, SIMPLIFY = FALSE), names(L))
  }
  
  if(!dependent){  
    
    j <- eval(substitute(metain(data = data, by = by, r = r)))
    
    ds <-  sapply(seq_along(j), function(i) j[[i]][1])
    sds <- sapply(seq_along(j), function(i) j[[i]][2])
    X <- cbind(1, sds)
    
    m <- rma.uni(ds, sei = sds, mods = X, intercept = FALSE)
    
    h <- round(data.frame(b1 = m$b[[2]], z.value = m$zval[2], p.value = m$pval[2], b1.lower = m$ci.lb[2], b1.upper = m$ci.ub[2]), 4)
  }
  else {
    
    
    if(!missing(by)) { 
      
      s <- substitute(by) 
      data <- subset(data, eval(s))
    }
    
    sds <- data[[se.name]]
    X <- cbind(1, sds)
    
    f <- formula(paste(ef.name, "~ X - 1"))
    
    m <- robu(f, data = data, var = sds^2, study = data$study.name, rho = r, model = trimws(dep.model))  
    
    h <- round(data.frame(b1 = m$reg_table$b.r[2], t.value = m$reg_table$t[2], p.value = m$reg_table$p[2], b1.lower = m$reg_table$CI.L[2], b1.upper = m$reg_table$CI.U[2]), 4)
  }
  
  result <- symnum(h$p.value, cut = c(0, .001, .01, .05, .1, 1), na = FALSE, symbols = c("***", "**", "*", ":-)", ":-))"), corr = FALSE)
  h <- cbind(h, result)
  
  attr(h, "rclab") <- c("", "(H0: funnel is symmetric)\nEgger symmetry test:")
  class(h) <- c("labdf", class(h)) 
  return(h)
}
                
#================================================================================================================================================================
                           
rma.robu <- function(f, var, id, data, w.model = "CORR", rho = .8, small = TRUE, group, ...){
  
  f <- formula(f)  
  data <- roundi(data)
  
  if(!missing(group)) { 
    
    s <- substitute(group) 
  data <- subset(data, eval(s)) 
  f <- formula(bquote(.(f[[2]]) ~ 1))
  }
  
  m <- eval(substitute(robu(f, data = data, var = var, study = id, model = w.model, rho = rho, small = small, ...)))
  
  w <- m$data.full$r.weights
  
  res <- eval(substitute(rma.uni(f, vi = var, data = data, slab = id, weights = w, subset = NULL, ...)))
  
  res$se <- m$reg_table$SE
  res$zval <- m$reg_table$t
  res$pval <- m$reg_table$prob
  res$ci.lb <- m$reg_table$CI.L
  res$ci.ub <- m$reg_table$CI.U
  res$tau2 <- m$mod_info$tau.sq
  res$I2 <- m$mod_info$I.2
  
  return(res)
}                           
     
#================================================================================================================================================================                           

infl <- function(x, add = "id"){

s <- as.character(getCall(x)$slab)  
  
cols <- c(as.character(getCall(x)$yi[[2]]), if(length(s) != 0) s else NULL, add)
  
out <- influence(x, progbar = TRUE)

plot.infl.rma.uni(out, plotinf = 8, bg.infl = "cyan")

 data <- eval(getCall(x)$data)
 
 izeh <- data[cols]
 
 h <- izeh[out$is.infl,]
 
 message("\nNote: Don't remove before careful examination.")
 attr(h, "rclab") <- c("", "\nInfluential Effect Size(s):")
 class(h) <- c("labdf", class(h)) 
 return(h)
}                           
       
#================================================================================================================================================================
          
rrbinom <- function(n, size, p1, p2, rho = 0){
  
  UseMethod("rrbinom")
  
}
               
rrbinom.default <- function(n, size, p1, p2, rho = 0){

p <- p1
q <- p2

a <- function(rho, p, q) {
  rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
}

a.0 <- a(rho, p, q)

prob <- c(`(0,0)`= a.0, `(1,0)`= 1-q-a.0, `(0,1)`= 1-p-a.0, `(1,1)`= a.0+p+q-1)
if(min(prob) < 0) {
  print(prob)
  stop("A probability is negative.", call. = FALSE)
}

u <- sample.int(4, n * size, replace = TRUE, prob = prob)
y <- floor((u-1)/2)
x <- 1 - u %% 2
x <- colSums(matrix(x, nrow = size)) 
y <- colSums(matrix(y, nrow = size)) 

list(x = x, y = y)
}                        

#================================================================================================================================================================
          
ddint.plot <- function(dppc, dppt, nc, nt, rev.sign = FALSE, ...){
  
  a <- dppc
  b <- dppt
  
  din <- b - a  
  
  test <- if(!rev.sign || rev.sign & b < 0 & a < 0 & abs(b) < abs(a)) FALSE else TRUE
  
  like1 <- function(x) dt(dppc*sqrt(nc), df = nc - 1, ncp = x*sqrt(nc))
  like2 <- function(x) dt(dppt*sqrt(nt), df = nt - 1, ncp = x*sqrt(nt))
  
  d1 <- AbscontDistribution(d = like1, low1 = -5e1, up1 = 5e1, withStand = TRUE)
  d2 <- AbscontDistribution(d = like2, low1 = -5e1, up1 = 5e1, withStand = TRUE)
  
  like.dif <- function(x) distr::d(if(test) -(d2 - d1) else d2 - d1)(x)
  
  Mean <- integrate(function(x) x*like.dif(x), -Inf, Inf)[[1]]
  SD <- sqrt(integrate(function(x) x^2*like.dif(x), -Inf, Inf)[[1]] - Mean^2)
  
  din <- if(test) -din else din
  
  curve(like.dif, din-(5*SD), din+(5*SD), n = 1e3, panel.f = abline(v = din, lty = 3, col = 2), lwd = 2, xlab = "dint (dpost - dpre)",
        panel.l = text(din, .6, round(din, 4), pos = 3, font = 2, col = 2, srt = 90), ylab = "Density", ...)
  
  return(c(MEAN = din, SD = SD))
}
 
#================================================================================================================================================================                       
                       
rob.fig <- function(n = 2e3, cluster = 5, adj = 1, col = 4, cex = .7, pch = 19, ann = FALSE, alpha = .4, seed = NULL, ...)
{
  
  size <- cluster - 1  
  
  set.seed(seed)
  a <- rrbinom(n, size, .5, .5)
  
  plot(jitter(a$x, adj), jitter(a$y, adj), xaxt = "n", yaxt = "n",
       pch = pch, cex = cex, col = adjustcolor(col, alpha), ann = ann, ...)
}
 
#================================================================================================================================================================
                       
rdint <- function(n, dppc, dppt, nc, nt, rev.sign = FALSE){
  
  a <- dppc
  b <- dppt
  
  test <- if(!rev.sign || rev.sign & b < 0 & a < 0 & abs(b) < abs(a)) FALSE else TRUE
  
  like1 <- function(x) dt(dppc*sqrt(nc), df = nc - 1, ncp = x*sqrt(nc))
  like2 <- function(x) dt(dppt*sqrt(nt), df = nt - 1, ncp = x*sqrt(nt))
  
  d1 <- AbscontDistribution(d = like1, low1 = -15, up1 = 15, withStand = TRUE)
  d2 <- AbscontDistribution(d = like2, low1 = -15, up1 = 15, withStand = TRUE)
  
  distr::r(if(test) -(d2 - d1) else d2 - d1)(n)
}

#================================================================================================================================================================

qdint <- function(p, dppc, dppt, nc, nt, rev.sign = FALSE, lower.tail = TRUE){
  
  a <- dppc
  b <- dppt 
  
  test <- if(!rev.sign || rev.sign & b < 0 & a < 0 & abs(b) < abs(a)) FALSE else TRUE
  
  like1 <- function(x) dt(dppc*sqrt(nc), df = nc - 1, ncp = x*sqrt(nc))
  like2 <- function(x) dt(dppt*sqrt(nt), df = nt - 1, ncp = x*sqrt(nt))
  
  d1 <- AbscontDistribution(d = like1, low1 = -15, up1 = 15, withStand = TRUE)
  d2 <- AbscontDistribution(d = like2, low1 = -15, up1 = 15, withStand = TRUE)
  
  res <- distr::q.l(if(test) -(d2 - d1) else d2 - d1)(p)
  if(lower.tail) res else 1 - res
}

#================================================================================================================================================================

pdint <- function(q, dppc, dppt, nc, nt, rev.sign = FALSE, lower.tail = TRUE){
  
  a <- dppc
  b <- dppt
  
  test <- if(!rev.sign || rev.sign & b < 0 & a < 0 & abs(b) < abs(a)) FALSE else TRUE
  
  like1 <- function(x) dt(dppc*sqrt(nc), df = nc - 1, ncp = x*sqrt(nc))
  like2 <- function(x) dt(dppt*sqrt(nt), df = nt - 1, ncp = x*sqrt(nt))
  
  d1 <- AbscontDistribution(d = like1, low1 = -15, up1 = 15, withStand = TRUE)
  d2 <- AbscontDistribution(d = like2, low1 = -15, up1 = 15, withStand = TRUE)
  
  res <- distr::p.r(if(test) -(d2 - d1) else d2 - d1)(q)
  if(lower.tail) res else 1 - res
}

#================================================================================================================================================================


ddint <- function(x, dppc, dppt, nc, nt, rev.sign = FALSE){
  
  a <- dppc
  b <- dppt
  
  test <- if(!rev.sign || rev.sign & b < 0 & a < 0 & abs(b) < abs(a)) FALSE else TRUE
  
  like1 <- function(y) dt(dppc*sqrt(nc), df = nc - 1, ncp = y*sqrt(nc))
  like2 <- function(y) dt(dppt*sqrt(nt), df = nt - 1, ncp = y*sqrt(nt))
  
  d1 <- AbscontDistribution(d = like1, low1 = -5e1, up1 = 5e1, withStand = TRUE)
  d2 <- AbscontDistribution(d = like2, low1 = -5e1, up1 = 5e1, withStand = TRUE)
  
  like.dif <- function(y) distr::d(if(test) -(d2 - d1) else d2 - d1)(y)
  
  like.dif(x)
}                       

#================================================================================================================================================================                                   
                                   
rev.dint.norm <- function(percent){ 
  
a <- as.numeric(substr(percent, 1, nchar(percent)-1)) / 1e2
round(qnorm(a + pnorm(0)), 4)

}
                                   
#================================================================================================================================================================
                                                                     
impute <- function(D, FUN = median){
  y <- sapply(D, is.numeric)
  D[y] <- lapply(D[y], function(x){x[is.na(x)] <- FUN(x, na.rm = TRUE); x})
  return(D)
}                                   

                                   
#================================================================================================================================================================
                                   
hist.pt <- function(x, breaks = "Sturges", xlab = "x", ylab = if(freq) "Frequency" else "Density", freq = TRUE, ...){
  
  h <- hist(x, breaks = breaks, plot = FALSE)

  x <- rep(h$mids, h$counts)
  y <- if(freq) sequence(h$counts) else unlist(lapply(seq_len(length(h$mids)), function(i) seq(0, h$density[i], length.out = h$counts[i])))

plot(x, y, xlab = xlab, ylab = ylab, ...)

invisible(h)
}                                  
 
#================================================================================================================================================================                                                      
                                                                                                            
pt.curve <- function(X, adjust = 1, compact = NULL, pch = 16, col = 2, cex = .7, seed = 0, reset = TRUE, add = FALSE, na.rm = TRUE, ...) {
  
  if(na.rm) X <- na.omit(X)  
  n.target <- length(X)
  
  d <- density(X, adjust = adjust, n = n.target)
  
  n <- if(!is.null(compact)) { 
    
    auc <- sum(d$y*median(diff(d$x)))/(diff(range(d$x))*max(d$y))
    
    compact*ceiling(n.target/auc)
    
  } else { n.target }
  
  set.seed(seed)
  pts <- data.frame(x = runif(n, min(d$x), max(d$x)), y = runif(n, 0, max(d$y)))
  
  pts <- pts[pts$y < approx(d$x, d$y, xout = pts$x)$y, ]
  
  if(nrow(pts) == 0) stop("Increase the size of sample 'X' OR use 'compact = NULL'.", call. = FALSE)
  
  pts <- pts[sample(seq_len(nrow(pts)), n, replace = TRUE), ]
  
  if(!add){
  
  if(reset) graphics.off()    
  plot(pts, pch = pch, col = col, cex = cex, ...)
    
  } else {
    
  points(pts, pch = pch, col = col, cex = cex, ...)
    
  }
}             

#================================================================================================================================================================                                                      

qplot <- function(..., col.pt = 1, pch = 21, cex = 1, col.line = 2, main = NULL){

M <- list(...)
n <- substitute(...())

foo <- function(x, col.pt, pch, cex, col.line, main){ 
  
qqnorm(x, col = col.pt, cex = cex, pch = pch, main = main)
qqline(x, col = col.line)
}

graphics.off()             
org.par <- par(no.readonly = TRUE)
on.exit(par(org.par))

N <- length(M)

if(N > 1L) { par(mfrow = n2mfrow(N)) ; set.margin() }

invisible(lapply(1:N, function(i) foo(M[[i]], col.pt = col.pt, pch = pch, cex = cex, col.line = col.line, main = if(is.null(main)) n[[i]] else if(is.na(main)) NA else main[i])))
}                                                      
       
#================================================================================================================================================================
                                      
plan.item <- function(margin = .2, S.index = NA){
  
  margin[margin > .99] <- .99  
  margin[margin < .01] <- .01
  
  data.frame(n.item = ceiling(1/margin^2), margin = paste0("+-", margin), S.index = S.index, lower = if(!is.na(S.index)) S.index - margin else NA, upper = if(!is.na(S.index)) S.index + margin else NA)
}

#================================================================================================================================================================
                            
plan.coder <- function(se2irr = .2, S.index = NA){
  
  se2irr[se2irr > .99] <- .99  
  se2irr[se2irr < .01] <- .01 
  
  n <- ceiling(2/se2irr)
  se <- if(!is.na(S.index)) se2irr*S.index else NA
  
  data.frame(n.coder = n, se2irr = se2irr, S.index = S.index, lower = if(!is.na(S.index)) S.index - 2*se else NA, upper = if(!is.na(S.index)) S.index + 2*se else NA, conf.lev = .95)
}
      
#================================================================================================================================================================
             
irr.diag <- function(X, useNA = "ifany"){
  
  a <- detail2(X, useNA = useNA)
  b <- detail(X, useNA = useNA)
  
  round(data.frame(KAPPA = a, SA = b), 3)
}

#==========================================================================================================================================

find.irr <- function(X, what){
  
  s <- substitute(what)
  X <- trim(X)
  res <- try(rm.colrowNA(do.call("subset", list(X, s))[c("study.name", "group.name")]), silent = TRUE)
  if(inherits(res, "try-error")) NULL else res
}
             
#============================================================================================================================================
                         
exam.efa <- function(x, factors, data = NULL, covmat = NULL, n.obs = NA,
                     subset = NULL, na.action = "na.omit", start = NULL,
                     scores = c("none", "regression", "Bartlett"),
                     rotation = "varimax", control = NULL, cutoff = .5, digits = 6, plot = TRUE, file.name = NULL, ...){
  
  
  cc <- match.call(expand.dots = FALSE)
  cc[[1]] <- quote(factanal)
  fit <- eval.parent(cc)
  fit$call <- match.call(expand.dots = FALSE)
  
  
  res <- round(transform(subset(as.data.frame.table(fit[[2]]), Freq >= cutoff),
                         Var2 = match(Var2, unique(Var2)), Var1 = as.numeric(Var1)), digits = digits)
  
  names(res) <- c("Item", "Factor", "Loading")
  
  rownames(res) <- NULL
  
  file.name <- trimws(file.name)  
  
  if(length(file.name) != 0){
    
    nm <- paste0(file.name, ".csv")
    ur <- try(write.csv(res, nm, row.names = FALSE), silent = TRUE)
    if(inherits(ur, "try-error")) stop(paste0("\nClose the Excel file '", nm, "' and try again OR pick another file name."), call. = FALSE)
    message(paste0("\nNote: Check folder '", basename(getwd()),"' for the Excel file '", nm, "'.\n"))
  } 
  
  f <- res$Factor
  i <- res$Item
  u <- unique(f)
    
  if(plot){
    
    graphics.off()
    org.par <- par(no.readonly = TRUE)
    on.exit(par(org.par))
    par(mar = c(3.8, 1.1, 1.5, 4.1), mgp = c(1.7, .5, 0))
    
    y <- unlist(lapply(u, function(j) seq_len(sum(f == j))))
    
    plot(f, y, las = 1, pch = 22, cex = 1.2, xlim = c(-.1, max(f)+.1), axes = FALSE, xlab = NA, main = NA, font.lab = 2, ylab = NA)
    
    at <- mean(u)
    
    mtext("FACTORS", 1, line = 2.5, font = 2, at = at)
    
    text(f, 0, f, pos = 1, xpd = NA, font = 2, cex = 1.75)
    
    rect(u-.5, 0, u+.5, tapply(y, f, FUN = max)+.5, col = adjustcolor(1:8, .25), xpd = NA, border = NA)
    
    dup <- duplicated(i) | duplicated(i, fromLast = TRUE)
    
    text(f, y, i, pos = 4, cex = .7, xpd = NA, font = 2, col = ifelse(dup, 2, 1))
    
    legend(at, par('usr')[4], legend = "ITEMS", pch = 22, horiz = TRUE, bty = "n", text.font = 2, xpd = NA, pt.cex = 1.4, yjust = .5)
  }
  return(res)
}                  

#============================================================================================================================================                      

setRowNames <- function (object = nm, nm) 
{
  rownames(object) <- nm
  object
}

#============================================================================================================================================
                       
setColNames <- function (object = nm, nm) 
{
  colnames(object) <- nm
  object
}

#============================================================================================================================================                       
                       
ave.effect <- function(data = NULL, by, r = .5, ef.name = "dint", se.name = "SD"){  
  
  data <- rm.colrowNA(trim(data))
  check <- "study.name" %in% names(data)
  if(!check) stop("Add a new column called 'study.name'.", call. = FALSE)
  
    m <- split(data, data$study.name)
    
    L <- if(missing(by)) { 
      
      m  
    
    } else { 
    
      s <- substitute(by) 
      res <- Filter(NROW, lapply(m, function(x) do.call("subset", list(x, s)))) 
      if(length(res) == 0) {
        
        message("Note: Moderator(s) not found; returning unmoderated averages.")
        m 
        
        } else { res }
    }
    
    ds <- lapply(seq_along(L), function(i) L[[i]][[trimws(ef.name)]])
    sds <-lapply(seq_along(L), function(i) L[[i]][[trimws(se.name)]])
    
    j <- setNames(mapply(option1, ds = ds, sds = sds, r = r, SIMPLIFY = FALSE), names(L))
  
  setRowNames(setNames(do.call(rbind.data.frame, j), c("dint", "SD")), names(j))
}

#=============================================================================================================================
                                          
data.str <- function(dat, drop = NULL){
  
  dat <- if(is.null(drop)) dat else drop.col(dat, vec = drop)
  
  setNames(lapply(names(dat), function(i) sort(unique(dat[[i]]))), names(dat))
}                                          

#=============================================================================================================================
                  
long <- function(data, one.cols, multi.cols = NULL, multi.colnames = NULL, time.name = "time", 
                 time.val = NULL, replace = "#N/A", with = NA, drop = NULL, order.by = one.cols[1]){
  
  data <- rm.colrowNA(trim(data)) 
  
  data[sapply(data, `%in%`, replace)] <- with
  
  varying <- if(is.null(multi.cols)) setdiff(names(data), one.cols) else multi.cols
  
  time.val <- if(is.null(time.val) & !is.null(multi.cols)) seq_along(varying[[1L]]) else time.val
  
  if(is.null(time.val) || length(time.val) == 1) stop("Please provide unique values for 'time.val'.", call. = FALSE)
  
  res <- tryCatch({na.omit(reshape(data, direction = 'long', timevar = time.name, times = time.val,
                                   idvar = one.cols, v.names = multi.colnames,
                                   varying = varying, drop = drop))}, 
                  error = function(e) {
                    if(as.character(e) == "Error in guess(varying): failed to guess time-varying variables from their names\n") stop("Please provide the 'multi.cols'.", call. = FALSE) 
                    else e
                  })
  
  if(inherits(res, "simpleError")) stop(paste0(res[[2]]), call. = FALSE)
  res <- res[order(res[[order.by]]), ]
  
  rownames(res) <- NULL
  
  res[] <- lapply(res, function(x) type.convert(as.character(x), as.is = TRUE))
  return(res)
}                

#=============================================================================================================================
                  
mask <- function(data, what, full = FALSE){
  
  data[] <- lapply(data, function(x) type.convert(as.character(x), as.is = TRUE))

  f1 <- function(x) as.numeric(factor(x, levels = unique(x)))
  f2 <- function(x) {
    temp <- substr(x, 1, 1)
    paste0(temp, ave(x, temp, FUN = function(y) match(y, unique(y))))
  }
  cols <- names(data)[sapply(data, is.numeric)]
  num.cols <- cols[cols %in% what]
  cols <- names(data)[sapply(data, is.character)]
  char.cols <- cols[cols %in% what]

  if(!full){
    
    if(length(num.cols))  data[num.cols] <- lapply(data[num.cols], f1)
    if(length(char.cols)) data[char.cols] <- lapply(data[char.cols], f2)
    
  }else{
    
    data[what] <- lapply(data[what], f1)
  }
  return(data)
}

#=============================================================================================================================
                     
make.dummy <- function(data, what){

ff <- function(data, what){
  
  if(!inherits(data, 'data.frame')) stop("'data' must be a data.frame.", call. = FALSE)
  if(!is.character(what)) stop("'what' must be a character.", call. = FALSE)
  what <- trimws(what)
  data <- trim(data)
  if(!is.character(as.vector(data[[what]]))) stop('Not a character variable.', call. = FALSE)
  formula <- as.formula(paste('~', what))
  data.frame(model.matrix(formula, data = data))[-1L]
}

cbind(data, do.call(cbind, lapply(what, ff, data = data)))
}                     

#=============================================================================================================================
                     
plot.cor <- function (corr, outline = FALSE, col = colorRampPalette(c(4, 2))(choose(ncol(corr), 2)), 
                      upper.panel = c("ellipse", "number", "none"), lower.panel = c("ellipse", "number", "none"), diag = c("none", "ellipse", "number"), digits = 2, bty = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, cex.lab = par("cex.lab"), cex = 0.75 * par("cex"), mar = c(2, 2, 4, 2)-.2, ...)
{
  savepar <- par(pty = "s", mar = mar)
  on.exit(par(savepar))
  if (is.null(corr))
    return(invisible())
  if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1))
    stop("Please input a correlation matrix")
  plot.new()
  par(new = TRUE)
  rowdim <- dim(corr)[1]
  coldim <- dim(corr)[2]
  rowlabs <- dimnames(corr)[[1]]
  collabs <- dimnames(corr)[[2]]
  if (is.null(rowlabs))
    rowlabs <- 1:rowdim
  if (is.null(collabs))
    collabs <- 1:coldim
  rowlabs <- as.character(rowlabs)
  collabs <- as.character(collabs)
  col <- rep(col, length = length(corr))
  dim(col) <- dim(corr)
  upper.panel <- match.arg(upper.panel)
  lower.panel <- match.arg(lower.panel)
  diag <- match.arg(diag)
  cols <- 1:coldim
  rows <- 1:rowdim
  maxdim <- max(length(rows), length(cols))
  plt <- par("plt")
  xlabwidth <- max(strwidth(rowlabs[rows], units = "figure", cex = cex.lab))/(plt[2] - plt[1])
  xlabwidth <- xlabwidth * maxdim/(1 - xlabwidth)
  ylabwidth <- max(strwidth(collabs[cols], units = "figure", cex = cex.lab))/(plt[4] - plt[3])
  ylabwidth <- ylabwidth * maxdim/(1 - ylabwidth)
  plot(c(-xlabwidth - 0.5, maxdim + 0.5), c(0.5, maxdim + 1 + ylabwidth), type = "n", bty = bty, axes = axes, xlab = "", ylab = "", asp = asp, cex.lab = cex.lab, ...)
  text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], adj = 1, cex = cex.lab)
  text(cols, rep(length(rows) + 1, length(cols)), labels = collabs[cols], srt = 90, adj = 0, cex = cex.lab)
  mtext(xlab, 1, 0)
  mtext(ylab, 2, 0)
  mat <- diag(c(1, 1))
  plotcorrInternal <- function() {
    if (i == j){ 
      if (diag == 'none'){
        return()
      } else if (diag == 'number'){
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else if (diag == 'ellipse') {
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      }
    } else if (i >= j){ 
      if (lower.panel == 'ellipse') { 
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (lower.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    } else { 
      if (upper.panel == 'ellipse') { 
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (upper.panel == 'number') { 
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    }
  }
  for (i in 1:dim(corr)[1]) {
    for (j in 1:dim(corr)[2]) {
      plotcorrInternal()
    }
  }
  invisible()
}                     

#=============================================================================================================================
                     
rbetab <- function(n, mu.p, disp){
  
if(mu.p < 0 || mu.p > 1) { message("Warning: 'mu.p' is 'average probability' of a 'beta dist.' bound between '0' & '1'.") ;
mu.p[mu.p < 0] <- 0 ;
mu.p[mu.p > 1] <- 1 }
    
  shape1 <- mu.p * disp
  shape2 <- (1 - mu.p) * disp
  rbeta(n, shape1 = shape1, shape2 = shape2)
}                     
                     
                     
#=============================================================================================================================
           
dint.sem <- function(n.c = 45, mpre.c = 10, sdpre.c = 1.5, mpos.c = 13, 
                   sdpos.c = 1.3, n.t = 39, mpre.t = 11, sdpre.t = 1.3, 
                   mpos.t = 18, sdpos.t = 1.2, correct = FALSE,
                   r.c = .6, r.t = .6, plot = FALSE, ...){

cor.C <- cor.mat(r.c, 2)

Cov.C <- lavaan::cor2cov(cor.C, sds = c(sdpre.c, sdpos.c), names = c("x_pre.C","x_post.C"))

Mean.C <- c(mpre.c, mpos.c)

ct.C <- if(correct) cfactor(n.c-1) else 1

#--------- Effect size for Control group (everything ends in ".C"):

model.C <- paste('

eta_pre.C = ~ sd_pre.C*x_pre.C
eta_post.C = ~ sd_post.C*x_post.C

eta_pre.C ~~ r*eta_post.C

x_pre.C ~~ 0*x_pre.C
x_post.C ~~ 0*x_post.C

x_pre.C ~ m_pre.C*1
x_post.C ~ m_post.C*1

SMD.dif.C := (m_post.C - m_pre.C) / sqrt(sd_pre.C^2+sd_post.C^2-2*sd_pre.C*sd_post.C*r) *', ct.C)


#--------- Effect size for Treatment group (everything ends in ".T"):


cor.T <- cor.mat(r.t, 2)

Cov.T <- lavaan::cor2cov(cor.T, sds = c(sdpre.t, sdpos.t), names = c("x_pre.T","x_post.T"))

Mean.T <- c(mpre.t, mpos.t)

ct.T <- if(correct) cfactor(n.t-1) else 1

model.T <- paste('
eta_pre.T = ~ sd_pre.T*x_pre.T
eta_post.T = ~ sd_post.T*x_post.T

eta_pre.T ~~ r*eta_post.T

x_pre.T ~~ 0*x_pre.T
x_post.T ~~ 0*x_post.T

x_pre.T ~ m_pre.T*1
x_post.T ~ m_post.T*1

SMD.dif.T := (m_post.T - m_pre.T) / sqrt(sd_pre.T^2+sd_post.T^2 -2*sd_pre.T*sd_post.T*r) *', ct.T)

# / ( (sqrt(sd_pre.C^2+sd_post.C^2-2*sd_pre.C*sd_post.C*r) + sqrt(sd_pre.T^2+sd_post.T^2-2*sd_pre.T*sd_post.T*r)) / 2 )

mod <- c(' group: Control ', model.C, ' group: Treatment ', model.T,
         ' dint.sem := SMD.dif.T - SMD.dif.C ' )

fit <- lavaan::sem(mod, std.lv = TRUE,
           sample.cov = list(Control = Cov.C, Treatment = Cov.T),
           sample.mean = list(Mean.C, Mean.T),
           sample.nobs = list(n.c, n.t),
           sample.cov.rescale = FALSE)

est.er <- parameterEstimates(fit)[23:25, 6:8]   #[25, 7:8]

#semPlot::semPaths(fit, edge.label.cex = 1.5, sizeLat = 9, sizeInt = 3, combineGroups = TRUE, ...)

if(plot){
  
library(semPlot)

fit.C <- lavaan::sem(model.C, sample.cov = Cov.C, sample.mean = Mean.C,
             sample.nobs= n.c, std.lv = TRUE,
             sample.cov.rescale = FALSE)

fit.T <- lavaan::sem(model.T, sample.cov = Cov.T, sample.mean = Mean.T,
             sample.nobs= n.c, std.lv = TRUE,
             sample.cov.rescale = FALSE)

graphics.off()
org.par <- par(no.readonly = TRUE)
on.exit(par(org.par))
par(mfrow = 2:1)

semPlot::semPaths(fit.C, label.prop=.9, sizeInt = 4, sizeLat = 9, curvePivot = TRUE,
                  equalizeManifests = FALSE, optimizeLatRes = TRUE, node.width = 1.5, 
                  edge.width = 0.5, shapeMan = "rectangle", shapeLat = "ellipse", edge.label.cex = 1.8, 
                  shapeInt = "triangle", sizeMan = 9, asize = 6, unCol = "red4", label.cex = 1.4, ...)

semPlot::semPaths(fit.T, label.prop=.9, sizeInt = 4, sizeLat = 9, curvePivot = TRUE, 
                  equalizeManifests = FALSE, optimizeLatRes = TRUE, node.width = 1.5, 
                  edge.width = 0.5, shapeMan = "rectangle", shapeLat = "ellipse", edge.label.cex = 1.8,
                  shapeInt = "triangle", sizeMan = 9, asize = 6, unCol = "blue4", label.cex = 1.4, ...)
}

return(est.er)

}                                  
   
#============================================================================================================================                     
             
                     
samp.dist <- function(n, pop.dist = c('nor','exp','uni','poi','bin','gam','chi','tds', 'bet'), param1 = NULL, param2 = NULL, times = 1e4, others = FALSE, xlab = NA, ylab = "Density", obs.mean = NULL, rev.page = FALSE, seed.pop = 0, digits = 3, xaxt = "s", reset = FALSE, unknown = FALSE, obs.sd = NULL, ...){
  
  pop.dist <- match.arg(pop.dist)
  
  if(unknown & is.null(obs.sd)) stop("Provide 'obs.sd' to see the place of your observed test statistic.", call. = FALSE)
  
  samples <- switch(pop.dist,
                    "exp" = replicate(times, rexp(n, param1)),
                    "nor" = replicate(times, rnorm(n, param1, param2)),
                    "uni" = replicate(times, runif(n, param1, param2)),
                    "poi" = replicate(times, rpois(n, param1)),
                    "bin" = replicate(times, rbinom(n, param1, param2)),
                    "gam" = replicate(times, rgamma(n, param1, param2)),
                    "chi" = replicate(times, rchisq(n, param1)),
                    "tds" = replicate(times, rt(n, param1)),
                    "bet" = replicate(times, rbetab(n, param1, param2)))
  
  set.seed(seed.pop)
  pop <- switch(pop.dist,
                "exp" = rexp(1e4, param1),
                "nor" = rnorm(1e4, param1, param2),
                "uni" = runif(1e4, param1, param2),
                "poi" = rpois(1e4, param1),
                "bin" = rbinom(1e4, param1, param2),
                "gam" = rgamma(1e4, param1, param2),
                "chi" = rchisq(1e4, param1),
                "tds" = rt(1e4, param1),
                "bet" = rbetab(1e4, param1, param2))
  
  
  all.sample.means <- colMeans(samples, na.rm = TRUE)   
  all.sample.sd <- apply(samples, 2, sd, na.rm = TRUE)
  
  if(reset){
  graphics.off()
  org.par <- par(no.readonly = TRUE)
  on.exit(par(org.par))
  }
  
  h <- if(others) 5 else 3
  dev <- if(!rev.page) n2mfrow(h) else rev(n2mfrow(h))
  par(mfcol = dev, mar = c(2.5, 2.6, 1.8, .5), mgp = c(1.65, .4, 0))
  
  
  pt.curve(pop, main = "Population", cex.main = 1, col = 4, ylab = ylab, xlab = xlab, pch = '.', xaxt = xaxt, ...)
  sd.pop <- sd(pop, na.rm = TRUE)
  m <- mean(pop, na.rm = TRUE)
  abline(v = c(m, m-sd.pop, m+sd.pop), col = 3, lty = c(2, 1,1))
  at1 <- axTicks(1)
  
  pt.curve(all.sample.means,main = "Sampling Distribution of Means (Normal)", cex.main = 1, xlim = range(pop, na.rm = TRUE), xlab = xlab, ylab = ylab, pch = '.', xaxt = "n", ...)
  at2 <- axTicks(1)
  if(xaxt == "s") axis(1, at = union(at1, at2), ...)
  
  se <- sd(all.sample.means, na.rm = TRUE)
  m <- mean(all.sample.means, na.rm = TRUE)
  if(!is.null(obs.mean)) {  
    points(obs.mean, 0, pch = 23, bg = "cyan", col = 'magenta', cex = 2, xpd = NA)
    obs.z <-  if(!unknown) (obs.mean-m) / (sd.pop/sqrt(n)) else (obs.mean-m) / (obs.sd/sqrt(n))
  }
  
  abline(v = c(m, m-se, m+se), col = 3, lty = c(2, 1,1))
  
  z <- if(!unknown)(all.sample.means - m) / (sd.pop/sqrt(n)) else (all.sample.means - m) / (all.sample.sd/sqrt(n))
  
  pt.curve(z, main = paste("Sampling Distribution of",if(unknown) "'t'" else "'z'","(Std. Normal)"), cex.main = 1, xlab = xlab, ylab = ylab, col = 'purple', pch = '.', xaxt = xaxt, xlim = if(!is.null(obs.mean)) range(z, obs.z, na.rm = TRUE) else range(z, na.rm = TRUE), ...)
  se.z <- sd(z, na.rm = TRUE)
  m.z <- mean(z, na.rm = TRUE)
  if(!is.null(obs.mean)) points(obs.z, 0, pch = 23, bg = "cyan", col = 'magenta', cex = 2, xpd = NA)
  qs <- quantile(z, names = F, probs = c(.025, .975), na.rm = TRUE)
  arrows(qs[1], .1, -3, .1, code = 2, length = .12, angle = 20)
  arrows(qs[2], .1, 3, .1, code = 2, length = .12, angle = 20)
  abline(v = c(m.z, qs), col = 3, lty = c(2, 1,1))
  
  if(others){
    
    all.sample.sums <- colSums(samples, na.rm = TRUE)
    #all.sample.vars <- apply(samples,2,var, na.rm = TRUE) 
    pt.curve(all.sample.sums, main = "Sampling Distribution of Sum", cex.main = 1, xlab = xlab, pch = ".", ylab = ylab, xaxt = xaxt, ...)
    pt.curve(all.sample.sd, main = "Sampling Distribution of Sd", cex.main = 1, xlab = xlab, pch = ".", ylab = ylab, xaxt = xaxt, ...) 
    #m.sd <- mean(all.sample.sd) ; med.sd <- median(all.sample.sd)
    abline(v = mean(all.sample.sd, na.rm = TRUE), col = 3)
  }
  
  obj <- round(c(sd.pop = sd.pop, se = se, clt.se = sd.pop / sqrt(n), obs.z = if(!is.null(obs.mean)) obs.z else NA), digits)
  names(obj)[4] <- if(unknown) "obs.t" else "obs.z"
  setNames(obj, names(obj))
}                     
      
#=============================================================================================================================
                 
meta_bayes <- function(data = NULL, by, tau.prior = function(x){dhalfnormal(x)}, mu.prior = mu.norm(-6, 6), r = .5, ef.name = "dint", se.name = "SD"){
  
  j <- eval(substitute(ave.effect(data = data, by = by, r = r, ef.name = ef.name, se.name = se.name)))
  
  res <- bayesmeta(                y = j[[1]],
                                   sigma = j[[2]],
                                   labels = rownames(j), 
                                   mu.prior = mu.prior,
                                   tau.prior = tau.prior)
  res$call <- match.call(expand.dots = FALSE)
  
  return(res)
}                 


#=============================================================================================================================                 
                 
metal.dint2 <- function(data = NULL, by, over = time, mu.prior = mu.norm(-6, 6), tau.prior = function(x){dhalfnormal(x)}, 
         r = .5, ef.name = "dint", se.name = "SD", method = c("robust", "bayes")){
  
  over <- deparse(substitute(over))
  
  method <- match.arg(method)
  
  chep <- sort(unique(na.omit(data[[over]])))
  
  
  f1 <- function(data, zy, r = .5){
    
    m <- split(data, data$study.name)
    
    L <- if(missing(zy)) { 
      
      m  
      
    } else { 
      
      s <- substitute(zy) 
      res <- Filter(NROW, lapply(m, function(x) do.call("subset", list(x, s)))) 
      if(length(res) == 0) {
        
        message("Note: Moderator(s) not found; returning unmoderated averages.")
        m 
        
      } else { res }
    }
    
    ds <- lapply(seq_along(L), function(i) L[[i]][[trimws(ef.name)]])
    sds <-lapply(seq_along(L), function(i) L[[i]][[trimws(se.name)]])
    
    setNames(mapply(option1, ds = ds, sds = sds, r = r, SIMPLIFY = FALSE), names(L))
    
  }
  
  
  f2 <- function(j, tau.prior, mu.prior){  
    
    ds <- sapply(seq_along(j), function(i) j[[i]][1])
    sds <- sapply(seq_along(j), function(i) j[[i]][2])
    
    test <- length(ds) >= 2
    
    if(!test) return(NA)
    
    res <- bayesmeta(        y = ds,
                             sigma = sds,
                             labels = names(j), 
                             tau.prior = tau.prior,
                             mu.prior = mu.prior)
    res$call <- match.call(expand.dots = FALSE)
    
    return(res)
  }  
  
  G <- if(missing(by)) { lapply(chep, function(y) bquote(.(as.name(noquote(over))) == .(y))) 
    
  } else {
    
    s <- substitute(by)
    lapply(chep, function(x) bquote(.(s) & .(as.name(noquote(over))) == .(x)))
  }
  
  if(method == "robust"){
    
    dd <- Filter(length, lapply(seq_along(G), function(j) subset(data, eval(G[[j]]))))
    
    if(length(dd) == 0) return(NA)
    
    f <- try(setNames(lapply(seq_along(dd), function(x) robu(dint~1, data = dd[[x]], studynum = as.vector(study.name), var = SD^2)), paste(over, chep)), silent = TRUE)
    if(inherits(f, "try-error")) NA else f[sapply(seq_along(f), function(i) is.finite(f[[i]]$reg_table$dfs[1]) & f[[i]]$reg_table$dfs[1] > 1)]
  } 
  
  else 
    
  {
    
    go <- length(G)
    
    k <- vector("list", go)
    
    for(w in seq_len(go)) k[[w]] <- f1(data = data, zy = eval(G[[w]]), r = r)
    
    so <- length(k)
    
    z <- vector("list", so)
    
    for(a in seq_len(so)) z[[a]] <- f2(j = k[[a]], tau.prior = tau.prior, mu.prior = mu.prior)
    
    setNames(z, paste(over, chep))
  }
}    
                    
#=============================================================================================================================
                    
dint.plot3 <- function(..., main = NA, ylab = "Effect Size (dint)", labels = NULL, file = NULL,
                      percent = FALSE, lwd = 1, reset = TRUE, cex.txt = .9, cex.pt = 6.3){
  
  
  m <- Filter(NROW, lapply(list(...), function(x) x[!is.na(x)]))
  L <- length(m)
  n <- substitute(...())
  
  if(reset){
    graphics.off()
    org.par <- par(no.readonly = TRUE)
    on.exit(par(org.par))
  }                          
  
  
  if(L > 1L) { par(mfrow = n2mfrow(L)) ; set.margin() ; if(percent) par(mar = c(1.5, 2.6, 1.8, 1.6)) }
  
  G <- function(fit, main, labels){  
    
    L <- length(fit)  
    
    bs <- all(sapply(fit, inherits, "bayesmeta"))
    
    if(bs){
      
      mu <- sapply(1:L, function(i) fit[[i]]$summary["mean","mu"])
      lo <- sapply(1:L, function(i) fit[[i]]$summary["95% lower","mu"])
      hi <- sapply(1:L, function(i) fit[[i]]$summary["95% upper","mu"])
      k <- sapply(1:L, function(i) fit[[i]]$k)
      
    } else {
      
      mu <- sapply(1:L, function(i) fit[[i]]$reg_table$b.r[[1]])
      lo <- sapply(1:L, function(i) fit[[i]]$reg_table$CI.L[[1]])            
      hi <- sapply(1:L, function(i) fit[[i]]$reg_table$CI.U[[1]])
      k <- sapply(1:L, function(i) fit[[i]]$N)
      dfs <- sapply(1:L, function(i) fit[[i]]$reg_table$dfs[[1]] < 4)
    }
    
    x <- 0:(L-1)
    
    plot(x, mu, type = "l", xlim = range(x)+c(-.05, .05), ylim = range(lo, hi), ylab = ylab, lwd = lwd, lty = 2, lend = 1, font.lab = 2,
         xaxt = "n", xlab = NA, panel.last = axis(1, at = x, labels = labels), main = main, las = 1, cex.axis = .9, padj = .3)
    
    invisible(lapply(seq_len(L), function(i) if(!is.na(mu[i])) lines(c(i-1, i-1), c(lo[i], hi[i]), lwd = 4, lend = 1, col = if(bs) 2 else if(!bs & dfs[i]) 8 else 4)))
    
    axis(3, at = x, labels = paste0("k=", k), cex.axis = .65, font = 2, xpd = NA, las = 3, mgp = c(1,.6,0))
    
    points(x, mu, pch = 22, cex = cex.pt, bg = "cyan", col = "magenta", xpd = NA)
    
    text(x, mu, round(mu, 3), cex = cex.txt, font = 2, xpd = NA)
    
    text(x, lo, round(lo, 3), cex = cex.txt, font = 2, xpd = NA, pos = 1)
    
    text(x, hi, round(hi, 3), cex = cex.txt, font = 2, xpd = NA, pos = 3)
    
    if(percent){
      
      text(x*1.04, c(lo, mu, hi),
           paste0("[", dint.norm(c(lo, mu, hi)),"]"), cex = .7, font = 2, xpd = NA, pos = 4, col = "magenta")
    }
    mu
  }
  
 invisible(lapply(seq_len(L), function(i) G(m[[i]], main = if(is.null(main)) n[[i]] else if(is.na(main)) NA else main[i], labels = if(is.null(labels)) names(m[[i]]) else labels[[i]])))
  

  out <- data.frame(lapply(meta.stats(...), unlist))
  ( out2 <- cbind(code = rownames(out), out)[-c(6,9)] )
  
  if(!is.null(file)){   
  file <- paste0(file, ".doc")
  tab_df(out2,
         file=file) }
}                     

#=======
     
dint.plot3 <- function(..., main = NA, ylab = "Effect Size (dint)", labels = NULL, file = NULL,
                       percent = FALSE, lwd = 1, reset = TRUE, cex.txt = .9, cex.pt = 6.3, digits = 2){
  
  
  m <- Filter(NROW, lapply(list(...), function(x) x[!is.na(x)]))
  L <- length(m)
  n <- substitute(...())
  
  if(reset){
    graphics.off()
    org.par <- par(no.readonly = TRUE)
    on.exit(par(org.par))
  }                          
  
  
  if(L > 1L) { par(mfrow = n2mfrow(L)) ; set.margin() ; if(percent) par(mar = c(1.5, 2.6, 1.8, 1.6)) }
  
  G <- function(fit, main, labels){  
    
    L <- length(fit)  
    
    bs <- all(sapply(fit, inherits, "bayesmeta"))
    
    if(bs){
      
      mu <- sapply(1:L, function(i) fit[[i]]$summary["mean","mu"])
      lo <- sapply(1:L, function(i) fit[[i]]$summary["95% lower","mu"])
      hi <- sapply(1:L, function(i) fit[[i]]$summary["95% upper","mu"])
      k <- sapply(1:L, function(i) fit[[i]]$k)
      
    } else {
      
      mu <- sapply(1:L, function(i) fit[[i]]$reg_table$b.r[[1]])
      lo <- sapply(1:L, function(i) fit[[i]]$reg_table$CI.L[[1]])            
      hi <- sapply(1:L, function(i) fit[[i]]$reg_table$CI.U[[1]])
      k <- sapply(1:L, function(i) fit[[i]]$N)
      dfs <- sapply(1:L, function(i) fit[[i]]$reg_table$dfs[[1]] < 4)
    }
    
    x <- 0:(L-1)
    
    plot(x, mu, type = "l", xlim = range(x)+c(-.05, .05), ylim = range(lo, hi), ylab = ylab, lwd = lwd, lty = 2, lend = 1, font.lab = 2,
         xaxt = "n", xlab = NA, panel.last = axis(1, at = x, labels = labels), main = main, las = 1, cex.axis = .9, padj = .3)
    
    invisible(lapply(seq_len(L), function(i) if(!is.na(mu[i])) lines(c(i-1, i-1), c(lo[i], hi[i]), lwd = 4, lend = 1, col = if(bs) 2 else if(!bs & dfs[i]) 8 else 4)))
    
    axis(3, at = x, labels = paste0("k=", k), cex.axis = .65, font = 2, xpd = NA, las = 3, mgp = c(1,.6,0))
    
    points(x, mu, pch = 22, cex = cex.pt, bg = "cyan", col = "magenta", xpd = NA)
    
    text(x, mu, round(mu, 3), cex = cex.txt, font = 2, xpd = NA)
    
    text(x, lo, round(lo, 3), cex = cex.txt, font = 2, xpd = NA, pos = 1)
    
    text(x, hi, round(hi, 3), cex = cex.txt, font = 2, xpd = NA, pos = 3)
    
    if(percent){
      
      text(x*1.04, c(lo, mu, hi),
           paste0("[", dint.norm(c(lo, mu, hi)),"]"), cex = .7, font = 2, xpd = NA, pos = 4, col = "magenta")
    }
    mu
  }
  
  lapply(seq_len(L), function(i) G(m[[i]], main = if(is.null(main)) n[[i]] else if(is.na(main)) NA else main[i], labels = if(is.null(labels)) names(m[[i]]) else labels[[i]]))
  
  
  out <- data.frame(lapply(meta.stats(..., digits = digits), unlist))
  out2 <- cbind(code = if(is.null(labels))rownames(out) else unlist(labels), out[-1]) 
  rownames(out2) <- NULL
  
   
    file <- paste0(file, ".doc")
    tab_df(out2,
           file=file)
  out2
}                                        
#=============================================================================================================================
                                   
                                   
meta_bayes_dif <- function(list_fit, labels = NULL){
  
  list_fit <- list_fit[!is.na(list_fit)]    
  
  if(!all(sapply(list_fit, inherits, "bayesmeta"))) stop("Non-bayesmeta model detected.", call. = FALSE)
  
  comboMatrix <- t(combn(length(list_fit),2))
  
  conList <- lapply(1:nrow(comboMatrix),function(a){
    con <-    convolve(dens1 = function(x){list_fit[[comboMatrix[a,1]]]$dposterior(mu = x)},
                       dens2 = function(x){list_fit[[comboMatrix[a,2]]]$dposterior(mu = -x)},
                       cdf1 = function(x){list_fit[[comboMatrix[a,1]]]$pposterior(mu = x)},
                       cdf2 = function(x){1 - list_fit[[comboMatrix[a,2]]]$pposterior(mu = -x)})
    con$quantile(c(0.025, 0.975))
  })
  
  bnms <- if(is.null(labels)) names(list_fit) else labels
  
  names(conList) <- unlist(lapply(1:nrow(comboMatrix),function(x){
    paste(bnms[comboMatrix[x,1]],"vs",bnms[comboMatrix[x,2]])
  }))
  
  res <- data.frame(t(sapply(conList, c)))
  names(res) <- c("lower", "upper")
  
  transform(res, diff. = lower > 0 | upper < 0)
}    
                                   
#=============================================================================================================================
       
                                   
trim_fill <- function(x, side = NULL, xlab = "Effect Size (dint)", mgp = c(2, .5, 0), legend = FALSE, main = NA, col.main = 1, ...){

if(!inherits(x, "bayesmeta")) stop("Non-bayesmeta model detected.", call. = FALSE)  
  
aa <- metafor::rma.uni(x$y, sei = x$sigma, slab = x$labels)
taf <- trimfill.rma.uni(aa, side = side)
ok <- isFALSE(any(taf$fill))

metafor::funnel.rma(taf, xlab = xlab, mgp = mgp, main = if(ok) "No Missing \nStudies Found" else main, col.main = if(ok) 2 else col.main, ...)

if(legend) graphics::legend("topright", c("Actual Studies", "Filled Studies"), pch=21, pt.bg = 1:0)
}                                   
                                   
#=============================================================================================================================
                    
veva_hedge <- function(x, steps = c(0.025, 1), mods = NULL, weights = NULL, 
                       fe = FALSE, table = FALSE, pval = NULL, digits = 4){
  
  if(!inherits(x, "bayesmeta")) stop("Non-bayesmeta model detected.", call. = FALSE) 
  
  x <- weightfunct(x$y, x$sigma^2, steps = steps, mods = if(!is.null(mods)) as.formula(mods) else mods,
                   weights = weights, fe = fe, table = table, pval = pval)
  
  df <- length(x[[2]]$par) - length(x[[1]]$par)
  lrchisq <- 2 * (abs(x[[1]]$value - x[[2]]$value))
  pvalue <- 1 - pchisq(lrchisq, df)
  
  .W <- round(data.frame(df, lrchisq, pvalue, row.names = "Likelihood Ratio Test:"), digits)
  names(.W) <- c("Df", "X^2", "p.value")
  return(.W)
}                    

#=============================================================================================================================
                   

z_outlier <- function(fit, by){
    
    if(!inherits(fit, "bayesmeta")) stop("Non-bayesmeta model detected.", call. = FALSE) 
    
    s <- substitute(by)
    
    ss <- lapply(meta.stats(fit)[-9], as.numeric)
    
    z <- (fit$y - ss$mu)/ss$tau
    
    d <- data.frame(study.name = fit$labels, z = z)
    
    if(!missing(by)) d <- subset(d, eval(s))
    
    if(nrow(d) == 0) return(NA) else d
  }
                   
#=============================================================================================================================
                   
                   
robu2 <- function(formula, data, studynum, var.eff.size, userweights, 
                   modelweights = c("CORR", "HIER"), rho = 0.8, small = TRUE, 
                   ...) 
{
  
  
  tmp <- lm(as.formula(formula), data = data)
  
  coef.na <- is.na(coef(tmp))
  
  dropped <- colnames(model.matrix(tmp)[,coef.na, drop = FALSE])
  
  if(any(coef.na)) message("\n", toString(dQuote(dropped), width = 50), " dropped due to lack of data.\n")
  
  modelweights <- match.arg(modelweights)
  if (modelweights == "CORR" && rho > 1 | rho < 0) 
    stop("Rho must be a value between 0 and 1.")
  if (missing(userweights)) {
    user_weighting = FALSE
  }
  else {
    user_weighting = TRUE
  }
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  ml <- mf[[2]]
  m <- match(c("formula", "data", "studynum", "var.eff.size", 
               "userweights"), names(mf))
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  if (!user_weighting) {
    dframe <- data.frame(effect.size = mf[, 1], stats::model.matrix(formula, 
                                                                    mf)[, !coef.na, drop = FALSE], studynum = mf[["(studynum)"]], var.eff.size = mf[["(var.eff.size)"]])
    X.full.names <- names(dframe)[-match(c("effect.size", 
                                           "studynum", "var.eff.size"), names(dframe))]
  }
  else {
    dframe <- data.frame(effect.size = mf[, 1], stats::model.matrix(formula, 
                                                                    mf)[, !coef.na, drop = FALSE], studynum = mf[["(studynum)"]], var.eff.size = mf[["(var.eff.size)"]], 
                         userweights = mf[["(userweights)"]])
    X.full.names <- names(dframe)[-match(c("effect.size", 
                                           "studynum", "userweights", "var.eff.size"), names(dframe))]
  }
  study_orig_id <- dframe$studynum
  dframe$study <- as.factor(dframe$studynum)
  dframe$study <- as.numeric(dframe$study)
  dframe <- dframe[order(dframe$study), ]
  k_temp <- as.data.frame(unclass(rle(sort(dframe$study))))
  dframe$k <- k_temp[[1]][match(dframe$study, k_temp[[2]])]
  dframe$avg.var.eff.size <- stats::ave(dframe$var.eff.size, 
                                        dframe$study)
  dframe$sd.eff.size <- sqrt(dframe$var.eff.size)
  switch(modelweights, HIER = {
    dframe$weights <- 1/dframe$var.eff.size
  }, CORR = {
    dframe$weights <- 1/(dframe$k * dframe$avg.var.eff.size)
  })
  X.full <- dframe[c("study", X.full.names)]
  data.full.names <- names(dframe)[-match(c("studynum", X.full.names), 
                                          names(dframe))]
  data.full <- dframe[c(data.full.names)]
  k <- data.full[!duplicated(data.full$study), ]$k
  k_list <- as.list(k)
  M <- nrow(data.full)
  p <- ncol(X.full) - 2
  N <- max(data.full$study)
  W <- as.matrix(by(data.full$weights, data.full$study, function(x) diag(x, 
                                                                         nrow = length(x)), simplify = FALSE))
  X <- data.matrix(X.full)
  X <- lapply(split(X[, 2:(p + 2)], X[, 1]), matrix, ncol = p + 
                1)
  y <- by(data.full$effect.size, data.full$study, function(x) matrix(x))
  J <- by(rep(1, nrow(X.full)), X.full$study, function(x) matrix(x, 
                                                                 nrow = length(x), ncol = length(x)))
  sigma <- by(data.full$sd.eff.size, data.full$study, function(x) tcrossprod(x))
  vee <- by(data.full$var.eff.size, data.full$study, function(x) diag(x, 
                                                                      nrow = length(x)))
  SigmV <- Map(function(sigma, V) sigma - V, sigma, vee)
  sumXWX <- Reduce("+", Map(function(X, W) t(X) %*% W %*% 
                              X, X, W))
  sumXWy <- Reduce("+", Map(function(X, W, y) t(X) %*% W %*% 
                              y, X, W, y))
  sumXWJWX <- Reduce("+", Map(function(X, W, J) t(X) %*% W %*% 
                                J %*% W %*% X, X, W, J))
  Matrx_WKXX <- Reduce("+", Map(function(X, W, k) {
    t(X) %*% (W/k) %*% X
  }, X, W, k_list))
  Matrx_wk_XJX_XX <- Reduce("+", Map(function(X, W, J, k) {
    (W/k)[1, 1] * (t(X) %*% J %*% X - t(X) %*% X)
  }, X, W, J, k_list))
  switch(modelweights, HIER = {
    tr.sumJJ <- Reduce("+", Map(function(J) sum(diag(J %*% 
                                                       J)), J))
    sumXJX <- Reduce("+", Map(function(X, J) t(X) %*% J %*% 
                                X, X, J))
    sumXWJJX <- Reduce("+", Map(function(X, W, J) t(X) %*% 
                                  W %*% J %*% J %*% X, X, W, J))
    sumXJJWX <- Reduce("+", Map(function(X, W, J) t(X) %*% 
                                  J %*% J %*% W %*% X, X, W, J))
    sumXWWX <- Reduce("+", Map(function(X, W) t(X) %*% W %*% 
                                 W %*% X, X, W))
    sumXJWX <- Reduce("+", Map(function(X, W, J) t(X) %*% 
                                 J %*% W %*% X, X, W, J))
    sumXWJX <- Reduce("+", Map(function(X, W, J) t(X) %*% 
                                 W %*% J %*% X, X, W, J))
  })
  b <- solve(sumXWX) %*% sumXWy
  Xreg <- as.matrix(X.full[-c(1)], dimnames = NULL)
  data.full$pred <- Xreg %*% b
  data.full$e <- data.full$effect.size - data.full$pred
  if (!user_weighting) {
    switch(modelweights, HIER = {
      sumV <- sum(data.full$var.eff.size)
      W <- diag(1/data.full$var.eff.size)
      sumW <- sum(W)
      Qe <- t(data.full$e) %*% W %*% data.full$e
      e <- by(data.full$e, data.full$study, function(x) matrix(x))
      sumEJE <- Reduce("+", Map(function(e, J) t(e) %*% 
                                  J %*% e, e, J))
      Qa <- sumEJE
      V.i <- solve(sumXWX)
      A1 <- tr.sumJJ - sum(diag(V.i %*% sumXJJWX)) - sum(diag(V.i %*% 
                                                                sumXWJJX)) + sum(diag(V.i %*% sumXJX %*% V.i %*% 
                                                                                        sumXWJWX))
      B1 <- length(data.full$study) - sum(diag(V.i %*% 
                                                 sumXWJX)) - sum(diag(V.i %*% sumXJWX)) + sum(diag(V.i %*% 
                                                                                                     sumXJX %*% V.i %*% sumXWWX))
      C1 <- sumV - sum(diag(V.i %*% sumXJX))
      A2 <- sumW - sum(diag(V.i %*% sumXWJWX))
      B2 <- sumW - sum(diag(V.i %*% sumXWWX))
      C2 <- length(data.full$study) - (p + 1)
      omega.sq1 <- ((Qa - C1) * A2 - (Qe - C2) * A1)/(B1 * 
                                                        A2 - B2 * A1)
      omega.sq <- ifelse(omega.sq1 < 0, 0, omega.sq1)
      tau.sq1 <- ((Qe - C2)/A2) - omega.sq * (B2/A2)
      tau.sq <- ifelse(tau.sq1 < 0, 0, tau.sq1)
      data.full$r.weights <- (1/(as.vector(data.full$var.eff.size) + 
                                   as.vector(tau.sq) + as.vector(omega.sq)))
      mod_info <- list(omega.sq = omega.sq, tau.sq = tau.sq)
    }, CORR = {
      W <- diag(data.full$weights)
      sumW <- sum(data.full$weights)
      Qe <- t(data.full$e) %*% W %*% data.full$e
      denom <- sumW - sum(diag(solve(sumXWX) %*% sumXWJWX))
      termA <- sum(diag(solve(sumXWX) %*% Matrx_WKXX))
      termB <- sum(diag(solve(sumXWX) %*% Matrx_wk_XJX_XX))
      term1 <- (Qe - N + termA)/denom
      term2 <- termB/denom
      tau.sq1 <- term1 + rho * term2
      tau.sq <- ifelse(tau.sq1 < 0, 0, tau.sq1)
      df <- N - termA - rho * (termB)
      I.2.1 <- ((Qe - df)/Qe) * 100
      I.2 <- ifelse(I.2.1 < 0, 0, I.2.1)
      data.full$r.weights <- 1/(as.vector(data.full$k) * 
                                  (as.vector(data.full$avg.var.eff.size) + as.vector(tau.sq)))
      mod_info <- list(rho = rho, I.2 = I.2, tau.sq = tau.sq, 
                       term1 = term1, term2 = term2)
    })
  }
  else {
    data.full$r.weights <- data.full$userweights
    mod_info <- list(k = k, N = N, p = p, M = M)
  }
  W.r.big <- diag(data.full$r.weights)
  W.r <- by(data.full$r.weights, data.full$study, function(x) diag(x, 
                                                                   nrow = length(x)))
  sumXWX.r <- Reduce("+", Map(function(X, W) t(X) %*% W %*% 
                                X, X, W.r))
  sumXWy.r <- Reduce("+", Map(function(X, W, y) t(X) %*% W %*% 
                                y, X, W.r, y))
  b.r <- solve(sumXWX.r) %*% sumXWy.r
  data.full$pred.r <- Xreg %*% b.r
  data.full$e.r <- cbind(data.full$effect.size) - data.full$pred.r
  data.full$e.r <- as.numeric(data.full$e.r)
  sigma.hat.r <- by(data.full$e.r, data.full$study, function(x) tcrossprod(x))
  if (!small) {
    sumXWeeWX.r <- Reduce("+", Map(function(X, W, V) t(X) %*% 
                                     W %*% V %*% W %*% X, X, W.r, sigma.hat.r))
    VR.r <- solve(sumXWX.r) %*% sumXWeeWX.r %*% solve(sumXWX.r)
    SE <- sqrt(diag(VR.r)) * sqrt(N/(N - (p + 1)))
    t <- b.r/SE
    dfs <- N - (p + 1)
    prob <- 2 * (1 - stats::pt(abs(t), dfs))
    CI.L <- b.r - stats::qt(0.975, dfs) * SE
    CI.U <- b.r + stats::qt(0.975, dfs) * SE
  }
  else {
    Q <- solve(sumXWX.r)
    Q.list <- rep(list(Q), N)
    H <- Xreg %*% Q %*% t(Xreg) %*% W.r.big
    ImH <- diag(c(1), dim(Xreg)[1], dim(Xreg)[1]) - H
    data.full$ImH <- cbind(ImH)
    ImHj <- lapply(split(x = ImH, f = as.factor(data.full$study)), 
                   function(x) {
                     matrix(x, ncol = M)
                   })
    diag_one <- by(rep(1, M), X.full$study, function(x) diag(x, 
                                                             nrow = length(x)))
    ImHii <- Map(function(X, Q, W, D) D - X %*% Q %*% t(X) %*% 
                   W, X, Q.list, W.r, diag_one)
    if (!user_weighting) {
      Working_Matrx_E <- diag(1/data.full$r.weights)
      Working_Matrx_E_j <- by(data.full$r.weights, data.full$study, 
                              function(x) diag(1/x, nrow = length(x)))
      switch(modelweights, HIER = {
        InsideMatrx_list <- Map(function(W_E_j, ImH_jj) {
          sqrt(W_E_j) %*% ImH_jj %*% (W_E_j^1.5)
        }, Working_Matrx_E_j, ImHii)
        eigenres_list <- lapply(InsideMatrx_list, function(x) eigen(x))
        eigenval_list <- lapply(eigenres_list, function(x) x$values)
        eigenvec_list <- lapply(eigenres_list, function(x) x$vectors)
        A.MBB <- Map(function(eigenvec, eigenval, k, 
                              W_E_j) {
          eigenval_InvSqrt <- ifelse(eigenval < 10^-10, 
                                     0, 1/sqrt(eigenval))
          sqrt(W_E_j) %*% eigenvec %*% diag(eigenval_InvSqrt, 
                                            k, k) %*% t(eigenvec) %*% sqrt(W_E_j)
        }, eigenvec_list, eigenval_list, k_list, Working_Matrx_E_j)
      }, CORR = {
        eigenres_list <- lapply(ImHii, function(x) eigen(x))
        eigenval_list <- lapply(eigenres_list, function(x) x$values)
        eigenvec_list <- lapply(eigenres_list, function(x) x$vectors)
        A.MBB <- Map(function(eigenvec, eigenval, k, 
                              W_E_j) {
          eigenval_InvSqrt <- ifelse(eigenval < 10^-10, 
                                     0, 1/sqrt(eigenval))
          eigenvec %*% diag(eigenval_InvSqrt, k, k) %*% 
            t(eigenvec)
        }, eigenvec_list, eigenval_list, k_list, Working_Matrx_E_j)
      })
    }
    else {
      V.big <- diag(c(1), dim(Xreg)[1], dim(Xreg)[1]) %*% 
        diag(data.full$avg.var.eff.size)
      v.j <- by(data.full$avg.var.eff.size, data.full$study, 
                function(x) diag(x, nrow = length(x)))
      v.j.sqrt_list <- lapply(v.j, function(x) sqrt(x))
      Working_Matrx_E_j <- v.j
      Working_Matrx_E <- V.big
      InsideMatrx_list <- Map(function(ImH_j) {
        ImH_j %*% Working_Matrx_E %*% t(ImH_j)
      }, ImHj)
      eigenres_list <- lapply(InsideMatrx_list, function(x) eigen(x))
      eigenval_list <- lapply(eigenres_list, function(x) x$values)
      eigenvec_list <- lapply(eigenres_list, function(x) x$vectors)
      A.MBB <- Map(function(eigenvec, eigenval, k, v.j.sqrt) {
        eigenval_InvSqrt <- ifelse(eigenval < 10^-10, 
                                   0, 1/sqrt(eigenval))
        v.j.sqrt %*% eigenvec %*% diag(eigenval_InvSqrt, 
                                       k, k) %*% t(eigenvec)
      }, eigenvec_list, eigenval_list, k_list, v.j.sqrt_list)
    }
    sumXWA.MBBeeA.MBBWX.r <- Map(function(X, W, A, S) t(X) %*% 
                                   W %*% A %*% S %*% A %*% W %*% X, X, W.r, A.MBB, 
                                 sigma.hat.r)
    sumXWA.MBBeeA.MBBWX.r <- Reduce("+", sumXWA.MBBeeA.MBBWX.r)
    giTemp <- Map(function(I, A, W, X, Q) t(I) %*% A %*% 
                    W %*% X %*% Q, ImHj, A.MBB, W.r, X, Q.list)
    giTemp <- do.call(rbind, giTemp)
    gi_matrix <- lapply(X = 1:(p + 1), FUN = function(i) {
      matrix(giTemp[, i], nrow = M)
    })
    if (!user_weighting) {
      W.mat <- matrix(rep(1/sqrt(data.full$r.weights), 
                          times = N), nrow = M)
      B_matrix_half <- lapply(X = gi_matrix, FUN = function(gi_mat) {
        W.mat * gi_mat
      })
    }
    else {
      B_matrix_half <- gi_matrix
    }
    B_mat <- lapply(X = B_matrix_half, FUN = tcrossprod)
    B_trace_square <- sapply(X = B_mat, FUN = function(B) {
      (sum(diag(B)))^2
    })
    B_square_trace <- sapply(X = B_mat, FUN = function(B) {
      sum(B * B)
    })
    dfs <- B_trace_square/B_square_trace
    VR.MBB1 <- solve(sumXWX.r) %*% sumXWA.MBBeeA.MBBWX.r %*% 
      solve(sumXWX.r)
    VR.r <- VR.MBB1
    SE <- sqrt(diag(VR.r))
    t <- b.r/SE
    prob <- 2 * (1 - stats::pt(abs(t), df = dfs))
    CI.L <- b.r - stats::qt(0.975, dfs) * SE
    CI.U <- b.r + stats::qt(0.975, dfs) * SE
  }
  reg_table <- data.frame(cbind(b.r, SE, t, dfs, prob, CI.L, 
                                CI.U))
  labels <- c(colnames(X.full[2:length(X.full)]))
  sig <- ifelse(prob < 0.01, "***", ifelse(prob > 0.01 & prob < 
                                             0.05, "**", ifelse(prob > 0.05 & prob < 0.1, "*", "")))
  reg_table <- cbind(labels, reg_table, sig)
  colnames(reg_table) <- c("labels", "b.r", "SE", "t", "dfs", 
                           "prob", "CI.L", "CI.U", "sig")
  if (!small) {
    mod_label_sm <- ""
    mod_notice <- ""
  }
  else {
    mod_label_sm <- "with Small-Sample Corrections"
    mod_notice <- "Note: If df < 4, do not trust the results"
  }
  if (!user_weighting) {
    switch(modelweights, HIER = {
      mod_label <- c("RVE: Hierarchical Effects Model", 
                     mod_label_sm)
    }, CORR = {
      mod_label <- c("RVE: Correlated Effects Model", 
                     mod_label_sm)
    })
  }
  else {
    mod_label <- c("RVE: User Specified Weights", mod_label_sm)
  }
  res <- list(data.full = data.full, X.full = X.full, reg_table = reg_table, 
              mod_label = mod_label, mod_notice = mod_notice, modelweights = modelweights, 
              mod_info = mod_info, user_weighting = user_weighting, 
              ml = ml, cl = cl, N = N, M = M, k = k, k_list = k_list, 
              p = p, X = X, y = y, Xreg = Xreg, b.r = b.r, VR.r = VR.r, 
              dfs = dfs, small = small, data = data, labels = labels, 
              study_orig_id = study_orig_id)
  class(res) <- "robu"
  res
}                   
                   
                   
#===========================# Datasets # ===================================================================================== 
   
table1 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/irr1.csv", row.names = 1)
table2 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/irr2.csv", row.names = 1)          
table3 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/X.csv", row.names = 1)
table5 <- read.csv('https://raw.githubusercontent.com/hkil/m/master/t5.csv', row.names = 1)                 
c1 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/c1.csv")
c2 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/c2.csv")
c3 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/c3.csv")
c4 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/c4.csv")             
             
#================================================================================================================================================================
                                      
                                      
need <- c("bayesmeta", "distr", "robumeta", "ellipse", "zoo", "lavaan", "semPlot", "tidyverse", "weightr", "nlme", "lme4","effects") # "sjPlot", "ggeffects"
not.have <- need[!(need %in% installed.packages()[,"Package"])]
if(length(not.have)) install.packages(not.have)
 
#options(warn = -1)
                    
suppressWarnings(                                         
suppressMessages({ 
  
for(i in need){
  library(i, character.only = TRUE)
}
}))                                        
