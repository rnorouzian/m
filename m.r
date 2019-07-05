d.prepos <- function(d = NA, study.name = NA, group.name = NA, n = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, autoreg = FALSE, t = NA, sdif = NA, sdp = NA, F1 = NA, df2 = NA, post, control, outcome, ...) 
{
  
  if(missing(control) || missing(post) || missing(outcome)) stop("'post', 'outcome' and/or 'control' missing.", call. = FALSE)  
  
  r <- ifelse(autoreg == TRUE & !is.na(r), autoreg(max(post, na.rm = TRUE), r)[,1][-1][post], r)
  
  d <- ifelse(!is.na(d), d, ifelse(!is.na(t) & !missing(n), t2d(t, n), ifelse(!is.na(F1) & !missing(n), t2d(sqrt(F1), n), ifelse(!is.na(F1) & missing(n) & !is.na(df2), t2d(sqrt(F1), df2+2), NA))))
  mdif <- ifelse(!is.na(mpre) & !is.na(mpre), mpos - mpre, NA)
  sdif <- ifelse(is.na(sdif), sdif(sdpre = sdpre, sdpos = sdpos, t = t, r = r, n = n, mpos = mpos, mpre = mpre, F1 = F1, sdp = sdp), sdif)
  cor. <- ifelse(is.na(r), rdif(n = n, mpre = mpre, mpos = mpos, sdpre = sdpre, sdpos = sdpos, sdif = sdif, sdp = sdp), r)
  d <- ifelse(!is.na(mdif) & is.na(d) & !is.na(sdif), mdif/sdif, d)
  
  se <- se.d(d, n1 = n, g = TRUE)
  
  d <- d*cfactor(n-1)
  
  out <- data.frame(d = d, n = n, sdif = sdif, rpr.po = cor., post, control, outcome, ...)
  
  #if(!anyNA(group.name) & length(group.name) == nrow(out)) row.names(out) <- as.character(group.name) else if(!anyNA(group.name) & length(group.name) != nrow(out)) stop("'group.name' incorrectly specified.", call. = FALSE)
  
  if(all(is.na(out$d))) stop("\ninsufficient info. to calculate effect size(s).", call. = FALSE)
  
  return(out)
  
}

#================================================================================================================================


dint <- function(..., per.study = NULL, study.name = NA, group.name = NA, n.sim = 1e5, digits = 6, by, data = NULL)
{
  
  L <- if(!is.null(data)){
    
    m <- split(data, data$study.name)        
    
    m[[1]] <- NULL                          
    
    if(is.null(reget(m, control, F))) stop("Required 'control' group not found.", call. = FALSE)
    
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
  
  G <- function(m, study.name, group.name, n.sim, digits)
  {
    
    cdel1 <- reget(m, control & post == 2 & outcome == 1, F)
    cdel2 <- reget(m, control & post == 3 & outcome == 1, F)
    cs <- reget(m, control & post == 1 & outcome == 1, F)
    
    tdel1 <- reget(m, !control & post == 2 & outcome == 1, F)
    tdel2 <- reget(m, !control & post == 3 & outcome == 1, F)
    ts <- reget(m, !control & post == 1 & outcome == 1, F) 
    
    if(all(sapply(list(cdel1, cdel2, tdel1, tdel2, ts, cs), is.null))) stop("Either 'control' or 'post' incorrectly coded.", call. = FALSE)
    
    short <- all(sapply(list(cs, ts), function(x) !is.null(x)))
    
    del1 <- all(sapply(list(cdel1, tdel1), function(x) !is.null(x)))
    
    del2 <- all(sapply(list(cdel2, tdel2), function(x) !is.null(x)))
    
    
    cdel1..2 <- reget(m, control & post == 2 & outcome == 2, F)
    cdel2..2 <- reget(m, control & post == 3 & outcome == 2, F)
    cs..2 <- reget(m, control & post == 1 & outcome == 2, F)
    
    tdel1..2 <- reget(m, !control & post == 2 & outcome == 2, F)
    tdel2..2 <- reget(m, !control & post == 3 & outcome == 2, F)
    ts..2 <- reget(m, !control & post == 1 & outcome == 2, F)


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
      SHORT <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim, digits = digits)))
     # row.names(SHORT) <- group.name1
    }
    
    
    if(short..2){
      nc1 <- m$n[m$control & m$post == 1 & m$outcome == 2]
      nt1 <- m$n[m$control == FALSE & m$post == 1 & m$outcome == 2]
      dps <- pair(cs..2, ts..2)  
      dppc1 <- sapply(1:length(dps), function(i) dps[[i]][[1]][1])
      dppt1 <- sapply(1:length(dps), function(i) dps[[i]][[1]][2])
      #group.name1 <- unlist(lapply(1:length(dps), function(i) names(dps[[i]])))
      SHORT..2 <- data.frame(t(dit(dppc = dppc1, dppt = dppt1, nc = nc1, nt = nt1, n.sim = n.sim, digits = digits)))
     # row.names(SHORT..2) <- group.name1
    }
    
    
    if(del1){
      nc2 <- m$n[m$control & m$post == 2 & m$outcome == 1]
      nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 1]
      dpdel1 <- pair(cdel1, tdel1)
      dppc2 <- dppcdel1 <- sapply(1:length(dpdel1), function(i) dpdel1[[i]][[1]][1])
      dppt2 <- dpptdel1 <- sapply(1:length(dpdel1), function(i) dpdel1[[i]][[1]][2])
      #group.name2 <- unlist(lapply(1:length(dpdel1), function(i) names(dpdel1[[i]])))
      DEL1 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim, digits = digits)))
      #row.names(DEL1) <- group.name2
    }
    
    
    if(del1..2){
      nc2 <- m$n[m$control & m$post == 2 & m$outcome == 2]
      nt2 <- m$n[m$control == FALSE & m$post == 2 & m$outcome == 2]
      dpdel1 <- pair(cdel1..2, tdel1..2)
      dppc2 <- sapply(1:length(dpdel1), function(i) dpdel1[[i]][[1]][1])
      dppt2 <- sapply(1:length(dpdel1), function(i) dpdel1[[i]][[1]][2])
     # group.name2 <- unlist(lapply(1:length(dpdel1), function(i) names(dpdel1[[i]])))
      DEL1..2 <- data.frame(t(dit(dppc = dppc2, dppt = dppt2, nc = nc2, nt = nt2, n.sim = n.sim, digits = digits)))
      #row.names(DEL1..2) <- group.name2
    }
    
    
    if(del2){
      nc3 <- m$n[m$control & m$post == 3 & m$outcome == 1]
      nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 1]
      dpdel2 <- pair(cdel2, tdel2)
      dppc3 <- dppcdel2 <- sapply(1:length(dpdel2), function(i) dpdel2[[i]][[1]][1])
      dppt3 <- dpptdel2 <- sapply(1:length(dpdel2), function(i) dpdel2[[i]][[1]][2])
    #  group.name3 <- unlist(lapply(1:length(dpdel2), function(i) names(dpdel2[[i]])))
      DEL2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, digits = digits)))
      #row.names(DEL2) <- group.name3
    }
    
    if(del2..2){
      nc3 <- m$n[m$control & m$post == 3 & m$outcome == 2]
      nt3 <- m$n[m$control == FALSE & m$post == 3 & m$outcome == 2]
      dpdel2 <- pair(cdel2..2, tdel2..2)
      dppc3 <- sapply(1:length(dpdel2), function(i) dpdel2[[i]][[1]][1])
      dppt3 <- sapply(1:length(dpdel2), function(i) dpdel2[[i]][[1]][2])
      #group.name3 <- unlist(lapply(1:length(dpdel2), function(i) names(dpdel2[[i]])))
      DEL2..2 <- data.frame(t(dit(dppc = dppc3, dppt = dppt3, nc = nc3, nt = nt3, n.sim = n.sim, digits = digits)))
      #row.names(DEL2..2) <- group.name3
    }
    
    
    list(SHORT = if(short) SHORT else NULL, SHORT..2 = if(short..2) SHORT..2 else NULL, DEL1 = if(del1) DEL1 else NULL, DEL1..2 = if(del1..2) DEL1..2 else NULL, DEL2 = if(del2) DEL2 else NULL, DEL2..2 = if(del2..2) DEL2..2 else NULL) 
  }
  
  h <- lapply(1:length(L), function(i) G(m = L[[i]], study.name = study.name, group.name = group.name, n.sim = n.sim, digits = digits))
  
  if(!is.null(data)) study.name <- names(L)
  
  names(h) <- if(anyNA(study.name)) paste0("Study", seq_along(h)) else if(!anyNA(study.name) & length(study.name) == length(h)) study.name else if(!anyNA(study.name) & length(study.name) != length(h)) stop("'study.name' incorrectly specified.", call. = FALSE)
  
  return(h)
  
}
              

#================================================================================================================================
              
              
              
meta.within <- function(..., per.study = NULL, study.name = NA, group.name = NA, tau.prior = function(x){dhnorm(x)}, by, data = NULL){
  
  L <- eval(substitute(dint(... = ..., per.study = per.study, group.name = group.name, study.name = study.name, by = by, data = data)))
  
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
      result1 <- bayesmeta(     y = d1,
                                sigma = sd1,
                                labels = NULL, tau.prior = tau.prior)
      result1$call <- match.call(expand.dots = FALSE)
      
      short <- c(result1$summary["mean","mu"], result1$summary["sd","mu"])
      
    }
    
    
    if(Short & length(d1) > 1 & Short..2 & length(d1..2) > 1) { 
      
      
      res1 <- bayesmeta(         y = d1..2,
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
      
      
      resi3 <- bayesmeta(         y = ds,
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
              

meta.bayes <- function(..., per.study = NULL, group.name = NA, study.name = NA, tau.prior = function(x){dhnorm(x)}, by, long = FALSE, data = NULL)
{
  

j <- eval(substitute(meta.within(... = ..., per.study = per.study, group.name = group.name, study.name = study.name, tau.prior = tau.prior, by = by, data = data)))
  
if(!is.null(data)) study.name <- names(j)

if(anyNA(study.name)) study.name <- NULL


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
               
               
               
              
