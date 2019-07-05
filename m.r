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