pt.curve <- function(X, adjust = 1, compact = .3, pch = 16, col = 2, cex = .7, ...) {
  
  n.target <- length(X)
  
  d <- density.default(X, adjust = adjust, n = n.target, na.rm = TRUE)
  
  auc <- sum(d$y*median(diff(d$x)))/(diff(range(d$x))*max(d$y))
  
  n <- compact*ceiling(n.target/auc)
  
  pts <- data.frame(x = runif(n, min(d$x), max(d$x)), y = runif(n, 0, max(d$y)))
  
  pts <- pts[pts$y < approx(d$x, d$y, xout = pts$x)$y, ]
  
  pts <- pts[sample(seq_len(nrow(pts)), n, replace = TRUE), ]
  
  plot(pts, pch = pch, col = col, cex = cex, ...)
}

#==============================================================================

#Question 1: Suppose a human skill can be measured by 15 questions that are each worth up to 2 points (i.e., partial credit allowed), 
#what could be a likely distribution of overall scores (i.e., sum of 15 questions) of 5000 randomly selected test-takers?

# USE R to find an answer:
add.norm <- function(n.test.taker = 5e3, n.question = 15, pt.worth = 2, compact = .1){
  
  dis <- replicate(n.test.taker, sum(runif(n.question, 0, pt.worth)))
  
  pt.curve(dis, ylab = "Density", main = "Score Distribution", xlab = "Test Scores", compact = compact)
  abline(v = mean(dis))
}


#Question 2: Suppose the growth rate of a human skill (e.g., driving) is enhanced by 5 other interacting subskills (i.e., their effects multiply), each of which can stimulate the growth by a small percentage (e.g., 1% maximum), 
# what can be a likely distribution for that human skill’s growth rate among 5000 randomly selected subjects?

# USE R to find an answer:
mult.norm <- function(n.subject = 5e3, n.subskill = 5, max.small.growth = .01, compact = .1){
  
  growth <- replicate(n.subject, prod(1 + runif(n.subskill, 0, max.small.growth)))
  
  pt.curve(growth, ylab = "Density", main = "Growth Rate Distribution", xlab = "Growth Rate", compact = compact)
  abline(v = mean(growth))
}


#Question 3: Go back to question 2. Now, suppose that each interacting subskill can stimulate the growth rate of driving skill by a much higher percentage (e.g., 100% max). 
# (A) Is the distribution of driving skill’s growth rate among 5000 randomly selected subjects the same as in question 2?
# (B) Suppose we are open to re-scale growth rate of driving skill to logrithmic scale (`log = TRUE`) as a way to reduce the size of the growth rates,
#     Now what is the likely distribution of the log-scaled growth rates?


mult.norm.free <- function(n.subject = 5e3, n.subskill = 5, max.big.growth = 1, compact = .1, log = FALSE){
  
  growth <- replicate(n.subject, if(log)log(prod(1 + runif(n.subskill, 0, max.big.growth))) else prod(1 + runif(n.subskill, 0, max.big.growth)))
  
  pt.curve(growth, ylab = "Density", main = if(log)"Log Growth Rate Distribution" else "Growth Rate Distribution", xlab = if(log)"Log Growth Rate" else "Growth Rate", compact = compact)
  abline(v = mean(growth))
}
           
