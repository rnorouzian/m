
#Question 1: Suppose a human skill can be measured by asking 15 questions that are each worth 2 points, 
# what can be a likely distribution of overall scores (i.e., sum of 15 questions) for 5000 randomly selected test takers?

# USE R to find an answer:
add.norm <- function(n.test.taker = 5e3, n.question = 15, pt.worth = 2){
  
  dis <- replicate(n.test.taker, sum(round(runif(n.question, 0, pt.worth))))
  
  plot(density(dis, adjust = 1.2), main = "Score Distribution", xlab = "Test Scores")
  mean(dis)
}


#Question 2: Suppose the growth rate of a human skill is enhanced by 5 other interacting subskills (i.e., their effects multiply), each of which can stimulate the growth by a small percentage (e.g., 1% maximum), 
#what can be a likely distribution for that human skill’s growth rate among 5000 randomly selected subjects?

# USE R to find an answer:
mult.norm <- function(n.subject = 5e3, n.subskill = 5, max.small.growth = .01){
  
  growth <- replicate(n.subject, prod(1 + runif(n.subskill, 0, max.small.growth)))
  
  plot(density(growth), main = "Growth Distribution", xlab = "Growth Rate")
  mean(growth)
}



mult.norm.free <- function(n.subject = 5e3, n.subskill = 5, max.big.growth = 2){
  
  growth <- replicate(n.subject, log(prod(1 + runif(n.subskill, 0, max.big.growth))))
  
  plot(density(growth), main = "Growth Distribution", xlab = "Growth Rate")
  mean(growth)
}
