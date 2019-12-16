
# Suppose a human skill can be measured by asking 15 questions that are each worth 2 points, 
# what can be a likely distribution of overall scores of 5000 randomly selected test-takers?

# USE R to find an answer:
add.norm <- function(n.test.taker = 5e3, n.question = 15, pt.worth = 2){
  
  dis <- replicate(n.test.taker, sum(round(runif(n.question, 0, pt.worth))))
  
  plot(density(dis, adjust = 1.2), main = "Score Distribution", xlab = "Test Scores")
  mean(dis)
}


# Suppose growth of a human skill depends on small growth in 5 other interacting subskills (their effects multiply), 
# what can be that human skill's growth distribution among 10000 randomly selected subjects? 

# USE R to find an answer:
mult.norm <- function(n.subject = 5e3, n.subskill = 5, max.small.growth = .1){
  
  growth <- replicate(n.subject, prod(1 + runif(n.subskill, 0, max.small.growth)))
  
  plot(density(growth), main = "Growth Distribution", xlab = "Growth Rate")
}