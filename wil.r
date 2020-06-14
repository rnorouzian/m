
#-- Williamson's test statistic and sample size:
t.wil <- 4.75 
n.wil <- 166


#-- p-value obtained correctly & incorrectly:
wrong <- 2*pnorm(t.wil, lower.tail = FALSE)
right <- 2*pt(t.wil, df = n.wil-1, lower.tail = FALSE)


#-- Ratio of right to wrong p-value:
right/wrong