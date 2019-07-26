
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#**** Run each line of code below by placing cursor next to the line, and clicking "Run" (top-right corner) **** 
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>



source("https://raw.githubusercontent.com/rnorouzian/m/master/m.r")                 # the software

D <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/k.csv", h = T) # The data



#=================General Meta-Analysis======================= 

p <- meta.bayes(D)  # Bayesian meta-analysis: is WCF effective?


forestplot(p$SHORT) # distribution of short-term effects


forestplot(p$DEL1)  # distribution of long-term effects


funnel(p$SHORT)     # discover publication bias (short-term effects)


funnel(p$DEL1)      # discover publication bias (long-term effects)


dint.plot(p)        # assess the short- and long-term meta-analystic results


#==================Specific Meta-Analysis===========================


s <- meta.bayes(D, ESL == 2) # what is the effct of WCF in 'ESL' settings?


forestplot(s$SHORT)         # distribution of short-term effects


forestplot(s$DEL1)          # distribution of long-term effects


funnel(s$SHORT)             # discover publication bias (short-term effects)


funnel(s$DEL1)              # discover publication bias (long-term effects)


dint.plot(p, s)             # compare the short- and long-term meta-analystic results with general meta-nanalytic results above


#==================Assess theoretical Positions===========================


u <- meta.bayes(D, type == 2 & prof == 2)  # Ferris(2003): 'Direct' WCF has a postive effect on 'beggining' learners.


forestplot(u$SHORT)                        # distribution of short-term effects


forestplot(u$DEL1)                         # distribution of long-term effects


funnel(u$SHORT)                            # discover publication bias (short-term effects)


funnel(u$DEL1)                             # discover publication bias (long-term effects)


dint.plot(p, s, u)                         # compare the short- and long-term meta-analystic results with general and specific meta-nanalytic results above
