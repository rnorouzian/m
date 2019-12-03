

#=========== SLI PAPER DATA-ANALYSIS ==============#

# By Component 5 team members: Reza & Rudy

# Access Reza's suites of new programs:
source("https://raw.githubusercontent.com/rnorouzian/m/master/m.r")
source("https://raw.githubusercontent.com/rnorouzian/i/master/i.r")


# Access the SLI paper Data:
D <- read.csv("https://raw.githubusercontent.com/rnorouzian/i/master/SLI.csv", h = T)


# total participants:
total <- nrow(D)         

# non-responses:
none <- nrow(na.omit(D)) 


#Percetange of data missing:
mis <- noquote(paste0(round(1 - none/total, 4)*1e2, "%")) # 12.68% of responses missing


# Use Reza's newly developed program to recover the missing responses:

D <- impute(D)   # impute the missing using median responses

pre <- D[,1:8]
pos <- D[,9:16]

# Use Reza's newly developed program to measure latent constructs (pre-post):

a <- (efa(pre, factors = 1, scores = "reg")$score)$Fa  # measure pre latent construct (educational leadership)
b <- (efa(pos, factors = 1, scores = "reg")$score)$Fa  # measure pos latent construct (educational leadership)


# run a paired t-test on the latent contstructs before & after SLI rather than items:
test <- t.test(a, b, paired = TRUE)  

t.val <- unname(test$statistic)

N <- total

# Use Reza's newly developed program to measure effect size (pre-post):

cohen.d  <- t2d(t.val, N) # measure effect size of change after SLI
# a d of "0.5896104"


# Use Reza's newly developed program to measure 95% confidence interval for effect size (pre-post):
(ci <- d.ci(cohen.d, n1 = N))

#     Cohen.d     lower     upper     conf.level      ncp
#    0.5896104 0.3356008  0.8400179       0.95      4.968146

# Use Reza's newly developed program to interpret effect size in percentages (pre-post):
dint.norm(c(ci$lo, cohen.d, ci$up)) # 'd' shows SLI has made 22.23% improvement.


# Show the conceptual model:

m2 <- " 
LBpr = ~Q1_a+Q6_a
BAVpr= ~Q2_a+Q3_a
OSpr = ~Q4_a+Q5_a
EFpr = ~Q7_a+Q8_a
LBps = ~Q1_b+Q6_b
BAVps= ~Q2_b+Q3_b
OSps = ~Q4_b+Q5_b
EFps = ~Q7_b+Q8_b
EdLpr =~ LBpr + BAVpr + OSpr+ EFpr
EdLps =~ LBps + BAVps + OSps+ EFps
"


library(semPlot)
library(lavaan)
library(ReporteRs)


# fit the conceptual model:

fit2 <- cfa(m2, data = D)

# Create a figure of the conceptual model
G2 <- function() semPaths(fit2, layout = "spring", residuals = F)

G2()

# print the conceptual model for publication:

doc2 <- addPlot(docx(), fun = G2, vector.graphic = TRUE, width = 3.7, height = 4.1,  
                par.properties = parCenter(), editable = TRUE)


writeDoc(doc2, file = "CFA2.docx")


# Use Reza's program to generalize the results to a wider pool of participants:

nor <- function(di = cohen.d, ri = c(1e-6, .999999)){

xx <- range(c(qnorm(ri), qnorm(ri, di)), finite = TRUE)

a <- curve(dnorm(x), xx[1], xx[2], n = 1e4, xlab = "Effect Size", ylab = "Density")

b <- curve(dnorm(x, di), add = TRUE, col = 2, n = 1e4, lty = 2)
axis(1, at = di, labels = round(di, 3), col = 2, col.axis = 2)

u <- par('usr')[3]

segments(x <- c(0, di), c(u, u), x, c(max(a$y, na.rm = T), max(b$y, na.rm = T)), col = c(1, 2), lwd = 2, lend = 1)

xy <- a$x >= 0 & a$x <= di

x <- c(0, a$x[xy], di)
y <- c(u, pmin.int(a$y, b$y)[xy], u)

polygon(x, y, col = 4, border = NA, density = 15)
}


nor()
