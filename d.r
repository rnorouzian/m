
#================= Duolingo English Test (CAT Section):

### Response (xi) is on continous [0,1] scale.
### Two possible transformations include logit transformation & dyadic expansion.

### logit transformation and its linearizing effect on responses are well known.
### logit transformation also allows highly flexible SEM measurement models.

#============= raw dataset:


# > R version 3.6.1 (2019-07-05)

# Source my utility functions (e.g., data trimming, masking id variables etc.):

source('https://raw.githubusercontent.com/rnorouzian/m/master/m.r')

# R packages used:
library(lavaan)
library(semPlot)
library(lattice)
library(latticeExtra)

# 'trim' all potential 'white spaces' before or after the string variables + column names
dat1 <- trim(read.csv('https://raw.githubusercontent.com/hkil/m/master/task.csv', stringsAsFactors = F))


# Drop columns not working with for this task only:
dat2 <- drop.col(dat1, c('country', 'region', 'item_level', 'test_score', 'native_language', 'age'))


# Make easy id variables for tests and items:
dat <- mask(dat2, c('test_id', 'item_id'), full = T)


# rename 'test_id' & 'item_grade':
names(dat)[c(1,4)] <- c('person_id', 'item_score')

# see how sparse the item-person matrix is:
tab <- table(dat$item_id)
summary(as.data.frame(tab))


# plot items use frequency to see skewness (natural in CAT delievery):
# - Find the number of most "exposed items" (given to >= 100 test-takers):
expo.items <- 1:3*1e2
tx <- paste('items >=', expo.items)
ns <- sapply(expo.items, function(i) length(tab[tab >= i]))

# plot item exposure frequency             
plot(tab, xlab = "Items", ylab = "Exposure Frequency", lwd = 1)
legend('topright', paste(tx, '=', ns), bty = 'n', text.col = 4, cex = .8)

# So, due to extreme sparsity, item-specific difficulty estimation is not possible.
# As one option, first logit transform responses,
# Then, make composite (synthetic) variables based on item_type.
# Finally, try measurement models with robust standard errors uing SEM framework

#==============================================================================
## Widen the dataset for SEM package 'lavaan' and expose sparsity concentration:
#==============================================================================
             
# sort the data by person_id, item_type and item_id (change dataset name to 'w2'):
w2 <- dat[order(dat$person_id, dat$item_type, dat$item_id),]

# create item_type sequence (1st ctest, 2nd ctest etc.):
w2$item_seq <- with(w2, ave(seq_along(person_id), person_id,
                            item_type, FUN = seq_along))

# Add numeric suffix for item sequences:             
w2$item_type <- sprintf("%s_%d",w2$item_type,w2$item_seq)

# Reshape and drop temprary index variables:             
w2 <- reshape(w2, dir = 'wide', drop = c("item_id","item_seq"), 
              idvar = c('person_id','gender'),
              timevar = 'item_type')

# clean up column names:
names(w2) <- gsub("item_score.","", names(w2))


# rename dataset to op1:
op1 <- w2

# exclude non-item columns
unwant <- 1:2

# tranform to logit:
op1[-unwant] <- lapply(op1[-unwant], stats::qlogis)


# Visualize values after transformation:
boxplot(op1[-unwant], cex.axis = .7)

# SEM computes robust Standard Errors to account for non-normalities

## We are planning to control for 'gender'.
# SEM packages require a dummy variable for 'gender', let's add a dummy to the dataset:

gender <- op1$gender
dummy <- data.frame(model.matrix(~gender))[-1]

op1 <- cbind(op1, dummy)

#======================================================================
### Make Measurement Models Using SEM (Structural Equation Modeling)
#======================================================================
             
############################
## ASSUMPTIONS AND DECISIONS   
############################
             
# Exclude extremely sparse 7th sequence of item types to avoid convergence issues (over 80% of test takers didn't take the 7th ctest etc.)
# Robust estimators accomadete the extreme observations.

## In these SEM models, 'intecepts' function as difficulty parameters and, 
## 'path coeffiecients' function as discrimination parameters,
## 'latent factor scores' represent ability parameters.

## Test a first measurement model and evaluate the fit:
# This model assumes DET CAT section consists of uncorrelated item types that
# over a sequence of presentations to test takers measure their subskills i.e., first-order composite variables.
# Yet at a higher level, these first-order composites themselves in combination form the DET representing 
# "Proficiency" as measured in its CAT section.


m1 <- '
au_v=~ audio_vocab_1+audio_vocab_2+audio_vocab_3+audio_vocab_4+audio_vocab_5+audio_vocab_6
ctst=~ ctest_1+ctest_2+ctest_3+ctest_4+ctest_5+ctest_6
dict=~ dictation_1+dictation_2+dictation_3+dictation_4+dictation_5+dictation_6
el_s=~ elicited_speech_1+elicited_speech_2+elicited_speech_3+elicited_speech_4+elicited_speech_5+elicited_speech_6
tx_v=~ text_vocab_1+text_vocab_2+text_vocab_3+text_vocab_4+text_vocab_5+text_vocab_6

DET = ~ NA*au_v+ctst+dict+el_s+tx_v
DET ~~ 1*DET
 au_v ~~ 1*au_v
 ctst ~~ 1*ctst
 dict ~~ 1*dict
 el_s ~~ 1*el_s
 tx_v ~~ 1*tx_v
'

# Fit the first model:
m1 <- sem(m1, data = op1, meanstructure = TRUE, missing = 'ML', estimator = "MLF")

## Visualize the model:
semPlot::semPaths(m1)

# Evaluate the fit
summary(m1, fit.measures=TRUE) 

# CFI = .804     # Not very good fit comapred to 'null' model (model with no covariation among any variables)
# RMSEA = 0.044  # Acceptable fit compared to 'saturated' model (model with all variables correlated with one another)

## Conclusion: The model needs further improvement to obtain very good fit (e.g., CFI > .95)

## How? Consult modification idicies and then consider them theoritically:

## Modification indices:
mi1 <- modindices(m1, sort. = T, minimum.value = 3) 

# dication_1 and text_vocab_1 are suggested to relate with auto_vocab first-order composite variable 
# Similarly, text_vocab_1 is suggested to relate with dictation composite variable 
# Finally, dictation_1 is suggested to relate with text_vocab composite variable.
# These changes are implemented in model 2:

#==== Model 2:

m2 <- '

au_v=~ NA*audio_vocab_1+audio_vocab_2+audio_vocab_3+audio_vocab_4+audio_vocab_5+audio_vocab_6+text_vocab_1+dictation_1 
ctst=~ ctest_1+ctest_2+ctest_3+ctest_4+ctest_5+ctest_6
dict=~ dictation_1+dictation_2+dictation_3+dictation_4+dictation_5+dictation_6+text_vocab_1
el_s=~ elicited_speech_1+elicited_speech_2+elicited_speech_3+elicited_speech_4+elicited_speech_5+elicited_speech_6
tx_v=~ NA*text_vocab_1+text_vocab_2+text_vocab_3+text_vocab_4+text_vocab_5+text_vocab_6+dictation_1

DET = ~ NA*au_v+ctst+dict+el_s+tx_v
DET ~~ 1*DET
 au_v ~~ 1*au_v
 ctst ~~ 1*ctst
 dict ~~ 1*dict
 el_s ~~ 1*el_s
 tx_v ~~ 1*tx_v
 '

# Fit the model:             
m2 <- sem(m2, data = op1, meanstructure = TRUE, missing = "ML", estimator = "MLF")

# The fit has drastically improved:
summary(m2, fit.measures=TRUE)  ## CFI 0.964
## RMSEA 0.019

## Modification indicies suggest possible small improvements but because of large sample size these small change maybe sig.
mi2 <- modindices(m2, sort. = TRUE)

## Visualize the new model:
semPlot::semPaths(m2)


# To honor model 2's suggestion, we relate 'ctest_1' to 'au_vocab' composite variable 
# leading to model 3:

#==== Model 3:

m3 <- '

au_v=~ NA*audio_vocab_1+audio_vocab_2+audio_vocab_3+audio_vocab_4+audio_vocab_5+audio_vocab_6+text_vocab_1+dictation_1+ctest_1 
ctst=~ ctest_1+ctest_2+ctest_3+ctest_4+ctest_5+ctest_6
dict=~ dictation_1+dictation_2+dictation_3+dictation_4+dictation_5+dictation_6+text_vocab_1
el_s=~ elicited_speech_1+elicited_speech_2+elicited_speech_3+elicited_speech_4+elicited_speech_5+elicited_speech_6
tx_v=~ NA*text_vocab_1+text_vocab_2+text_vocab_3+text_vocab_4+text_vocab_5+text_vocab_6+dictation_1

DET = ~ NA*au_v+ctst+dict+el_s+tx_v
DET ~~ 1*DET
 au_v ~~ 1*au_v
 ctst ~~ 1*ctst
 dict ~~ 1*dict
 el_s ~~ 1*el_s
 tx_v ~~ 1*tx_v
 '

# Fit the model:               
m3 <- sem(m3, data = op1, meanstructure = TRUE, missing = "ML", estimator = "MLF")

# The fit has slightly improved:
summary(m3, fit.measures = TRUE)  ## CFI 0.968
## RMSEA 0.018

## Visualize the new model:
semPlot::semPaths(m3)


# Since models 2 and 3 are nested, a chi-sq difference test can determine their advantage:
anova(m2, m3)

# Model 3 seems to be a very good-fitting measurement model with sig. advantage over model 2!


#=== But given excellent fit of model 3, how DET CAT section, let's add a covariate
## - In model 4, let's control for 'gender' (the 'other' category is very small so perhaps non-significant)

#=== Model 4:
             
m4 <- '

au_v=~ NA*audio_vocab_1+audio_vocab_2+audio_vocab_3+audio_vocab_4+audio_vocab_5+audio_vocab_6+text_vocab_1+dictation_1+ctest_1 
ctst=~ ctest_1+ctest_2+ctest_3+ctest_4+ctest_5+ctest_6
dict=~ dictation_1+dictation_2+dictation_3+dictation_4+dictation_5+dictation_6+text_vocab_1
el_s=~ elicited_speech_1+elicited_speech_2+elicited_speech_3+elicited_speech_4+elicited_speech_5+elicited_speech_6
tx_v=~ NA*text_vocab_1+text_vocab_2+text_vocab_3+text_vocab_4+text_vocab_5+text_vocab_6+dictation_1

DET = ~ NA*au_v+ctst+dict+el_s+tx_v
DET ~~ 1*DET
 au_v ~~ 1*au_v
 ctst ~~ 1*ctst
 dict ~~ 1*dict
 el_s ~~ 1*el_s
 tx_v ~~ 1*tx_v
 
 DET ~ genderMALE
 DET ~ genderOTHER
 
 '
# Fit the model
m4 <- sem(m4, data = op1, meanstructure = TRUE, missing = "ML", estimator = "MLF")

# Summary of results:
summary(m4, fit.measures = TRUE)

# visualize the gender model
semPlot::semPaths(m4, thresholds = F)

# based on these simulated data:
# Females are shown to have overperformed males and other (other's n = 14) 


#==== Deriving the estimates from model 4:
             
parm <- parameterEstimates(m4)
item.diff <- setNames(parm$est[82:111], parm[82:111,1])    
x <- seq_len(length(item.diff))
nms <- names(item.diff)
it.dif <- item.diff

# clean name suffixes for a final data.frame of item diff. estimates:             
names(it.dif) <- sub("_\\d+$", "", names(it.dif))


# Estimated person abilities:
overall <- data.frame(lavPredict(m4))

# Ability empirical distribution:             
ability <- density(overall$DET)

         
# Plot difficulty parameters against test-takers ability:
a <- lattice::xyplot(x*.015~it.dif| names(it.dif), col = 'magenta', type = 'h', xlab = bquote(theta), ylab = 'Density')
a + layer(panel.polygon(ability, col = adjustcolor('orange', .3)))


#Make a final data.frame of difficulty parameters:
dif <- data.frame(split(unname(item.diff), sub("_\\d+$", "", names(item.diff))))














