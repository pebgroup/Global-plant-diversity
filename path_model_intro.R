# R script testing the lavaan tutorial available on
# https://lavaan.ugent.be/tutorial/tutorial.pdf

library(lavaan)

# these variables are all latent variables
HS.model <- "visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9"

fit <- cfa(HS.model, data=HolzingerSwineford1939)

# cfa() s used for confirmatory factor analysis models
## It is used to test whether measures of a construct are consistent with a researcher's understanding of the nature of that construct (or factor). As such, the objective of confirmatory factor analysis is to test whether the data fit a hypothesized measurement model. 

summary(fit)
summary(fit, fit.measures=TRUE)

# EVALUATION: a good model fit means only that the model is plausible, not that it fits well.
#  The RMSEA ranges from 0 to 1, with smaller values indicating better model fit. A value of .06 or less is indicative of acceptable model fit (doi:10.1080/10705519909540118. , Brown, Timothy (2015). Confirmatory factor analysis for applied research. New York London: The Guilford Press. p. 72. ISBN 978-1-4625-1779-4.)
# CFI values range from 0 to 1, with larger values indicating better fit. Previously, a CFI value of .90 or larger was considered to indicate acceptable model fit.[40] However, recent studies have indicated that a value greater than .90 is needed to ensure that misspecified models are not deemed acceptable (Hu & Bentler, 1999). Thus, a CFI value of .95 or higher is presently accepted as an indicator of good fit (Hu & Bentler, 1999).



model <- "
  # measurement model
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8  
    y6 ~~ y8"

fit2 <- sem(model, data=PoliticalDemocracy)
summary(fit, standardized=TRUE)


# SEM plot
library(semPlot)
semPaths(fit)
semPaths(fit, whatLabels = c("est"))


semPaths(fit, whatLabels = "est",
         edge.label.cex = 1.2,
         layout = "circle")
