#### Script information #### 

# What factors influence the activity of invasive mammalian 
# predators and other native mammals in the eastern Otway Ranges?

# Authors: Darcy Watchorn et al.

# Collaborators: Tim Doherty, Don Driscoll

# ____________________________________________________________________________________


#### Packages ####

library(corrplot)
library(lme4)
library(effects)
library(AICcmodavg)
library(MuMIn)
library(DHARMa)
library(rstanarm)
library(sjPlot)
library(sjmisc)
library(MASS)

setwd("C:/Users/dwatchorn/OneDrive - Deakin University/PhD/PhD/Chapters/Chapter 1/Statistical analysis/GLMMs_June2023")


# ____________________________________________________________________________________

#### Collinearity check ####

# Check for collinearity among predictor variables
Covs = read.csv("CovsForCollinTest.csv", stringsAsFactors = T)
predictor_vars <- Covs[, c("NumberOfFires100Yrs", 
                           "YrsSincePreviousFire", 
                           "HabitatType", 
                           "DistanceNearestTown",
                           "DistanceNearestRoad",
                           "DistanceNearestFarm",
                           "Elevation",
                           "LocalPlacement",
                           "LandscapePlacement",
                           "TPI100",
                           "TPI1000",
                           "BurnTreatment",
                           "PercentAreaBurnt",
                           "DistanceBurnEdge",
                           "NDVI100",
                           "PreyActivityLargeMammal",
                           "PreyActivityMediumBird",
                           "PreyActivityMediumMammal",
                           "PreyActivitySmallBird",
                           "PreyActivitySmallMammal",
                           "AvgNightsSinceBaited")]

# Calculate the correlation matrix
cor_matrix <- cor(predictor_vars)

# Print the correlation matrix
print(cor_matrix)

# Visualize the correlation matrix as a heatmap
corrplot(cor_matrix, method = "color")

# Set a threshold for identifying strong correlations
threshold <- 0.65

# Find highly correlated variable pairs
high_correlation <- which(abs(cor_matrix) > threshold & cor_matrix < 1 & !is.na(cor_matrix), arr.ind = TRUE)

# Print the pairs of highly correlated variables
if (nrow(high_correlation) > 0) {
  print("Highly correlated variables:")
  for (i in 1:nrow(high_correlation)) {
    var1 <- rownames(cor_matrix)[high_correlation[i, 1]]
    var2 <- colnames(cor_matrix)[high_correlation[i, 2]]
    print(paste(var1, "and", var2))
  }
} else {
  print("No highly correlated variables found.")
}

# Only PercentAreaBurnt and BurnTreatment are highly correlated (>0.7). 


# ____________________________________________________________________________________


#### Red fox ####

# Create the dataframe
FoxRR = read.csv("SiteSpecificCovariates.csv", stringsAsFactors = T)
str(FoxRR)

hist(FoxRR$FoxDet, 20) # Lots of zeroes
nrow(subset(FoxRR, FoxDet == 0))

# Do LocalPlacement and AvgNightsSinceBaited influence detections? 
# If yes, they will be included as fixed effects in the global model.

FoxModel1<-list()

FoxModel1[[1]] <- glmer(cbind(FoxDet,FoxNotDet) ~ 
                            1 + 
                            (1|Cam), 
                          family = binomial, 
                          data = FoxRR, 
                          control = glmerControl(optimizer="bobyqa"))
simulateResiduals(FoxModel1[[1]], plot = TRUE)
testDispersion(FoxModel1[[1]]) # slightly overdispersed

# what about fitting a model with a zero inflation term in glmmTMB
ziFoxModel <- glmmTMB(cbind(FoxDet,FoxNotDet) ~ 
                          1 + 
                          (1|Cam), 
                        family = binomial, 
                        zi = ~1,
                        data = FoxRR)
simulateResiduals(ziFoxModel, plot = TRUE) # looks nicer
testDispersion(ziFoxModel)


FoxModel1[[2]] <- glmer(cbind(FoxDet,FoxNotDet) ~ 
                            LocalPlacement + 
                            (1|Cam), 
                          family = binomial, 
                          data = FoxRR, 
                          control = glmerControl
                          (optimizer = "bobyqa")); summary(FoxModel1[[2]])

FoxModel1[[3]] <- glmer(cbind(FoxDet,FoxNotDet) ~ 
                            scale(AvgNightsSinceBaited) + 
                            (1|Cam), 
                          family = binomial,
                          data = FoxRR, 
                          control = glmerControl
                          (optimizer = "bobyqa")); summary(FoxModel1[[3]])

FoxModel1[[4]] <- glmer(cbind(FoxDet,FoxNotDet) ~ 
                            LocalPlacement + 
                            scale(AvgNightsSinceBaited) + 
                            (1|Cam), 
                          family = binomial, 
                          data = FoxRR, 
                          control = glmerControl
                          (optimizer = "bobyqa")); summary(FoxModel1[[4]])

FoxModel1names <- paste(c("1. Intercept",
                            "2. Local placement",
                            "3. Nights since baited", 
                            "4. Local placement + Nights since baited"))

(FoxTable1 <- aictab(cand.set = FoxModel1, 
                      modnames = FoxModel1names, 
                      sort = T))

FoxTable1

# How much better is the top model vs the 2nd top model?
evidence(FoxTable1, 
         model.high = "2. Local placement",
         model.low = "4. Local placement + Nights since baited")

# Yes, LocalPlacement is important. This model is x2.8 better than the next model. 
# Local placment will be included as a fixed effect in the global model dredge.

# Fox global model
FoxGlobalTD<- glmer(cbind(FoxDet,FoxNotDet) ~ 
                      LocalPlacement +
                      HabitatType * BeforeAfterFire +
                      scale(PreyActivitySmallMammal) * BeforeAfterFire +
                      scale(PreyActivityMediumMammal) * BeforeAfterFire +
                      scale(PreyActivitySmallBird) * BeforeAfterFire +
                      scale(PreyActivityLargeMammal) * BeforeAfterFire + 
                      scale(DistanceNearestFarm) * BeforeAfterFire +
                      scale(DistanceNearestTown) * BeforeAfterFire + 
                      scale(DistanceNearestRoad) * BeforeAfterFire + 
                      scale(Elevation) * BeforeAfterFire + 
                      #scale(TPI1000) * BeforeAfterFire +
                      scale(NDVI100) * BeforeAfterFire +
                      #BurnTreatment * BeforeAfterFire +
                      scale(PercentAreaBurnt) * BeforeAfterFire +
                      #scale(DistanceBurnEdge) * BeforeAfterFire +
                      scale(YrsSincePreviousFire) * BeforeAfterFire +
                      scale(NumberOfFires100Yrs) * BeforeAfterFire +
                      (1|Cam) +
                      (1|Session),
                    family = binomial, 
                    data = FoxRR,  
                    control = glmerControl(optimizer = "bobyqa")) 

# Perhaps I should combine the prey activity predictors to reduce the number of predictors in the model?

summary(FoxGlobalTD)
simulateResiduals(FoxGlobalTD, plot = T)
testDispersion(FoxGlobalTD) # slightly overdispersed

plot(predictorEffects(FoxGlobalTD))

# Check for singularity
# The definition of singularity is that some of the 
# constrained parameters of the random effects theta 
# parameters are on the boundary (equal to zero, or 
# very very close to zero, say <10-6)

tt <- getME(FoxGlobalTD,"theta")
ll <- getME(FoxGlobalTD,"lower")
min(tt[ll==0])

# Dredge the global model. This does model selection by systematically fitting 
# and comparing different model combinations, considering all possible 
# predictor variable subsets and their interactions.

options(na.action = "na.fail") 
FoxDredge2 <- dredge(FoxGlobalTD, 
                     m.lim = c(0, 4))

FoxDredge2
subset(FoxDredge2,delta <= 2,recalc.weights = FALSE)

# The summed Akaike weights for each predictor.
# Rule of thumb, once you get below ~0.7, the 
# effect of the predictor isn't very strong.
sw(FoxDredge2)


# Summary of 'best' model from dredge
summary(get.models(FoxDredge2, 1)[[1]])
BestFoxModel <- (get.models(FoxDredge2, 1)[[1]])
confint(get.models(FoxDredge2, 1)[[1]])
plot(predictorEffects(get.models(FoxDredge2, 2)[[1]]))
simulateResiduals(get.models(FoxDredge2, 1)[[1]], plot = T)

# Rank the top models based on AICc
ranked_models <- model.sel(FoxDredge2)

# Plot the model selection table
plot(ranked_models)

summary(BestFoxModel)

confint(BestFoxModel)

# ____________________________________________________________________________________
