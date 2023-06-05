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
library(reshape2)
library(ggplot2)
library(cowplot)
library(tictoc)

# ____________________________________________________________________________________
# Read in the data and scale the numerical predictor variables
cam_df = read.csv("SiteSpecificCovariates v2.csv", stringsAsFactors = T)


vars_to_scale = c("NumberOfFires100Yrs", 
                  "YrsSincePreviousFire", 
                  "DistanceNearestTown",
                  "DistanceNearestRoad",
                  "DistanceNearestFarm",
                  "Elevation",
                  "TPI100",
                  "TPI1000",
                  "PercentAreaBurnt",
                  "DistanceBurnEdge",
                  "NDVI100",
                  "PreyActivityLargeMammal",
                  "PreyActivityMediumBird",
                  "PreyActivityMediumMammal",
                  "PreyActivitySmallBird",
                  "PreyActivitySmallMammal",
                  "AvgNightsSinceBaited")

# Apply scaling to the selected variables
scaled_vars = scale(cam_df[, vars_to_scale])

# Create new variable names with 'scaled_' prefix
new_var_names = paste0("scaled_", vars_to_scale)

# Combine the scaled variables with the original dataframe
cam_df2 <- cbind(cam_df, setNames(as.data.frame(scaled_vars), new_var_names))
str(cam_df2)


# ____________________________________________________________________________________
#### Collinearity check ####

# Here I create a matrix for subsetting dredge later where the cutoff is r = 0.5
cor.matrix <- abs(cor(cam_df2[, c("scaled_NumberOfFires100Yrs", 
                            "scaled_YrsSincePreviousFire", 
                            "scaled_DistanceNearestTown",
                            "scaled_DistanceNearestRoad",
                            "scaled_DistanceNearestFarm",
                            "scaled_Elevation",
                            "scaled_TPI100",
                            "scaled_TPI1000",
                            "scaled_PercentAreaBurnt",
                            "scaled_DistanceBurnEdge",
                            "scaled_NDVI100",
                            "scaled_PreyActivityLargeMammal",
                            "scaled_PreyActivityMediumBird",
                            "scaled_PreyActivityMediumMammal",
                            "scaled_PreyActivitySmallBird",
                            "scaled_PreyActivitySmallMammal",
                            "scaled_AvgNightsSinceBaited")])) <= .5
cor.matrix[!lower.tri(cor.matrix)] <- NA

# glmmTMB names the model terms like cond(LocalPlacement), so we need to rename this matrix
# Define a function to add "cond()" prefix to column and row names
add_cond_prefix <- function(name) {
  return(paste0("cond(", name, ")"))
}

# Rename the column names
colnames(cor.matrix) <- sapply(colnames(cor.matrix), add_cond_prefix)

# Rename the row names
rownames(cor.matrix) <- sapply(rownames(cor.matrix), add_cond_prefix)

cor.matrix 

# ____________________________________________________________________________________
#### Red fox ####

# Check the distribution of the data
hist(cam_df2$FoxDet, 20) # Lots of zeroes
nrow(subset(cam_df2, FoxDet == 0))

# Do LocalPlacement and AvgNightsSinceBaited influence detections? 
# If yes, they will be included as fixed effects in the global model.

FoxModel1<-list()

FoxModel1[[1]] <- glmmTMB(cbind(FoxDet,FoxNotDet) ~ 
                          1 + 
                          (1|Cam) + (1|Session), 
                        family = binomial, 
                        data = cam_df2)
simulateResiduals(FoxModel1[[1]], plot = TRUE)
testDispersion(FoxModel1[[1]]) # slightly overdispersed

# try the betabinomial distribution instead
bbFoxModel <- glmmTMB(cbind(FoxDet,FoxNotDet) ~ 
                        1 + 
                        (1|Cam)  + (1|Session), 
                      family = betabinomial, 
                      data = cam_df2)
simulateResiduals(ziFoxModel, plot = TRUE) # looks much nicer
testDispersion(bbFoxModel)

FoxModel1[[1]] = bbFoxModel

FoxModel1[[2]] <- glmmTMB(cbind(FoxDet,FoxNotDet) ~ 
                          LocalPlacement + 
                          (1|Cam)  + (1|Session), 
                        family = betabinomial, 
                        data = cam_df2); summary(FoxModel1[[2]])

FoxModel1[[3]] <- glmmTMB(cbind(FoxDet,FoxNotDet) ~ 
                          scaled_AvgNightsSinceBaited + 
                          (1|Cam)  + (1|Session), 
                        family = betabinomial,
                        data = cam_df2); summary(FoxModel1[[3]])

FoxModel1[[4]] <- glmmTMB(cbind(FoxDet,FoxNotDet) ~ 
                          LocalPlacement + 
                          scaled_AvgNightsSinceBaited + 
                          (1|Cam)  + (1|Session), 
                        family = betabinomial, 
                        data = cam_df2); summary(FoxModel1[[4]])

FoxModel1names <- paste(c("1. Intercept",
                          "2. Local placement",
                          "3. Nights since baited", 
                          "4. Local placement + Nights since baited"))

(FoxTable1 <- aictab(cand.set = FoxModel1, 
                     modnames = FoxModel1names, 
                     sort = T))

# How much better is the top model vs the 2nd top model?
evidence(FoxTable1, 
         model.high = "2. Local placement",
         model.low = "4. Local placement + Nights since baited")

# Yes, LocalPlacement is important. This model is x2.8 better than the next model. 
# Local placment will be included as a fixed effect in the global model dredge.


## Fox make some exploratory plots
str(cam_df2)
which(colnames(cam_df2) == "FoxNotDet")
which(colnames(cam_df2) == "HabitatType")
which(colnames(cam_df2) == "LocalPlacement")
which(colnames(cam_df2) == "LandscapePlacement")
which(colnames(cam_df2) == "BurnTreatment")
which(colnames(cam_df2) == "BurnSession")
which(colnames(cam_df2) == "FireSeverity")

cam_df2_long = reshape2::melt(cam_df2, id.vars = c(1:26,31, 36, 37, 40, 41, 43))
head(cam_df2_long)
str(cam_df2_long)
levels(cam_df2_long$variable)

ggplot(subset(cam_df2_long, variable %in% 
                c("NumberOfFires100Yrs",
                  "YrsSincePreviousFire",
                  "DistanceNearestTown",
                  "DistanceNearestRoad",
                  "DistanceNearestFarm",
                  "Elevation",
                  "DistanceBurnEdge",
                  "NDVI100",
                  "PercentAreaBurnt",
                  "PreyActivityLargeMammal",
                  "PreyActivityMediumMammal",
                  "PreyActivitySmallMammal",
                  "PreyActivitySmallBird")),
  aes(x = value, y = FoxDet/FoxNotDet, group = BeforeAfterFire)) +
  geom_point(aes(colour=BeforeAfterFire),shape = 21) +
  geom_smooth(aes(colour = BeforeAfterFire)) +
  facet_wrap(~variable, scales = "free") +
  theme_cowplot()

## Fox global model
FoxGlobal = glmmTMB(cbind(FoxDet,FoxNotDet) ~ 
                      LocalPlacement +
                      HabitatType * BeforeAfterFire +
                      scaled_PreyActivitySmallMammal * BeforeAfterFire +
                      scaled_PreyActivityMediumMammal * BeforeAfterFire +
                      scaled_PreyActivitySmallBird * BeforeAfterFire +
                      scaled_PreyActivityLargeMammal * BeforeAfterFire + 
                      scaled_DistanceNearestFarm * BeforeAfterFire +
                      scaled_DistanceNearestTown * BeforeAfterFire + 
                      scaled_DistanceNearestRoad * BeforeAfterFire + 
                      scaled_Elevation * BeforeAfterFire + 
                      scaled_NDVI100 * BeforeAfterFire +
                      BurnTreatment * BeforeAfterFire +
                      scaled_PercentAreaBurnt * BeforeAfterFire +
                      #scaled_DistanceBurnEdge * BeforeAfterFire +
                      scaled_YrsSincePreviousFire * BeforeAfterFire +
                      scaled_NumberOfFires100Yrs * BeforeAfterFire +
                      (1|Cam) +
                      (1|Session),
                    family = betabinomial, # compare the residual tests for binomial vs betabinomial
                    data = cam_df2) 

# Perhaps I should combine the prey activity predictors to reduce the number of predictors in the model?

summary(FoxGlobal)
simulateResiduals(FoxGlobal, plot = T) 
# don't worry too much if this huge global model looks a bit funny. we're not using it for inference
testDispersion(FoxGlobal) 


# Dredge the global model. This does model selection by systematically fitting 
# and comparing different model combinations, considering all possible 
# predictor variable subsets and their interactions.

options(na.action = "na.fail") 
tic() # run these two either side of a function to find out how long it takes. must run at the same time
FoxDredge <- dredge(FoxGlobal, 
                     m.lim = c(0, 8),
                    fixed = c("cond(LocalPlacement)"),
                     subset = cor.matrix)
toc()

FoxDredge
subset(FoxDredge,delta <= 2,recalc.weights = FALSE)

subset(FoxDredge, is.na(delta),recalc.weights = FALSE) # find models that didn't converge

summary(get.models(FoxDredge, 1)[[1]])
plot(predictorEffects(get.models(FoxDredge, 1)[[1]]))
simulateResiduals(get.models(FoxDredge, 1)[[1]], plot = T)


summary(get.models(FoxDredge, 2)[[1]])
summary(get.models(FoxDredge, 3)[[1]])
summary(get.models(FoxDredge, 4)[[1]])
summary(get.models(FoxDredge, 5)[[1]])
summary(get.models(FoxDredge, 6)[[1]])
summary(get.models(FoxDredge, 7)[[1]])
summary(get.models(FoxDredge, 8)[[1]])
summary(get.models(FoxDredge, 9)[[1]])




# ____________________________________________________________________________________
