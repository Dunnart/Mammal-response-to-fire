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
library(snow)
library(glmmTMB)
library(parallel)   

# ____________________________________________________________________________________

# Read in the data and scale the numerical predictor variables
cam_df = read.csv("SiteSpecificCovariates v3.csv", stringsAsFactors = T)


vars_to_scale = c("NumberOfFires100Yrs", 
                  "YrsSincePreviousFire", 
                  "DistanceNearestTown",
                  "DistanceNearestRoad",
                  "DistanceNearestFarm",
                  "TPI100",
                  "PercentAreaBurnt",
                  "NDVI100",
                  "PreyActivityLargeMammal",
                  "PreyActivityBird",
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
                            "scaled_TPI100",
                            "scaled_NDVI100",
                            "scaled_PreyActivityLargeMammal",
                            "scaled_PreyActivityBird",
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
simulateResiduals(bbFoxModel, plot = TRUE) # looks much nicer
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
which(colnames(cam_df2) == "Treatment_BA")

cam_df2_long = reshape2::melt(cam_df2, id.vars = c(1:27,32, 37, 38, 41, 42, 44, 6))
head(cam_df2_long)
str(cam_df2_long)
levels(cam_df2_long$variable)

ggplot(subset(cam_df2_long, variable %in% 
                c("NumberOfFires100Yrs",
                  "YrsSincePreviousFire",
                  "DistanceNearestTown",
                  "DistanceNearestRoad",
                  "DistanceNearestFarm",
                  "NDVI100",
                  "PercentAreaBurnt",
                  "PreyActivityLargeMammal",
                  "PreyActivitySmallMammal",
                  "PreyActivityBird")),
  aes(x = value, y = FoxDet/FoxNotDet, group = BeforeAfterFire)) +
  geom_point(aes(colour=BeforeAfterFire),shape = 21) +
  geom_smooth(aes(colour = BeforeAfterFire)) +
  facet_wrap(~variable, scales = "free") +
  theme_cowplot()

## Fox global model
FoxGlobal = glmmTMB(cbind(FoxDet,FoxNotDet) ~ 
                      LocalPlacement +
                      HabitatType +
                      scaled_PreyActivitySmallMammal * Treatment_BA +
                      scaled_PreyActivityBird * Treatment_BA +
                      scaled_PreyActivityLargeMammal * Treatment_BA + 
                      scaled_DistanceNearestFarm * Treatment_BA +
                      scaled_DistanceNearestTown * Treatment_BA + 
                      scaled_DistanceNearestRoad * Treatment_BA + 
                      scaled_TPI100 * Treatment_BA + 
                      scaled_NDVI100 * Treatment_BA +
                      scaled_PercentAreaBurnt * Treatment_BA +
                      scaled_YrsSincePreviousFire +
                      scaled_NumberOfFires100Yrs * Treatment_BA +
                      (1|Cam) +
                      (1|Session),
                    family = betabinomial, # compare the residual tests for binomial vs betabinomial
                    data = cam_df2) 

# Proceed without fitting the interaction for HabitatType and YrsSincePreviousFire.
# Just write in the methods that the data wasn't suitable for testing that interaction

summary(FoxGlobal)
simulateResiduals(FoxGlobal, plot = T) 
# don't worry too much if this huge global model looks a bit funny. we're not using it for inference
testDispersion(FoxGlobal) 


## Dredge the global fox model. 
# The dredge function does model selection by systematically fitting 
# and comparing different model combinations, considering all possible 
# predictor variable subsets and their interactions.
# First, run code for parallel processing, to designate more CPU cores to the dredge task, to make it faster

options(na.action = "na.fail") 

# Detect number of cores and create cluster (leave one out to not overwhelm pc)
nCores <- detectCores() - 1
my_clust <- makeCluster(nCores, type = "SOCK")

# Export all objects to be used to all the cores in the cluster
clusterExport(my_clust, list("FoxGlobal","cam_df2","cor.matrix"))

# Load packages to be used
clusterEvalQ(my_clust,library(MuMIn,logical.return =T))
clusterEvalQ(my_clust,library(glmmTMB,logical.return =T))

# Run tic() and toc() either side of a function to find out how long it takes. must run at the same time
tic()
FoxDredge = MuMIn::dredge(FoxGlobal,
                          m.lim = c(0, 6),
                          fixed = c("cond(LocalPlacement)"),
                          subset = cor.matrix,
                          trace = 2,
                          cluster = my_clust)

toc()

# Took 12 min to run

FoxDredge # Doesn't give me the delta AIC value for each model
subset(FoxDredge,delta <= 2,recalc.weights = FALSE)  # Can't get this bit of code to work (prob because there's no delta AIC column)
subset(FoxDredge, is.na(delta),recalc.weights = FALSE) # find models that didn't converge

# Top model
summary(get.models(FoxDredge, 1)[[1]])
plot(predictorEffects(get.models(FoxDredge, 1)[[1]]))
simulateResiduals(get.models(FoxDredge, 1)[[1]], plot = T)

# Other top models
summary(get.models(FoxDredge, 2)[[1]])
summary(get.models(FoxDredge, 3)[[1]])
summary(get.models(FoxDredge, 4)[[1]])
summary(get.models(FoxDredge, 5)[[1]])
summary(get.models(FoxDredge, 6)[[1]])
summary(get.models(FoxDredge, 7)[[1]])
summary(get.models(FoxDredge, 8)[[1]])
summary(get.models(FoxDredge, 9)[[1]])


# ____________________________________________________________________________________



#### Feral cat ####

# Check the distribution of the data
hist(cam_df2$CatDet, 20) # Lots of zeroes
nrow(subset(cam_df2, CatDet == 0))

# Do LocalPlacement and AvgNightsSinceBaited influence detections? 
# If yes, they will be included as fixed effects in the global model.

CatModel1<-list()

CatModel1[[1]] <- glmmTMB(cbind(CatDet,CatNotDet) ~ 
                            1 + 
                            (1|Cam) + (1|Session), 
                          family = binomial, 
                          data = cam_df2)
simulateResiduals(CatModel1[[1]], plot = TRUE)
testDispersion(CatModel1[[1]]) # Looks pretty good

CatModel1[[1]]

CatModel1[[2]] <- glmmTMB(cbind(CatDet,CatNotDet) ~ 
                            LocalPlacement + 
                            (1|Cam)  + (1|Session), 
                          family = binomial, 
                          data = cam_df2); summary(CatModel1[[2]])

CatModel1[[3]] <- glmmTMB(cbind(CatDet,CatNotDet) ~ 
                            scaled_AvgNightsSinceBaited + 
                            (1|Cam)  + (1|Session), 
                          family = binomial,
                          data = cam_df2); summary(CatModel1[[3]])

CatModel1[[4]] <- glmmTMB(cbind(CatDet,CatNotDet) ~ 
                            LocalPlacement + 
                            scaled_AvgNightsSinceBaited + 
                            (1|Cam)  + (1|Session), 
                          family = binomial, 
                          data = cam_df2); summary(CatModel1[[4]])

CatModel1names <- paste(c("1. Intercept",
                          "2. Local placement",
                          "3. Nights since baited", 
                          "4. Local placement + Nights since baited"))

(CatTable1 <- aictab(cand.set = CatModel1, 
                     modnames = CatModel1names, 
                     sort = T))

# How much better is the top model vs the 2nd top model?
evidence(CatTable1, 
         model.high = "2. Local placement",
         model.low = "4. Local placement + Nights since baited")

# Yes, LocalPlacement is important. This model is x2.8 better than the next model. 
# Local placment will be included as a fixed effect in the global model dredge.


## Cat make some exploratory plots
str(cam_df2)
which(colnames(cam_df2) == "CatNotDet")
which(colnames(cam_df2) == "HabitatType")
which(colnames(cam_df2) == "LocalPlacement")
which(colnames(cam_df2) == "LandscapePlacement")
which(colnames(cam_df2) == "BurnTreatment")
which(colnames(cam_df2) == "BurnSession")
which(colnames(cam_df2) == "FireSeverity")
which(colnames(cam_df2) == "Treatment_BA")

cam_df2_long = reshape2::melt(cam_df2, id.vars = c(1:27,32, 37, 38, 41, 42, 44, 6))
head(cam_df2_long)
str(cam_df2_long)
levels(cam_df2_long$variable)

ggplot(subset(cam_df2_long, variable %in% 
                c("NumberOfFires100Yrs",
                  "YrsSincePreviousFire",
                  "DistanceNearestTown",
                  "DistanceNearestRoad",
                  "DistanceNearestFarm",
                  "NDVI100",
                  "PercentAreaBurnt",
                  "PreyActivityLargeMammal",
                  "PreyActivitySmallMammal",
                  "PreyActivityBird")),
       aes(x = value, y = CatDet/CatNotDet, group = BeforeAfterFire)) +
  geom_point(aes(colour=BeforeAfterFire),shape = 21) +
  geom_smooth(aes(colour = BeforeAfterFire)) +
  facet_wrap(~variable, scales = "free") +
  theme_cowplot()

## Cat global model
CatGlobal = glmmTMB(cbind(CatDet,CatNotDet) ~ 
                      LocalPlacement +
                      HabitatType +
                      scaled_PreyActivitySmallMammal * Treatment_BA +
                      scaled_PreyActivityBird * Treatment_BA +
                      scaled_PreyActivityLargeMammal * Treatment_BA + 
                      scaled_DistanceNearestFarm * Treatment_BA +
                      scaled_DistanceNearestTown * Treatment_BA + 
                      scaled_DistanceNearestRoad * Treatment_BA + 
                      scaled_TPI100 * Treatment_BA + 
                      scaled_NDVI100 * Treatment_BA +
                      scaled_PercentAreaBurnt * Treatment_BA +
                      scaled_YrsSincePreviousFire +
                      scaled_NumberOfFires100Yrs * Treatment_BA +
                      (1|Cam) +
                      (1|Session),
                    family = binomial, 
                    data = cam_df2) 


summary(CatGlobal)
simulateResiduals(CatGlobal, plot = T) 
# don't worry too much if this huge global model looks a bit funny. we're not using it for inference
testDispersion(CatGlobal) 

## Dredge the global cat model. 

# The dredge function does model selection by systematically fitting 
# and comparing different model combinations, considering all possible 
# predictor variable subsets and their interactions.
# First, run code for parallel processing, to designate more CPU cores to the dredge task, to make it faster

options(na.action = "na.fail") 

# Detect number of cores and create cluster (leave one out to not overwhelm pc)
nCores <- detectCores() - 1
my_clust <- makeCluster(nCores, type = "SOCK")

# Export all objects to be used to all the cores in the cluster
clusterExport(my_clust, list("CatGlobal","cam_df2","cor.matrix"))

# Load packages to be used
clusterEvalQ(my_clust,library(MuMIn,logical.return =T))
clusterEvalQ(my_clust,library(glmmTMB,logical.return =T))

# Run tic() and toc() either side of a function to find out how long it takes. must run at the same time
tic()
CatDredge = MuMIn::dredge(CatGlobal,
                          m.lim = c(0, 6),
                          fixed = c("cond(LocalPlacement)"),
                          subset = cor.matrix,
                          trace = 2,
                          cluster = my_clust)

toc()

# Took 3 min to run

CatDredge # Doesn't give me the delta AIC value for each model
subset(CatDredge,delta <= 2,recalc.weights = FALSE) # Can't get this bit of code to work (prob because there's no delta AIC column)
subset(CatDredge, is.na(delta),recalc.weights = FALSE) # Find models that didn't converge
summary(get.models(CatDredge, 1)[[1]])
plot(predictorEffects(get.models(CatDredge, 1)[[1]]))
simulateResiduals(get.models(CatDredge, 1)[[1]], plot = T)


summary(get.models(CatDredge, 2)[[1]])
summary(get.models(CatDredge, 3)[[1]])
summary(get.models(CatDredge, 4)[[1]])
summary(get.models(CatDredge, 5)[[1]])
summary(get.models(CatDredge, 6)[[1]])
summary(get.models(CatDredge, 7)[[1]])
summary(get.models(CatDredge, 8)[[1]])
summary(get.models(CatDredge, 9)[[1]])


# ____________________________________________________________________________________



#### Swamp wallaby ####

# Check the distribution of the data
hist(cam_df2$SWDet, 20) # Not that many zeros
nrow(subset(cam_df2, SWDet == 0))

# Do LocalPlacement and AvgNightsSinceBaited influence detections? 
# If yes, they will be included as fixed effects in the global model.

SWModel1<-list()

SWModel1[[1]] <- glmmTMB(cbind(SWDet,SWNotDet) ~ 
                            1 + 
                            (1|Cam) + (1|Session), 
                          family = binomial, 
                          data = cam_df2)
simulateResiduals(SWModel1[[1]], plot = TRUE)
testDispersion(SWModel1[[1]]) # Looks pretty good

SWModel1[[1]]

SWModel1[[2]] <- glmmTMB(cbind(SWDet,SWNotDet) ~ 
                            LocalPlacement + 
                            (1|Cam)  + (1|Session), 
                          family = binomial, 
                          data = cam_df2); summary(SWModel1[[2]])

SWModel1[[3]] <- glmmTMB(cbind(SWDet,SWNotDet) ~ 
                            scaled_AvgNightsSinceBaited + 
                            (1|Cam)  + (1|Session), 
                          family = binomial,
                          data = cam_df2); summary(SWModel1[[3]])

SWModel1[[4]] <- glmmTMB(cbind(SWDet,SWNotDet) ~ 
                            LocalPlacement + 
                            scaled_AvgNightsSinceBaited + 
                            (1|Cam)  + (1|Session), 
                          family = binomial, 
                          data = cam_df2); summary(SWModel1[[4]])

SWModel1names <- paste(c("1. Intercept",
                          "2. Local placement",
                          "3. Nights since baited", 
                          "4. Local placement + Nights since baited"))

(SWTable1 <- aictab(cand.set = SWModel1, 
                     modnames = SWModel1names, 
                     sort = T))

# How much better is the top model vs the 2nd top model?
evidence(SWTable1, 
         model.high = "2. Local placement",
         model.low = "4. Local placement + Nights since baited")

# Yes, LocalPlacement is important. This model is x2.8 better than the next model. 
# Local placment will be included as a fixed effect in the global model dredge.


## SW make some exploratory plots
str(cam_df2)
which(colnames(cam_df2) == "SWNotDet")
which(colnames(cam_df2) == "HabitatType")
which(colnames(cam_df2) == "LocalPlacement")
which(colnames(cam_df2) == "LandscapePlacement")
which(colnames(cam_df2) == "BurnTreatment")
which(colnames(cam_df2) == "BurnSession")
which(colnames(cam_df2) == "FireSeverity")
which(colnames(cam_df2) == "Treatment_BA")

cam_df2_long = reshape2::melt(cam_df2, id.vars = c(1:27,32, 37, 38, 41, 42, 44, 6))
head(cam_df2_long)
str(cam_df2_long)
levels(cam_df2_long$variable)

ggplot(subset(cam_df2_long, variable %in% 
                c("NumberOfFires100Yrs",
                  "YrsSincePreviousFire",
                  "DistanceNearestTown",
                  "DistanceNearestRoad",
                  "DistanceNearestFarm",
                  "NDVI100",
                  "PercentAreaBurnt",
                  "PreyActivityLargeMammal",
                  "PreyActivitySmallMammal",
                  "PreyActivityBird")),
       aes(x = value, y = SWDet/SWNotDet, group = BeforeAfterFire)) +
  geom_point(aes(colour=BeforeAfterFire),shape = 21) +
  geom_smooth(aes(colour = BeforeAfterFire)) +
  facet_wrap(~variable, scales = "free") +
  theme_cowplot()

## SW global model
SWGlobal = glmmTMB(cbind(SWDet,SWNotDet) ~ 
                      LocalPlacement +
                      HabitatType +
                      scaled_DistanceNearestFarm * Treatment_BA +
                      scaled_DistanceNearestTown * Treatment_BA + 
                      scaled_DistanceNearestRoad * Treatment_BA + 
                      scaled_TPI100 * Treatment_BA + 
                      scaled_NDVI100 * Treatment_BA +
                      scaled_PercentAreaBurnt * Treatment_BA +
                      scaled_YrsSincePreviousFire +
                      scaled_NumberOfFires100Yrs * Treatment_BA +
                      (1|Cam) +
                      (1|Session),
                    family = binomial, 
                    data = cam_df2) 


summary(SWGlobal)
simulateResiduals(SWGlobal, plot = T) 
# don't worry too much if this huge global model looks a bit funny. we're not using it for inference
testDispersion(SWGlobal) 

## Dredge the global swamp wallaby model. 
# The dredge function does model selection by systematically fitting 
# and comparing different model combinations, considering all possible 
# predictor variable subsets and their interactions.
# First, run code for parallel processing, to designate more CPU cores to the dredge task, to make it faster

options(na.action = "na.fail") 

# Detect number of cores and create cluster (leave one out to not overwhelm pc)
nCores <- detectCores() - 1
my_clust <- makeCluster(nCores, type = "SOCK")

# Export all objects to be used to all the cores in the cluster
clusterExport(my_clust, list("SWGlobal","cam_df2","cor.matrix"))

# Load packages to be used
clusterEvalQ(my_clust,library(MuMIn,logical.return =T))
clusterEvalQ(my_clust,library(glmmTMB,logical.return =T))

# Run tic() and toc() either side of a function to find out how long it takes. must run at the same time
tic()
SWDredge = MuMIn::dredge(SWGlobal,
                          m.lim = c(0, 6),
                          subset = cor.matrix,
                          trace = 2,
                          cluster = my_clust)

toc()

# Took 2.5 min to run

SWDredge
subset(SWDredge,delta <= 2,recalc.weights = FALSE)
summary(get.models(SWDredge, 1)[[1]])
plot(predictorEffects(get.models(SWDredge, 1)[[1]]))
simulateResiduals(get.models(SWDredge, 1)[[1]], plot = T)

# Only one good model for swamp wallabies


# ____________________________________________________________________________________



#### Eastern grey kangaroo ####

# Check the distribution of the data
hist(cam_df2$RooDet, 20) # Lots of zeroes
nrow(subset(cam_df2, RooDet == 0))

# Do LocalPlacement and AvgNightsSinceBaited influence detections? 
# If yes, they will be included as fixed effects in the global model.

RooModel1<-list()

RooModel1[[1]] <- glmmTMB(cbind(RooDet,RooNotDet) ~ 
                            1 + 
                            (1|Cam) + (1|Session), 
                          family = binomial, 
                          data = cam_df2)

simulateResiduals(RooModel1[[1]], plot = TRUE)
testDispersion(RooModel1[[1]]) # Looks pretty good

RooModel1[[1]]

RooModel1[[2]] <- glmmTMB(cbind(RooDet,RooNotDet) ~ 
                            LocalPlacement + 
                            (1|Cam)  + (1|Session), 
                          family = binomial, 
                          data = cam_df2); summary(RooModel1[[2]])

RooModel1[[3]] <- glmmTMB(cbind(RooDet,RooNotDet) ~ 
                            scaled_AvgNightsSinceBaited + 
                            (1|Cam)  + (1|Session), 
                          family = binomial,
                          data = cam_df2); summary(RooModel1[[3]])

RooModel1[[4]] <- glmmTMB(cbind(RooDet,RooNotDet) ~ 
                            LocalPlacement + 
                            scaled_AvgNightsSinceBaited + 
                            (1|Cam)  + (1|Session), 
                          family = binomial, 
                          data = cam_df2); summary(RooModel1[[4]])

RooModel1names <- paste(c("1. Intercept",
                          "2. Local placement",
                          "3. Nights since baited", 
                          "4. Local placement + Nights since baited"))

(RooTable1 <- aictab(cand.set = RooModel1, 
                     modnames = RooModel1names, 
                     sort = T))

# How much better is the top model vs the 2nd top model?
evidence(RooTable1, 
         model.high = "2. Local placement",
         model.low = "4. Local placement + Nights since baited")

# Yes, both LocalPlacement and AvgNightsSinceBaited are important 
# and will be included as a fixed effect in the global model dredge.


## Make some exploratory plots
str(cam_df2)
which(colnames(cam_df2) == "RooNotDet")
which(colnames(cam_df2) == "HabitatType")
which(colnames(cam_df2) == "LocalPlacement")
which(colnames(cam_df2) == "LandscapePlacement")
which(colnames(cam_df2) == "BurnTreatment")
which(colnames(cam_df2) == "BurnSession")
which(colnames(cam_df2) == "FireSeverity")
which(colnames(cam_df2) == "Treatment_BA")

cam_df2_long = reshape2::melt(cam_df2, id.vars = c(1:27,32, 37, 38, 41, 42, 44, 6))
head(cam_df2_long)
str(cam_df2_long)
levels(cam_df2_long$variable)

ggplot(subset(cam_df2_long, variable %in% 
                c("NumberOfFires100Yrs",
                  "YrsSincePreviousFire",
                  "DistanceNearestTown",
                  "DistanceNearestRoad",
                  "DistanceNearestFarm",
                  "NDVI100",
                  "PercentAreaBurnt",
                  "PreyActivityLargeMammal",
                  "PreyActivitySmallMammal",
                  "PreyActivityBird")),
       aes(x = value, y = RooDet/RooNotDet, group = BeforeAfterFire)) +
  geom_point(aes(colour=BeforeAfterFire),shape = 21) +
  geom_smooth(aes(colour = BeforeAfterFire)) +
  facet_wrap(~variable, scales = "free") +
  theme_cowplot()

## Roo global model
RooGlobal = glmmTMB(cbind(RooDet,RooNotDet) ~ 
                      LocalPlacement +
                      scaled_AvgNightsSinceBaited +
                      HabitatType +
                      scaled_DistanceNearestFarm * Treatment_BA +
                      scaled_DistanceNearestTown * Treatment_BA + 
                      scaled_DistanceNearestRoad * Treatment_BA + 
                      scaled_TPI100 * Treatment_BA + 
                      scaled_NDVI100 * Treatment_BA +
                      scaled_PercentAreaBurnt * Treatment_BA +
                      scaled_YrsSincePreviousFire +
                      scaled_NumberOfFires100Yrs * Treatment_BA +
                      (1|Cam) +
                      (1|Session),
                    family = binomial, 
                    data = cam_df2) 


summary(RooGlobal)
simulateResiduals(RooGlobal, plot = T) 
# don't worry too much if this huge global model looks a bit funny. we're not using it for inference
testDispersion(RooGlobal) 

## Dredge the global Roo model. 

# The dredge function does model selection by systematically fitting 
# and comparing different model combinations, considering all possible 
# predictor variable subsets and their interactions.
# First, run code for parallel processing, to designate more CPU cores to the dredge task, to make it faster

options(na.action = "na.fail") 

# Detect number of cores and create cluster (leave one out to not overwhelm pc)
nCores <- detectCores() - 1
my_clust <- makeCluster(nCores, type = "SOCK")

# Export all objects to be used to all the cores in the cluster
clusterExport(my_clust, list("RooGlobal","cam_df2","cor.matrix"))

# Load packages to be used
clusterEvalQ(my_clust,library(MuMIn,logical.return =T))
clusterEvalQ(my_clust,library(glmmTMB,logical.return =T))

# Run tic() and toc() either side of a function to find out how long it takes. must run at the same time
tic()
RooDredge = MuMIn::dredge(RooGlobal,
                          m.lim = c(0, 6),
                          fixed = c("cond(LocalPlacement, scaled_AvgNightsSinceBaited)"),
                          subset = cor.matrix,
                          trace = 2,
                          cluster = my_clust)

toc()

# Took 4 min to run

RooDredge
subset(RooDredge,delta <= 2,recalc.weights = FALSE)
subset(RooDredge, is.na(delta),recalc.weights = FALSE) # Find models that didn't converge
summary(get.models(RooDredge, 1)[[1]])
plot(predictorEffects(get.models(RooDredge, 1)[[1]]))
simulateResiduals(get.models(RooDredge, 1)[[1]], plot = T)


summary(get.models(RooDredge, 2)[[1]])


# Two models with delta AICs <2


# ____________________________________________________________________________________



#### Small mammals ####

# Check the distribution of the data
hist(cam_df2$SmMamDet, 20) # LOTS of zeroes
nrow(subset(cam_df2, SmMamDet == 0))

# Do LocalPlacement and AvgNightsSinceBaited influence detections? 
# If yes, they will be included as fixed effects in the global model.

SmMamModel1<-list()

SmMamModel1[[1]] <- glmmTMB(cbind(SmMamDet,SmMamNotDet) ~ 
                            1 + 
                            (1|Cam) + (1|Session), 
                          family = binomial, 
                          data = cam_df2)
simulateResiduals(SmMamModel1[[1]], plot = TRUE)
testDispersion(SmMamModel1[[1]]) 

SmMamModel1[[1]]

SmMamModel1[[2]] <- glmmTMB(cbind(SmMamDet,SmMamNotDet) ~ 
                            LocalPlacement + 
                            (1|Cam)  + (1|Session), 
                          family = binomial, 
                          data = cam_df2); summary(SmMamModel1[[2]])

SmMamModel1[[3]] <- glmmTMB(cbind(SmMamDet,SmMamNotDet) ~ 
                            scaled_AvgNightsSinceBaited + 
                            (1|Cam)  + (1|Session), 
                          family = binomial,
                          data = cam_df2); summary(SmMamModel1[[3]])

SmMamModel1[[4]] <- glmmTMB(cbind(SmMamDet,SmMamNotDet) ~ 
                            LocalPlacement + 
                            scaled_AvgNightsSinceBaited + 
                            (1|Cam)  + (1|Session), 
                          family = binomial, 
                          data = cam_df2); summary(SmMamModel1[[4]])

SmMamModel1names <- paste(c("1. Intercept",
                          "2. Local placement",
                          "3. Nights since baited", 
                          "4. Local placement + Nights since baited"))

(SmMamTable1 <- aictab(cand.set = SmMamModel1, 
                     modnames = SmMamModel1names, 
                     sort = T))

# How much better is the top model vs the 2nd top model?
evidence(SmMamTable1, 
         model.high = "2. Local placement",
         model.low = "4. Local placement + Nights since baited")

# Neither are <0.05, scaled_AvgNightsSinceBaited is close. 


## Make some exploratory plots
str(cam_df2)
which(colnames(cam_df2) == "SmMamNotDet")
which(colnames(cam_df2) == "HabitatType")
which(colnames(cam_df2) == "LocalPlacement")
which(colnames(cam_df2) == "LandscapePlacement")
which(colnames(cam_df2) == "BurnTreatment")
which(colnames(cam_df2) == "BurnSession")
which(colnames(cam_df2) == "FireSeverity")
which(colnames(cam_df2) == "Treatment_BA")

cam_df2_long = reshape2::melt(cam_df2, id.vars = c(1:27,32, 37, 38, 41, 42, 44, 6))
head(cam_df2_long)
str(cam_df2_long)
levels(cam_df2_long$variable)

ggplot(subset(cam_df2_long, variable %in% 
                c("NumberOfFires100Yrs",
                  "YrsSincePreviousFire",
                  "DistanceNearestTown",
                  "DistanceNearestRoad",
                  "DistanceNearestFarm",
                  "NDVI100",
                  "PercentAreaBurnt",
                  "PreyActivityLargeMammal",
                  "PreyActivitySmallMammal",
                  "PreyActivityBird")),
       aes(x = value, y = SmMamDet/SmMamNotDet, group = BeforeAfterFire)) +
  geom_point(aes(colour=BeforeAfterFire),shape = 21) +
  geom_smooth(aes(colour = BeforeAfterFire)) +
  facet_wrap(~variable, scales = "free") +
  theme_cowplot()

## SmMam global model
SmMamGlobal = glmmTMB(cbind(SmMamDet,SmMamNotDet) ~ 
                      LocalPlacement +
                      HabitatType +
                      scaled_DistanceNearestFarm * Treatment_BA +
                      scaled_DistanceNearestTown * Treatment_BA + 
                      scaled_DistanceNearestRoad * Treatment_BA + 
                      scaled_TPI100 * Treatment_BA + 
                      scaled_NDVI100 * Treatment_BA +
                      scaled_PercentAreaBurnt * Treatment_BA +
                      scaled_YrsSincePreviousFire +
                      scaled_NumberOfFires100Yrs * Treatment_BA +
                      (1|Cam) +
                      (1|Session),
                    family = binomial, 
                    data = cam_df2) 


summary(SmMamGlobal)
simulateResiduals(SmMamGlobal, plot = T) 
# don't worry too much if this huge global model looks a bit funny. we're not using it for inference
testDispersion(SmMamGlobal) 

## Dredge the global small mammal model. 
# The dredge function does model selection by systematically fitting 
# and comparing different model combinations, considering all possible 
# predictor variable subsets and their interactions.
# First, run code for parallel processing, to designate more CPU cores to the dredge task, to make it faster

options(na.action = "na.fail") 

# Detect number of cores and create cluster (leave one out to not overwhelm pc)
nCores <- detectCores() - 1
my_clust <- makeCluster(nCores, type = "SOCK")

# Export all objects to be used to all the cores in the cluster
clusterExport(my_clust, list("SmMamGlobal","cam_df2","cor.matrix"))

# Load packages to be used
clusterEvalQ(my_clust,library(MuMIn,logical.return =T))
clusterEvalQ(my_clust,library(glmmTMB,logical.return =T))

# Run tic() and toc() either side of a function to find out how long it takes. must run at the same time
tic()
SmMamDredge = MuMIn::dredge(SmMamGlobal,
                          m.lim = c(0, 6),
                          subset = cor.matrix,
                          trace = 2,
                          cluster = my_clust)

toc()

# Took 2 min to run

SmMamDredge
subset(SmMamDredge,delta <= 2,recalc.weights = FALSE) 
subset(SmMamDredge, is.na(delta),recalc.weights = FALSE) # Find models that didn't converge
summary(get.models(SmMamDredge, 1)[[1]])
plot(predictorEffects(get.models(SmMamDredge, 1)[[1]]))
simulateResiduals(get.models(SmMamDredge, 1)[[1]], plot = T)


summary(get.models(SmMamDredge, 2)[[1]])
summary(get.models(SmMamDredge, 3)[[1]])
summary(get.models(SmMamDredge, 4)[[1]])
summary(get.models(SmMamDredge, 5)[[1]])
summary(get.models(SmMamDredge, 6)[[1]])

# Six models with delta AIC <2


# ____________________________________________________________________________________



#### Medium marsupials ####

# Check the distribution of the data
hist(cam_df2$MedMarDet, 20) # LOTS of zeroes
nrow(subset(cam_df2, MedMarDet == 0))

# Do LocalPlacement and AvgNightsSinceBaited influence detections? 
# If yes, they will be included as fixed effects in the global model.

MedMarModel1<-list()

MedMarModel1[[1]] <- glmmTMB(cbind(MedMarDet,MedMarNotDet) ~ 
                              1 + 
                              (1|Cam) + (1|Session), 
                            family = binomial, 
                            data = cam_df2)

simulateResiduals(MedMarModel1[[1]], plot = TRUE)
testDispersion(MedMarModel1[[1]]) 

MedMarModel1[[1]]

MedMarModel1[[2]] <- glmmTMB(cbind(MedMarDet,MedMarNotDet) ~ 
                              LocalPlacement + 
                              (1|Cam)  + (1|Session), 
                            family = binomial, 
                            data = cam_df2); summary(MedMarModel1[[2]])

MedMarModel1[[3]] <- glmmTMB(cbind(MedMarDet,MedMarNotDet) ~ 
                              scaled_AvgNightsSinceBaited + 
                              (1|Cam)  + (1|Session), 
                            family = binomial,
                            data = cam_df2); summary(MedMarModel1[[3]])

MedMarModel1[[4]] <- glmmTMB(cbind(MedMarDet,MedMarNotDet) ~ 
                              LocalPlacement + 
                              scaled_AvgNightsSinceBaited + 
                              (1|Cam)  + (1|Session), 
                            family = binomial, 
                            data = cam_df2); summary(MedMarModel1[[4]])

MedMarModel1names <- paste(c("1. Intercept",
                            "2. Local placement",
                            "3. Nights since baited", 
                            "4. Local placement + Nights since baited"))

(MedMarTable1 <- aictab(cand.set = MedMarModel1, 
                       modnames = MedMarModel1names, 
                       sort = T))

# How much better is the top model vs the 2nd top model?
evidence(MedMarTable1, 
         model.high = "2. Local placement",
         model.low = "4. Local placement + Nights since baited")

# Neither are <0.05.

## Make some exploratory plots
str(cam_df2)
which(colnames(cam_df2) == "MedMarNotDet")
which(colnames(cam_df2) == "HabitatType")
which(colnames(cam_df2) == "LocalPlacement")
which(colnames(cam_df2) == "LandscapePlacement")
which(colnames(cam_df2) == "BurnTreatment")
which(colnames(cam_df2) == "BurnSession")
which(colnames(cam_df2) == "FireSeverity")
which(colnames(cam_df2) == "Treatment_BA")

cam_df2_long = reshape2::melt(cam_df2, id.vars = c(1:27,32, 37, 38, 41, 42, 44, 6))
head(cam_df2_long)
str(cam_df2_long)
levels(cam_df2_long$variable)

ggplot(subset(cam_df2_long, variable %in% 
                c("NumberOfFires100Yrs",
                  "YrsSincePreviousFire",
                  "DistanceNearestTown",
                  "DistanceNearestRoad",
                  "DistanceNearestFarm",
                  "NDVI100",
                  "PercentAreaBurnt",
                  "PreyActivityLargeMammal",
                  "PreyActivitySmallMammal",
                  "PreyActivityBird")),
       aes(x = value, y = MedMarDet/MedMarNotDet, group = BeforeAfterFire)) +
  geom_point(aes(colour=BeforeAfterFire),shape = 21) +
  geom_smooth(aes(colour = BeforeAfterFire)) +
  facet_wrap(~variable, scales = "free") +
  theme_cowplot()

## Medium marsupial global model
MedMarGlobal = glmmTMB(cbind(MedMarDet,MedMarNotDet) ~ 
                        #LocalPlacement +
                        #HabitatType +
                        #scaled_DistanceNearestFarm * Treatment_BA +
                        scaled_DistanceNearestTown * Treatment_BA + 
                        scaled_DistanceNearestRoad * Treatment_BA + 
                        scaled_TPI100 * Treatment_BA + 
                        scaled_NDVI100 * Treatment_BA +
                        scaled_PercentAreaBurnt * Treatment_BA +
                        #scaled_YrsSincePreviousFire +
                        #scaled_NumberOfFires100Yrs * Treatment_BA +
                        (1|Cam) +
                        (1|Session),
                      family = binomial, 
                      data = cam_df2) 

# Model only runs without convergence issues when considerably simplified
# and binomial distribution (as opposed to betabinomial)

summary(MedMarGlobal)
simulateResiduals(MedMarGlobal, plot = T) 
# don't worry too much if this huge global model looks a bit funny. we're not using it for inference
testDispersion(MedMarGlobal) 

## Dredge the global small mammal model. 
# The dredge function does model selection by systematically fitting 
# and comparing different model combinations, considering all possible 
# predictor variable subsets and their interactions.
# First, run code for parallel processing, to designate more CPU cores to the dredge task, to make it faster

options(na.action = "na.fail") 

# Detect number of cores and create cluster (leave one out to not overwhelm pc)
nCores <- detectCores() - 1
my_clust <- makeCluster(nCores, type = "SOCK")

# Export all objects to be used to all the cores in the cluster
clusterExport(my_clust, list("MedMarGlobal","cam_df2","cor.matrix"))

# Load packages to be used
clusterEvalQ(my_clust,library(MuMIn,logical.return =T))
clusterEvalQ(my_clust,library(glmmTMB,logical.return =T))

# Run tic() and toc() either side of a function to find out how long it takes. must run at the same time
tic()
MedMarDredge = MuMIn::dredge(MedMarGlobal,
                            m.lim = c(0, 6),
                            subset = cor.matrix,
                            trace = 2,
                            cluster = my_clust)

toc()

# Took <1 min to run

MedMarDredge # Doesn't give me the delta AIC value for each model
subset(MedMarDredge,delta <= 2,recalc.weights = FALSE) 
subset(MedMarDredge, is.na(delta),recalc.weights = FALSE) # Find models that didn't converge
summary(get.models(MedMarDredge, 1)[[1]])
plot(predictorEffects(get.models(MedMarDredge, 1)[[1]]))
simulateResiduals(get.models(MedMarDredge, 1)[[1]], plot = T)


summary(get.models(MedMarDredge, 2)[[1]])
summary(get.models(MedMarDredge, 3)[[1]])
summary(get.models(MedMarDredge, 4)[[1]])
summary(get.models(MedMarDredge, 5)[[1]])
summary(get.models(MedMarDredge, 6)[[1]])

# Six models with delta AIC <2

# ____________________________________________________________________________________




