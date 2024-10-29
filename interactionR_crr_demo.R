### check the "interactionR_crr.R"
# Cleaning environment 
rm(list = ls()) 
source("C:/Users/NN/Desktop/Code/interactionR_crr.R")
# Loading up libraries and setting directory 
library("tidyverse") 
#library("gmodels") 
#library("broom") 
library("survival") 
#library("prodlim") 
library("cmprsk") 
#library("sandwich") 
#library("tableone") 
#library("ltmle") 
#library("lmtest") 

###### Created data ############
# Dataset simulation
# This follows a DAG with baseline confounding and selection/mediation by a competing event 
# Definition of variables: 
# La: baseline confounder, P(La) = 0.5 (no parents) 
# Lb: predictor of C and Y, P(Lb) = 0.5 (no parents) 
# A: exposure, P(A | La=0) = 0.45 (one parent, La) 
# C: competing event, P(C | A=0, Lb=0) = 0.20 

# Y: outcome, P(Y | A=0, La=0, Lb=0, C=0) = 0.15 
# T: follow-up time 
# Associations between variables in the risk difference scale using rbinom 
# Generating dataset 
set.seed(1) 
n <- 10000 # Number of repetitions 
La <- rbinom(n, 1, 0.5) # Baseline confounder 
Lb <- rbinom(n, 1, 0.5) # Predictor of competing event 
A <- rbinom(n, 1, 0.45 + 0.05 * La) # Exposure (non-invasive ventilation) 
C <- rbinom(n, 1, 0.2 + 0.2 * Lb - 0.1 * A) # Competing event all-cause death 
Y <- ifelse(C==1, NA, rbinom(n, 1, 0.15 + 0.05*Lb - 0.04 * A + 0.05 * La)) # Outcome 
# Time; for simplicity only depends on the competing event, exposure and fixed values 
fup.days <- ifelse(C == 1, 
                   pmin((rgeom(n, plogis(-4.3 - 0.10 * A)) + 1), 730), 
                   pmin((rgeom(n, plogis(-5.0 - 0.18 * A)) + 1), 730)) 
aggregate(fup.days ~ C, FUN = median) # Simple description of time by censoring group 
aggregate(fup.days ~ A, FUN = median) # Simple description of time by treatment group 
max(fup.days) # Maximum follow-up ~ 2 years 
# Creating dataset 
data.comp <- data.frame(A, La, Lb, C, Y, fup.days) # Variables to use below 
data.comp <- mutate(data.comp, id = row_number()) # Adding id variable 
head(data.comp) # Glimpse of dataset 
# Creating matrix for Fine and Gray regression model (see below) 
covariates <- data.frame(data.comp$A, 
                         data.comp$La, 
                         data.comp$Lb) 
head(covariates)
# Data preparation----- 
# Generating variables for analysis 
# Competing where 0 is censored; 1 is re-hospitalization; 2 is the competing event 
data.comp$competing <- ifelse(data.comp$C == 1, 2, 
                              ifelse(data.comp$Y == 1, 1, 
                                     0) 
) 
# Combined outcome where 1 is either rehospitalization or the competing event 
data.comp$combined <- ifelse(data.comp$C == 1,1, 
                             ifelse(data.comp$Y == 1, 1, 
                                    0) 
) 
# Selection as the inverse of competing event 
data.comp$selection <- ifelse(data.comp$C == 1, 0, 
                              1) 
## Simple description---- 
description <- print(CreateTableOne( 
  data = data.comp, 
  vars = c("La", "Lb", "C", "Y", "fup.days"), 
  strata = c("A"), 
  factorVars = c("La","Lb", "C", "Y"), 
  addOverall = TRUE), 
  test = FALSE, 
  smd = TRUE) 
table_one <- CreateTableOne( 
  data = data.comp, 
  vars = c("La", "Lb", "C", "Y", "fup.days"), 
  strata = c("A"), 
  factorVars = c("La","Lb", "C", "Y"), 
  addOverall = TRUE) 
summary(table_one) 


data.comp$`A:La` =  data.comp$A * data.comp$La
# Creating matrix for Fine and Gray regression model (see below) 
covariates_inter <- data.frame(data.comp$A, 
                               data.comp$La,
                               data.comp$Lb,
                               data.comp$`A:La`) 
fg <- crr(ftime = data.comp$fup.days, # Time variable 
          fstatus = data.comp$competing, # Variable including all levels as above 
          cov1 = covariates_inter, # Covariate matrix 
          failcode = 1, # Outcome of interest 
          cencode = 0) # Censored observations 
summary(fg) 
inter_out <- interactionR_crr(fg,c(1,2,4), ci.type = "delta", ci.level = 0.95, 
                              em = F)
inter_out$dframe
