###### LSM pathway analysis #####
rm(list=ls())
setwd("~/Kutz Postdoc/Caribou body condition/Data analysis")

library(dplyr)
# install "devtools"
#install.packages("devtools") 
library(devtools)
# install "plspm"
#install_github("gastonstat/plspm")
library(plspm)

# read data
data <- read.csv("caribougreenland.csv")
nqmes(data)
str(data)

# make sure vars are numeric

data$TB_total <- as.numeric(data$TB_total)
data$MM_total <- as.numeric(data$MM_total)
data$fecundity <- as.numeric(data$fecundity)
data$fecundity_state <- as.factor(data$fecundity_state)
data$fecundity_state <- factor(data$fecundity_state, levels = c("none", "pregnant", "calf", "both"))
levels(data$fecundity_state)

## Change pregnant and calf columns to 0,1 
data <- data %>%
  mutate(pregnant = ifelse(Pregnant == "N",0,1),
         calfatheel = ifelse(calf_present == "N", 0,1))

data$both <- ifelse(data$fecundity_state == "both", 1, 0)
#data$pregnant <- ifelse(data$fecundity_state == "both", 0,  data$pregnant)
#data$calfatheel <- ifelse(data$fecundity_state == "both", 0, data$calfatheel)

#### subset for KS ####
KSdata <- subset(data, site == "KS")
#KSdata2 <- na.omit(KSdata)
KSdata <- KSdata[!is.na(KSdata$MM_total),]

ggplot(KSdata, aes(x = kidney_fat))+
  geom_histogram(position="identity")
mean(KSdata$kidney_fat)
sd(KSdata$kidney_fat)
mean(KSdata$kidney_fat)+(sd(KSdata$kidney_fat)*3)


ggplot(KSdata, aes(x = carcass_weight))+
  geom_histogram(position="identity")
mean(KSdata$carcass_weight)
sd(KSdata$carcass_weight)
mean(KSdata$carcass_weight)+(sd(KSdata$carcass_weight)*3)

ggplot(KSdata, aes(x = protein_mass))+
  geom_histogram(position="identity")
mean(KSdata$protein_mass)
sd(KSdata$protein_mass)
mean(KSdata$protein_mass)+(sd(KSdata$protein_mass)*3)

ggplot(KSdata, aes(x = back_fat))+
  geom_histogram(position="identity")
mean(KSdata$back_fat)
sd(KSdata$back_fat)
mean(KSdata$back_fat)+(sd(KSdata$back_fat)*3)

ggplot(KSdata, aes(x = marrow_fat))+
  geom_histogram(position="identity")
mean(KSdata$marrow_fat)
sd(KSdata$marrow_fat)
mean(KSdata$marrow_fat)+(sd(KSdata$marrow_fat)*3)

### remove outlier for protein mass. Fecundity is no longer retained in the model. 
#KSdata <- subset(KSdata, protein_mass < 14)

### remove outlier for protein mass. Fecundity is no longer retained in the model. 
KSdata <- subset(KSdata, protein_mass < 14)

###### comparison with pregnant/not pregnant females ###
head(KSdata)
str(KSdata)
summary(aov(MM_total~Pregnant, data= KSdata))

mean(KSdata$MM_total[KSdata$Pregnant == "Y"])
mean(KSdata$MM_total[KSdata$Pregnant == "N"])
var(KSdata$MM_total[KSdata$Pregnant == "Y"])
var(KSdata$MM_total[KSdata$Pregnant == "N"])

shapiro.test(KSdata$MM_total)

wilcox.test(KSdata$MM_total ~ as.factor(KSdata$Pregnant) )

###### comparison with pregnant/not pregnant females ###
head(KSdata)
str(KSdata)
summary(aov(TB_total~Pregnant, data= KSdata))

mean(KSdata$TB_total[KSdata$Pregnant == "Y"])
mean(KSdata$TB_total[KSdata$Pregnant == "N"])
var(KSdata$TB_total[KSdata$Pregnant == "Y"])
var(KSdata$TB_total[KSdata$Pregnant == "N"])

shapiro.test(KSdata$TB_total)

wilcox.test(KSdata$TB_total ~ as.factor(KSdata$Pregnant) , alternative = "two.sided", var.equal = FALSE)

###
###
###### comparison with pregnant/not pregnant females ##

var(KSdata$carcass_weight[KSdata$Pregnant == "Y"])
var(KSdata$carcass_weight[KSdata$Pregnant == "N"])

shapiro.test(KSdata$carcass_weight)

t.test(KSdata$carcass_weight ~ as.factor(KSdata$Pregnant) , alternative = "two.sided", var.equal = FALSE)

###
###### comparison with pregnant/not pregnant females ###

var(KSdata$protein_mass[KSdata$Pregnant == "Y"])
var(KSdata$protein_mass[KSdata$Pregnant == "N"])

shapiro.test(KSdata$protein_mass)

t.test(KSdata$protein_mass ~ as.factor(KSdata$Pregnant) , alternative = "two.sided", var.equal = TRUE)

####
###### comparison with pregnant/not pregnant females ###
mean(KSdata$back_fat[KSdata$Pregnant == "Y"])
mean(KSdata$back_fat[KSdata$Pregnant == "N"])
var(KSdata$back_fat[KSdata$Pregnant == "Y"])
var(KSdata$back_fat[KSdata$Pregnant == "N"])

shapiro.test(KSdata$back_fat)

wilcox.test(KSdata$back_fat ~ as.factor(KSdata$Pregnant) , alternative = "two.sided", var.equal = FALSE)

####
###### comparison with pregnant/not pregnant females ###
mean(KSdata$kidney_fat[KSdata$Pregnant == "Y"])
mean(KSdata$kidney_fat[KSdata$Pregnant == "N"])
var(KSdata$kidney_fat[KSdata$Pregnant == "Y"])
var(KSdata$kidney_fat[KSdata$Pregnant == "N"])

shapiro.test(KSdata$kidney_fat)

t.test(KSdata$kidney_fat ~ as.factor(KSdata$Pregnant) , alternative = "two.sided", var.equal = FALSE)

####
###### comparison with pregnant/not pregnant females ###
mean(KSdata$marrow_fat[KSdata$Pregnant == "Y"])
mean(KSdata$marrow_fat[KSdata$Pregnant == "N"])
var(KSdata$marrow_fat[KSdata$Pregnant == "Y"])
var(KSdata$marrow_fat[KSdata$Pregnant == "N"])

shapiro.test(KSdata$marrow_fat)

wilcox.test(KSdata$marrow_fat ~ as.factor(KSdata$Pregnant) , alternative = "two.sided", var.equal = FALSE)



## subset
nKSes(KSdata)
KSdata <- subset(KSdata[c(5, 9, 36:37, 38, 12:17, 24, 28)])
str(KSdata) 



#### Create the inner matrix. ie. the model paths #######
# rows of theinner model matrix
Cow_age <- c(0, 0, 0, 0, 0)
Calf <- c(1, 0, 0, 0, 0)
Parasite_intensity <- c(1, 1, 0, 0, 0)
Body_condition <- c(0, 1, 1, 0, 0)
Pregnant <- c(1, 0, 1, 1, 0)

# create path matrix and add column nKSes
foot_path <- rbind(Cow_age, Calf, Parasite_intensity, Body_condition, Pregnant)
colnKSes(foot_path) <-  rownKSes(foot_path)

# check what has been made
foot_path
innerplot(foot_path)

##### Create the outer model list ####
#### This defines the columns which will be used #####
#### to inform each aspect of the inner model ####
nKSes(KSdata)
foot_blocks <- list(1, 4, c(12:13), 6:10, 3)

# set the foot modes, to assume measurement of latent 
# variables is in reflective mode (I dont understand)
foot_modes <- c("A", "A", "A", "A", "B")

#### Run the pls_pm ####
foot_pls <- plspm(KSdata, foot_path, foot_blocks, modes = foot_modes,
                  boot.val = TRUE, br = 200)

# summarized results
summary(foot_pls) # larval intensity does not contribute

## Remove larval intensity
foot_blocks <- list(1, 4, 12, 6:10, 3)
foot_pls <- plspm(KSdata, foot_path, foot_blocks, modes = foot_modes,
                  boot.val = TRUE, br = 200)
summary(foot_pls)

plot(foot_pls)
# plotting loadings of the outer model
plot(foot_pls, what = "weights", arr.width = 0.5)

#### Assessing the model
# cronbach's alpha - to assess whether the blocks are good
# should be > 0.7
foot_pls$unidim[, 3, drop = FALSE]

# dillon-goldstein rho - assesses unidimentionality. 
# > 0.7 = unidimensional (which is good)
foot_pls$unidim[, 4, drop = FALSE]

# eigenvalues
foot_pls$unidim[, 5:6]

## model assessment
foot_pls$inner_summary

foot_pls$gof

#### vaildation using bootstrapping #####
# running bootstrap validation
foot_val = plspm(KSdata, foot_path, foot_blocks, modes = foot_modes,
                 boot.val = TRUE, br = 200)
summary(foot_val)

foot_val$boot

foot_val$boot$weights
foot_val$boot$rsq
foot_val$boot$path
foot_val$boot$total.efs

