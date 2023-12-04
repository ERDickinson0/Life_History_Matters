###### LSM pathway analysis #####
rm(list=ls())
setwd("~/Kutz Postdoc/Caribou body condition/Data analysis")

library(ggplot2)
library(reshape)
library(dplyr)
# install "devtools"
#install.packages("devtools") 
library(devtools)
# install "plspm"
#install_github("gastonstat/plspm")
library(plspm)

# read data
data <- read.csv("caribougreenland.csv")
names(data)
str(data)

# make sure vars are numeric
data$MM_total <- as.numeric(data$MM_total)
data$larval_intensity <- as.numeric(data$larval_intensity)
data$fecundity <- as.numeric(data$fecundity)

plot(as.numeric(data$adult_intensity)~data$age_cemenetum)#
ggplot(data, aes(x=age_cemenetum, y=adult_intensity))+
  geom_point()+
  geom_smooth()

## Change pregnant and calf columns to 0,1 
data <- data %>%
  mutate(Pregnant = ifelse(Pregnant == "N",0,1),
         calf_present = ifelse(calf_present == "N", 0, 1))

#### subset for AM ####
AMdata <- subset(data, site == "AM")
AMdata <- na.omit(AMdata)

###### comparison with pregnant/not pregnant females ###
head(AMdata)
str(AMdata)
summary(aov(OG_total~Pregnant, data= AMdata))

var(AMdata$OG_total[AMdata$Pregnant == "Y"])
var(AMdata$OG_total[AMdata$Pregnant == "N"])

shapiro.test(AMdata$OG_total[AMdata$Pregnant == "Y"])
shapiro.test(AMdata$OG_total[AMdata$Pregnant == "N"])

wilcox.test(AMdata$OG_total ~ as.factor(AMdata$Pregnant) , alternative = "two.sided", var.equal = FALSE)

###
###### comparison with pregnant/not pregnant females ###
head(AMdata)
str(AMdata)
summary(aov(OG_total~Pregnant, data= AMdata))

var(AMdata$OG_total[AMdata$Pregnant == "Y"])
var(AMdata$OG_total[AMdata$Pregnant == "N"])

shapiro.test(AMdata$OG_total[AMdata$Pregnant == "Y"])
shapiro.test(AMdata$OG_total[AMdata$Pregnant == "N"])

wilcox.test(AMdata$OG_total ~ as.factor(AMdata$Pregnant) , alternative = "two.sided", var.equal = FALSE)

###
###### comparison with pregnant/not pregnant females ##

var(AMdata$carcass_weight[AMdata$Pregnant == "Y"])
var(AMdata$carcass_weight[AMdata$Pregnant == "N"])

shapiro.test(AMdata$carcass_weight[AMdata$Pregnant == "Y"])
shapiro.test(AMdata$carcass_weight[AMdata$Pregnant == "N"])

t.test(AMdata$carcass_weight ~ as.factor(AMdata$Pregnant) , alternative = "two.sided", var.equal = FALSE)

###
###### comparison with pregnant/not pregnant females ###

var(AMdata$protein_mass[AMdata$Pregnant == "Y"])
var(AMdata$protein_mass[AMdata$Pregnant == "N"])

shapiro.test(AMdata$protein_mass[AMdata$Pregnant == "Y"])
shapiro.test(AMdata$protein_mass[AMdata$Pregnant == "N"])

t.test(AMdata$protein_mass ~ as.factor(AMdata$Pregnant) , alternative = "two.sided", var.equal = FALSE)

####
###### comparison with pregnant/not pregnant females ###
mean(AMdata$back_fat[AMdata$Pregnant == "Y"])
mean(AMdata$back_fat[AMdata$Pregnant == "N"])
var(AMdata$back_fat[AMdata$Pregnant == "Y"])
var(AMdata$back_fat[AMdata$Pregnant == "N"])

shapiro.test(AMdata$back_fat[AMdata$Pregnant == "Y"])
shapiro.test(AMdata$back_fat[AMdata$Pregnant == "N"])

wilcox.test(AMdata$back_fat ~ as.factor(AMdata$Pregnant) , alternative = "two.sided", var.equal = FALSE)

####
###### comparison with pregnant/not pregnant females ###
mean(AMdata$kidney_fat[AMdata$Pregnant == "Y"])
mean(AMdata$kidney_fat[AMdata$Pregnant == "N"])
var(AMdata$kidney_fat[AMdata$Pregnant == "Y"])
var(AMdata$kidney_fat[AMdata$Pregnant == "N"])

shapiro.test(AMdata$kidney_fat[AMdata$Pregnant == "Y"])
shapiro.test(AMdata$kidney_fat[AMdata$Pregnant == "N"])

t.test(AMdata$kidney_fat ~ as.factor(AMdata$Pregnant) , alternative = "two.sided", var.equal = FALSE)

####
###### comparison with pregnant/not pregnant females ###
mean(AMdata$marrow_fat[AMdata$Pregnant == "Y"])
mean(AMdata$marrow_fat[AMdata$Pregnant == "N"])
var(AMdata$marrow_fat[AMdata$Pregnant == "Y"])
var(AMdata$marrow_fat[AMdata$Pregnant == "N"])

shapiro.test(AMdata$marrow_fat)

wilcox.test(AMdata$marrow_fat ~ as.factor(AMdata$Pregnant) , alternative = "two.sided", var.equal = FALSE)

### Exploratory plots
AMdata$fecundity_state <- factor(AMdata$fecundity_state, levels = c("none", "pregnant", "calf", "both"))

#png(file="AM_larval_fecundity.png",width=500,height=500,res=120)
ggplot(AMdata, aes(y = larval_intensity, x = fecundity_state))+
  geom_boxplot()+theme_classic()
#dev.off()

ggplot(AMdata, aes(x = kidney_fat))+
  geom_histogram(position="identity")
mean(AMdata$kidney_fat)
sd(AMdata$kidney_fat)

ggplot(AMdata, aes(x = carcass_weight))+
  geom_histogram(position="identity")
mean(AMdata$carcass_weight)
sd(AMdata$carcass_weight)
mean(AMdata$carcass_weight)+(sd(AMdata$carcass_weight)*3)

ggplot(AMdata, aes(x = protein_mass))+
  geom_histogram(position="identity")
mean(AMdata$protein_mass)
sd(AMdata$protein_mass)
mean(AMdata$protein_mass)+(sd(AMdata$protein_mass)*3)

ggplot(AMdata, aes(x = back_fat))+
  geom_histogram(position="identity")
mean(AMdata$back_fat)
sd(AMdata$back_fat)
mean(AMdata$back_fat)+(sd(AMdata$back_fat)*3)

ggplot(AMdata, aes(x = marrow_fat))+
  geom_histogram(position="identity")
mean(AMdata$marrow_fat)
sd(AMdata$marrow_fat)
mean(AMdata$marrow_fat)+(sd(AMdata$marrow_fat)*3)

# subset
AMdata <- subset(AMdata[c(5:7, 12:18, 20)])
str(AMdata)
names(AMdata)

#### Create the inner matrix. ie. the model paths #######
# rows of theinner model matrix
Cow_age <- c(0, 0, 0, 0, 0)
Calf <- c(1, 0, 0, 0, 0)
Parasite_intensity <- c(1, 1, 0, 0, 0)
Body_condition <- c(0, 1, 1, 0, 0)
Pregnant <- c(1, 0, 1, 1, 0)

# create path matrix and add column names
foot_path <- rbind(Cow_age, Calf, Parasite_intensity, Body_condition, Pregnant)
colnames(foot_path) <-  rownames(foot_path)

# check what has been made
foot_path
innerplot(foot_path)

##### Create the outer model list ####
#### This defines the columns which will be used #####
#### to inform each aspect of the inner model ####
names(AMdata)
foot_blocks <- list(1, 3, c(11), 4:8, 2)

# set the foot modes, to assume measurement of latent 
# variables is in reflective mode (I dont understand)
foot_modes <- c("A", "A", "A", "A", "B")

#### Run the pls_pm ####
foot_pls <- plspm(AMdata, foot_path, foot_blocks, modes = foot_modes,
                  boot.val = TRUE, br = 200)

# summarized results
summary(foot_pls) # OG adults arent good, marrow fat isnt good and calf present isnt good

# re run the model without marrow fat
foot_blocks <- list(1, 3, c(9,11), 4:7, 2)
foot_pls <- plspm(AMdata, foot_path, foot_blocks, modes = foot_modes,
                  boot.val = TRUE, br = 200)


summary(foot_pls)

plot(foot_pls)

# re run the model without
foot_blocks <- list(1, 3, c(11), 4:7, 2)
foot_pls <- plspm(AMdata, foot_path, foot_blocks, modes = foot_modes,
                  boot.val = TRUE, br = 200)


summary(foot_pls)

plot(foot_pls)

# plotting loadings of the outer model
plot(foot_pls, what = "weights", arr.width = 0.1)

##### Interpreting the results #####
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
foot_val = plspm(AMdata, foot_path, foot_blocks, modes = foot_modes,
                 boot.val = TRUE, br = 200)

foot_val$boot


foot_val$boot$rsq
foot_val$boot$path
foot_val$boot$total.efs
foot_val$boot$loadings

weight <- foot_val$boot$weights


write.csv(data, file = "AM_weights.csv")