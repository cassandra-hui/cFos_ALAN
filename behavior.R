#####Setup#####

#Set working Directory
setwd("~/UNR/cFos/R/cFos")


#Install libraries 

library(tidyverse) 
library(ggforce) 
library(ggsci)
library(patchwork)
library(Hmisc)
library(ggplot2)
library(moments)
library(plyr)
library(plotrix)
library(lme4)



##########################################
#Behavior
##########################################


#load data

behav_df <- read.csv('Behavior.csv')
#confirm data
head(behav_df)
str(behav_df)

###data cleanup###
#Change variables to factors for analysis

behav_df$bird=as.factor(behav_df$bird)
behav_df$treatment=as.factor(behav_df$treatment)
behav_df$behavior=as.factor(behav_df$behavior)
behav_df$total=as.numeric(behav_df$total)

require(pscl)
require(MASS)
require(boot)






inputData <- read.csv('Dir.csv')

kruskal.test(inputData$Grooming~inputData$Treatment)
kruskal.test(inputData$Feeding~inputData$Treatment)
kruskal.test(inputData$Hoping~inputData$Treatment)
kruskal.test(inputData$Inactive~inputData$Treatment)


pairwise.wilcox.test(inputData$Feeding, inputData$Treatment, p.adjust.method="BH")

#number of hops by treatment group, with kruskal test 

hop <- read.csv('hop_data.csv')

kruskal.test(hop$hops~hop$treatment)

pairwise.wilcox.test(hop$hops, hop$treatment, p.adjust.method="BH")


hopsf <- ggplot(data=hop, aes(x = treatment, y = hops, fill = treatment)) +
  #facet_wrap(~treatment) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  geom_point(aes(fill = treatment), pch = 21, color = "black", size = 5) + 
  stat_summary(fun = mean, geom = "point", size = 5, #color = c("yellow", "grey", "black")
  ) + 
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) +
  scale_fill_manual(name = "group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  ggtitle("Perch Hops 75 to 105 Minutes before Perfusion") +
  xlab('Treatment') + 
  ylab('Number of Hops') 


hopsf







############
library (DirichletReg)


library (DirichletReg)
inputData <- ArcticLake  # plug-in your data here.
set.seed(100)
train <- sample (1:nrow (inputData), round (0.7*nrow (inputData)))  # 70% training sample
inputData_train <- inputData [train, ] # training Data
inputData_test <- inputData [-train, ] # test Data
inputData$Y <- DR_data (inputData[,1:3])  # prepare the Y's
inputData_train$Y <- DR_data (inputData_train[,1:3])
inputData_test$Y <- DR_data (inputData_test[,1:3])

res1 <- DirichReg(Y ~ depth + I(depth^2), inputData_train)  # modify the predictors and input data here
res2 <- DirichReg(Y ~ depth + I(depth^2) | depth, inputData_train, model="alternative") 
summary(res1)

inputData <- read.csv('Dir.csv')

kruskal.test(inputData$Grooming~inputData$Treatment)
kruskal.test(inputData$Feeding~inputData$Treatment)
kruskal.test(inputData$Hoping~inputData$Treatment)
kruskal.test(inputData$Inactive~inputData$Treatment)


pairwise.wilcox.test(inputData$Feeding, inputData$Treatment, p.adjust.method="BH")

#number of hops by treatment group, with kruskal test 

hop <- read.csv('hop_data.csv')

kruskal.test(hop$hops~hop$treatment)

pairwise.wilcox.test(hop$hops, hop$treatment, p.adjust.method="BH")


hopsf <- ggplot(data=hop, aes(x = treatment, y = hops, fill = treatment)) +
  #facet_wrap(~treatment) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  geom_point(aes(fill = treatment), pch = 21, color = "black", size = 5) + 
  stat_summary(fun = mean, geom = "point", size = 5, #color = c("yellow", "grey", "black")
  ) + 
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) +
  scale_fill_manual(name = "group", values = c("Yellow", "grey", "black")) +
  ggtitle("Perch Hops 75 to 105 Minutes before Perfusion") +
  xlab('Treatment') + 
  ylab('Number of Hops') 


hopsf

################################


inputdat <-read.csv('Dir.csv')
inputdat$Y=DR_data (inputdat[,1:4])
inputdat_train$Y

r1= DirichReg(Y~Treatment, inputdat)
summary(r1)

#x <- as.vector(inputData[, 5])

set.seed(100)
train <- sample (1:nrow (inputData), round (0.7*nrow (inputData)))  # 70% training sample
inputData_train <- inputData [train, ] # training Data
inputData_test <- inputData [-train, ] # test Data
inputData$Y <- DR_data (inputData[,1:4])  # prepare the Y's
inputData_train$Y <- DR_data (inputData_train[,1:4])
inputData_test$Y <- DR_data (inputData_test[,1:4])



res1 <- DirichReg(Y ~ Treatment + I(Treatment^2), inputData)  # modify the predictors and input data here
res2 <- DirichReg(Y ~ Treatment + I(Treatment^2) | Treatment, inputData_train, model="alternative") 

summary(res1)


#####################
library(DirichletReg)
library(Compositional)

Data <- read.csv('Dir.csv')


x <- as.vector(Data[, 5])
y <- as.matrix(Data[, 1:4])
y <- y / rowSums(y)
y <- na.omit(y)

mod1 <- diri.reg(y, x)
mod2 <-diri.reg2(y, x)
mod3 <- comp.reg(y, x)


diri.reg(y, x, plot = FALSE, xnew = NULL)

diri.reg2(y, x, xnew = NULL)






####################


str(behav_df)


b1=glmer(cbind(time, total - time) ~ treatment*behavior +(1|bird) , family = binomial,
        data = behav_df)

summary(b1)

behav_df$intime=behav_df$time*100

m1 <- zeroinfl(intime ~ treatment*behavior,
               data = behav_df, dist = "negbin", EM = TRUE)


res1 <- DirichReg(time ~ treatment*behavior, data=behav_df) 
b2=glmer(cbind(time, total - time) ~ treatment +(1| bird), family = binomial,
         data=subset(behav_df, behavior=="Inactive"))

summary(b2)


b3=glmer(cbind(time, total - time) ~ treatment +(1| bird), family = binomial,
         data=subset(behav_df))

summary(b3)



#b4=glmer(cbind(time, total - time) ~ behavior +(1| bird), family = binomial,
#         data=subset(behav_df))

#summary(b4)



##############################


activity <- ggplot(data=behav_df, aes(x = behavior, y = percent, fill = treatment)) +
  facet_wrap(~treatment) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  geom_point(aes(fill = treatment), position = position_jitter(width = 0.24), pch = 21, color = "black", size = 3.5) + 
  stat_summary(fun = mean, geom = "point", size = 5, #color = c("yellow", "grey", "dark grey")
              show.legend = F) + 
  theme(legend.position = "right", title = element_text(size = 15, colour = "black", face = "bold")) +
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  xlab('Behavior') + 
  ylab('Time (Percentage)') 
  


activity + theme_classic()







ggplot( data=behav_df, aes(x = behavior, y = percent, fill = treatment))+ #manually set your variables
  scale_fill_manual(values = c("yellow", "grey", "dark grey")) + #manually set colors
  geom_point(
    alpha = 0.2,
    position = position_jitter(width = 0.1),
    show.legend = F
  ) + # the raw data points w/ jitter adjust the width of the jitter to your taste. 
  stat_summary(geom = "point", fun = mean, size = 3.25, alpha = 0.94, show.legend = F) + #mean of raw data
  stat_summary(geom = "errorbar", fun.data = mean_cl_normal,
               fun.args=list(conf.int=0.95), width = 0.2, size = 0.8, alpha = 0.9, show.legend = F ) + #error bars of raw data
  theme_classic( ) +
  facet_wrap(~treatment, row_number(1))  +
  labs(
    x = "Behavior",
    y = "Percent of Time",
    title = "whatever"
  )




#act=lm(total_active_minutes ~treatment, 
#       data=behav_df)

#summary(act)


#TukeyHSD(aov(act), conf.level=.95)







activity <- ggplot(data=behav_df, aes(x = treatment, y = total_active_minutes, fill = treatment)) +
  #facet_wrap(~area, scales = "free") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  geom_point(aes(fill = treatment), pch = 21, color = "black") + 
  stat_summary(fun = mean, geom = "point", size = 5, color = c("yellow", "grey", "black")) + 
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) +
  scale_fill_manual(name = "group", values = c("yellow", "grey", "black")) +
  
  #geom_signif(y_position = c(7.3,25.3), xmin = c(0.8,1.8), 
  #           xmax = c(1.2,2.2), annotation = c("NS","**"),
  #          tip_length = 0)+
  #geom_signif(comparisons = list(c("ALAN","control day")), y_position = 28,
  #            tip_length = .2, vjust = .1)+
  
  
  #geom_signif(comparisons = list(c("ALAN", "control night")), annotations="ns", 
  #           map_signif_level=TRUE) +
  #geom_signif(comparisons = list(c("ALAN", "control day")), annotations="ns", 
#           map_signif_level=TRUE) +
#ggtitle("75 to 105 Minutes before Perfusion") +
#  xlab('Treatment') + 
#  ylab('Total Minutes Active') 

#activity 

#one.wayb <- aov(total_active_minutes ~ treatment, data = behav_df)

#summary(one.wayb)


##########
# INFO
##########

# simple frequency table
freq_major <- table(behav_df$total_active_minutes)
freq_major

summary(behav_df$treatment)

#for mode
getmode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}
# And you can use this function
getmode(behav_df$total_active_minutes, na.rm = TRUE)


#############################
# something a little bit advanced
# you can write your own function to calcuate IQR
############################
getIQR <- function(x, na.rm = FALSE){
  
  # Remove missing values if specified
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  # Calculate IQR
  x_iqr <- as.vector(diff(quantile(x, c(0.25, 0.75))))
  
  # return to the IQR value
  return(x_iqr)
}
# to use your function
getIQR(behav_df$total_active_minutes, na.rm = T)
############################




tab_activity <- ddply(behav_df, .(treatment), summarise,
                      mean = mean(total_active_minutes, na.rm = T),
                      median = median(total_active_minutes, na.rm = T),
                      mode = getmode(total_active_minutes, na.rm = T),   
                      # make sure that you have the getmode function and have run it in R!!!
                      # getmode function is not a built-in function in R
                      # it is created by you in this script
                      
                      range = diff(range(total_active_minutes, na.rm = T)),
                      IQR = getIQR(total_active_minutes, na.rm = T),
                      # same for the getIQR function
                      
                      var = var(total_active_minutes, na.rm = T),
                      sd = sd(total_active_minutes, na.rm = T),
                      SE = std.error(total_active_minutes, na.rm = T),
                      
                      skew = skewness(total_active_minutes, na.rm = T),
                      kur = kurtosis(total_active_minutes, na.rm = T)
)

# Print out the table
tab_activity

