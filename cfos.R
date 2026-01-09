#######Immediate early gene expression in day, night, and ALAN groups#########

#####Setup#####

#Set working Directory
setwd("~/UNR/cFos/R/cFos")


#Install libraries 

library(tidyverse) 
library(ggforce) 
library(ggsci)
library(patchwork)
library(Hmisc)
library(multcomp)

#install.packages("multcomp")

library(ggplot2)
library(ggsignif)
library(lattice)

library(lme4)

library(nlme)

library(pbkrtest)

library(pbkrtest)

library(AICcmodavg)

library(pwr)
library(moments)
library(plyr)
library(plotrix)
#install.packages("plotrix")


############################
#Codes for Table Data


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





#load data

cfos_df <- read.csv('cfos_df.csv')
#confirm data
head(cfos_df)

###data cleanup###
#Change variables to factors for analysis

cfos_df$bird=as.factor(cfos_df$bird)
cfos_df$pathway=as.factor(cfos_df$pathway)
cfos_df$area=as.factor(cfos_df$area)
cfos_df$group=as.factor(cfos_df$group)
cfos_df$active=as.factor(cfos_df$active)
#cfos_df$bird=as.numeric(cfos_df$bird) If you wanted to set something as numeric
#See setup
str(cfos_df)


###exploratory plots####
plot(cfos_df$group, cfos_df$cFos_exact, 
     ylab="cFos (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))


plot(cfos_df$group, cfos_df$ZENK_exact, 
     ylab="ZENK (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))
    




plot(cfos_df$group, cfos_df$cFos_exact, col=factor(cfos_df$pathway))



####Overall analysis

#Look at the effect of group
m01=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, data=cfos_df)

summary(m01)

#Look at the effect of group
m02=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, data=cfos_df)

summary(m02)



#Total cFos expression 

cfos_tot <- ggplot(data = cfos_df, aes(x = group, y = cFos_percent, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  geom_point(aes(fill = group), position = position_jitter(width = 0.2), alpha = 0.5, pch = 21, color = "black", size = 3.5) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  geom_signif(comparisons = list(c("ALAN", "night")), annotations="*", 
              map_signif_level=TRUE) +
  xlab('Treatment') + 
  ylab('cFos Expression (Percentage)')

cfos_tot + theme_classic()


cfos_total <- ggplot(data=cfos_df, aes(x = group, y = cFos_percent, fill = group)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  geom_point(aes(fill = group), position = position_jitter(width = 0.2), alpha = 0.5, pch = 21, color = "black", size = 3.5) + 
  stat_summary(fun = mean, geom = "point", size = 5, #color = c("yellow", "grey", "dark grey")
               show.legend = F) + 
  #theme(legend.position = "none", axis.text = "none") +
        #, title = element_text(size = 15, colour = "black", face = "bold")
       # ) +
 
  scale_fill_manual(name = "Treatment Group", values = c("Yellow", "grey", "black")) 
  #geom_signif(comparisons = list(c("ALAN", "night")), annotations="*", 
  #            map_signif_level=TRUE) +
  #xlab('Treatment') + 
  #ylab('cFos Expression (Percentage)') 

cfos_total + theme_classic()

cfos_total + theme(axis.text = element_text(size = 20))

cfos_total + theme(axis.title.x = element_text(face="bold", size=15),
                   axis.text.x  = element_text(angle=90, vjust=0.5, size=16))


#####
#Violin
#####

cfos_t <- ggplot(data=cfos_df, aes(x = group, y = cFos_percent, fill = group)) +
  geom_violin(aes(fill = group)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
 
  #theme(legend.position = "none", axis.text = "none") +
  #, title = element_text(size = 15, colour = "black", face = "bold")
  # ) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) 


cfos_t + theme_classic()


ZENK_t <- ggplot(data=cfos_df, aes(x = group, y = ZENK_percent, fill = group)) + 

  geom_violin(aes(fill = group)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  #theme(legend.position = "none", axis.text = "none") +
  #, title = element_text(size = 15, colour = "black", face = "bold")
  # ) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) 


ZENK_t + theme_classic()

###################


ZENK_total <- ggplot(data=cfos_df, aes(x = group, y = ZENK_percent, fill = group)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  geom_point(aes(fill = group), position = position_jitter(width = 0.2), alpha = 0.5, pch = 21, color = "black", size = 3.5) + 
  stat_summary(fun = mean, geom = "point", size = 5, #color = c("yellow", "grey", "dark grey")
               show.legend = F) + 
  #theme(legend.position = "right", title = element_text(size = 15, colour = "black", face = "bold")) +
  scale_fill_manual(name = "Treatment Group", values = c("Yellow", "grey", "black")) +
  #geom_signif(comparisons = list(c("ALAN", "night")), annotations="*", 
  #            map_signif_level=TRUE) +
  xlab('Treatment') + 
  ylab('ZENK Expression (Percentage)') 

ZENK_total + theme_classic()


cfos_total + theme_classic() + ZENK_total + theme_classic()


#Total ZENK expression

zenk_tot <- ggplot(data = cfos_df, aes(x = group, y = ZENK_percent, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  geom_signif(comparisons = list(c("ALAN", "night")), annotations="*", 
              map_signif_level=TRUE) +
  xlab('Treatment') + 
  ylab('ZENK Expression (Percentage)')


zenk_tot

cfos_tot + zenk_tot


#explore plots by area or pathway






#m3=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~pathway+(1| bird), family = poisson, data=cfos_df)

#summary(m3)



#linear model = lm
#Y ~ X
#glm = account for skew, (allow for different dist.)
#binomial is for 0,1

# to deal with 0-inflated data, first turn data into 0s and 1s, run a binomial distribution model
#then, remove 0s and run the rest of the data (still right skewed) using poisson distribution

####ZENK analyses####

#Look at the interaction between pathway and group
mz1=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~pathway*group+(1| bird), family = binomial, data=cfos_df)

summary(mz1)


#Look at the effect of group
mz2=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, data=cfos_df)

summary(mz2)

#Anova
a1=aov(cfos_df$ZENK_Int~cfos_df$group)
#Tukeys test
TukeyHSD(a1)



####cFos analyses####

#Look at interaction of pathway and group
mc1=glmer(cbind(cFos_Int,Total_Int-cFos) ~pathway*group+(1| bird), family = binomial, data=cfos_df)

summary(mc1)


mc2=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, data=cfos_df)

summary(mc2)


ac1=aov(cfos_df$cFos_Int~cfos_df$group)
TukeyHSD(ac1)



#looking at effect of pathway
mc3=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~pathway+(1| bird), family = binomial, data=cfos_df)

summary(mc3)



####subset####

###Subset Pathways###
###cFos#####
#cFos Extra
mpc1=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
          data=subset(cfos_df, pathway=="Extra"))

summary(mpc1)

Extra <- subset(cfos_df, pathway == c("Extra"))

#plot(Extra$group, Extra$cFos_exact, 
#     ylab="cFos in Extra (percentage)", xlab="Treatment", las=1,
#     col=c('yellow', 'grey', 'black'))



ggplot(data = Extra, aes(x = group, y = cFos_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("Miscellaneous Areas") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('cFos Expression (Percentage)')




#cFos Visual 
mpc2=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, pathway=="Visual"))

summary(mpc2)


# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -4.0854     0.2895 -14.113   <2e-16 ***
#   groupday     -0.7542     0.4366  -1.727   0.0841 .  
# groupnight   -0.7838     0.4363  -1.797   0.0724 .  


Visual <- subset(cfos_df, pathway == c("Visual"))


#plot(Visual$group, Visual$cFos_exact, 
#     ylab="cFos in the Visual Pathway (percentage)", xlab="Treatment", las=1,
#     col=c('yellow', 'grey', 'black'))


ggplot(data = Visual, aes(x = group, y = cFos_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("Visual Pathway Expression") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('cFos Expression (Percentage)')




#cFos Motor 
mpc3=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, pathway=="Movement"))

summary(mpc3)

# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -4.3754     0.2746 -15.934   <2e-16 ***
#   groupday     -0.3513     0.4131  -0.850    0.395    
# groupnight   -0.6138     0.4129  -1.487    0.137    


Movement <- subset(cfos_df, pathway == c("Movement"))



#plot(Movement$group, Movement$cFos_exact, 
#     ylab="cFos in the Motor Pathway (percentage)", xlab="Treatment", las=1,
#     col=c('yellow', 'grey', 'black'))




ggplot(data = Movement, aes(x = group, y = cFos_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("Motor Pathway") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('cFos Expression (Percentage)')




#####ZENK#####
#ZENK Extra
mpz1=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, pathway=="Extra"))

summary(mpz1)

#Extra <- subset(cfos_df, pathway == c("Extra"))

#plot(Extra$group, Extra$ZENK_exact, 
#     ylab="ZENK in Extra (percentage)", xlab="Treatment", las=1,
#     col=c('yellow', 'grey', 'black'))




ggplot(data = Extra, aes(x = group, y = ZENK_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("Miscellaneous Areas") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('ZENK Expression (Percentage)')




#ZENK Visual 
mpz2=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, pathway=="Visual"))

summary(mpz2)




# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.9383     0.2384  -8.132 4.22e-16 ***
#   groupday     -0.1840     0.3576  -0.514    0.607    
# groupnight   -0.5713     0.3579  -1.596    0.110    


#Visual <- subset(cfos_df, pathway == c("Visual"))


#plot(Visual$group, Visual$ZENK_exact, 
#     ylab="ZENK in the Visual Pathway (percentage)", xlab="Treatment", las=1,
#     col=c('yellow', 'grey', 'black'))



ggplot(data = Visual, aes(x = group, y = ZENK_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("Visual Pathway Expression") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('ZENK Expression (Percentage)')





#ZENK Motor 
mpz3=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, pathway=="Movement"))

summary(mpz3)


# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -1.89622    0.17835 -10.632   <2e-16 ***
#   groupday     0.09539    0.26752   0.357    0.721    
# groupnight  -0.40002    0.26762  -1.495    0.135    


#Movement <- subset(cfos_df, pathway == c("Movement"))



#plot(Movement$group, Movement$ZENK_exact, 
#     ylab="ZENK in the Motor Pathway (percentage)", xlab="Treatment", las=1,
#     col=c('yellow', 'grey', 'black'))




ggplot(data = Movement, aes(x = group, y = ZENK_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("Motor Pathway") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('ZENK Expression (Percentage)')







###########################################

#######Subset by Area######

#Subset to see cFos effect in group of only one area

####cFos#####

#APH
mac1=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
          data=subset(cfos_df, area=="APH"))

summary(mac1)



APH <- subset(cfos_df, area == c("APH"))




#aac1=aov(APH$cFos_Int~APH$group, p.adjust.method = "bonferroni")

#pairwise.t.test(APH$cFos_exact, APH$group, p.adjust.method="bonferroni")


plot(APH$group, APH$cFos_exact, 
     ylab="cFos in the APH (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))



ggplot(data = APH, aes(x = group, y = cFos_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("Area Parahippocampalis") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('cFos Expression (Percentage)')


APH_Fig <- ggplot(data=APH, aes(x = group, y = cFos_percent, fill = group)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  geom_point(aes(fill = group), position = position_jitter(width = 0.2), alpha = 0.5, pch = 21, color = "black", size = 3.5) + 
  stat_summary(fun = mean, geom = "point", size = 5, #color = c("yellow", "grey", "dark grey")
               show.legend = F) + 
  #theme(legend.position = "right", title = element_text(size = 15, colour = "black", face = "bold")) +
  scale_fill_manual(name = "Treatment Group", values = c("Yellow", "grey", "black")) +
  #geom_signif(comparisons = list(c("ALAN", "night")), annotations="*", 
  #            map_signif_level=TRUE) +
  #ggtitle("Area Parahippocampalis") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Treatment') + 
  ylab('cFos Expression (Percentage)') 

APH_Fig + theme_classic()


############



#MMD
mac2=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="MMD"))

summary(mac2)



MMD <- subset(cfos_df, area == c("MMD"))


plot(MMD$group, MMD$cFos_exact, 
     ylab="cFos in the MMD (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))




#E
mac3=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="E"))

summary(mac3)



# 2) adjusted models => regression table# Function to adjust p-values
adjMC <- function( mac3 ) {
  model_glht <- glht(mac3)
  model_MCadj <- summary(model_glht, test = adjusted('holm')) # Bonferroni-Holm 
  return(model_MCadj)}# Apply function to models

#model_lm0_adj <- adjMC( mac3 = model_lm0 )  Unsure what this is for
#model_lm_adj <- adjMC( mac3 = model_lm ) 




E <- subset(cfos_df, area == c("E"))


#plot(E$group, E$cFos_exact, 
#     ylab="cFos in the E (percentage)", xlab="Treatment", las=1,
#     col=c('yellow', 'grey', 'black'))


#Need to figure out how to change error bars

ggplot(data = E, aes(x = group, y = cFos_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("Entopallium") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('cFos Expression (Percentage)')



E_Fig <- ggplot(data=E, aes(x = group, y = cFos_percent, fill = group)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  geom_point(aes(fill = group), position = position_jitter(width = 0.2), alpha = 0.5, pch = 21, color = "black", size = 3.5) + 
  stat_summary(fun = mean, geom = "point", size = 5, #color = c("yellow", "grey", "dark grey")
               show.legend = F) + 
  #theme(legend.position = "right", title = element_text(size = 15, colour = "black", face = "bold")) +
  scale_fill_manual(name = "Treatment Group", values = c("Yellow", "grey", "black")) +
  #geom_signif(comparisons = list(c("ALAN", "night")), annotations="***", 
  #            map_signif_level=TRUE) +
  #geom_signif(comparisons = list(c("ALAN", "day")), annotations="***", y_position = 2.7,
  #            map_signif_level=TRUE) +
  #ggtitle("Entopallium") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Treatment') + 
  ylab('cFos Expression (Percentage)') 

E_Fig + theme_classic()




#StP
mac4=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="StP"))

summary(mac4)



StP <- subset(cfos_df, area == c("StP"))


plot(StP$group, StP$cFos_exact, 
     ylab="cFos in the StP (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))





#MLV
mac5=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="MLV"))

summary(mac5)



MLV <- subset(cfos_df, area == c("MLV"))


plot(MLV$group, MLV$cFos_exact, 
     ylab="cFos in the MLV (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))



ggplot(data = MLV, aes(x = group, y = cFos_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("Lateral Ventral Mesopallium") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('cFos Expression (Percentage)')



MLV_Fig <- ggplot(data=MLV, aes(x = group, y = cFos_percent, fill = group)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  geom_point(aes(fill = group), position = position_jitter(width = 0.2), alpha = 0.5, pch = 21, color = "black", size = 3.5) + 
  stat_summary(fun = mean, geom = "point", size = 5, #color = c("yellow", "grey", "dark grey")
               show.legend = F) + 
  #theme(legend.position = "right", title = element_text(size = 15, colour = "black", face = "bold")) +
  scale_fill_manual(name = "Treatment Group", values = c("Yellow", "grey", "black")) +
  #geom_signif(comparisons = list(c("ALAN", "night")), annotations="**", 
  #            map_signif_level=TRUE) +
  #ggtitle("Lateral Ventral Mesopallium") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Treatment') + 
  ylab('cFos Expression (Percentage)') 

MLV_Fig + theme_classic()










#Eco
mac6=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="Eco"))

summary(mac6)



Eco <- subset(cfos_df, area == c("Eco"))


plot(Eco$group, Eco$cFos_exact, 
     ylab="cFos in the Eco (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#NE
mac7=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="NE"))

summary(mac7)



NE <- subset(cfos_df, area == c("NE"))


plot(NE$group, NE$cFos_exact, 
     ylab="cFos in the NE (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))





#StE
mac8=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="StE"))

summary(mac8)



StE <- subset(cfos_df, area == c("StE"))


plot(StE$group, StE$cFos_exact, 
     ylab="cFos in the StE (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#MVe
mac9=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="MVe"))

summary(mac9)



MVe <- subset(cfos_df, area == c("MVe"))


plot(MVe$group, MVe$cFos_exact, 
     ylab="cFos in the MVe (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))





#PMD
mac10=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="PMD"))

summary(mac10)



PMD <- subset(cfos_df, area == c("PMD"))


plot(PMD$group, PMD$cFos_exact, 
     ylab="cFos in the PMD (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))




#PH
mac11=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="PH"))

summary(mac11)



PH <- subset(cfos_df, area == c("PH"))


plot(PH$group, PH$cFos_exact, 
     ylab="cFos in the PH (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))







#AN
mac12=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="AN"))

summary(mac12)



AN <- subset(cfos_df, area == c("AN"))


plot(AN$group, AN$cFos_exact, 
     ylab="cFos in the AN (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#Ast
mac13=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="Ast"))

summary(mac13)



Ast <- subset(cfos_df, area == c("Ast"))


plot(Ast$group, Ast$cFos_exact, 
     ylab="cFos in the Ast (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#AH
mac14=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="AH"))

summary(mac14)



AH <- subset(cfos_df, area == c("AH"))


plot(AH$group, AH$cFos_exact, 
     ylab="cFos in the AH (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#AMD
mac15=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="AMD"))

summary(mac15)



AMD <- subset(cfos_df, area == c("AMD"))


plot(AMD$group, AMD$cFos_exact, 
     ylab="cFos in the AMD (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))





#AMV
mac16=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="AMV"))

summary(mac16)



AMV <- subset(cfos_df, area == c("AMV"))


plot(AMV$group, AMV$cFos_exact, 
     ylab="cFos in the AMV (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#PLN
mac17=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="PLN"))

summary(mac17)



PLN <- subset(cfos_df, area == c("PLN"))


plot(PLN$group, PLN$cFos_exact, 
     ylab="cFos in the PLN (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#PLMV
mac18=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="PLMV"))

summary(mac18)



PLMV <- subset(cfos_df, area == c("PLMV"))


plot(PLMV$group, PLMV$cFos_exact, 
     ylab="cFos in the PLMV (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))




#DLN
mac19=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="DLN"))

summary(mac19)



DLN <- subset(cfos_df, area == c("DLN"))


plot(DLN$group, DLN$cFos_exact, 
     ylab="cFos in the DLN (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))




#LAI
mac20=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="LAI"))

summary(mac20)



LAI <- subset(cfos_df, area == c("LAI"))


plot(LAI$group, LAI$cFos_exact, 
     ylab="cFos in the LAI (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#MVb
mac21=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="MVb"))

summary(mac21)



MVb <- subset(cfos_df, area == c("MVb"))


plot(MVb$group, MVb$cFos_exact, 
     ylab="cFos in the MVb (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))





#Nb
mac22=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="Nb"))

summary(mac22)



Nb <- subset(cfos_df, area == c("Nb"))


plot(Nb$group, Nb$cFos_exact, 
     ylab="cFos in the Nb (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))




#S
mac23=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="S"))

summary(mac23)



S <- subset(cfos_df, area == c("S"))


plot(S$group, S$cFos_exact, 
     ylab="cFos in the S (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))


ggplot(data = S, aes(x = group, y = cFos_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  #stat_summary(fun.data = mean_cl_boot, 
  #geom = "errorbar", 
  #color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("cFos Expression") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('Septum (Percentage)')





#HP
mac24=glmer(cbind(cFos_Int,Total_Int-cFos_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="HP"))

summary(mac24)



HP <- subset(cfos_df, area == c("HP"))


plot(HP$group, HP$cFos_exact, 
     ylab="cFos in the S (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))


ggplot(data = HP, aes(x = group, y = cFos_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  #stat_summary(fun.data = mean_cl_boot, 
  #geom = "errorbar", 
  #color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("cFos Expression") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('Hippocampus (Percentage)')





########################################


####ZENK#####

#APH
maz1=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="APH"))

summary(maz1)



#APH <- subset(cfos_df, area == c("APH"))


plot(APH$group, APH$ZENK_exact, 
     ylab="ZENK in the APH (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))



#MMD
maz2=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="MMD"))

summary(maz2)



#MMD <- subset(cfos_df, area == c("MMD"))


plot(MMD$group, MMD$ZENK_exact, 
     ylab="ZENK in the MMD (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))

ggplot(data = MMD, aes(x = group, y = ZENK_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("Medial Dorsal Mesopallium") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('ZENK Expression (Percentage)')





MMD_Z_Fig <- ggplot(data=MMD, aes(x = group, y = ZENK_percent, fill = group)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  geom_point(aes(fill = group), position = position_jitter(width = 0.2), alpha = 0.5, pch = 21, color = "black", size = 3.5) + 
  stat_summary(fun = mean, geom = "point", size = 5, #color = c("yellow", "grey", "dark grey")
               show.legend = F) + 
  #theme(legend.position = "right", title = element_text(size = 15, colour = "black", face = "bold")) +
  scale_fill_manual(name = "Treatment Group", values = c("Yellow", "grey", "black")) +
  #geom_signif(comparisons = list(c("ALAN", "night")), annotations="***", 
  #            map_signif_level=TRUE) +
  #geom_signif(comparisons = list(c("ALAN", "day")), annotations="***", y_position = 43,
  #            map_signif_level=TRUE) +
  #ggtitle("Medial Dorsal Mesopallium") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Treatment') + 
  ylab('ZENK Expression (Percentage)') 

MMD_Z_Fig + theme_classic()






#E
maz3=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="E"))

summary(maz3)



#E <- subset(cfos_df, area == c("E"))

ggplot(data = E, aes(x = group, y = ZENK_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("Entopallium") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('ZENK Expression (Percentage)')



E_Z_Fig <- ggplot(data=E, aes(x = group, y = ZENK_percent, fill = group)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  geom_point(aes(fill = group), position = position_jitter(width = 0.2), alpha = 0.5, pch = 21, color = "black", size = 3.5) + 
  stat_summary(fun = mean, geom = "point", size = 5, #color = c("yellow", "grey", "dark grey")
               show.legend = F) + 
  #theme(legend.position = "right", title = element_text(size = 15, colour = "black", face = "bold")) +
  scale_fill_manual(name = "Treatment Group", values = c("Yellow", "grey", "black")) +
  #geom_signif(comparisons = list(c("ALAN", "night")), annotations="***", 
  #            map_signif_level=TRUE) +
  #geom_signif(comparisons = list(c("ALAN", "day")), annotations="*", y_position = 21,
  #            map_signif_level=TRUE) +
  #ggtitle("Entopallium") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Treatment') + 
  ylab('ZENK Expression (Percentage)') 

E_Z_Fig + theme_classic()





plot(E$group, E$ZENK_exact, 
     ylab="ZENK in the E (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))




#StP (CSt)
maz4=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="StP"))

summary(maz4)



#StP <- subset(cfos_df, area == c("StP"))


plot(StP$group, StP$ZENK_exact, 
     ylab="ZENK in the StP (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))





#MLV
maz5=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="MLV"))

summary(maz5)



#MLV <- subset(cfos_df, area == c("MLV"))


plot(MLV$group, MLV$ZENK_exact, 
     ylab="ZENK in the MLV (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#Eco
maz6=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="Eco"))

summary(maz6)



#Eco <- subset(cfos_df, area == c("Eco"))


plot(Eco$group, Eco$ZENK_exact, 
     ylab="ZENK in the Eco (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))





#NE
maz7=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="NE"))

summary(maz7)



#NE <- subset(cfos_df, area == c("NE"))


plot(NE$group, NE$ZENK_exact, 
     ylab="ZENK in the NE (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))





#StE
maz8=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="StE"))

summary(maz8)



#StE <- subset(cfos_df, area == c("StE"))


plot(StE$group, StE$ZENK_exact, 
     ylab="ZENK in the StE (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#MVe
maz9=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="MVe"))

summary(maz9)



#MVe <- subset(cfos_df, area == c("MVe"))


plot(MVe$group, MVe$ZENK_exact, 
     ylab="ZENK in the MVe (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))





#PMD
maz10=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="PMD"))

summary(maz10)



#PMD <- subset(cfos_df, area == c("PMD"))


plot(PMD$group, PMD$ZENK_exact, 
     ylab="ZENK in the PMD (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))




#PH
maz11=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="PH"))

summary(maz11)



#PH <- subset(cfos_df, area == c("PH"))


plot(PH$group, PH$ZENK_exact, 
     ylab="ZENK in the PH (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))







#AN
maz12=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="AN"))

summary(maz12)



#AN <- subset(cfos_df, area == c("AN"))


plot(AN$group, AN$ZENK_exact, 
     ylab="ZENK in the AN (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#Ast
maz13=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="Ast"))

summary(maz13)



#Ast <- subset(cfos_df, area == c("Ast"))


plot(Ast$group, Ast$ZENK_exact, 
     ylab="ZENK in the Ast (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#AH
maz14=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="AH"))

summary(maz14)



#AH <- subset(cfos_df, area == c("AH"))


plot(AH$group, AH$ZENK_exact, 
     ylab="ZENK in the AH (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#AMD
maz15=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="AMD"))

summary(maz15)



#AMD <- subset(cfos_df, area == c("AMD"))


plot(AMD$group, AMD$ZENK_exact, 
     ylab="ZENK in the AMD (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))





#AMV
maz16=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="AMV"))

summary(maz16)



#AMV <- subset(cfos_df, area == c("AMV"))


plot(AMV$group, AMV$ZENK_exact, 
     ylab="ZENK in the AMV (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#PLN
maz17=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="PLN"))

summary(maz17)



#PLN <- subset(cfos_df, area == c("PLN"))


plot(PLN$group, PLN$ZENK_exact, 
     ylab="ZENK in the PLN (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#PLMV
maz18=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="PLMV"))

summary(maz18)



#PLMV <- subset(cfos_df, area == c("PLMV"))


plot(PLMV$group, PLMV$ZENK_exact, 
     ylab="ZENK in the PLMV (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))




#DLN
maz19=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="DLN"))

summary(maz19)



#DLN <- subset(cfos_df, area == c("DLN"))


plot(DLN$group, DLN$ZENK_exact, 
     ylab="ZENK in the DLN (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))




#LAI
maz20=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="LAI"))

summary(maz20)



#LAI <- subset(cfos_df, area == c("LAI"))


plot(LAI$group, LAI$ZENK_exact, 
     ylab="ZENK in the LAI (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))






#MVb
maz21=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="MVb"))

summary(maz21)



#MVb <- subset(cfos_df, area == c("MVb"))


plot(MVb$group, MVb$ZENK_exact, 
     ylab="ZENK in the MVb (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))





#Nb
maz22=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="Nb"))

summary(maz22)



#Nb <- subset(cfos_df, area == c("Nb"))


plot(Nb$group, Nb$ZENK_exact, 
     ylab="ZENK in the Nb (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))



#S
maz23=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="S"))

summary(maz23)



S <- subset(cfos_df, area == c("S"))


plot(S$group, S$ZENK_exact, 
     ylab="ZENK in the S (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))


ggplot(data = S, aes(x = group, y = ZENK_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  #stat_summary(fun.data = mean_cl_boot, 
               #geom = "errorbar", 
               #color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("ZENK Expression") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('Septum (Percentage)')





#HP
maz24=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
            data=subset(cfos_df, area=="HP"))

summary(maz24)



HP <- subset(cfos_df, area == c("HP"))


plot(HP$group, HP$ZENK_exact, 
     ylab="ZENK in the HP (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))


ggplot(data = HP, aes(x = group, y = ZENK_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("Hippocampus") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('ZENK Expression (Percentage)')





HP_Z_Fig <- ggplot(data=HP, aes(x = group, y = ZENK_percent, fill = group)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  geom_point(aes(fill = group), position = position_jitter(width = 0.2), alpha = 0.5, pch = 21, color = "black", size = 3.5) + 
  stat_summary(fun = mean, geom = "point", size = 5, #color = c("yellow", "grey", "dark grey")
               show.legend = F) + 
  #theme(legend.position = "right", title = element_text(size = 15, colour = "black", face = "bold")) +
  scale_fill_manual(name = "Treatment Group", values = c("Yellow", "grey", "black")) +
  #geom_signif(comparisons = list(c("ALAN", "night")), annotations="*", 
  #            map_signif_level=TRUE) +
  #ggtitle("Hippocampus") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Treatment') + 
  ylab('ZENK Expression (Percentage)') 

HP_Z_Fig + theme_classic()







####

#Check Distribution 
hist(cfos_df$cFos_exact)
hist(cfos_df$ZENK_exact)


#Expression ~ Group*Area

#g=glmer(cFos_exact  ~group*area + (1|bird), family=poisson, data=cfos_df)
#g
#summary(g)

#areaPH<-subset(cfos_df, area=PH)
#glm_ph <- glm(cFos_Int ~ group, family=poisson, data = areaPH)
#summary(glm_ph)

#glmer_ph <- glmer(cFos_Int ~ group + (1|bird), family=poisson, data = areaPH)
#summary(glmer_ph)


#Use binomial code, cbind(y-y1)


#####Basic Plots of all data#####


ggplot(data=cfos_df, aes(x = group, y = cFos, color = group)) +
  facet_wrap(~area, scales = "free") +
  geom_point(aes(fill = group), pch = 21, color = "black") + 
  scale_fill_manual(name = "group", values = c("yellow", "grey", "black")) +
  stat_summary(fun = mean, geom = "point", color = "black") + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black")


ggplot(data=cfos_df, aes(x = group, y = ZENK, color = group)) +
  facet_wrap(~area, scales = "free") +
  geom_point(aes(fill = group), pch = 21, color = "black") + 
  scale_fill_manual(name = "group", values = c("yellow", "grey", "black")) +
  stat_summary(fun = mean, geom = "point", color = "black") + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black")


ggplot(data=cfos_df, aes(x = group, y = ZENK, color = group)) +
  facet_wrap(~pathway, scales = "free") +
  geom_point(aes(fill = group), pch = 21, color = "black") + 
  scale_fill_manual(name = "group", values = c("yellow", "grey", "black")) +
  stat_summary(fun = mean, geom = "point", color = "black") + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black")



ggplot(data=cfos_df, aes(x = group, y = cFos, color = group)) +
  facet_wrap(~pathway, scales = "free") +
  geom_point(aes(fill = group), pch = 21, color = "black") + 
  scale_fill_manual(name = "group", values = c("yellow", "grey", "black")) +
  stat_summary(fun = mean, geom = "point", color = "black") + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black")


