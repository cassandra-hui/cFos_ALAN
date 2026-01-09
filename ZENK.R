#######Immediate early gene expression in day, night, and ALAN groups#########

#####Setup#####

#Set working Directory
setwd("~/UNR/cFos/R")


#Install libraries 

library(tidyverse) 
library(ggforce) 
library(ggsci)
library(patchwork)
library(Hmisc)
library(multcomp)

#install.packages("multcomp")

library(ggplot2)
library(lattice)

library(lme4)

library(nlme)

library(pbkrtest)

library(pbkrtest)

library(AICcmodavg)

library(pwr)

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




#Total ZENK expression

ggplot(data = cfos_df, aes(x = group, y = ZENK_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("All Areas") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('ZENK Expression (Percentage)')


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




#E
maz3=glmer(cbind(ZENK_Int,Total_Int-ZENK_Int) ~group+(1| bird), family = binomial, 
           data=subset(cfos_df, area=="E"))

summary(maz3)



#E <- subset(cfos_df, area == c("E"))


plot(E$group, E$ZENK_exact, 
     ylab="ZENK in the E (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))




#StP
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



#S <- subset(cfos_df, area == c("S"))


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



#HP <- subset(cfos_df, area == c("HP"))


plot(HP$group, HP$ZENK_exact, 
     ylab="ZENK in the HP (percentage)", xlab="Treatment", las=1,
     col=c('yellow', 'grey', 'black'))


ggplot(data = HP, aes(x = group, y = ZENK_exact, fill = group)) +
  geom_sina(size = 4, pch = 21) +
  #stat_summary(fun.data = mean_cl_boot, 
  #geom = "errorbar", 
  #color = "black") +
  stat_summary(fun = mean, geom = "point", size = 5, colour = "black") +
  ggtitle("ZENK Expression") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", title = element_text(size = 20, colour = "black", face = "bold")) + 
  scale_fill_manual(name = "group:", values = c("yellow", "grey", "black")) + 
  xlab('Treatment') + 
  ylab('Hippocampus (Percentage)')


