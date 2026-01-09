
#Set working Directory
setwd("~/UNR/cFos/R/cFos")


#Install libraries 

library(tidyverse) 
library(ggforce) 
library(ggsci)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(Hmisc)




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



#####
#Violin Total IEG
#####

cfos_t <- ggplot(data=cfos_df, aes(x = group, y = cFos_percent, fill = group)) +
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black", size = 1, width = 0.75) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  #theme(legend.position = "none", axis.text = "none") +
  #, title = element_text(size = 15, colour = "black", face = "bold")
  # ) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) 


cfos_t + theme_classic()


ZENK_t <- ggplot(data=cfos_df, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black", size = 1, width = 0.75) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  #theme(legend.position = "none", axis.text = "none") +
  #, title = element_text(size = 15, colour = "black", face = "bold")
  # ) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) 


ZENK_t + theme_classic()

###################


##cFos##

E <- subset(cfos_df, area == c("E"))


E_t <- ggplot(data=E, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) 


E_t + theme_classic()


MLV <- subset(cfos_df, area == c("MLV"))


MLV_t <- ggplot(data=MLV, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black", size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) 


MLV_t + theme_classic()


AMD <- subset(cfos_df, area == c("AMD"))

AMD_t <- ggplot(data=AMD, aes(x = group, y = cFos_percent, fill = group)) +
  geom_violin(size = 1) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black", size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F"))

AMD_t + theme_classic()


StE_t <- ggplot(data=StE, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) 


StE_t + theme_classic()



PH_t <- ggplot(data=PH, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) 


PH_t + theme_classic()



####ZENK


E_t <- ggplot(data=E, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black", size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) 

E_t + theme_classic()




MMD <- subset(cfos_df, area == c("MMD"))


MMD_t <- ggplot(data=MMD, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black", size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) 

MMD_t + theme_classic()


HP <- subset(cfos_df, area == c("HP"))


HP_t <- ggplot(data=HP, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black", size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) 

HP_t + theme_classic()

############
#Behavior
##########


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



activity <- ggplot(data=behav_df, aes(x = behavior, y = percent, fill = treatment)) +
  facet_wrap(~treatment) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black", size = 1, width = 0.75) +
  geom_point(aes(fill = treatment), position = position_jitter(width = 0.24), pch = 21, color = "black", size = 3.5) + 
  stat_summary(fun = mean, geom = "point", size = 5, #color = c("yellow", "grey", "dark grey")
               show.legend = F) + 
  theme(legend.position = "right", title = element_text(size = 15, colour = "black", face = "bold")) +
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  xlab('Behavior') + 
  ylab('Time (Percentage)') 



activity + theme_classic()



###########################
#Supplemental 
###########################


##cFos

# AH = Anterior Hyperpallium

AH_t <- ggplot(data=AH, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


AH_t


# AMV = Anterior Mesopallium Ventral

AMV_t <- ggplot(data=AMV, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


AMV_t


# AN = Anterior Nidopallium


AN_t <- ggplot(data=AN, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


AN_t

# APH = Area Parahippocampalis
APH_t <- ggplot(data=APH, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


APH_t



# ASt = Anterior Striatum

ASt_t <- ggplot(data=Ast, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


ASt_t

# CSt = Caudal Striatum (StP)

StP_t <- ggplot(data=StP, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


StP_t


# DLN = Dorsal Lateral Nidopallium

DLN_t <- ggplot(data=DLN, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


DLN_t

# Eco = Core of the Entopallium

Eco_t <- ggplot(data=Eco, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


Eco_t

# HP = Hippocampus

HP_t <- ggplot(data=HP, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


HP_t

# LAi = Lateral Intermediate Arcopallium\

LAi_t <- ggplot(data=LAI, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


LAi_t


# MLV = Lateral Ventral Mesopallium

MLV_t <- ggplot(data=MLV, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


MLV_t

# MMD = Medial Dorsal Mesopallium

MMD_t <- ggplot(data=MMD, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


MMD_t

# MVb = Ventral Mesopallium adjacent to Basorostral Nucleus

MVb_t <- ggplot(data=MVb, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


MVb_t

# MVe = Ventral Mesopallium near Eco

MVe_t <- ggplot(data=MVe, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


MVe_t

# Nb = Nidopallium adjacent to Basorostral Nucleus

Nb_t <- ggplot(data=Nb, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


Nb_t

# NE = Nidopallium adjacent to Eco

NE_t <- ggplot(data=NE, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


NE_t

# PLMV = Posterior Lateral Ventral Mesopallium

PLMV_t <- ggplot(data=PLMV, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


PLMV_t

# PLN = Posterior Lateral Nidopallium

PLN_t <- ggplot(data=PLN, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


PLN_t

# PMD = Posterior Dorsal Mesopallium

PMD_t <- ggplot(data=PMD, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


PMD_t

# S = Septum

S_t <- ggplot(data=S, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


S_t

# StE = Striatum adjacent to Eco

StE_t <- ggplot(data=StE, aes(x = group, y = cFos_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


StE_t






########

##ZENK


# AH = Anterior Hyperpallium

AH_t <- ggplot(data=AH, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


AH_t

# AMD = Anterior Mesopallium Dorsal

AMD_t <- ggplot(data=AMD, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


AMD_t


# AMV = Anterior Mesopallium Ventral

AMV_t <- ggplot(data=AMV, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


AMV_t


# AN = Anterior Nidopallium


AN_t <- ggplot(data=AN, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


AN_t

# APH = Area Parahippocampalis
APH_t <- ggplot(data=APH, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


APH_t



# ASt = Anterior Striatum

ASt_t <- ggplot(data=Ast, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


ASt_t

# CSt = Caudal Striatum (StP)

StP_t <- ggplot(data=StP, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


StP_t


# DLN = Dorsal Lateral Nidopallium

DLN_t <- ggplot(data=DLN, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


DLN_t

# Eco = Core of the Entopallium

Eco_t <- ggplot(data=Eco, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


Eco_t


# LAi = Lateral Intermediate Arcopallium

LAi_t <- ggplot(data=LAI, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


LAi_t


# MLV = Lateral Ventral Mesopallium

MLV_t <- ggplot(data=MLV, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


MLV_t


# MVb = Ventral Mesopallium adjacent to Basorostral Nucleus

MVb_t <- ggplot(data=MVb, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


MVb_t

# MVe = Ventral Mesopallium near Eco

MVe_t <- ggplot(data=MVe, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


MVe_t

# Nb = Nidopallium adjacent to Basorostral Nucleus

Nb_t <- ggplot(data=Nb, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


Nb_t

# NE = Nidopallium adjacent to Eco

NE_t <- ggplot(data=NE, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


NE_t

# PH = Posterior Hyperpallium

PH_t <- ggplot(data=PH, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


PH_t

# PLMV = Posterior Lateral Ventral Mesopallium

PLMV_t <- ggplot(data=PLMV, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


PLMV_t

# PLN = Posterior Lateral Nidopallium (NCL)

PLN_t <- ggplot(data=PLN, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


PLN_t

# PMD = Posterior Dorsal Mesopallium

PMD_t <- ggplot(data=PMD, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


PMD_t

# S = Septum

S_t <- ggplot(data=S, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


S_t

# StE = Striatum adjacent to Eco

StE_t <- ggplot(data=StE, aes(x = group, y = ZENK_percent, fill = group)) + 
  geom_violin(size = 1) + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black",size = 1, width = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 4, color = c("black"),
               show.legend = F) +
  
  scale_fill_manual(name = "Treatment Group", values = c("#F9D731", "#C4C4C3", "#50504F")) +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=15, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=12, family= "Times"),
    axis.title.y = element_text( size=15, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=19, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


StE_t


