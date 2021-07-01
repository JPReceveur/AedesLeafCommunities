
library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggpubr)
library(plyr)
library(dplyr)
library(rstatix)
Trial1=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Leaf oviposition\\AllEggDataMergeLossRemoved.txt",header=TRUE)
MixTrial=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Leaf oviposition\\MixedtrialLossRemoved.txt",header=TRUE)
Survivalship=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Leaf oviposition\\LabSurvivalR.txt",header=TRUE)
head(Survivalship)
theme_set(theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


###################
#Trial1
##################
sum(Trial1$Eggs)
head(Trial1)
av=aov( Eggs~ Leaf+as.character(Block), data=Trial1)
summary(av)
hist(resid(av))
posthoc <- TukeyHSD(av, 'Leaf', conf.level=0.95)
posthoc

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = Eggs ~ Leaf + as.character(Block), data = Trial1)
# 
# $Leaf
# diff       lwr        upr     p adj
# Oak-Honey       263.6026  121.1940  406.01112 0.0000896
# Sycamore-Honey   15.7619 -124.1841  155.70788 0.9901417
# Water-Honey    -133.3000 -271.0761    4.47614 0.0609571
# Sycamore-Oak   -247.8407 -384.8577 -110.82362 0.0001263
# Water-Oak      -396.9026 -531.7026 -262.10250 0.0000000
# Water-Sycamore -149.0619 -281.2577  -16.86608 0.0219736



Trial1 %>% select(Leaf, Eggs) %>% group_by(Leaf) %>% 
  summarise(n = n(), 
            mean = mean(Eggs, na.rm = TRUE), 
            sd = sd(Eggs, na.rm = TRUE),
            stderr = sd/sqrt(n), 
            LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
            UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
            median=median(Eggs, na.rm = TRUE),
            min=min(Eggs, na.rm = TRUE), 
            max=max(Eggs, na.rm = TRUE),
            IQR=IQR(Eggs, na.rm = TRUE))


Trtdata <- ddply(Trial1, c("Leaf2"), summarise,
                 N    = length(Eggs),
                 mean = mean(Eggs),
                 sd   = sd(Eggs),
                 se   = sd / sqrt(N)
)
Trtdata
Trial1$Leaf2 = factor(Trial1$Leaf2, levels = c("OAK","SYC","HON","WAT")) #fixes x-axis labels
geom_text(aes(x=DecompStage, y=mean+se+10,label=vec))
Trial1Plot<-ggplot(Trial1, aes(x=Leaf2,y=Eggs,color=Leaf2))+geom_boxplot()+ylab("Eggs laid")+xlab("Leachate Type")+scale_color_brewer(palette="Dark2")+theme(legend.position = "none")
Trial1Plot#NEED TO ADD MULTICOMP LETTERS, OAK=A SYC=B HON=BC WAT=C

dev.off()
tiff("Figures/LeachateOviPlot.tiff", width = 84, height = 84,units = 'mm', res = 600)
Trial1Plot

dev.off()

#######################3
#MixTrial
#########################
LeafData=MixTrial
#LeafData
Block=LeafData$Block
Eggs=LeafData$Eggs
LogEggs=log(Eggs)
Leaf=LeafData$Leaf
#Removed=LeafData$Removed
#Signif=LeafData$Signif
hist(MixTrial$Eggs)
hist(MixTrial$LogEggs)
#Untransformed egg numbers used for Oak trial based on hist(resid(av))
MixTrialOak<-subset(MixTrial,MixTrial$Type=="Oak"|MixTrial$Type=="Oak_Sterile")
av=aov( Eggs~ Type+as.character(Block), data=MixTrialOak)
summary(av)
hist(resid(av))

#plot(av)

posthoc <- TukeyHSD(av, 'Type', conf.level=0.95)
posthoc


#Log transformation used for mix trial based on residual plot 
MixTrialMix<-subset(MixTrial,MixTrial$Type!="Oak"&MixTrial$Type!="Oak_Sterile")
av=aov( LogEggs~ Type+as.character(Block), data=MixTrialMix)
summary(av)
hist(resid(av))

#plot(av)

posthoc <- TukeyHSD(av, 'Type', conf.level=0.95)
posthoc
# diff        lwr         upr    p adj
# Mix_Sterile-Mix -0.2739241 -0.5360013 -0.01184687 0.042492


Trtdata <- ddply(MixTrial, c("Type"), summarise,
                 N    = length(Eggs),
                 mean = mean(Eggs),
                 sd   = sd(Eggs),
                 se   = sd / sqrt(N)
)
Trtdata

MixTrial$Type = factor(MixTrial$Type, levels = c("Oak","Oak_Sterile","Mix","Mix_Sterile")) #fixes x-axis labels

MixFig<-ggplot(MixTrial,aes(x=Type,y=Eggs))+
  geom_boxplot()+ylab("Eggs Laid")+xlab("Leachate Type")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank())
MixFig

dev.off()
tiff("Figures/MixedPlot.tiff", width = 84, height = 84,units = 'mm', res = 600)
MixFig

dev.off()
#< 0.05= *
#< 0.01= **
#< 0.001 = ***

#####################################
###############Survivalship
########################################

Survivalship=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Leaf oviposition\\LabSurvivalR.txt",header=TRUE)
head(Survivalship)
plot=ggplot(Survivalship, aes(Time,Alive))+facet_wrap(~LeafType)+stat_summary(fun=mean,geom="point", size=2)
plot
plot+stat_smooth()+geom_point()+ylab("Number Alive")+xlab("Time (hrs)")+ ylim(0, 21)
av=aov( Alive~ Time*LeafType, data=Survivalship)
summary(av)
posthoc <- TukeyHSD(av, 'LeafType', conf.level=0.95)
posthoc

###################Suvivorship Trial 2
Survivalship=read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Leaf oviposition\\Survivorship2018Trial2.csv",header=TRUE)
head(Survivalship)
TrialDay7=subset(Survivalship, LeachatePrep=="7Day")
head(TrialDay7)


TrialDay7$Concentration<-factor(TrialDay7$Concentration, levels = c("1X","1/2X","1/4X","Sterile"))

plot=ggplot(TrialDay7, aes(Time,Alive,color=Concentration,shape=Concentration,group=Concentration))+facet_wrap(~LeafType)+geom_point(size=3.5)
Alive<-plot+geom_line()+ylab("Number Alive")+xlab("Time (hrs)")+ theme(legend.justification=c(1,0), legend.position=c(.98,0.05),legend.title = element_blank())+
  scale_color_brewer(palette = "Dark2")
Alive
#ggsave("Trial2Survival9.20.2018.tiff")


dev.off()
tiff("Figures/SurvivalByLeachate.tiff", width = 120, height = 120, units = 'mm', res = 600)
Alive
dev.off()

mod <- glm(Alive ~ LeafType + Time+Concentration, family=poisson,data=TrialDay7)
summary(mod)
plot(mod)

#Day 3 graph
TrialDay3=subset(Survivalship,LeachatePrep=="3Day")
head(TrialDay3)
plot=ggplot(TrialDay3, aes(Time,Alive,color=LeafType))+facet_wrap(~LeafType)+geom_point(size=3)
plot+geom_line()+ylab("Number Alive")+xlab("Time (hrs)")+ ylim(0, 21)+xlim(0,80)+ theme(legend.justification=c(1,0), legend.position=c(1,0))
ggsave("Trial2Survival9.20.2018.tiff")



#############Kaplin-Meyer
#install.packages("survminer")
library(survival)
library(survminer)
library(dplyr)
KMCurve=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Leaf oviposition\\LabSurvivalKMC.txt",header=TRUE)
head(KMCurve)
head(KMCurve)
KMCurve$SurvObj <- with(KMCurve, Surv(Time, Status == 0))
lung$SurvObj <- with(lung, Surv(time, status == 2))

km.by.leaf <- survfit(SurvObj ~ LeafType, data = KMCurve, conf.type = "log-log")
theme_set(theme_bw(base_size = 25))

#fit= survfit(Surv(Time,SurvObj)~LeafType,data=KMCurve)
#summary(fit)
SurvivalPlot1<-ggsurvplot(km.by.leaf, data = KMCurve, risk.table = FALSE, pval=TRUE,label.curves=TRUE, conf.int = TRUE,palette="Dark2")+xlab("Time (hrs)")
SurvivalPlot1

head(KMCurve)
km.as.one <- survfit(SurvObj ~ 1, data = KMCurve, conf.type = "log-log")
km.by.leaf <- survfit(SurvObj ~ LeafType, data = KMCurve, conf.type = "log-log")
km.by.leaf
plot(km.by.leaf)
plot= plot(km.by.leaf)


#install.packages('survminer')
#library(survminer)

source("https://bioconductor.org/biocLite.R")
biocLite("RTCGA.clinical")



##############
# Growth Expt
#############
theme_set(theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

Growthdata=read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Leaf oviposition\\GrowthExptDataMay2021.csv",header=TRUE)
head(Growthdata)
Growthdata<-subset(Growthdata,Mass_mg!="REDO"|Mass_mg!="NA") #Remove Redo and NA values
GrowthdataDay7<-subset(Growthdata,Day=="Seven")
GrowthdataDay7$Mass_mg<-as.numeric(GrowthdataDay7$Mass_mg)
GrowthdataDay7$Concentration<-factor(GrowthdataDay7$Concentration, levels = c("1X","1/2X","1/4X","Sterile"))

hist(GrowthdataDay7$HeadCapsule_mm,breaks=25,xlab= "Head Capsule Width (mm)",main="")


GrowthPlot<-ggplot(GrowthdataDay7, aes(x=Concentration, y=Mass_mg,color=Concentration))+geom_boxplot()+ylab("Larval mass (mg) after 7 days")+theme(legend.position = "none",axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~Lechate)+ scale_color_brewer(palette = "Dark2")

dev.off()
tiff("Figures/GrowthPlot.tiff", width = 174, height = 174, units = 'mm', res = 600)
GrowthPlot
dev.off()


GrowthPlotV2<-ggplot(GrowthdataDay7, aes(x=Concentration, y=Mass_mg,color=Concentration))+geom_boxplot()+ylab("Mass (mg) at day 7")+theme(legend.position = "none",axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+facet_grid(~Lechate)+ scale_color_brewer(palette = "Dark2")

dev.off()
tiff("Figures/GrowthPlotV2.tiff", width = 174, height = 84, units = 'mm', res = 600)
GrowthPlotV2
dev.off()



#< 0.05= *
#< 0.01= **
#< 0.001 = ***


#Honeysuckle growth
HoneysuckleGrowth<-subset(GrowthdataDay7,Lechate=="Honeysuckle")
kruskal.test(Mass_mg~Concentration,data=HoneysuckleGrowth)
# 
# kruskal.test(Mass_mg~Concentration,data=HoneysuckleGrowth)
# 
# Kruskal-Wallis rank sum test
# 
# data:  Mass_mg by Concentration
# Kruskal-Wallis chi-squared = 10.589, df = 2, p-value = 0.005019

compare_means(Mass_mg~Concentration, data=HoneysuckleGrowth,method = "wilcox.test",p.adjust.method="fdr")
# # A tibble: 3 x 8
# .y.     group1 group2       p  p.adj p.format p.signif method  
# <chr>   <chr>  <chr>    <dbl>  <dbl> <chr>    <chr>    <chr>   
#   1 Mass_mg 1X     1/2X   0.130   0.13   0.1304   ns       Wilcoxon
# 2 Mass_mg 1X     1/4X   0.0375  0.056  0.0375   *        Wilcoxon
# 3 Mass_mg 1/2X   1/4X   0.00326 0.0098 0.0033   **       Wilcoxon

#Mix Growth
MixGrowth<-subset(GrowthdataDay7,Lechate=="Mix")
kruskal.test(Mass_mg~Concentration,data=MixGrowth)
# 
# Kruskal-Wallis rank sum test
# 
# data:  Mass_mg by Concentration
# Kruskal-Wallis chi-squared = 23.314, df = 3, p-value = 3.474e-05
compare_means(Mass_mg~Concentration, data=MixGrowth,method = "wilcox.test",p.adjust.method="fdr")
# # A tibble: 6 x 8
# .y.     group1 group2         p  p.adj p.format p.signif method  
# <chr>   <chr>  <chr>      <dbl>  <dbl> <chr>    <chr>    <chr>   
#   1 Mass_mg 1X     1/2X    0.000275 0.0013 0.00027  ***      Wilcoxon
# 2 Mass_mg 1X     1/4X    0.00267  0.0053 0.00267  **       Wilcoxon
# 3 Mass_mg 1X     Sterile 0.000440 0.0013 0.00044  ***      Wilcoxon
# 4 Mass_mg 1/2X   1/4X    0.789    0.79   0.78859  ns       Wilcoxon
# 5 Mass_mg 1/2X   Sterile 0.0115   0.017  0.01151  *        Wilcoxon
# 6 Mass_mg 1/4X   Sterile 0.106    0.13   0.10590  ns       Wilcoxon

#Oak
OakGrowth<-subset(GrowthdataDay7,Lechate=="Oak")
kruskal.test(Mass_mg~Concentration,data=OakGrowth)

# Kruskal-Wallis rank sum test
# 
# data:  Mass_mg by Concentration
# Kruskal-Wallis chi-squared = 10.366, df = 2, p-value = 0.005612

compare_means(Mass_mg~Concentration, data=OakGrowth,method = "wilcox.test",p.adjust.method="fdr")

# .y.     group1 group2        p p.adj p.format p.signif method  
# <chr>   <chr>  <chr>     <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 Mass_mg 1/2X   1/4X    0.0207  0.031 0.0207   *        Wilcoxon
# 2 Mass_mg 1/2X   Sterile 0.00451 0.014 0.0045   **       Wilcoxon
# 3 Mass_mg 1/4X   Sterile 0.245   0.25  0.2450   ns       Wilcoxon

#Sycamore
SycamoreGrowth<-subset(GrowthdataDay7,Lechate=="Sycamore")
kruskal.test(Mass_mg~Concentration,data=SycamoreGrowth)

# Kruskal-Wallis rank sum test
# 
# data:  Mass_mg by Concentration
# Kruskal-Wallis chi-squared = 26.548, df = 3, p-value = 7.321e-06

compare_means(Mass_mg~Concentration, data=SycamoreGrowth,method = "wilcox.test",p.adjust.method="fdr")

# .y.     group1 group2         p   p.adj p.format p.signif method  
# <chr>   <chr>  <chr>      <dbl>   <dbl> <chr>    <chr>    <chr>   
#   1 Mass_mg 1X     1/2X    0.00118  0.0023  0.00118  **       Wilcoxon
# 2 Mass_mg 1X     1/4X    0.000291 0.0017  0.00029  ***      Wilcoxon
# 3 Mass_mg 1X     Sterile 0.00140  0.0023  0.00140  **       Wilcoxon
# 4 Mass_mg 1/2X   1/4X    0.0464   0.046   0.04640  *        Wilcoxon
# 5 Mass_mg 1/2X   Sterile 0.00155  0.0023  0.00155  **       Wilcoxon
# 6 Mass_mg 1/4X   Sterile 0.00400  0.00480 0.00400  **       Wilcoxon


WaterGrowth<-subset(GrowthdataDay7,Lechate=="Water")
kruskal.test(Mass_mg~Concentration,data=WaterGrowth)
# 
# Kruskal-Wallis rank sum test
# 
# data:  Mass_mg by Concentration
# Kruskal-Wallis chi-squared = 0.24407, df = 1, p-value = 0.6213
compare_means(Mass_mg~Concentration, data=WaterGrowth,method = "wilcox.test",p.adjust.method="fdr")
# 
# .y.     group1 group2      p p.adj p.format p.signif method  
# <chr>   <chr>  <chr>   <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 Mass_mg 1X     Sterile 0.711  0.71 0.71     ns       Wilcoxon

head(GrowthdataDay7)


SterileDay7<-subset(GrowthdataDay7,Concentration=="Sterile")
kruskal.test(Mass_mg~Lechate,data=SterileDay7)

# Kruskal-Wallis rank sum test
# 
# data:  Mass_mg by Lechate
# Kruskal-Wallis chi-squared = 9.303, df = 3, p-value = 0.02552



Trtdata <- ddply(SterileDay7, c("Lechate"), summarise,
                 N    = length(Mass_mg),
                 mean = mean(Mass_mg),
                 sd   = sd(Mass_mg),
                 se   = sd / sqrt(N)
)
Trtdata

compare_means(Mass_mg~Lechate, data=SterileDay7,method = "wilcox.test",p.adjust.method="fdr")


head(GrowthdataDay7)


Instar <- ddply(GrowthdataDay7, c( "Lechate","Concentration","Instar"), summarise,
                 N    = length(Mass_mg)
)

Instar
PercentInstar<-read.csv("PercentageInstarByConcentration.csv",header=T)
head(PercentInstar)

PercentInstar$Concentration<-factor(PercentInstar$Concentration, levels = c("1X","1/2X","1/4X","Sterile"))

InstarPlot<-ggplot(PercentInstar, aes(x=Concentration,y=Percentage,fill=Instar))+geom_bar(stat="identity")+facet_wrap(~Leachate)+scale_fill_viridis_d()+
  theme(legend.justification=c(1,0), legend.position=c(.9,.03))+ylab("Relative abundance of each instar (%)")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
InstarPlot
dev.off()
tiff("Figures/InstarPlot.tiff", width = 174, height = 174, units = 'mm', res = 600)
InstarPlot
dev.off()
##############
#Antibiotic trial
############
#file.choose()
theme_set(theme_bw(base_size = 16)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

AntiTrial<-read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Leaf oviposition\\AntibioticTrial9.17.19.csv",header=T)
head(AntiTrial)
AntiTrial$PercentVsOak<-AntiTrial$PercentVsOak*100
Trtdata <- ddply(AntiTrial, c( "Treatment"), summarise,
                     N    = length(Eggs),
                     mean = mean(Eggs),
                     sd   = sd(Eggs),
                     se   = sd / sqrt(N)
)
Trtdata


av=aov( Eggs~ Treatment+as.character(Block), data=AntiTrial)
summary(av)
hist(resid(av))
hsd <- TukeyHSD(av, which = 'Treatment') 
hsd$Treatment



hist(AntiTrial$PercentVsOak)
AntiTrial$Treatment = factor(AntiTrial$Treatment, levels = c("AMP ","AMPHO","AMP+AMPHO","OAK")) #fixes x-axis labels

Antiplot<-ggplot(AntiTrial, aes(x=Treatment, y=Eggs))+geom_boxplot()+ theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank())+ylab("Eggs laid")
Antiplot


###############
#Water chemistry
##############

#Data labeled as trial 1 in file not used due to issues with probe (see Do or SPC/ORP values)
#Trials 2 +3 conducted at same time as survivorship experiment and averaged since from same preperation of leachate which was then split into seperate containers
#Values 



WaterChem<-read.table("LeafOviWaterChemOct2018.txt",header=T)
WaterChem<-subset(WaterChem, Trial!="One")
head(WaterChem)
WaterChem$Concentration<-factor(WaterChem$Concentration, levels = c("1X","1/2X","1/4X","Sterile"))

HoneysuckleWaterChem<-subset(WaterChem,Type=="HON")
OakWaterChem<-subset(WaterChem,Type=="OAK")
SycamoreWaterChem<-subset(WaterChem,Type=="SYC")
MixWaterChem<-subset(WaterChem,Type=="MIX")

ggplot(WaterChem, aes(x=Time, y=DO.,color=Concentration,shape=Trial))+geom_point()+facet_wrap(~Type)+geom_line()


#Water chem 2+3 combined

WaterChemCombined <- ddply(WaterChem, c( "Type","Time","Concentration"), summarise,
                 N    = length(pH),
                 pH = mean(pH),
                 DO = mean(DO.),
                 DO.mg.ul = mean(DO.MG.L.),
                 SPC.Us.cm = mean(SPC.Us.cm),
                 C.uS.cm = mean(C.uS.cm),
                 ORP.mV = mean(ORP.mV))
WaterChemCombined
ORPPlot<-ggplot(WaterChemCombined, aes(x=as.character(Time),shape=Concentration, y=ORP.mV,color=Concentration,group=Concentration))+
  geom_point(size=3)+facet_wrap(~Type)+geom_line()+ylab("Oxidative Reduction Potential (ORP, mV)")+xlab("Time (hrs)")+
  scale_fill_viridis_d()+ scale_color_brewer(palette = "Dark2")



dev.off()
tiff("Figures/ORPFIG.tiff", width = 174, height = 174, units = 'mm', res = 600)
ORPPlot
dev.off()

DOPlot<-ggplot(WaterChemCombined, aes(x=(Time),shape=Concentration, y=DO,color=Concentration,group=Concentration))+
  geom_point(size=3)+facet_grid(~Type)+geom_line()+ylab("Dissolved Oxygen (%)")+xlab("Time (hrs)")+
  theme(legend.justification=c(1,0), legend.position=c(.95,.15),legend.title = element_blank())+ scale_color_brewer(palette = "Dark2")

dev.off()
tiff("Figures/DOPlot.tiff", width = 174, height = 84, units = 'mm', res = 600)
DOPlot
dev.off()



pHPlot<-ggplot(WaterChemCombined, aes(x=as.character(Time),shape=Concentration, y=pH,color=Concentration,group=Concentration))+
  geom_point(size=3)+facet_wrap(~Type)+geom_line()+ylab("pH")+xlab("Time (hrs)")+
  scale_color_brewer(palette = "Dark2")

dev.off()
tiff("Figures/pHPlot.tiff", width = 174, height = 174, units = 'mm', res = 600)
pHPlot
dev.off()

#Conductivity
#Mostlystable over time ()

WaterChemAllHrsCombined <- ddply(WaterChem, c( "Type","Concentration"), summarise,
                           N    = length(pH),
                           pH = mean(pH),
                           DO = mean(DO.),
                           DO.mg.ul = mean(DO.MG.L.),
                           SPC.Us.cm = mean(SPC.Us.cm),
                           C.uS.cm = mean(C.uS.cm),
                           ORP.mV = mean(ORP.mV))
WaterChemCombined

ConductivityAllHours <- ddply(WaterChemCombined, c( "Type","Concentration"), summarise,
                 N    = length(C.uS.cm),
                 meanConductivity = mean(C.uS.cm),
                 sd   = sd(C.uS.cm),
                 se   = sd / sqrt(N)
)
ConductivityAllHours

WaterChemCombined$TypeConc<-paste0(WaterChemCombined$Type,":",WaterChemCombined$Concentration)
res.aov <- anova_test(C.uS.cm~ Type*Concentration,data = WaterChemCombined, within = c(TypeConc))
get_anova_table(res.aov)

av=aov( C.uS.cm~ Type,data = WaterChemCombined)
summary(av)
hist(resid(av))
posthoc <- TukeyHSD(av, 'Type', conf.level=0.95)
posthoc



Conc1X<-subset(WaterChemCombined)

ConductivityPlot<-ggplot(WaterChemCombined, aes(x=(Time),shape=Concentration, y=C.uS.cm,color=Concentration,group=Concentration))+
  geom_point(size=3)+facet_grid(~Type)+geom_line()+ylab(expression("Conductivity ("~mu~"S/cm)"))+xlab("Time (hrs)")+
  theme(legend.justification=c(1,0), legend.position=c(.95,.15),legend.title = element_blank())+ scale_color_brewer(palette = "Dark2")
ConductivityPlot
dev.off()
tiff("Figures/ConductivityPlot.tiff", width = 174, height = 84, units = 'mm', res = 600)
ConductivityPlot
dev.off()




pHAllHours <- ddply(WaterChemCombined, c( "Type"), summarise,
                              N    = length(pH),
                              meanpH = mean(pH),
                              sd   = sd(pH),
                              se   = sd / sqrt(N)
)
pHAllHours

