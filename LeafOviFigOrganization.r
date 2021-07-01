#Figure 1
#a Leachate Type
Trial1Plot
#b Mixed trial

#c Effect of sterilization
MixFig
#d antibiotics
Antiplot

dev.off()
tiff("Figures/AntiPlot1.tiff", width = 84, height = 84,units = 'mm', res = 600)
Antiplot

dev.off()


#Figure2

#a Kaplin-Meyer
SurvivalPlot1
dev.off()
tiff("Figures/SurvivalPlot1.tiff", width = 120, height = 120,units = 'mm', res = 600)
SurvivalPlot1

dev.off()


#b suvivor2


#c Growth
GrowthPlot

dev.off()
tiff("Figures/GrowthPlot.tiff", width = 120, height = 120,units = 'mm', res = 600)
GrowthPlot

dev.off()

#d Relative abu of instar


InstarPlot<-ggplot(PercentInstar, aes(x=Concentration,y=Percentage,fill=Instar))+geom_bar(stat="identity")+facet_wrap(~Leachate)+scale_fill_viridis_d()+
  theme(legend.justification=c(1,0), legend.position=c(0.99,-.08))+ylab("Relative abundance of each instar (%)")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
InstarPlot
InstarPlot

dev.off()
tiff("Figures/InstarPlot2.tiff", width = 125, height = 125,units = 'mm', res = 600)
InstarPlot

dev.off()


#GGARRANGE doesnt work with survival1 plot
# dev.off()
# tiff("Figures/Fig2.tiff", width = 174, height = 174, units = 'mm', res = 600)
# 
# ggarrange(SurvivalPlot1,GrowthPlot,InstarPlot,
#           labels = c("a", "b","c"),
#           ncol = 2, nrow = 2)
# 
# dev.off()

#Water Chemistry figure
#a
ConductivityPlot
#b
DOPlot


dev.off()
tiff("Figures/WaterChem.tiff", width = 174, height = 174, units = 'mm', res = 600)

ggarrange(ConductivityPlot,DOPlot,
                     labels = c("A", "B"),
                     ncol = 1, nrow = 2)
dev.off()

#Bacterial Figure
#a Simpson's D
SimpsonPlot


dev.off()
tiff("Figures/SimpsonPlot.tiff", width = 84, height = 84,units = 'mm', res = 600)
SimpsonPlot
dev.off()


#b Shannon
ShannonPlot

dev.off()
tiff("Figures/ShannonPlot.tiff", width = 84, height = 84,units = 'mm', res = 600)
ShannonPlot
dev.off()
#c#Variance plot

#d PCoA
ordplot


set.seed(1245)
theme_set(theme_bw(base_size = 15)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))#Set the general plotting theme see gggplot2 package for more info on other themes

ord=ordinate(physeq,"PCoA", "wunifrac")
ordplot=plot_ordination(physeq, ord,"samples", color="Treatment",shape="Treatment")+geom_point(size=5)+ 
  scale_color_brewer(palette = "Dark2")+facet_wrap(~Trial)+theme(legend.position = "bottom",legend.title = element_blank())#+theme(legend.justification=c(1,0), legend.position=c(.9,.05))#+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot


dev.off()
tiff("Figures/OrdPlotBacteriaWunifrac.tiff", width = 120, height = 120,units = 'mm', res = 600)
ordplot
dev.off()

theme_set(theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))#Set the general plotting theme see gggplot2 package for more info on other themes

#a Phylum level 
PhylumPlot

theme_set(theme_bw(base_size = 13)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))#Set the general plotting theme see gggplot2 package for more info on other themes

dev.off()
tiff("Figures/BacterialPhylumPlot.tiff", width = 84, height = 84,units = 'mm', res = 600)
PhylumPlot
dev.off()

#b family level?

dev.off()
tiff("Figures/BacterialFamilyLevel.tiff", width = 174, height = 174,units = 'mm', res = 600)
FamilyPlot
dev.off()

#c #genus level indicators long
theme_set(theme_bw(base_size = 15)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))#Set the general plotting theme see gggplot2 package for more info on other themes

cdataplot

dev.off()
tiff("Figures/BacterialGenusLevelIndicators.tiff", width = 174, height = 90,units = 'mm', res = 600)
cdataplot
dev.off()


#Fungal fig
#phyla plot
PhylumTrial
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))#Set the general plotting theme see gggplot2 package for more info on other themes

dev.off()
tiff("Figures/FungalPhylumTrial.tiff", width = 100, height = 100,units = 'mm', res = 600)
PhylumTrial
dev.off()

#familyplot
dev.off()
tiff("Figures/FungalFamily.tiff", width = 115, height = 115,units = 'mm', res = 600)
FamilyPlot
dev.off()




######
#Orplot
#####
set.seed(1245)

ord=ordinate(physeq,"PCoA", "wunifrac")
ordplot=plot_ordination(physeq, ord,"samples", color="Treatment",shape="Treatment")+geom_point(size=5)+ scale_color_brewer(palette = "Dark2")+facet_wrap(~Trial)#+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot

dev.off()
tiff("Figures/FungalOrdplot.tiff", width = 110, height = 110,units = 'mm', res = 600)
ordplot
dev.off()


#phylum
#Top 4-5 indicators


dev.off()
tiff("Figures/FungalTop6indicators.tiff", width = 110, height = 110,units = 'mm', res = 600)
cdataplot
dev.off()
