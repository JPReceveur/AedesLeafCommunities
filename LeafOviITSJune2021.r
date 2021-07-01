#Leaf Ovi ITS June 2021

set.seed(3246)
#Parse Silva function

parse_taxonomy_silva_138 <- function(char.vec){
  # Use default to assign names to elements in case problem with greengenes prefix
  char.vec = parse_taxonomy_default(char.vec)
  # Check for unassigned taxa
  if (char.vec["Rank1"] == "Unassigned") {
    char.vec <- c(Rank1="d__Unassigned", Rank2="d__Unassigned", Rank3="d__Unassigned", Rank4="d__Unassigned",
                  Rank5="d__Unassigned", Rank6="d__Unassigned", Rank7="d__Unassigned")
  }
  # Define the meaning of each prefix according to SILVA taxonomy
  Tranks = c(Rank1="Kingdom", Rank2="Phylum", Rank3="Class", Rank4="Order", Rank5="Family", Rank6="Genus", Rank7="Species")
  # Check for prefix using regexp, warn if there were none. trim indices, ti
  ti = grep("[[:alpha:]]\\_\\_", char.vec)
  if( length(ti) == 0L ){
    warning(
      "No silva prefixes were found. \n",
      "Consider using parse_taxonomy_delfault() instead if true for all OTUs. \n",
      "Dummy ranks may be included among taxonomic ranks now."
    )
    # Will want to return without further modifying char.vec
    taxvec = char.vec
    # Replace names of taxvec according to prefix, if any present...
  } else {
    # Format character vectors for Ambiguous taxa
    if( length(ti) < 7 ){
      for (key in names(char.vec)) {
        if ( char.vec[key] == "Ambiguous_taxa" ) {
          tax_no <- (as.numeric(substr(key, 5, 5)) - 1)
          char.vec[key] = sprintf("d__Ambiguous_taxa", tax_no)
        }
      }
      # Reset the trimmed indicies if Ambiguous taxa
      ti = grep("[[:alpha:]]\\_\\_", char.vec)
    }
    # Remove prefix using sub-"" regexp, call result taxvec
    taxvec = gsub("[[:alpha:]]\\_\\_", "", char.vec)
    # Define the ranks that will be replaced
    repranks = Tranks[substr(char.vec[ti], 1, 3)]
    # Replace, being sure to avoid prefixes notK present in Tranks
    names(taxvec)[ti[!is.na(repranks)]] = repranks[!is.na(repranks)]
  }
  return(taxvec)
}


#Load the packages needed
#To install these packages, run the command install.packages("package name here"), you can also install multiple packages at a time (see R documentation online for more info)
library(vegan)
library(MASS)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(ape)
library(ggpubr)
library(viridis)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7") #Create a user defined color palette 
theme_set(theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))#Set the general plotting theme see gggplot2 package for more info on other themes
biom=import_biom("LeafOviITS23Jan2021_4k.biom",parseFunction=parse_taxonomy_silva_138)
colnames(tax_table(biom))<-c(Rank1="Kingdom", Rank2="Phylum", Rank3="Class", Rank4="Order", Rank5="Family", Rank6="Genus", Rank7="Species")
#Import biom file and use the parse function built above
#If you put the file inside the folder of the R project, you can call it by just typing the file name as in the line above. IF its somewhere else and you need to find the path use the command file.choose() to bring up a browser to find the file path.
#biom
#tax_table(biom)
#tax_table(biom) <- tax_table(biom)[,-c(5:10,14)]#remove dummy ranks #Not needed

metadata=read.table("LeafOviMetadataJan2021.txt",header = TRUE) #Pull in your metadata and name the R object as metadata
metadata$Treatment = factor(metadata$Treatment, levels = c("OAK","SYC","HON","MIX","WAT","NC")) #fixes x-axis labels

tree=read_tree("LOITSTree.nwk") #Read in the phylogenetic tree of the bacteria that was built in earlier processing


#This next few lines takes the individual files and combines them into a single object to make it easier to work with
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$id
physeq=merge_phyloseq(biom,sampdat,tree)

physeq<-filter_taxa(physeq, function (x) {sum(x > 0) > 1}, prune=TRUE)#Remove singletons and zeros
physeq
physeq<-subset_samples(physeq,Treatment!="WAT") #Removal of water only samples (rarefaction already removed extraction and sequencing controls, none above 1,400 reads)
physeq


#############
#Alpha diversity
#############
ObservedSpecies<-plot_richness(physeq, x="Treatment",color="Treatment", measures=c("Observed"))+ylab("Observed Species")+geom_boxplot()+scale_fill_manual(cbPalette)+ scale_color_brewer(palette = "Dark2")
ObservedSpecies
ObservedRichness<-ObservedSpecies$data
head(ObservedRichness)

kruskal.test(value~Treatment, data=ObservedRichness)
ddply(ObservedRichness, c("Treatment"), summarise,
      N    = length(value),
      mean = mean(value),
      sd   = sd(value),
      se   = sd / sqrt(N)
)
ggplot(ObservedRichness, aes(x=Treatment,y= value))+geom_boxplot()

#Simpson
InvSimpson<-plot_richness(physeq, x="Treatment",color="Treatment", measures=c("invsimpson"))+ylab("Inverse Simpson ")+scale_fill_manual(cbPalette)+geom_boxplot()
InvSimpson
SimpsonData<-InvSimpson$data

kruskal.test(value~Treatment, data=SimpsonData)
compare_means(value~Treatment, data=SimpsonData,method = "wilcox.test",p.adjust.method="fdr")


ddply(ObservedRichness, c("Treatment"), summarise,
      N    = length(value),
      mean = mean(value),
      sd   = sd(value),
      se   = sd / sqrt(N)
)
SimpsonPlot<-ggplot(SimpsonData, aes(x=Treatment,y= value,color=Treatment))+geom_boxplot()+ylab("Inverse Simpson (1/D)")+ scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")+xlab("Leachate Type")

SimpsonPlot

#Shannon
Shannon<-plot_richness(physeq, x="Treatment",color="Treatment", measures=c("shannon"))+ylab("Shannon Diversity")+scale_fill_manual(cbPalette)+facet_wrap(~Trial)
Shannon
ShannonData<-Shannon$data
kruskal.test(value~Treatment, data=ShannonData)

############3
#Fungal taxa plots
############

GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
PhylumAll=tax_glom(GPr, "Phylum")
PhylumLevel = filter_taxa(PhylumAll, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%
FamilyAll=tax_glom(GPr,"Family")
FamilyLevel = filter_taxa(FamilyAll, function(x) mean(x) > 3e-2, TRUE) #filter out any taxa lower than 3%
GenusAll=tax_glom(GPr,"Genus")
GenusLevel = filter_taxa(GenusAll, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%

df <- psmelt(PhylumAll)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
#ddply from the ddplyr package is being used here to calculate the mean, standard deviation and standard error of your samples, see below for other examples using it to do the same thing but also between different sample groups
Trtdata


df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100

compare_means(Abundance~Trial, data=df,group.by="Phylum",method = "kruskal.test",p.adjust.method="fdr")
compare_means(Abundance~Treatment, data=df,group.by="Phylum",method = "wilcox.test",p.adjust.method="fdr")

Trtdata <- ddply(df, c("Phylum","Treatment","Trial"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd=sd(Abundance),
                 se=sd/sqrt(N)
)
Trtdata


PhylumTrial=ggplot(Trtdata, aes(x=Treatment,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+xlab("Leachate")+
ylab("Relative Abundance (%, > 1%)") + theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank())+scale_fill_viridis_d()+facet_wrap(~Trial)
#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette)
PhylumTrial



############
#Family plot
############

df <- psmelt(FamilyLevel)
df$Abundance=df$Abundance*100


compare_means(Abundance~Treatment, data=df,group.by="Family",method = "kruskal.test",p.adjust.method="fdr")
compare_means(Abundance~Treatment, data=df,group.by="Family",method = "wilcox.test",p.adjust.method="fdr")


summary(aov(Abundance~Treatment,group.by="Family",data=df))
#Sig families Erwiniaceae, Oxalobacteraceae,
Means<-compare_means(Abundance~Treatment, data=df,group.by="Family",method = "wilcox.test",p.adjust.method="fdr")

Trtdata <- ddply(df, c("Family","Treatment"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd=sd(Abundance),
                 se=sd / sqrt(N)
)
Trtdata
Trtdata[17,1]<-"Microbotryomycetes\n inc.sed."
Trtdata[18,1]<-"Microbotryomycetes\n inc.sed."
Trtdata[19,1]<-"Microbotryomycetes\n inc.sed."
Trtdata[20,1]<-"Microbotryomycetes\n inc.sed."
Trtdata[21,1]<-"Rhyncho-\n gastremataceae"
Trtdata[22,1]<-"Rhyncho-\n gastremataceae"
Trtdata[23,1]<-"Rhyncho-\n gastremataceae"
Trtdata[24,1]<-"Rhyncho-\n gastremataceae"


FamilyPlot=ggplot(Trtdata, aes(x=Treatment,y=mean))+geom_bar(aes(fill = Treatment),colour="black", stat="identity")+
  xlab("Leachate")+ylab("Relative Abundance (%, > 1%)") + theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none",axis.title.x=element_blank())+
  scale_fill_brewer(palette="Dark2")+facet_wrap(~Family)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))
#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette)
FamilyPlot





#################
#Fungal beta diversity
###################

set.seed(1245)
GPdist=phyloseq::distance(physeq, "jaccard")

adonis(GPdist ~ Treatment*as.factor(Trial), as(sample_data(physeq), "data.frame"))


set.seed(1245)
GPdist=phyloseq::distance(physeq, "jaccard")
beta=betadisper(GPdist, sample_data(physeq)$Treatment)
permutest(beta)
boxplot(beta,xlab="Leachate Type",col=brewer.pal(n = 4, name = "Dark2"))



################3
#Random Forest indicators
################

set.seed(1245)
ForestData=GenusLevel#Change this one so you dont have to rewrite all variables
predictors=t(otu_table(ForestData))
response <- as.factor(sample_data(ForestData)$Treatment)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000,importance=T)

print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseAccuracy))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]#


destroyX = function(es) {
  f = es
  for (row in c(1:nrow(f))){ #for each column in dataframe
    if (startsWith(row.names(f)[row], "X") == TRUE)  { #if starts with 'X' ..
      row.names(f)[row] <- substr(row.names(f)[row], 2, 100) #get rid of it
    }
    assign(deparse(substitute(es)), f, inherits = TRUE)
  }
}
destroyX(imp.20)




###Top 6 genus level indicator taxa for differences between Urban V Rural (MDA)

imp.10 <- imp.sort[1:6, ]#
destroyX(imp.10)
otunames <- row.names(imp.10)
r <- rownames(tax_table(ForestData)) %in% otunames
GenusRandomForestSubset = subset_taxa(GenusAll, row.names(tax_table(GenusAll))%in% otunames)
#GenusRandomForestSubset

df <- psmelt(GenusRandomForestSubset)
df$Abundance=df$Abundance*100
compare_means(Abundance ~ Treatment, data = df, group.by = "Genus", p.adjust.method = "fdr",method="kruskal.test")
Means<-compare_means(Abundance ~ Treatment, data = df, group.by = "Genus", p.adjust.method = "fdr",method="wilcox.test")


SigList<-list()
for (i in (Means$Genus)){
  Tax<-i
  TaxAbundance<-subset(Means,Genus==i )
  Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
  difference<-TaxAbundance$p.adj
  names(difference)<-Hyphenated
  Letters<-multcompLetters(difference)
  #print(Letters)
  SigList[i]<-Letters
  
}
SigList
vec<-unlist(SigList)
vec<-vec[-1]



Trtdata <- ddply(df, c("Genus","Treatment"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd=sd(Abundance),
                 se=sd/sqrt(N)
)
Trtdata


cdataplot=ggplot(Trtdata, aes(x=Treatment,y=mean))+geom_bar(aes(fill = Treatment),colour="black", stat="identity")+ facet_wrap(~Genus)+xlab("Leachate")+ylab("Relative Abundance (%)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_fill_brewer(palette = "Dark2")+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+
  theme(legend.position = "none",axis.title.x=element_blank())
cdataplot


