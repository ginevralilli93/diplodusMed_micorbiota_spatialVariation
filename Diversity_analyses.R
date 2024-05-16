# Packages ----------------------------------------------------------------------
#Loading the packages that we will use
pkgs <- c("tidyverse", "phyloseq", "ggpubr", "ggplot2", 
          "vegan", "DESeq2","microbiome","colorspace",
          "picante","RColorBrewer","zCompositions","ggord","forecast",
          "eulerr","philr","pairwiseAdonis", "olsrr")

lapply(pkgs, require, character.only = TRUE)
#set the graphic enviroment
theme_set(theme_minimal())

#############################################################################################
# Functions ---

# Delete the ASVs with only 0s (common after subsetting a phyloseq object)
nozero_taxa<-function(phy_rare){
  colsum.nozero<-(colSums(otu_table(phy_rare), na.rm=T) !=0)
  nozero1<-otu_table(phy_rare)[,colsum.nozero] #matrix of asv with only the OTU that do not have all zeros 
  
  #recreate the same phyloseq object including only the txa without 0 as rowsum
  phy_rare2<-phyloseq(otu_table(nozero1, taxa_are_rows = FALSE),
                      tax_table(tax_table(phy_rare)), taxa_names(taxa_names(phy_rare)),
                      sample_data(sample_data(phy_rare)),phy_tree(phy_tree(phy_rare)))
  
}

nozero_taxaNOtree<-function(phy_rare){
  colsum.nozero<-(colSums(otu_table(phy_rare), na.rm=T) !=0)
  nozero1<-otu_table(phy_rare)[,colsum.nozero] #matrix of asv with only the OTU that do not have all zeros 
  
  #recreate the same phyloseq object including only the txa without 0 as rowsum
  phy_rare2<-phyloseq(otu_table(nozero1, taxa_are_rows = FALSE),
                      tax_table(tax_table(phy_rare)), taxa_names(taxa_names(phy_rare)),
                      sample_data(sample_data(phy_rare)))
  
}

# Rarefy 
tran.rare<-function(phy.obj, sample.size){
  rare<-rarefy_even_depth(phy.obj,rngseed = 123,sample.size=sample.size, replace=FALSE)
  #select only the asv that have colSums !=0. In this way the next steps can be performed
  colsum.nozero<-(colSums(otu_table(rare), na.rm=T) !=0)
  nozero<-otu_table(phy.obj)[,colsum.nozero] #matrix of asv with only the OTU that do not have all zeros 
  #recreate the same phyloseq object including only the txa without 0 as rowsum
  filt.obj<-phyloseq(otu_table(nozero, taxa_are_rows = FALSE),
                     tax_table(tax_table(rare)), taxa_names(taxa_names(rare)),
                     sample_data(sample_data(rare)),phy_tree(phy_tree(rare)))
}

#Centered log-ratio trasformation - without agglomeration of the ASV to higher taxonomical levels
tran.clr.not.glom<-function(phy.obj){
  
  #agglomerate at the genus level
  #phy.obj.glom<-tax_glom(phy.obj,taxrank = "Genus")
  #select only the asv that have colSums !=0. In this way the next steps can be performed
  colsum.nozero<-(colSums(otu_table(phy.obj), na.rm=T) !=0)
  nozero<-otu_table(phy.obj)[,colsum.nozero] #matrix of asv with only the OTU that do not have all zeros 
  
  #transform the 0s in probability vectors
  asv.prob<-zCompositions::cmultRepl(nozero, method = 'CZM', delta = 0.5, output = 'p-counts')
  
  #create an alternative ps object with the asv table with probabilities 
  prob<-phyloseq::phyloseq(otu_table(asv.prob, taxa_are_rows = FALSE),
                           tax_table(tax_table(phy.obj)), taxa_names(taxa_names(phy.obj)),
                           sample_data(sample_data(phy.obj)),phy_tree(phy_tree(phy.obj)))
  
  clr <- microbiome::transform(prob, "clr")  
  
}

#Compute prevalence of each feature, store as data.frame
prevalence<-function(phy_object,taxrank){
  phy.glom<-tax_glom(phy_object,taxrank = {{taxrank}}, NArm = FALSE)
  prevdf = apply(X = otu_table(phy.glom),
                 MARGIN = ifelse(taxa_are_rows(phy.glom), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  
  #Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(phy.glom),
                      tax_table(phy.glom))
  
  #add colum with prevalence/number of samples
  prevdf$Prevalence_percent<-prevdf$Prevalence / nsamples(phy.glom)
  
  prevdf
  
}

#create the alpha div table 
alpha.div<-function(phy,metadata){ #phy: phyloseq object rarefied
  alpha_div<- phy %>% 
    estimate_richness(split= T, measures = c("Shannon","Simpson","Observed")) %>%
    mutate(X.SampleID = sample_names(phy)) # add sample IDs
  colnames(alpha_div)[3]<-"FishID"
  #estimate richness by phylogeny: Faith's index - calculate alpha diversity
  tree<-phy_tree(phy)
  asv.table<-otu_table(phy)
  faith.index<-pd(asv.table,tree,include.root = FALSE)
  faith.index$FishID<-rownames(faith.index)
  
  #add the next column with the alpha diversity indexes in the metadata table.
  alpha_div <- cbind(metadata, # add metadata
                     alpha_div,
                     faith.index,
                     by = "FishID")
  #  alpha_div<-alpha_div[,-c(63:64)] #remove "by" and "FishID" since there are several columns called like that
}


stacked.barplot<-function(ps.object,taxrank,x.axis.value,grid){ 
  p=plot_bar(ps.object,x=x.axis.value,fill=taxrank)
  p=p+geom_bar(stat = "identity")
  p=p+facet_wrap(grid,scale="free_x")
  p=p+theme(legend.key.size = unit(3, "mm"), 
            legend.position = c("bottom"), 
            legend.box = c("horizontal"),
            legend.text = element_text(size=18,face = "italic"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  return(p)
}

# Ordination plots - PCA - MPA
pcaMPA<-function(phy.clr,palette){
  pca.clr <- ordinate(phy.clr, "RDA") #Performs redundancy analysis, or optionally principal components analysis, via rda
  sample_data(phy.clr)$MPA<-as.character(sample_data(phy.clr)$MPA)
  plot_ordination(phy.clr,pca.clr, color="MPA")+
    geom_point(size=4, alpha=1, aes(color=as.character(MPA) ))+
    scale_color_manual(values = c("pink2","purple","orange"))+    
    ggtitle("Taxonomic dissimilarity of the gut microbiota",
            "Beta diversity, Aitchinson's distance")+
    #  facet_wrap(~Species)+
    #stat_ellipse(geom = "polygon",  alpha=0.09,size=1)+
    theme(axis.text.x = element_text(size=15))+
    theme(axis.text.y = element_text(size=15))+
    theme(axis.title.x  = element_text(size=15))+
    theme(axis.title.y  = element_text(size=15))+
    theme(legend.text = element_text(face = "italic", size= 15))+
    theme(legend.title = element_text( size= 15))
}

# Welch MANOVA
Welch.Manova<-function(clr,var){ 
  #Robust distance-based multivariate analysis of variance (https://doi.org/10.1186/s40168-019-0659-9)
  matrix<-phyloseq::distance(clr,method="euclidean") #generate a distance matrix
  df<- sample_data(clr) %>% data.frame() # extract the dataframe
  df[,var]<-as.factor(df[,var]) 
  set.seed(1000)
  WDS<-WdS.test(dm = matrix, f = df[,var], nrep =999) 
}

#PERMANOVA
permanova<-function(clr,var){ 
  #Robust distance-based multivariate analysis of variance (https://doi.org/10.1186/s40168-019-0659-9)
  matrix<-phyloseq::distance(clr,method="euclidean") #generate a distance matrix
  df<- sample_data(clr) %>% data.frame() # extract the dataframe
  df[,var]<-as.factor(df[,var]) 
  set.seed(1000)
  permanova<-adonis2(matrix~df[,var])
}

permanovaLength<-function(clr,var){ 
  clr.sub<-subset_samples(clr,!is.na(Length))
  #Robust distance-based multivariate analysis of variance (https://doi.org/10.1186/s40168-019-0659-9)
  matrix<-phyloseq::distance(clr.sub,method="euclidean") #generate a distance matrix
  df<- sample_data(clr.sub) %>% data.frame() # extract the dataframe
  set.seed(1000)
  permanova<-adonis2(matrix~df[,var])
}

# PERMDISP - betadisperser

permdisp<-function(clr,var){ 
  matrix<-phyloseq::distance(clr,method="euclidean") #generate a distance matrix
  df<- sample_data(clr) %>% data.frame() # extract the dataframe
  df$MPA<-as.factor(df$MPA) 
  dispr<-betadisper(matrix,df[,var]) #betadisperser test 
  permutest(dispr) 
}

# PAIRWISE ADONIS 
pairwise.adonis<-function(clr,var){ 
  matrix<-phyloseq::distance(clr,method="euclidean") #generate a distance matrix
  df<- sample_data(clr) %>% data.frame() # extract the dataframe
  df[,var]<-as.factor(df[,var]) 
  set.seed(1000)
  pairwise<-pairwiseAdonis::pairwise.adonis2(matrix, factors=df[,var]) 
}

#create the alpha div table 
alpha.div.Kegg<-function(phy,metadata){ #phy: phyloseq object rarefied
  alpha_div<- phy %>% 
    estimate_richness(split= T, measures = c("Shannon","Observed")) %>%
    mutate(X.SampleID = sample_names(phy)) # add sample IDs
  colnames(alpha_div)[3]<-"FishID"
  # #estimate richness by phylogeny: Faith's index - calculate alpha diversity 
  # tree<-phy_tree(phy)
  # asv.table<-otu_table(phy)
  # faith.index<-pd(asv.table,tree,include.root = FALSE)
  # faith.index$FishID<-rownames(faith.index)
  # 
  #add the next column with the alpha diversity indexes in the metadata table. 
  alpha_div <- cbind(metadata, # add metadata
                     alpha_div,
                     #                     faith.index,
                     by = "FishID")
  #alpha_div<-alpha_div[,-c(60:61)] #remove "by" and "FishID" since there are several columns called like that
}

########################################################################################

# Loading object
library(readxl)

dv.med2<-get(load("GitHub/Phyloseq_139_dv.Rdata"))

#Filtering object ---------------------------------------------------------------------------
# Each ASV was assigned at the taxonomical level of kingdom, phylum, class, order, family and genus.
# A few filtration steps were performed on the ASV table in order to remove all the ASVs belonging to mitochondria, 
# chloroplasts and archaea and those occurring in only 1 sample. 
# Additionally, samples with less than 8000 reads were also excluded from further analyses. 


#REMOVE TAXA with ALL ZEROs 
#Delete the ASVs with colsums = 0
dv.nozero<-prune_taxa(taxa_sums(dv.med2)>0, dv.med2) 

# Filter by taxa and prevalence

### Filter by taxa: Mitocondria, chloroplasts and Archea
dv.flt.tax<-subset_taxa(dv.nozero, Family!="Mitochondria") #9162   taxa -> 20 % Euk
dv.flt.tax<-subset_taxa(dv.flt.tax, Family!="Chloroplasts") # 9162  taxa -> No Cloroplasts
dv.flt.tax<- subset_taxa(dv.flt.tax, Kingdom != "Archaea") #9040   taxa  ->  1,4 %


#CHECK-FILTER away the ASV < 2 samples

### Filter by prevalence: delete ASVs found in less than 2 samples

#Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(dv.flt.tax),
               MARGIN = ifelse(taxa_are_rows(dv.flt.tax), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

#Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(dv.flt.tax),
                    tax_table(dv.flt.tax))



#add colum with prevalence/number of samples
prevdf$Prevalence_percent<-prevdf$Prevalence / nsamples(dv.flt.tax)

# Execute prevalence filter, using prune_taxa()
keepTaxa = rownames(prevdf)[(prevdf$Prevalence > 1)]
dv.flt.tax.prev = prune_taxa(keepTaxa, dv.flt.tax) # 1743  ASVs, many are lost

#Taxonomic filtering of the ASVs unassigned at the Class level 

# Depth Filtering of samples with less than 8000 reads 

dv.meta<-sample_data(dv.flt.tax.prev)%>% data.frame()
sequencing_sums<-as.data.frame(sample_sums(dv.flt.tax.prev))%>%rownames_to_column(var="FishID")
colnames(sequencing_sums)[2]<-"Sequencing_depth"

dv.meta.seq.sum<-left_join(dv.meta,sequencing_sums,by="FishID")

#Order the samples in ps_meta from those with highest Sequencing reads to those
dv.meta.seq.sum<-dv.meta.seq.sum %>% arrange(desc(Sequencing_depth))

## Barplot - Sequencing Depth
#generate a barplot with the samples ordered from those with highest to lowest coverage
b<-barplot(dv.meta.seq.sum$Sequencing_depth, ylab = "Sequencing depth", xlab="Samples")
b+abline(h=20000,col="red")
b+abline(h=8000,col="green")

dv.flt.low.reads<-prune_samples(sample_sums(dv.flt.tax.prev) >= 8000, dv.flt.tax.prev)

# 9 samples lost, who were they? Can I affort to loose them?

#check how many reads and asv in the filtered object
sum(sample_sums(dv.flt.low.reads)) # 4621925 reads
avg.final<-mean(as.vector(sample_sums(dv.flt.low.reads))) # 35553.27 reads per sample
length(colnames(otu_table(dv.flt.low.reads))) # 1743 taxa, present in at least 2 samples 

# General information about the filtered object

# How many ASVs are unclassified at the genus level? 
tax<-as.data.frame(tax_table(dv.flt.low.reads))

tax.na.genus<-subset(tax, is.na(Genus))
tax.na.genus$Genus<-paste(tax.na.genus$Phylum,tax.na.genus$Class,
                          tax.na.genus$Order,tax.na.genus$Family,tax.na.genus$Genus,sep="_")
length(unique(tax.na.genus$Genus)) #92 unclassified genera

# How many classified genura?
tax.no.NA.genus<-subset(tax, !is.na(Genus))
tax.no.NA.genus$Genus<-paste(tax.no.NA.genus$Phylum,tax.no.NA.genus$Class,
                             tax.no.NA.genus$Order,tax.no.NA.genus$Family,tax.no.NA.genus$Genus,sep="_")
length(unique(tax.no.NA.genus$Genus)) # 435 classified genera


# Load the phyloseq object of the microbiota (GDP samples only)
# Results section 3.3 ----------------------------------------------------------------

# Codes for the article 	"The taxonomical and functional composition of the D. vulgaris gut mucosal microbiota vary spatially at large scale"

## Alpha diveristy -------------------------------------------------------------------
# Results in Supplmentary Table 3
# subset for the large scale dataset (81 samples)

dv.fsm<-subset_samples(dv.flt.low.reads,Method.of.sampling %in% c("FISHERS"))
dv.fsm<-nozero_taxa(dv.fsm)
dv.meta.fsm<-sample_data(dv.fsm)%>% data.frame()

dv.fsm.rare<-tran.rare(dv.fsm,sample.size = 8000)

# Agglomeration at higher taxonomical levels

dv.fsm.rareP<-tax_glom(dv.fsm.rare, taxrank = "Phylum", NArm = FALSE)
dv.fsm.rareC<-tax_glom(dv.fsm.rare, taxrank = "Class", NArm = FALSE)
dv.fsm.rareO<-tax_glom(dv.fsm.rare, taxrank = "Order", NArm = FALSE)
dv.fsm.rareF<-tax_glom(dv.fsm.rare, taxrank = "Family", NArm = FALSE)
dv.fsm.rareG<-tax_glom(dv.fsm.rare, taxrank = "Genus", NArm = FALSE)

# LS (large scale)
LS.P.alpha<-alpha.div(dv.fsm.rareP,dv.meta.fsm)
LS.C.alpha<-alpha.div(dv.fsm.rareC,dv.meta.fsm)
LS.O.alpha<-alpha.div(dv.fsm.rareO,dv.meta.fsm)
LS.F.alpha<-alpha.div(dv.fsm.rareF,dv.meta.fsm)
LS.G.alpha<-alpha.div(dv.fsm.rareG,dv.meta.fsm)
LS.alpha<-alpha.div(dv.fsm.rare,dv.meta.fsm)

LS.alpha<-LS.alpha[,-c(63,64)] # remove the repeteade columns in the alpha diversity table 
my_comparisons<-list(c('1','2'),c('1','3'),c('2','3'))

theme_set(theme_minimal())
comparisons = my_comparisons,(aes(label = ..p.signif..))

Shannon.plot<-ggplot(LS.alpha,aes(as.factor(MPA),Shannon, fill=as.factor(MPA)))+
  geom_boxplot(data=LS.alpha, aes(alpha=0.7))+
  geom_point(size=4,alpha=0.5, aes(color=as.factor(MPA)))+
  scale_fill_manual(values = c("pink2","purple","orange"))+
  scale_color_manual(values = c("pink2","purple","orange"))+
  stat_compare_means(method = "anova", paired=TRUE)+
  ylab("Shannon's index")+ xlab("")+
  scale_x_discrete(labels=c("1" = "BA", "2" = "CR", "3"="BO"))+
  stat_summary(fun.y="mean", color ="black", shape= 20)+
  theme(legend.key.size = unit(4, "mm"), # resize legend so it fits in the figure
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.position = "none",
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20,face="bold"),
        axis.text.x = element_text(size = 18))




# Does alpha diversity vary at a large scale?
## Shannon 
# LS
shapiro.test(LS.P.alpha$Shannon )#S
shapiro.test(LS.C.alpha$Shannon )#S
shapiro.test(LS.O.alpha$Shannon )#S
shapiro.test(LS.F.alpha$Shannon ) #S
shapiro.test(LS.G.alpha$Shannon ) # NS
shapiro.test(LS.alpha$Shannon )#NS

# Kruskal test/ANOVA on Shannon test --
kruskal.test(Shannon~as.factor(MPA),data=LS.P.alpha)# NS
kruskal.test(Shannon~as.factor(MPA),data=LS.C.alpha)# NS
kruskal.test(Shannon~as.factor(MPA),data=LS.O.alpha)# NS
kruskal.test(Shannon~as.factor(MPA),data=LS.F.alpha)# NS
summary(aov<-aov(Shannon~as.factor(MPA),data=LS.G.alpha)) #NS
summary(aov<-aov(Shannon~as.factor(MPA),data=LS.alpha)) #S

# Pairwise
TukeyHSD(aov,method="bonferroni",kw=TRUE) # None

# Length 
cor.test(LS.P.alpha$Shannon,as.numeric(LS.P.alpha$Length),method = "spearman") # NS
cor.test(LS.C.alpha$Shannon,as.numeric(LS.C.alpha$Length),method = "spearman") # NS
cor.test(LS.O.alpha$Shannon,as.numeric(LS.O.alpha$Length),method = "spearman") # NS
cor.test(LS.F.alpha$Shannon,as.numeric(LS.F.alpha$Length),method = "spearman") # NS
cor.test(LS.G.alpha$Shannon,as.numeric(LS.G.alpha$Length),method = "spearman") # NS
cor.test(LS.alpha$Shannon,as.numeric(LS.alpha$Length),method = "spearman") # NS

# All not significant. 

## Beta Diversity ----------------------------------------------------------------

dv.fsm

# Add the dates
meta.fsm<-sample_data(dv.fsm)%>% data.frame()
meta.fsm$Sampling_date<-as.Date.character(meta.fsm$Sampling_date, "%d/%m/%Y")
meta.fsm$MPA<-as.factor(meta.fsm$MPA)
sample_data(dv.fsm)<-meta.fsm

dv.fsm.clr<-tran.clr.not.glom(dv.fsm) # centered log ratio transformation

dv.atchinson.matrix<-phyloseq::distance(dv.fsm.clr,method="euclidean") #generate a distance matrix

theme_set(theme_minimal())
Regions.PCA.ordinate<-ordinate(dv.fsm.clr,"RDA")
Regions.PCA<-pcaMPA(dv.fsm.clr,MPA)



# Test

# Test variance 
permdisp(dv.fsm.clr,"MPA") # pvalue 0.543

matrix<-phyloseq::distance(dv.fsm.clr,method="euclidean") #generate a distance matrix
df<- sample_data(dv.fsm.clr) %>% data.frame() # extract the dataframe
df$MPA<-as.factor(df$MPA) 
dispr<-betadisper(matrix,df[,"MPA"]) #betadisperser test 

boxplot(dispr)

#Test - PERMANOVA
set.seed(1234)
adonis2(matrix~MPA+!is.na(Length),data = df )

#Permanova.MPA<-permanova(dv.fsm.clr,"MPA") # p.value = 0.001  ; statistics = 4.751    

resMPA<-  pairwiseAdonis::pairwise.adonis(dv.atchinson.matrix, factors=as.factor(dv.meta.fsm$MPA))

# At the Genus level the statistics are mantained? 

dv.fsm.genus<-tax_glom(dv.fsm, taxrank="Genus",NArm = TRUE)

dv.fsm.genus.clr<-tran.clr.not.glom(dv.fsm.genus)

dv.gen.atchinson.matrix<-phyloseq::distance(dv.fsm.genus.clr,method="euclidean") #generate a distance matrix

pcaMPA(dv.fsm.genus.clr,MPA)

# Test

# Test variance 
permdisp(dv.fsm.genus.clr,"MPA") # pvalue 0.971

#Test - PERMANOVA -  Does the micorbial dissimilarity vary at a large scale?
#(Results reported in Supplemntary Table 1)
#Test - PERMANOVA

set.seed(1234)
adonis2(dv.gen.atchinson.matrix~MPA+!is.na(Length),data = df )

resMPA<-  pairwiseAdonis::pairwise.adonis(dv.gen.atchinson.matrix, factors=as.factor(dv.meta.fsm$MPA))

# At the family level ?

dv.fsm.family<-tax_glom(dv.fsm, taxrank="Family",NArm = TRUE)

dv.fsm.family.clr<-tran.clr.not.glom(dv.fsm.family)

dv.fam.atchinson.matrix<-phyloseq::distance(dv.fsm.family.clr,method="euclidean") #generate a distance matrix

# Test

# Test variance 
permdisp(dv.fsm.family.clr,"MPA")

#Test - PERMANOVA
set.seed(1234)
adonis2(dv.fam.atchinson.matrix~MPA+!is.na(Length),data = df )

resMPA<-pairwiseAdonis::pairwise.adonis(dv.fam.atchinson.matrix, factors=as.factor(dv.meta.fsm$MPA))
Permanova.LengthF<-permanovaLength(dv.fsm.family.clr,"Length")

# At the Order level ?

dv.fsm.O<-tax_glom(dv.fsm, taxrank="Order",NArm = TRUE)

dv.fsm.O.clr<-tran.clr.not.glom(dv.fsm.O)

dv.O.atchinson.matrix<-phyloseq::distance(dv.fsm.O.clr,method="euclidean") #generate a distance matrix

# Test

# Test variance 
permdisp(dv.fsm.O.clr,"MPA") 

#Test - PERMANOVA
set.seed(1234)
adonis2(dv.O.atchinson.matrix~MPA+!is.na(Length),data = df )

Permanova.MPA.O<-permanova(dv.fsm.O.clr,"MPA") # p.value = 0.001  ; statistics = 3.6866        

resMPA.O<-pairwiseAdonis::pairwise.adonis(dv.O.atchinson.matrix, factors=as.factor(dv.meta.fsm$MPA))
Permanova.LengthO<-permanovaLength(dv.fsm.O.clr,"Length")


# At the Class level ?
dv.fsm.C<-tax_glom(dv.fsm, taxrank="Class",NArm = TRUE)

dv.fsm.C.clr<-tran.clr.not.glom(dv.fsm.C)

dv.C.atchinson.matrix<-phyloseq::distance(dv.fsm.C.clr,method="euclidean") #generate a distance matrix

# Test

# Test variance 
permdisp(dv.fsm.C.clr,"MPA") 

#Test - PERMANOVA
set.seed(1234)
adonis2(dv.C.atchinson.matrix~MPA+!is.na(Length),data = df )

resMPA.C<-pairwiseAdonis::pairwise.adonis(dv.C.atchinson.matrix, factors=as.factor(dv.meta.fsm$MPA))

# At the Phylum level ?

dv.fsm.P<-tax_glom(dv.fsm, taxrank="Phylum",NArm = TRUE)

dv.fsm.P.clr<-tran.clr.not.glom(dv.fsm.P)

dv.P.atchinson.matrix<-phyloseq::distance(dv.fsm.P.clr,method="euclidean") #generate a distance matrix

# Test

# Test variance 
permdisp(dv.fsm.P.clr,"MPA") 
#Test - PERMANOVA
set.seed(1234)
adonis2(dv.P.atchinson.matrix~MPA+!is.na(Length),data = df )

Permanova.MPA.p<-permanova(dv.fsm.P.clr,"MPA") # p.value = 0.002   ; statistics = 2.5514            

resMPA.P<-pairwiseAdonis::pairwise.adonis(dv.P.atchinson.matrix, factors=as.factor(dv.meta.fsm$MPA))
Permanova.LengthP<-permanovaLength(dv.fsm.P.clr,"Length")

## Core microbiota  ----------------------------------------------------------------
# Results in Supplementry Table 2
## Genus Agglomerated
# use the rarefied object to relative abundances

dv.genus<-tax_glom(dv.fsm.rare,taxrank = "Genus", NArm = FALSE)
dv.ra.gen<-transform_sample_counts(dv.genus,function(x){x / sum(x)})
sample_data(dv.ra.gen)$MPA<-as.character(sample_data(dv.ra.gen)$MPA)

#Make a list of the species
mpa.med <- unique(as.character(meta(dv.ra.gen)$MPA))

#Write a for loop to go through each of the region one by one and combine identified core taxa into a list
list_core.75.01.gen <- c() # an empty object to store information

for (n in mpa.med){ # for each variable n in Scorp.sp
  #print(paste0("Identifying Core Taxa for ", n))
  
  rare.sub <- subset_samples(dv.ra.gen, MPA == n) # Choose sample from Species by n
  
  core_m <- core_members(rare.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.0001  
                         prevalence = 0.75)
  print(paste0("No. of core genera in ", n, " : ", length(core_m))) # print core taxa identified in each Species
  list_core.75.01.gen[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

#Venn plot

p75.001.gen<-plot(venn(list_core.75.01.gen),
                  fills = MPA ,
                  main = c("Core Genera"),
                  quantities=TRUE,
                  cex.main=2)


# [1] "No. of core genera in 3 : 10"
# [1] "No. of core genera in 1 : 15"
# [1] "No. of core genera in 2 : 12"

list_core.75.01.gen.BA<-as.data.frame(list_core.75.01.gen$'1')
colnames(list_core.75.01.gen.BA)<-"ASV"


list_core.75.01.gen.CR<-as.data.frame(list_core.75.01.gen$'2')
colnames(list_core.75.01.gen.CR)<-"ASV"


list_core.75.01.gen.BO<-as.data.frame(list_core.75.01.gen$'3')
colnames(list_core.75.01.gen.BO)<-"ASV"


# combine the taxonomy
tax<-as.data.frame(tax_table(dv.fsm.rare))%>% rownames_to_column(var="ASV")
list_core.75.01.gen.BA<-inner_join(list_core.75.01.gen.BA,tax,by="ASV")
list_core.75.01.gen.CR<-inner_join(list_core.75.01.gen.CR,tax,by="ASV")
list_core.75.01.gen.BO<-inner_join(list_core.75.01.gen.BO,tax,by="ASV")

# obtain the relative abundance of all the genera 

# subset for species 
BA.gen.ra<-subset_samples(dv.ra.gen, MPA == "1")
CR.gen.ra<-subset_samples(dv.ra.gen, MPA == "2")
BO.gen.ra<-subset_samples(dv.ra.gen, MPA == "3")

# extract the asv table 
BA.gen.asvtb<-as.data.frame(t(otu_table(BA.gen.ra))) %>% rownames_to_column(var="ASV")
CR.gen.asvtb<-as.data.frame(t(otu_table(CR.gen.ra))) %>% rownames_to_column(var="ASV")
BO.gen.asvtb<-as.data.frame(t(otu_table(BO.gen.ra))) %>% rownames_to_column(var="ASV")

# join the asv table with the core taxa tables by ASV 
BA.core.gen<-inner_join(BA.gen.asvtb,list_core.75.01.gen.BA,by="ASV")
CR.core.gen<-inner_join(CR.gen.asvtb,list_core.75.01.gen.CR,by="ASV")
BO.core.gen<-inner_join(BO.gen.asvtb,list_core.75.01.gen.BO,by="ASV")

# Is the Beta diversity determined mainly by a variation in the rel abundance of the 7 shared taxa? 

dv.genus.clr<-tran.clr.not.glom(dv.genus)
otu.genus<-as.data.frame(t(otu_table(dv.genus.clr)))%>% rownames_to_column("ASV")
tax.genus<-as.data.frame(tax_table(dv.genus.clr))%>% rownames_to_column("ASV")

#subset for the core taxa

# Find the core taxa names to subset the phy object for these taxa
ba.names<-as.data.frame(BA.core.gen$ASV)
colnames(ba.names)<-"ASV"

cr.names<-as.data.frame(CR.core.gen$ASV)
colnames(cr.names)<-"ASV"


bo.names<-as.data.frame(BO.core.gen$ASV)
colnames(bo.names)<-"ASV"

names<-rbind(ba.names,cr.names,bo.names)

unique.asvs<-unique(names$ASV) #These are the 17 unique core taxa found across the 3 MPAs

# trasform the"unclassified" genera to "Family_name_unclussified"

otu.genus.core<-subset(otu.genus, ASV %in% unique.asvs)
tax.genus.core<-subset(tax.genus, ASV %in% unique.asvs)
tax.genus.core[2,7]<-"Neisseriaceae_unclassified"
tax.genus.core[6,7]<-"Vibrionaceae_unclassified"
tax.genus.core[13,7]<-"Bacillaceae_unclassified"

tax.names<-tax.genus.core$ASV

rownames(otu.genus.core)<-NULL
otu.genus.core<-otu.genus.core %>% column_to_rownames(var="ASV")
otu.genus.core<-as.matrix(t(otu.genus.core))

rownames(tax.genus.core)<-NULL
tax.genus.core<-tax.genus.core %>% column_to_rownames(var="ASV")
tax.genus.core<-as.matrix(tax.genus.core)


# Generate a phyloseq object to run the tests

dv.genus.clr.sub.core<-phyloseq(otu_table(otu.genus.core,taxa_are_rows=FALSE),
                                tax_table(tax.genus.core),taxa_names(tax.names),
                                sample_data(sample_data(dv.genus.clr)))



pcaMPA(dv.genus.clr.sub.core,MPA) # visualize the differences among the samples based only on the core taxa

# Test

# Test variance 

permdisp(dv.genus.clr.sub.core,"MPA") # pvalue 0.637

#Test - PERMANOVA

Permanova.MPA<-permanova(dv.genus.clr.sub.core,"MPA") # p.value = 0.001  ; statistics = 8.2251  

dv.atchinson.matrix.core<-phyloseq::distance(dv.genus.clr.sub.core,method="euclidean") #generate a distance matrix

resMPA.core<-  pairwiseAdonis::pairwise.adonis(dv.atchinson.matrix.core, factors=as.factor(dv.meta.fsm$MPA))

# The taxa not included in the core microbiota are also determining a variation? 

# Create an object only with the "non core" taxa
otu.genus.NOcore<-subset(otu.genus, !ASV %in% unique.asvs)

tax.genus.NOcore<-subset(tax.genus, !ASV %in% unique.asvs)

tax.namesNOcore<-tax.genus.NOcore$ASV

rownames(otu.genus.NOcore)<-NULL

otu.genus.NOcore<-otu.genus.NOcore %>% column_to_rownames(var="ASV")

otu.genus.NOcore<-as.matrix(t(otu.genus.NOcore))

rownames(tax.genus.NOcore)<-NULL

tax.genus.NOcore<-tax.genus.NOcore %>% column_to_rownames(var="ASV")

tax.genus.NOcore<-as.matrix(tax.genus.NOcore)


dv.genus.clr.sub.NOcore<-phyloseq(otu_table(otu.genus.NOcore,taxa_are_rows=FALSE),
                                  tax_table(tax.genus.NOcore),taxa_names(tax.namesNOcore),
                                  sample_data(sample_data(dv.genus.clr)))


# Test variance 

permdisp(dv.genus.clr.sub.NOcore,"MPA") # pvalue 0.637

#Test - PERMANOVA

Permanova.MPA<-permanova(dv.genus.clr.sub.NOcore,"MPA") # p.value = 0.001  ; statistics = 2.8947    

dv.atchinson.matrix.NOcore<-phyloseq::distance(dv.genus.clr.sub.NOcore,method="euclidean") #generate a distance matrix

resMPA.NOcore<-  pairwiseAdonis::pairwise.adonis(dv.atchinson.matrix.NOcore, factors=as.factor(dv.meta.fsm$MPA))

## Potential functionality --------------------------
# Load the function object (KEGG)
kegg.phy<-get(load("GitHub/kegg_phyloseq_NEW.RData"))

kegg.phy.fsm<-subset_samples(kegg.phy,Method.of.sampling %in% c("FISHERS"))
kegg.fsm.meta<-sample_data(kegg.phy.fsm)%>% data.frame()

kegg.hell<-kegg.phy.fsm
otu_table(kegg.hell) <- otu_table(decostand(otu_table(kegg.hell), method = "hellinger"), taxa_are_rows=FALSE)

# Subset Function 

kegg.bray.tot<-phyloseq::distance(kegg.hell,method="bray") #generate a distance matrix

kegg.sub1<-subset_taxa(kegg.hell, pathways_Level3 %in% c(" Carbohydrate metabolism"))
kegg.bray1<-phyloseq::distance(kegg.sub1,method="bray") #generate a distance matrix

kegg.sub2<-subset_taxa(kegg.hell, pathways_Level3 %in% c(" Glycan biosynthesis and metabolism"))
kegg.bray2<-phyloseq::distance(kegg.sub2,method="bray") #generate a distance matrix

kegg.sub3<-subset_taxa(kegg.hell, pathways_Level3 %in% c(" Amino acid metabolism"))
kegg.bray3<-phyloseq::distance(kegg.sub3,method="bray") #generate a distance matrix

kegg.sub4<-subset_taxa(kegg.hell, pathways_Level3 %in% c(" Lipid metabolism"))
kegg.bray4<-phyloseq::distance(kegg.sub4,method="bray") #generate a distance matrix

kegg.sub5<-subset_taxa(kegg.hell, pathways_Level3 %in% c(" Energy metabolism"))
kegg.bray5<-phyloseq::distance(kegg.sub5,method="bray") #generate a distance matrix

kegg.sub6<-subset_taxa(kegg.hell, pathways_Level3 %in% c(" Metabolism of cofactors and vitamins"))
kegg.bray6<-phyloseq::distance(kegg.sub6,method="bray") #generate a distance matrix

kegg.sub7<-subset_taxa(kegg.hell, pathways_Level3 %in% c(" Metabolism of terpenoids and polyketides"))
kegg.bray7<-phyloseq::distance(kegg.sub7,method="bray") #generate a distance matrix

kegg.sub8<-subset_taxa(kegg.hell, pathways_Level3 %in% c(" Xenobiotics biodegradation and metabolism"))
kegg.bray8<-phyloseq::distance(kegg.sub8,method="bray") #generate a distance matrix

# Alpha diversity funcitons 

alpha.sub1<-alpha.divKegg(kegg.sub1,kegg.fsm.meta)
alpha.sub2<-alpha.divKegg(kegg.sub2,kegg.fsm.meta)
alpha.sub3<-alpha.divKegg(kegg.sub3,kegg.fsm.meta)
alpha.sub4<-alpha.divKegg(kegg.sub4,kegg.fsm.meta)
alpha.sub5<-alpha.divKegg(kegg.sub5,kegg.fsm.meta)
alpha.sub6<-alpha.divKegg(kegg.sub6,kegg.fsm.meta)
alpha.sub7<-alpha.divKegg(kegg.sub7,kegg.fsm.meta)
alpha.sub8<-alpha.divKegg(kegg.sub8,kegg.fsm.meta)


#Shapiro 
shapiro.test(alpha.sub1$Shannon )#S
shapiro.test(alpha.sub2$Shannon )#S
shapiro.test(alpha.sub3$Shannon )#S
shapiro.test(alpha.sub4$Shannon )#S
shapiro.test(alpha.sub5$Shannon )#S
shapiro.test(alpha.sub6$Shannon )#S
shapiro.test(alpha.sub7$Shannon )#S
shapiro.test(alpha.sub8$Shannon )#S

# Kruskal test/ANOVA on Shannon test --
kruskal.test(Shannon~as.factor(MPA),data=alpha.sub1)# NS
kruskal.test(Shannon~as.factor(MPA),data=alpha.sub2)# NS
kruskal.test(Shannon~as.factor(MPA),data=alpha.sub3)# S
kruskal.test(Shannon~as.factor(MPA),data=alpha.sub4)# NS
kruskal.test(Shannon~as.factor(MPA),data=alpha.sub5)# NS
kruskal.test(Shannon~as.factor(MPA),data=alpha.sub6)# NS
kruskal.test(Shannon~as.factor(MPA),data=alpha.sub7)# NS
kruskal.test(Shannon~as.factor(MPA),data=alpha.sub8)# NS

dunn.test::dunn.test(alpha.sub3$Shannon,as.factor(alpha.sub3$MPA), method = "bonferroni") #CR vs the other two

boxplot(Shannon~as.factor(MPA), data=alpha.sub3) #CR had the lowest diversity

# Beta diversity Functions 

## PERMANOVA on the FUCNTIONS 
# Permanova or Welch MANOVA on the functions by the Region 
set.seed(123)
disp.tot<-betadisper(kegg.bray.tot,as.factor(kegg.fsm.meta$MPA)) #betadisperser test
permutest(disp.tot) # Not sign
adonis2(kegg.bray.tot~as.factor(kegg.fsm.meta$MPA)) #  Significant
pairwise.adonis2(kegg.bray.tot~as.factor(MPA),data=kegg.fsm.meta) # All but 3vs 1 

pca.phyHell <- ordinate(kegg.hell, "RDA") #Performs redundancy analysis, or optionally principal components analysis, via rda

vegan::protest(Regions.PCA.ordinate,pca.phyHell, symmetric =TRUE )
vegan::procrustes (Regions.PCA.ordinate,pca.phyHell, symmetric =TRUE )

pcaMPAKegg(kegg.hell)


# Subsets
set.seed(123)
as.data.frame(tax_table(kegg.sub1))$pathways_Level3
disp.1<-betadisper(kegg.bray1,as.factor(kegg.fsm.meta$MPA)) #betadisperser test
permutest(disp.1) # Not sign
adonis2(kegg.bray1~as.factor(kegg.fsm.meta$MPA)) #  Significant
pairwise.adonis2(kegg.bray1~as.factor(MPA),data=kegg.fsm.meta) # All but 3vs 1 



set.seed(123)
as.data.frame(tax_table(kegg.sub2))$pathways_Level3
disp.2<-betadisper(kegg.bray2,as.factor(kegg.fsm.meta$MPA)) #betadisperser test
permutest(disp.2) # Not sign
adonis2(kegg.bray2~as.factor(kegg.fsm.meta$MPA)) #  Not Significant

set.seed(123)
as.data.frame(tax_table(kegg.sub3))$pathways_Level3
disp.3<-betadisper(kegg.bray3,as.factor(kegg.fsm.meta$MPA)) #betadisperser test
permutest(disp.3) # Not sign
adonis2(kegg.bray3~as.factor(kegg.fsm.meta$MPA)) #  Significant
pairwise.adonis2(kegg.bray3~as.factor(MPA),data=kegg.fsm.meta) #All but 3 vs 1

set.seed(123)
as.data.frame(tax_table(kegg.sub4))$pathways_Level3
disp.4<-betadisper(kegg.bray4,as.factor(kegg.fsm.meta$MPA)) #betadisperser test
permutest(disp.4) # Not sign
adonis2(kegg.bray4~as.factor(kegg.fsm.meta$MPA)) #  Significant
pairwise.adonis2(kegg.bray4~as.factor(MPA),data=kegg.fsm.meta) #All but 3 vs 1

set.seed(123)
as.data.frame(tax_table(kegg.sub5))$pathways_Level3
disp.5<-betadisper(kegg.bray5,as.factor(kegg.fsm.meta$MPA)) #betadisperser test
permutest(disp.5) # Not sign
adonis2(kegg.bray5~as.factor(kegg.fsm.meta$MPA)) # Not  Significant

set.seed(123)
as.data.frame(tax_table(kegg.sub6))$pathways_Level3
disp.6<-betadisper(kegg.bray6,as.factor(kegg.fsm.meta$MPA)) #betadisperser test
permutest(disp.6) # Not sign
adonis2(kegg.bray6~as.factor(kegg.fsm.meta$MPA)) #  Significant
pairwise.adonis2(kegg.bray6~as.factor(MPA),data=kegg.fsm.meta) #All but 3 vs 1 

set.seed(123)
as.data.frame(tax_table(kegg.sub7))$pathways_Level3
disp.7<-betadisper(kegg.bray7,as.factor(kegg.fsm.meta$MPA)) #betadisperser test
permutest(disp.7) # Not sign
adonis2(kegg.bray7~as.factor(kegg.fsm.meta$MPA)) #  Significant
pairwise.adonis2(kegg.bray7~as.factor(MPA),data=kegg.fsm.meta) #All but 3 vs 1

set.seed(123)
as.data.frame(tax_table(kegg.sub8))$pathways_Level3
disp.8<-betadisper(kegg.bray8,as.factor(kegg.fsm.meta$MPA)) #betadisperser test
permutest(disp.8) # Not sign
w<-MicEco::WdS.test (dm=kegg.bray8, f=as.factor(kegg.fsm.meta$MPA),nrep=999) #  Significant
pairwise.adonis2(kegg.bray8~as.factor(MPA),data=kegg.fsm.meta) # only 3 vs 2 

pcaMPAKegg(kegg.sub1)
pcaMPAKegg(kegg.sub2)
pcaMPAKegg(kegg.sub3)
pcaMPAKegg(kegg.sub4)
pcaMPAKegg(kegg.sub5)
pcaMPAKegg(kegg.sub6)
pcaMPAKegg(kegg.sub7)
pcaMPAKegg(kegg.sub8)


# Result section 3.4 ---------------------------------------------------------------
# Result section 3.4.1--------------------------------------------------------------
# "The gut microbiota of Diplodus vulgaris vary taxonomically and functionally across a small spatial range"

# palette for plots 
#Palette
gdp<-c("turquoise3","firebrick2","purple","khaki","magenta1","ivory3","grey9")

# Objects needed

dv.flt.low.reads # micorbiota phyloseq object

# Subset fo the samples from the small dataset

dv.BO.gdp<-subset_samples(dv.flt.low.reads,Dataset %in% c("Small-scale dataset"))
dv.BO.gdp<-nozero_taxa(dv.BO.gdp)

# remove the GDP-point with only one observation
micro.meta.full<-sample_data(dv.BO.gdp) %>% data.frame()
micro.meta.full%>% group_by(GDP_point)%>%summarise(freq=n()) # GDP-3 has only one observation ,so it is exluded 
dv.BO.gdp.noGDP3<-subset_samples(dv.BO.gdp, !GDP_point %in% "GDP-3")

## Alpha diversity -----------------------------------------------------------------

#Rarefaction
micro.rare<-dv.BO.gdp.noGDP3
micro.rare<-tran.rare(micro.rare,sample.size = 8000)
micro.meta<-sample_data(micro.rare)%>% data.frame()

# Agglomerate at higher tax levels 

# Micro 
micro.rare.G<-tax_glom(micro.rare,taxrank = "Genus",NArm = FALSE)
micro.rare.F<-tax_glom(micro.rare,taxrank = "Family",NArm = FALSE)
micro.rare.O<-tax_glom(micro.rare,taxrank = "Order",NArm = FALSE)
micro.rare.C<-tax_glom(micro.rare,taxrank = "Class",NArm = FALSE)
micro.rare.P<-tax_glom(micro.rare,taxrank = "Phylum",NArm = FALSE)

micro.gdp.alphaP<-alpha.div(micro.rare.P,micro.meta)
micro.gdp.alphaC<-alpha.div(micro.rare.C,micro.meta)
micro.gdp.alphaO<-alpha.div(micro.rare.O,micro.meta)
micro.gdp.alphaF<-alpha.div(micro.rare.F,micro.meta)
micro.gdp.alphaG<-alpha.div(micro.rare.G,micro.meta)
micro.gdp.alpha<-alpha.div(micro.rare,micro.meta)

# mean alpha diversity?
mean(micro.gdp.alpha$Shannon) # 2.99
sd(micro.gdp.alpha$Shannon) # 0.7 

mean(micro.gdp.alpha$Observed) # 82.7
sd(micro.gdp.alpha$Observed) # 32.7


shapiro.test(micro.gdp.alphaP$Shannon) # Not Sign
shapiro.test(micro.gdp.alphaC$Shannon) # Sign
shapiro.test(micro.gdp.alphaO$Shannon) # Sign
shapiro.test(micro.gdp.alphaF$Shannon) # Sign
shapiro.test(micro.gdp.alphaG$Shannon) # Sign
shapiro.test(micro.gdp.alpha$Shannon) # Sign


## Does alpha diversity vary at a small scale? 

summary(aov(Shannon~as.factor(GDP_point),data=micro.gdp.alphaP)) # NS
kruskal.test(Shannon~as.factor(GDP_point),data=micro.gdp.alphaC) # NS
kruskal.test(Shannon~as.factor(GDP_point),data=micro.gdp.alphaO) # Sign
kruskal.test(Shannon~as.factor(GDP_point),data=micro.gdp.alphaF) # Sign
kruskal.test(Shannon~as.factor(GDP_point),data=micro.gdp.alphaG) # Sign
kruskal.test(Shannon~as.factor(GDP_point),data=micro.gdp.alpha) # Sign

# Pairwise
set.seed(1234)
dunn.test::dunn.test(micro.gdp.alphaO$Shannon,as.factor(micro.gdp.alphaO$GDP_point),method="bonferroni",kw=TRUE)
dunn.test::dunn.test(micro.gdp.alphaF$Shannon,as.factor(micro.gdp.alphaF$GDP_point),method="bonferroni",kw=TRUE)
dunn.test::dunn.test(micro.gdp.alphaG$Shannon,as.factor(micro.gdp.alphaG$GDP_point),method="bonferroni",kw=TRUE)
dunn.test::dunn.test(micro.gdp.alpha$Shannon,as.factor(micro.gdp.alpha$GDP_point),method="bonferroni",kw=TRUE)

micro.gdp.alpha<-micro.gdp.alpha[,-64]

# Plot 
ShannonGDP<-ggplot(micro.gdp.alpha, aes(x=factor(GDP_point,levels=c("GDP-1",
                                                                    "GDP-4",
                                                                    "GDP-5",
                                                                    "GDP-6",
                                                                    "GDP-7.2",
                                                                    "GDP-8",
                                                                    "GDP-12")),y=Shannon, fill=GDP_point))+
  geom_boxplot()+geom_point(size=3,alpha=0.2)+
  scale_fill_manual(values=gdp)+
  ylab("Shannon's index of richness and eveness")+ xlab("")




## Beta diversity -------------------------------------------------------------------
# Agglomerate at different taxonomic levels 
dv.BO.gdp.noGDP3.P<-tax_glom(dv.BO.gdp.noGDP3, taxrank = "Phylum",NArm = FALSE)
dv.BO.gdp.noGDP3.C<-tax_glom(dv.BO.gdp.noGDP3, taxrank = "Class",NArm = FALSE)
dv.BO.gdp.noGDP3.O<-tax_glom(dv.BO.gdp.noGDP3, taxrank = "Order",NArm = FALSE)
dv.BO.gdp.noGDP3.F<-tax_glom(dv.BO.gdp.noGDP3, taxrank = "Family",NArm = FALSE)
dv.BO.gdp.noGDP3.G<-tax_glom(dv.BO.gdp.noGDP3, taxrank = "Genus",NArm = FALSE)


#Transform
dv.BO.gdp.noGDP3.clrP<-tran.clr.not.glom(dv.BO.gdp.noGDP3.P)
dv.BO.gdp.noGDP3.clrC<-tran.clr.not.glom(dv.BO.gdp.noGDP3.C)
dv.BO.gdp.noGDP3.clrO<-tran.clr.not.glom(dv.BO.gdp.noGDP3.O)
dv.BO.gdp.noGDP3.clrF<-tran.clr.not.glom(dv.BO.gdp.noGDP3.F)
dv.BO.gdp.noGDP3.clrG<-tran.clr.not.glom(dv.BO.gdp.noGDP3.G)
dv.BO.gdp.noGDP3.clr<-tran.clr.not.glom(dv.BO.gdp.noGDP3)

# Calculate distance

noGDP3.distP<-phyloseq::distance(dv.BO.gdp.noGDP3.clrP,method="euclidean") #generate a distance matrix
noGDP3.distC<-phyloseq::distance(dv.BO.gdp.noGDP3.clrC,method="euclidean") #generate a distance matrix
noGDP3.distO<-phyloseq::distance(dv.BO.gdp.noGDP3.clrO,method="euclidean") #generate a distance matrix
noGDP3.distF<-phyloseq::distance(dv.BO.gdp.noGDP3.clrF,method="euclidean") #generate a distance matrix
noGDP3.distG<-phyloseq::distance(dv.BO.gdp.noGDP3.clrG,method="euclidean") #generate a distance matrix
noGDP3.dist<-phyloseq::distance(dv.BO.gdp.noGDP3.clr,method="euclidean") #generate a distance matrix


# PCA of the microbial dissimilarity
#Performs redundancy analysis, or optionally principal components analysis, via rda
pca.phyl<-ordinate(dv.BO.gdp.noGDP3.clrP, "RDA")
pca.class<-ordinate(dv.BO.gdp.noGDP3.clrC, "RDA")
pca.ord<-ordinate(dv.BO.gdp.noGDP3.clrO, "RDA")
pca.fam<-ordinate(dv.BO.gdp.noGDP3.clrF, "RDA")
pca.gen<-ordinate(dv.BO.gdp.noGDP3.clrG, "RDA")
pca.asv<-ordinate(dv.BO.gdp.noGDP3.clr, "RDA")

plot.pca<-function(clr,pca){plot_ordination(clr,pca, color="GDP_point")+
    geom_point(size=4, alpha=1, aes(color=as.character(GDP_point) ))+
    scale_color_manual(values = gdp)
}

pcaP<-plot.pca(dv.BO.gdp.noGDP3.clrP,pca.phyl)+
  ggtitle("PCA of microbial dissimilarity - Phylum agglomerated")
pcaC<-plot.pca(dv.BO.gdp.noGDP3.clrC,pca.class)+
  ggtitle("PCA of microbial dissimilarity - Class agglomerated")
pcaO<-plot.pca(dv.BO.gdp.noGDP3.clrO,pca.ord)+
  ggtitle("PCA of microbial dissimilarity - Order agglomerated")
pcaF<-plot.pca(dv.BO.gdp.noGDP3.clrF,pca.fam)+
  ggtitle("PCA of microbial dissimilarity - Family agglomerated")
pcaG<-plot.pca(dv.BO.gdp.noGDP3.clrG,pca.gen)+
  ggtitle("PCA of microbial dissimilarity - Genus agglomerated")
pcasp<-plot.pca(dv.BO.gdp.noGDP3.clr,pca.asv)+
  ggtitle("Taxonomic dissimilarity of the gut microbiota", "Aitchinson's distances of CLR transformed data")


ggarrange(pcaP,pcaC,
          pcaO,pcaF,
          pcaG,pcasp,nrow=3,ncol=2, common.legend = TRUE)

# Test - Does the micorbial dissimilarity vary at a small scale?
## PERMDISP and PERMANOVA (Results reported in Supplemntary Table 1)
set.seed(1234)
disprP<-betadisper(noGDP3.distP,as.factor(micro.meta$GDP_point)) #betadisperser test 
permutest(disprP) #  Not sign
adonis2(noGDP3.distP~as.factor(micro.meta$GDP_point)+micro.meta$Length) # Not Sign 

set.seed(1234)
disprC<-betadisper(noGDP3.distC,as.factor(micro.meta$GDP_point)) #betadisperser test 
permutest(disprC) #  Not sign
adonis2(noGDP3.distC~as.factor(micro.meta$GDP_point)+micro.meta$Length) 


set.seed(1234)
disprO<-betadisper(noGDP3.distO,as.factor(micro.meta$GDP_point)) #betadisperser test 
permutest(disprO) #  Not sign
adonis2(noGDP3.distO~as.factor(micro.meta$GDP_point)+micro.meta$Length) # Sign 
adonis2(noGDP3.distO~micro.meta$Length) # Not Sign 

set.seed(1234)
disprF<-betadisper(noGDP3.distF,as.factor(micro.meta$GDP_point)) #betadisperser test 
permutest(disprF) #  Not sign
adonis2(noGDP3.distF~as.factor(micro.meta$GDP_point)+micro.meta$Length) # Sign 

set.seed(1234)
disprG<-betadisper(noGDP3.distG,as.factor(micro.meta$GDP_point)) #betadisperser test 
permutest(disprG) #  Not sign
adonis2(noGDP3.distG~as.factor(micro.meta$GDP_point)+micro.meta$Length) # Sign 

set.seed(1234)
disprS<-betadisper(noGDP3.dist,as.factor(micro.meta$GDP_point)) #betadisperser test 
permutest(disprS) #  Not sign
adonis2(noGDP3.dist~as.factor(micro.meta$GDP_point)+micro.meta$Length) # Not Sign 


# GDP location is significant at every taxonomical level except Phylum level
# Total Length of the fish is never significantly correlated with the dissimilarity of the micorbiota community

## PAIRWISE ADONIS (Results in Supplmetary Table 6)
set.seed(1234)
# pairwiseP<-pairwiseAdonis::pairwise.adonis2(noGDP3.distP~GDP_point, data=micro.meta.bentos) # do not run it at phylum level because PERMANOVA was not significant
pairwiseC<-pairwiseAdonis::pairwise.adonis2(noGDP3.distC~GDP_point, data=micro.meta.bentos)
pairwiseO<-pairwiseAdonis::pairwise.adonis2(noGDP3.distO~GDP_point, data=micro.meta.bentos)
pairwiseF<-pairwiseAdonis::pairwise.adonis2(noGDP3.distF~GDP_point, data=micro.meta.bentos)
pairwiseG<-pairwiseAdonis::pairwise.adonis2(noGDP3.distG~GDP_point, data=micro.meta.bentos)
pairwiseS<-pairwiseAdonis::pairwise.adonis2(noGDP3.dist~GDP_point, data=micro.meta.bentos)


## Potential funcitonality ---------------------------------------------------------
### Generate functions ##################################################################################"
betadisp_kegg<-function(kegg.sub){
  # Transform hellinger 
  kegg.hell<-kegg.sub
  otu_table(kegg.hell) <- otu_table(decostand(otu_table(kegg.hell), method = "hellinger"), taxa_are_rows=FALSE)
  kegg.bray<-phyloseq::distance(kegg.hell,method="bray") #generate a distance matrix
  meta.kegg<-sample_data(kegg.sub)%>% data.frame()
  # Permanova or Welch MANOVA on the functions by the GDP point $
  set.seed(123)
  disp.kegg<-betadisper(kegg.bray,meta.kegg$GDP_point) #betadisperser test 
  permutest(disp.kegg)
}


distance.matrix<-function(kegg.sub){
  # Transform hellinger 
  kegg.hell<-kegg.sub
  otu_table(kegg.hell) <- otu_table(decostand(otu_table(kegg.hell), method = "hellinger"), taxa_are_rows=FALSE)
  kegg.bray<-phyloseq::distance(kegg.hell,method="bray") #generate a distance matrix
}

permanova_kegg<-function(kegg.sub){
  # Transform hellinger 
  kegg.hell<-kegg.sub
  otu_table(kegg.hell) <- otu_table(decostand(otu_table(kegg.hell), method = "hellinger"), taxa_are_rows=FALSE)
  kegg.bray<-phyloseq::distance(kegg.hell,method="bray") #generate a distance matrix
  meta.kegg<-sample_data(kegg.sub)%>% data.frame()
  # Permanova or Welch MANOVA on the functions by the GDP point $
  set.seed(123)
  adonis2(kegg.bray~meta.kegg$GDP_point)
}

welch_kegg<-function(kegg.sub){
  # Transform hellinger 
  kegg.hell<-kegg.sub
  otu_table(kegg.hell) <- otu_table(decostand(otu_table(kegg.hell), method = "hellinger"), taxa_are_rows=FALSE)
  kegg.bray<-phyloseq::distance(kegg.hell,method="bray") #generate a distance matrix
  meta.kegg<-sample_data(kegg.sub)%>% data.frame()
  # Permanova or Welch MANOVA on the functions by the GDP point $
  set.seed(123)
  wdTest<-MicEco::WdS.test(dm = kegg.bray, f = as.factor(meta.kegg$GDP_point), nrep =999) # F 1.73, P 0.07
}

# test effect of distance from FPA on the funcitonality
permanova_keggDIST<-function(kegg.sub){
  # Transform hellinger 
  kegg.hell<-kegg.sub
  otu_table(kegg.hell) <- otu_table(decostand(otu_table(kegg.hell), method = "hellinger"), taxa_are_rows=FALSE)
  kegg.bray<-phyloseq::distance(kegg.hell,method="bray") #generate a distance matrix
  meta.kegg<-sample_data(kegg.sub)%>% data.frame()
  # Permanova or Welch MANOVA on the functions by the GDP point $
  set.seed(123)
  adonis2(kegg.bray~meta.kegg$distance_closest_fpa)
}

welch_keggDIST<-function(kegg.sub){
  # Transform hellinger 
  kegg.hell<-kegg.sub
  otu_table(kegg.hell) <- otu_table(decostand(otu_table(kegg.hell), method = "hellinger"), taxa_are_rows=FALSE)
  kegg.bray<-phyloseq::distance(kegg.hell,method="bray") #generate a distance matrix
  meta.kegg<-sample_data(kegg.sub)%>% data.frame()
  # Permanova or Welch MANOVA on the functions by the GDP point $
  set.seed(123)
  wdTest<-MicEco::WdS.test(dm = kegg.bray, f = as.factor(meta.kegg$distance_closest_fpa), nrep =999) # F 1.73, P 0.07
}
####################################################################################
# Subset the Kegg object for the samples from Bonifacio

kegg.gdp<-subset_samples(kegg.phy, Dataset %in% c("Small-scale dataset"))
kegg.gdp.no3<-subset_samples(kegg.gdp, !GDP_point %in% c("GDP-3"))

kegg.gdp<-nozero_taxaNOtree(kegg.gdp.no3)
meta.kegg.gdp<-sample_data(kegg.gdp)%>% data.frame()

### PERMANOVA 
set.seed(1234)
kegg.sub1<-subset_taxa(kegg.gdp, pathways_Level3 %in% c(" Carbohydrate metabolism"))
betadisp_kegg(kegg.sub1)
permanova_kegg(kegg.sub1) #NS
permanova_keggDIST(kegg.sub1)

set.seed(1234)
kegg.sub2<-subset_taxa(kegg.gdp, pathways_Level3 %in% c(" Glycan biosynthesis and metabolism"))
betadisp_kegg(kegg.sub2)
permanova_kegg(kegg.sub2) #NS
permanova_keggDIST(kegg.sub2)

set.seed(1234)
kegg.sub3<-subset_taxa(kegg.gdp, pathways_Level3 %in% c(" Amino acid metabolism"))
betadisp_kegg(kegg.sub3)
permanova_kegg(kegg.sub3) #NS
permanova_keggDIST(kegg.sub3) #NS

set.seed(1234)
kegg.sub4<-subset_taxa(kegg.gdp, pathways_Level3 %in% c(" Lipid metabolism"))
betadisp_kegg(kegg.sub4)
permanova_kegg(kegg.sub4) #NS
permanova_keggDIST(kegg.sub4) #NS

set.seed(1234)
kegg.sub5<-subset_taxa(kegg.gdp, pathways_Level3 %in% c(" Energy metabolism"))
betadisp_kegg(kegg.sub5)
permanova_kegg(kegg.sub5) #NS
w5.dist<-permanova_keggDIST(kegg.sub5) #NS

set.seed(1234)
kegg.sub6<-subset_taxa(kegg.gdp, pathways_Level3 %in% c(" Metabolism of cofactors and vitamins"))
betadisp_kegg(kegg.sub6)
permanova_kegg(kegg.sub6) #NS
permanova_keggDIST(kegg.sub6) #NS

set.seed(1234)
kegg.sub7<-subset_taxa(kegg.gdp, pathways_Level3 %in% c(" Metabolism of terpenoids and polyketides"))
betadisp_kegg(kegg.sub7)
permanova_kegg(kegg.sub7) #S
permanova_keggDIST(kegg.sub7) #NS

set.seed(1234)
kegg.sub8<-subset_taxa(kegg.gdp, pathways_Level3 %in% c(" Xenobiotics biodegradation and metabolism"))
betadisp_kegg(kegg.sub8)
permanova_kegg(kegg.sub8) #S
permanova_keggDIST(kegg.sub8) #NS

# Only  Metabolism of terpenoids and polyketides and Xenobiotics biodegradation and metabolism is significantly different between the GDP points
# So focus on these two metabolisms for defining which specific funcitons are enriched or underrepresented in the different GDP zones
# Estimate differently abundant functions 

# Subset only the metabolism you are interested in
kegg.hell<-kegg.gdp
otu_table(kegg.hell) <- otu_table(decostand(otu_table(kegg.hell), method = "hellinger"), taxa_are_rows=FALSE)

kegg.hell.sub<-subset_taxa(kegg.hell,pathways_Level3 %in% c(" Metabolism of terpenoids and polyketides",
                                                            " Xenobiotics biodegradation and metabolism"))
otu.K.sub<-as.data.frame(t(otu_table(kegg.hell.sub)))
otu.K.sub$ASV<-rownames(otu.K.sub)

tax.K.sub<-as.data.frame(tax_table(kegg.hell.sub))
tax.K.sub$ASV<-rownames(tax.K.sub)
tax.K.sub<-tax.K.fsm.sub[order(tax.K.sub$pathways_Level3),]
tax.K.sub2<-tax.K.sub[,c("ASV","pathway_map")]

otu.tax.kegg<-dplyr::left_join(tax.K.sub2,otu.K.sub,by="ASV")
otu.tax.kegg<-otu.tax.kegg[,-1]
otu.tax.kegg<-otu.tax.kegg %>% column_to_rownames(var="pathway_map")
otu.tax.kegg<-as.data.frame(t(otu.tax.kegg))
otu.tax.kegg$FishID<-rownames(otu.tax.kegg)

kegg.meta<-sample_data(kegg.hell.sub)%>% data.frame()
kegg.meta<-kegg.meta[,c("FishID","GDP_point")]

otu.tax.meta.kegg<-dplyr::left_join(kegg.meta,otu.tax.kegg,by="FishID")
otu.tax.meta.kegg$GDP_point<-as.factor(otu.tax.meta.kegg$GDP_point)

### Loops 
# generate a loop to have the abundance for each sample of the ko pathway in order to run the KW test afterwards 
# for the loop 
colnames(otu.tax.meta.kegg)
otu.tax.meta.kegg.noASV<-otu.tax.meta.kegg[,c(3:33)]

df<-data.frame(0,0,0,0)
colnames(df)<-c("FishID"    ,    "GDP_point"        ,   "hell_abundance", "pathway_map" )

for (i in 1:ncol(otu.tax.meta.kegg.noASV)) {
  sub<-otu.tax.meta.kegg[,c(1,2,i)]
  sub$pathway_map<-colnames(sub)[3]
  colnames(sub)[3]<-"hell_abundance"
  df<-rbind(df,sub)
}

# The first rows are only the names of the column (error of the loop)
df.kegg.hell.abundances<-df
df.kegg.hell.abundances<-df.kegg.hell.abundances[-c(1:101),]
df.kegg.hell.abundances$hell_abundance<-as.numeric(df.kegg.hell.abundances$hell_abundance)
df.kegg.hell.abundances$GDP_point<-as.character(df.kegg.hell.abundances$GDP_point)

# Two pathways were excluded, however I recover them

pathways<-unique(df.kegg.hell.abundances$pathway_map) #29 --> why?
excl<-tax.K.sub2 %>% subset(!pathway_map %in% pathways)
excl1<-otu.tax.meta.kegg.noASV$`Biosynthesis of type II polyketide backbone`
excl2<-otu.tax.meta.kegg.noASV$`Biosynthesis of type II polyketide products`

sub.exc1<-otu.tax.meta.kegg[,c("FishID","GDP_point","Biosynthesis of type II polyketide backbone")]
sub.exc1$pathway_map<-colnames(sub.exc1)[3]
colnames(sub.exc1)[3]<-"hell_abundance"

sub.exc2<-otu.tax.meta.kegg[,c("FishID","GDP_point","Biosynthesis of type II polyketide products")]
sub.exc2$pathway_map<-colnames(sub.exc2)[3]
colnames(sub.exc2)[3]<-"hell_abundance"

excluded.recovered<-rbind(sub.exc1,sub.exc2)

excluded.recovered$hell_abundance<-as.numeric(excluded.recovered$hell_abundance)
excluded.recovered$GDP_point<-as.character(excluded.recovered$GDP_point)

# COmbine with the rest of the keggs

df.kegg.hell.abundances<-rbind(df.kegg.hell.abundances,excluded.recovered)

# All the functions are recovered

# Generate a loop to run shapiro wilks on every ko 
#Shapiro
pathways.recovered<-unique(df.kegg.hell.abundances$pathway_map) #31

data.path.s<-data.frame()

for (i in pathways.recovered) {
  subset.k<-subset(df.kegg.hell.abundances, pathway_map == i)
  s<-shapiro.test(subset.k$hell_abundance)
  data.s<-data.frame(s$p.value)
  data.path.s<-rbind(data.path.s,data.s)
}

rownames(data.path.s)<-pathways.recovered

data.path.s.NotSign<-subset(data.path.s, s.p.value > 0.05)

for.aov<-rownames(data.path.s.NotSign)

# Run anova on those that have not sign KW
subset1<-subset(df.kegg.hell.abundances, pathway_map == "Nitrotoluene degradation"    )
summary(aov(hell_abundance~GDP_point,data=subset1)) # NS

subset2<-subset(df.kegg.hell.abundances, pathway_map == "Zeatin biosynthesis"   )
summary(aov(hell_abundance~GDP_point,data=subset2)) # NS

subset3<-subset(df.kegg.hell.abundances, pathway_map == "Biosynthesis of ansamycins"   )
summary(aov(hell_abundance~GDP_point,data=subset3)) #NS

# Nitrotoluene degradation 1.463   0.214
# Zeatin biosynthesis 1.867  0.109
# Biosynthesis of ansamycins 0.725  0.632

### Kruskal wallis on every pathway 

pathways.recovered.KW<-unique(df.kegg.hell.abundances$pathway_map) #31
pathways.recovered.KW<-pathways.recovered.KW[!pathways.recovered.KW %in% for.aov]

data.path.k<-data.frame()

for (i in pathways.recovered.KW) {
  subset.k<-subset(df.kegg.hell.abundances, pathway_map == i)
  k<-kruskal.test(subset.k$hell_abundance,subset.k$GDP_point)
  data.k<-data.frame(k$p.value,k$statistic)
  data.path.k<-rbind(data.path.k,data.k)
}

rownames(data.path.k)<-pathways.recovered.KW

data.path.k.sign<-subset(data.path.k, k.p.value < 0.05)
dim(data.path.k.sign) # 5 pathways were significantly different  + 5 sign different by ANOVA 

# Pairwise difference
subset.k.signif<-subset(df.kegg.hell.abundances, pathway_map == rownames(data.path.k.sign))
dunnTest<-function(sub){
  dunn.test::dunn.test(sub$hell_abundance,sub$GDP_point, method = "bonferroni")
}

sub.KW1<-subset(df.kegg.hell.abundances, pathway_map =="Tetracycline biosynthesis" )
sub.KW2<-subset(df.kegg.hell.abundances, pathway_map =="Dioxin degradation" )
sub.KW3<-subset(df.kegg.hell.abundances, pathway_map == "Aminobenzoate degradation")
sub.KW4<-subset(df.kegg.hell.abundances, pathway_map =="Carotenoid biosynthesis" )
sub.KW5<-subset(df.kegg.hell.abundances, pathway_map =="Caprolactam degradation" )

dunnTest(sub.KW1) # -> GDP-1 vs GDP-12  and GDP-6
dunnTest(sub.KW2) # None
dunnTest(sub.KW3) # GDP-1 vs GDP-8
dunnTest(sub.KW4) # GDP-1 vs GDP-6 and vs GDP-4
dunnTest(sub.KW5) # GDP-7.2 vs GDP-4 and vs GDP-6 


# Result secrtion 3.4.2 -------------------------------------------------------------
#"The type of benthic habitat and the distance from fully protected areas are predictors of gut mucosal microbiota features in Diplodus vulgaris"

# DBRDA with the different types of benthic substrata and distance from the FPA as covariates (Supplementary Table 8)
# Load the information about the benthic substrata 
benthic.habitat<-read.csv("benthic_habitat.csv")

# Add the information about the habitat type (KUD 95-double radius)
micro.meta.bentos<-dplyr::left_join(micro.meta,benthic.habitat,by="GDP_point")

## benthic habitat and distance FPA - DBRDA ---------------------------------------

set.seed(1234)
LOCUSdbrda.bentosP<-dbrda(noGDP3.distP~LOCUS_perc_Mosaique_Posidonia_oceanica+
                            LOCUS_perc_Biocenose_des_algues_infralittorales+
                            LOCUS_cumulative_perc_fine_sand+
                            LOCUS_Number_of_bentic_habitats+
                            LOCUS_perc_Fonds_meubles+
                            LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
                            LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond+
                            distance_closest_fpa+
                            LOCUS_perc_Biocenose_du_coralligene,
                          data=micro.meta.bentos)

LOCUSdbrda.bentosC<-dbrda(noGDP3.distC~LOCUS_perc_Mosaique_Posidonia_oceanica+
                            LOCUS_perc_Biocenose_des_algues_infralittorales+
                            LOCUS_cumulative_perc_fine_sand+
                            LOCUS_Number_of_bentic_habitats+
                            LOCUS_perc_Fonds_meubles+
                            LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
                            LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond+
                            distance_closest_fpa+LOCUS_perc_Biocenose_du_coralligene,
                          data=micro.meta.bentos)

LOCUSdbrda.bentosO<-dbrda(noGDP3.distO~LOCUS_perc_Mosaique_Posidonia_oceanica+
                            LOCUS_perc_Biocenose_des_algues_infralittorales+
                            LOCUS_cumulative_perc_fine_sand+
                            LOCUS_Number_of_bentic_habitats+
                            LOCUS_perc_Fonds_meubles+
                            LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
                            LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond+
                            distance_closest_fpa+LOCUS_perc_Biocenose_du_coralligene,
                          data=micro.meta.bentos)

LOCUSdbrda.bentosF<-dbrda(noGDP3.distF~LOCUS_perc_Mosaique_Posidonia_oceanica+
                            LOCUS_perc_Biocenose_des_algues_infralittorales+
                            LOCUS_cumulative_perc_fine_sand+
                            LOCUS_Number_of_bentic_habitats+
                            LOCUS_perc_Fonds_meubles+
                            LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
                            LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond+
                            distance_closest_fpa+LOCUS_perc_Biocenose_du_coralligene,
                          data=micro.meta.bentos)

LOCUSdbrda.bentosG<-dbrda(noGDP3.distG~LOCUS_perc_Mosaique_Posidonia_oceanica+
                            LOCUS_perc_Biocenose_des_algues_infralittorales+
                            LOCUS_cumulative_perc_fine_sand+
                            LOCUS_Number_of_bentic_habitats+
                            LOCUS_perc_Fonds_meubles+
                            LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
                            LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond+
                            distance_closest_fpa+LOCUS_perc_Biocenose_du_coralligene,
                          data=micro.meta.bentos)

LOCUSdbrda.bentos<-dbrda(noGDP3.dist~LOCUS_perc_Mosaique_Posidonia_oceanica+
                           LOCUS_perc_Biocenose_des_algues_infralittorales+
                           LOCUS_cumulative_perc_fine_sand+
                           LOCUS_Number_of_bentic_habitats+
                           LOCUS_perc_Fonds_meubles+
                           LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
                           LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond+
                           distance_closest_fpa+LOCUS_perc_Biocenose_du_coralligene,
                         data=micro.meta.bentos)

# VIF.CCA 
vif.cca(LOCUSdbrda.bentosP)
vif.cca(LOCUSdbrda.bentosC)
vif.cca(LOCUSdbrda.bentosO)
vif.cca(LOCUSdbrda.bentosF)
vif.cca(LOCUSdbrda.bentosG)

# ANOVA test the dbRDA
anova(LOCUSdbrda.bentosP) #NS
anova(LOCUSdbrda.bentosC) # Sign
anova(LOCUSdbrda.bentosO) # Sign
anova(LOCUSdbrda.bentosF) # Sign
anova(LOCUSdbrda.bentosG) # Sign
anova(LOCUSdbrda.bentos)  # Sign

# ANOVA test the dbRDA by "TERMS"
anova(LOCUSdbrda.bentosC, by="terms") # Sign
anova(LOCUSdbrda.bentosO, by="terms") # Sign
anova(LOCUSdbrda.bentosF, by="terms") # Sign
anova(LOCUSdbrda.bentosG, by="terms") # Sign
anova(LOCUSdbrda.bentos, by="terms")  # Sign

# Perform the stepwise selection on the Significant models
LOCUSmodC<-dbrda(noGDP3.distC~1,data=micro.meta.bentos)
LOCUSmodO<-dbrda(noGDP3.distO~1,data=micro.meta.bentos)
LOCUSmodF<-dbrda(noGDP3.distF~1,data=micro.meta.bentos)
LOCUSmodG<-dbrda(noGDP3.distG~1,data=micro.meta.bentos)
LOCUSmod<-dbrda(noGDP3.dist~1,data=micro.meta.bentos)

#### Class--------------------------------------------------------------------------------
set.seed(1234)
LOCUSselC <- ordistep (LOCUSmodC,
                       scope = formula (LOCUSdbrda.bentosC),
                       R2scope = TRUE, #adjR2.tbrda.bentos,
                       direction = 'forward', permutations = 999)
print(LOCUSselC$call) # LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond + LOCUS_perc_Fonds_meubles

LOCUSdbrda.bentos.selectedC<-dbrda(formula = noGDP3.distC ~ LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond + 
                                     LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica, data = micro.meta.bentos)

vif.cca(LOCUSdbrda.bentos.selectedC)
RsquareAdj(LOCUSdbrda.bentos.selectedC) # 5.18%
anova(LOCUSdbrda.bentos.selectedC,by="margin")
anova(LOCUSdbrda.bentos.selectedC,by="terms")

# LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond **
# LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica   *


#### Order--------------------------------------------------------------------------------
set.seed(1234)
LOCUSselO <- ordistep (LOCUSmodO,
                       scope = formula (LOCUSdbrda.bentosO),
                       R2scope = TRUE, #adjR2.tbrda.bentos,
                       direction = 'forward', permutations = 999)
print(LOCUSselO$call) # LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond + LOCUS_perc_Fonds_meubles

LOCUSdbrda.bentos.selectedO<-dbrda(formula = noGDP3.distO ~ LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond + 
                                     distance_closest_fpa + LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica, 
                                   data = micro.meta.bentos)
RsquareAdj(LOCUSdbrda.bentos.selectedO) #4.72 %
anova(LOCUSdbrda.bentos.selectedO,by="margin")
anova(LOCUSdbrda.bentos.selectedO,by="terms")

vif.cca(LOCUSdbrda.bentos.selectedO)

# LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond **
# distance_closest_fpa                                                            *
# LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica                              *

#### Family--------------------------------------------------------------------------------
set.seed(1234)
LOCUSselF <- ordistep (LOCUSmodF,
                       scope = formula (LOCUSdbrda.bentosF),
                       R2scope = TRUE, #adjR2.tbrda.bentos,
                       direction = 'forward', permutations = 999)
print(LOCUSselF$call) # LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond + LOCUS_perc_Fonds_meubles

LOCUSdbrda.bentos.selectedF<-dbrda(formula = noGDP3.distF ~ LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond + 
                                     LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica, data = micro.meta.bentos)

RsquareAdj(LOCUSdbrda.bentos.selectedF) #3.79
anova(LOCUSdbrda.bentos.selectedF,by="margin")
anova(LOCUSdbrda.bentos.selectedF,by="terms")

vif.cca(LOCUSdbrda.bentos.selectedF)

# LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond **
# LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica                              * 

#### Genus--------------------------------------------------------------------------------
set.seed(1234)
LOCUSselG <- ordistep (LOCUSmodG,
                       scope = formula (LOCUSdbrda.bentosG),
                       R2scope = TRUE, #adjR2.tbrda.bentos,
                       direction = 'forward', permutations = 999)
print(LOCUSselG$call) 

LOCUSdbrda.bentos.selectedG<-dbrda(formula = noGDP3.distG ~  LOCUS_Number_of_bentic_habitats+LOCUS_perc_Biocenose_des_algues_infralittorales + 
                                     distance_closest_fpa, data = micro.meta.bentos)

RsquareAdj(LOCUSdbrda.bentos.selectedG) #4.5%
anova(LOCUSdbrda.bentos.selectedG,by="margin")
anova(LOCUSdbrda.bentos.selectedG,by="terms")

vif.cca(LOCUSdbrda.bentos.selectedG)

# LOCUS_perc_Biocenose_des_algues_infralittorales  1    20.45 1.9983  0.003 **
# LOCUS_Number_of_bentic_habitats                  1    16.75 1.6363  0.017 * 
# distance_closest_fpa                             1    15.40 1.5049  0.017 * 

#### ASV--------------------------------------------------------------------------------

set.seed(1234)
LOCUSsel <- ordistep (LOCUSmod,
                      scope = formula (LOCUSdbrda.bentos),
                      R2scope = TRUE, #adjR2.tbrda.bentos,
                      direction = 'forward', permutations = 999)
print(LOCUSsel$call) 

LOCUSdbrda.bentos.selected<-dbrda(formula = noGDP3.dist ~ LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica + 
                                    LOCUS_Number_of_bentic_habitats + distance_closest_fpa + 
                                    LOCUS_perc_Fonds_meubles + LOCUS_perc_Biocenose_du_coralligene, data = micro.meta.bentos)

RsquareAdj(LOCUSdbrda.bentos.selected) # 11.4%
anova(LOCUSdbrda.bentos.selected,by="margin")
anova(LOCUSdbrda.bentos.selected,by="terms")

vif.cca(LOCUSdbrda.bentos.selected)

# remove the varible with the highest vif

LOCUSdbrda.bentos.selected2<-dbrda(formula = noGDP3.dist ~ LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica +
                                     LOCUS_Number_of_bentic_habitats + distance_closest_fpa +
                                     LOCUS_perc_Biocenose_du_coralligene, data = micro.meta.bentos)



set.seed(1234)
RsquareAdj(LOCUSdbrda.bentos.selected2) # 10.1 % 
anova(LOCUSdbrda.bentos.selected2,by="terms")
anova(LOCUSdbrda.bentos.selected2,by="margin")
vif.cca(LOCUSdbrda.bentos.selected2)


## Multiple linear regression-----
# Influence of benthic habitat, distance from closest FPA on the Shannon diversity 

micro.gdp.alphaP<-micro.gdp.alphaP[,-c(64,65)]

micro.gdp.alphaP.bent<-left_join(micro.gdp.alphaP, benthic.habitat, by="GDP_point")

micro.gdp.alphaC<-micro.gdp.alphaC[,-c(64,65)]
micro.gdp.alphaC.bent<-left_join(micro.gdp.alphaC, benthic.habitat, by="GDP_point")

micro.gdp.alphaO<-micro.gdp.alphaO[,-c(64,65)]
micro.gdp.alphaO.bent<-left_join(micro.gdp.alphaO, benthic.habitat, by="GDP_point")

micro.gdp.alphaF<-micro.gdp.alphaF[,-c(64,65)]
micro.gdp.alphaF.bent<-left_join(micro.gdp.alphaF, benthic.habitat, by="GDP_point")

micro.gdp.alphaG<-micro.gdp.alphaG[,-c(64,65)]
micro.gdp.alphaG.bent<-left_join(micro.gdp.alphaG, benthic.habitat, by="GDP_point")

micro.gdp.alpha<-micro.gdp.alpha[,-c(64,65)]
micro.gdp.alpha.bent<-left_join(micro.gdp.alpha, benthic.habitat, by="GDP_point")

# Check for collinearity
df<-data.frame(micro.gdp.alpha.bent$LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica,
               micro.gdp.alpha.bent$LOCUS_perc_Mosaique_Posidonia_oceanica,
               micro.gdp.alpha.bent$LOCUS_perc_Biocenose_des_algues_infralittorales,
               micro.gdp.alpha.bent$LOCUS_perc_Fonds_meubles,
               micro.gdp.alpha.bent$LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond,
               micro.gdp.alpha.bent$distance_closest_fpa,
               micro.gdp.alpha.bent$LOCUS_cumulative_perc_fine_sand,
               micro.gdp.alpha.bent$LOCUS_perc_Biocenose_du_coralligene,
               micro.gdp.alpha.bent$LOCUS_Number_of_bentic_habitats)
cor<-cor(df)

library(olsrr)

shapiro.test(micro.gdp.alphaP.bent$Shannon)
mirP<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
           LOCUS_perc_Mosaique_Posidonia_oceanica+
           LOCUS_perc_Biocenose_des_algues_infralittorales+
           LOCUS_perc_Fonds_meubles +
           LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond +
           distance_closest_fpa+
           LOCUS_cumulative_perc_fine_sand +
           LOCUS_perc_Biocenose_du_coralligene+
           LOCUS_Number_of_bentic_habitats, data=micro.gdp.alphaP.bent, singular.ok = TRUE)
mirP.all<-ols_step_all_possible(mirP)
best.mirP<-ols_step_best_subset(mirP)
plot (best.mirP)
# Best model is "3"
# LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica LOCUS_perc_Biocenose_des_algues_infralittorales LOCUS_cumulative_perc_fine_sand

mirP<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
           LOCUS_perc_Biocenose_des_algues_infralittorales +
           LOCUS_cumulative_perc_fine_sand  , data=micro.gdp.alphaP.bent)

summary(mirP)
vif(mirP)
shapiro.test(mirP$residuals)
hist(mirP$residuals)
ols_plot_resid_fit(mirP)


shapiro.test(micro.gdp.alphaC.bent$Shannon)
mirC<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
           LOCUS_perc_Mosaique_Posidonia_oceanica+
           LOCUS_perc_Biocenose_des_algues_infralittorales+
           LOCUS_perc_Fonds_meubles +
           LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond +
           distance_closest_fpa+
           LOCUS_cumulative_perc_fine_sand +
           LOCUS_perc_Biocenose_du_coralligene+
           LOCUS_Number_of_bentic_habitats, data=micro.gdp.alphaC.bent)

summary(mirC)
best.mirC<-ols_step_best_subset(mirC)
plot (best.mirC)

# Best model is "3"
# LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica LOCUS_perc_Biocenose_des_algues_infralittorales LOCUS_cumulative_perc_fine_sand   
mirC<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
           LOCUS_perc_Biocenose_des_algues_infralittorales +
           LOCUS_cumulative_perc_fine_sand  , data=micro.gdp.alphaC.bent)

summary(mirC)
vif(mirC)
shapiro.test(mirC$residuals)
hist(mirC$residuals)
ols_plot_resid_fit(mirC)

mirO<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
           LOCUS_perc_Mosaique_Posidonia_oceanica+
           LOCUS_perc_Biocenose_des_algues_infralittorales+
           LOCUS_perc_Fonds_meubles +
           LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond +
           distance_closest_fpa+
           LOCUS_cumulative_perc_fine_sand +
           LOCUS_perc_Biocenose_du_coralligene+
           LOCUS_Number_of_bentic_habitats, data=micro.gdp.alphaO.bent)

summary(mirO)
best.mirO<-ols_step_best_subset(mirO)
plot (best.mirO)

# Best model is "3"
#LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica LOCUS_perc_Biocenose_des_algues_infralittorales LOCUS_cumulative_perc_fine_sand   

mirO<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
           LOCUS_perc_Biocenose_des_algues_infralittorales +
           LOCUS_cumulative_perc_fine_sand  , data=micro.gdp.alphaO.bent)

summary(mirO)
vif(mirO)
shapiro.test(mirO$residuals)
hist(mirO$residuals)
ols_plot_resid_fit(mirO)

mirF<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
           LOCUS_perc_Mosaique_Posidonia_oceanica+
           LOCUS_perc_Biocenose_des_algues_infralittorales+
           LOCUS_perc_Fonds_meubles +
           LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond +
           distance_closest_fpa+
           LOCUS_cumulative_perc_fine_sand +
           LOCUS_perc_Biocenose_du_coralligene+
           LOCUS_Number_of_bentic_habitats, data=micro.gdp.alphaF.bent)

summary(mirF)
best.mirF<-ols_step_best_subset(mirF)
plot (best.mirF)
# Best model is "3"
#  LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica LOCUS_perc_Biocenose_des_algues_infralittorales LOCUS_cumulative_perc_fine_sand  

mirF<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
           LOCUS_perc_Biocenose_des_algues_infralittorales +
           LOCUS_cumulative_perc_fine_sand  , data=micro.gdp.alphaF.bent)

summary(mirF)
vif(mirF)
shapiro.test(mirF$residuals)
hist(mirF$residuals)
ols_plot_resid_fit(mirF)


mirG<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
           LOCUS_perc_Mosaique_Posidonia_oceanica+
           LOCUS_perc_Biocenose_des_algues_infralittorales+
           LOCUS_perc_Fonds_meubles +
           LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond +
           distance_closest_fpa+
           LOCUS_cumulative_perc_fine_sand +
           LOCUS_perc_Biocenose_du_coralligene+
           LOCUS_Number_of_bentic_habitats, data=micro.gdp.alphaG.bent)

summary(mirG)
best.mirG<-ols_step_best_subset(mirG)
plot (best.mirG)
# Best model is "3"
#  LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica LOCUS_perc_Biocenose_des_algues_infralittorales LOCUS_cumulative_perc_fine_sand  
mirG<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
          LOCUS_perc_Biocenose_des_algues_infralittorales +
          LOCUS_cumulative_perc_fine_sand  , data=micro.gdp.alphaG.bent)

summary(mirG)
vif(mirG)
shapiro.test(mirG$residuals)
ols_plot_resid_fit(mirG)

shapiro.test(micro.gdp.alpha.bent$Shannon)
mir<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
           LOCUS_perc_Mosaique_Posidonia_oceanica+
           LOCUS_perc_Biocenose_des_algues_infralittorales+
           LOCUS_perc_Fonds_meubles +
           LOCUS_perc_Biocenose_des_sables_et_graviers_sous_influence_des_courants_de_fond +
           distance_closest_fpa+
           LOCUS_cumulative_perc_fine_sand +
           LOCUS_perc_Biocenose_du_coralligene+
           LOCUS_Number_of_bentic_habitats, data=micro.gdp.alpha.bent)

summary(mir)
best.mir<-ols_step_best_subset(mir)
plot (best.mir)

mir<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica+
          LOCUS_perc_Biocenose_des_algues_infralittorales +
          LOCUS_cumulative_perc_fine_sand  , data=micro.gdp.alpha.bent)

summary(mir)
shapiro.test(mir$residuals)
ols_plot_resid_fit(mir)
mir$coefficients

# Best model is "3"
#  LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica LOCUS_perc_Biocenose_des_algues_infralittorales LOCUS_cumulative_perc_fine_sand  

### Perform the linear regression only with Posidonia oceanica at ASVs level ----
mir.pos<-lm(Shannon~LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica , data=micro.gdp.alpha.bent)

summary(mir.pos)
shapiro.test(mir$residuals)
ols_plot_resid_fit(mir.pos)

# Plot the linear regression with Posidonia oceanica

ggplot(micro.gdp.alpha.bent, aes(x = LOCUS_perc_Biocenose_de_herbier_Posidonia_oceanica, y = Shannon, color=GDP_point)) + 
  geom_point(size=4, alpha=0.4) +
  scale_color_manual(values = gdp)+
   geom_smooth(method = "lm", col = "black",formula = 'y ~ x')+
  xlab("Prop. of Posidonia oceanica meadows") + ylab("Shannon's index")


  
## Result section 3.4.3 --------------------------------------------------------------	
# "No significant relationship is found between the diet and the gut mucosal microbiota taxonomical structure"

# Load the object with the diet data (COI ASV table and taxonomy table in a phyloseq object)

coi.dv.flt3<-get(load("GitHub/phyloseq_filtered_COI_obj.RData"))

# This object was filtered for less the ASVs with less than 3 reads in the dataset, those unclassified at 
# the phylum level and those classified as Diplodus vulgaris and Homo sapiens were removed from the final dataset
## Diet at different taxonomic levels ------

# Change the ID names in the diet object
coi.otu.old<-as.data.frame(otu_table(coi.dv.flt3))%>% rownames_to_column(var="Allgenetics")
meta.coi.old<-sample_data(coi.dv.flt3)%>% data.frame()
meta.coi.new<-meta.coi.old

rownames(meta.coi.new)<-NULL
rownames(meta.coi.new)<-meta.coi.new$FishID

FishID<-meta.coi.new[,c("FishID","Allgenetics")]
coi.otu.new<-coi.otu.old

coi.otu.new<-dplyr::left_join(FishID,coi.otu.old,by="Allgenetics")
rownames(coi.otu.new)<-NULL
coi.otu.new<-coi.otu.new %>% column_to_rownames(var="FishID")
coi.otu.new<-coi.otu.new[,-1]
coi.otu.new<-as.matrix(coi.otu.new)

coi.dv.flt3.mod<-phyloseq(otu_table(coi.otu.new,taxa_are_rows = FALSE),
                          tax_table(tax_table(coi.dv.flt3)),taxa_names(taxa_names(coi.dv.flt3)),
                          sample_data(meta.coi.new),phy_tree(phy_tree(coi.dv.flt3)))


# However 2 samples have a very low number of reads, so I should remove them.
coi.dv.sub<-subset_samples(coi.dv.flt3.mod,sample_sums(coi.dv.flt3)>=7446)
coi.dv.sub.noGDP3<-subset_samples(coi.dv.sub, !GDP_point %in% "GDP-3")


# Reduce the number of samples in the microbial object to the number of samples in the diet object 
ids.diet<-sample_data(coi.dv.sub.noGDP3) %>% data.frame()
ids.diet<-ids.diet$FishID

dv.BO.gdp.noGDP3.diet<-subset_samples(dv.BO.gdp.noGDP3, FishID %in% ids.diet)
dv.BO.gdp.noGDP3.diet<-nozero_taxa(dv.BO.gdp.noGDP3.diet)


# Calculate the Aitchinson distance for the subset of samples used for the diet 
# Agglomerate at different taxonomic levels 
dv.BO.gdp.noGDP3.P.diet<-tax_glom(dv.BO.gdp.noGDP3.diet, taxrank = "Phylum",NArm = FALSE)
dv.BO.gdp.noGDP3.C.diet<-tax_glom(dv.BO.gdp.noGDP3.diet, taxrank = "Class",NArm = FALSE)
dv.BO.gdp.noGDP3.O.diet<-tax_glom(dv.BO.gdp.noGDP3.diet, taxrank = "Order",NArm = FALSE)
dv.BO.gdp.noGDP3.F.diet<-tax_glom(dv.BO.gdp.noGDP3.diet, taxrank = "Family",NArm = FALSE)
dv.BO.gdp.noGDP3.G.diet<-tax_glom(dv.BO.gdp.noGDP3.diet, taxrank = "Genus",NArm = FALSE)


#Transform
dv.BO.gdp.noGDP3.clrP.diet<-tran.clr.not.glom(dv.BO.gdp.noGDP3.P.diet)
dv.BO.gdp.noGDP3.clrC.diet<-tran.clr.not.glom(dv.BO.gdp.noGDP3.C.diet)
dv.BO.gdp.noGDP3.clrO.diet<-tran.clr.not.glom(dv.BO.gdp.noGDP3.O.diet)
dv.BO.gdp.noGDP3.clrF.diet<-tran.clr.not.glom(dv.BO.gdp.noGDP3.F.diet)
dv.BO.gdp.noGDP3.clrG.diet<-tran.clr.not.glom(dv.BO.gdp.noGDP3.G.diet)
dv.BO.gdp.noGDP3.clr.diet<-tran.clr.not.glom(dv.BO.gdp.noGDP3.diet)

# Calculate distance

noGDP3.diet.distP<-phyloseq::distance(dv.BO.gdp.noGDP3.clrP.diet,method="euclidean") #generate a distance matrix
noGDP3.diet.distC<-phyloseq::distance(dv.BO.gdp.noGDP3.clrC.diet,method="euclidean") #generate a distance matrix
noGDP3.diet.distO<-phyloseq::distance(dv.BO.gdp.noGDP3.clrO.diet,method="euclidean") #generate a distance matrix
noGDP3.diet.distF<-phyloseq::distance(dv.BO.gdp.noGDP3.clrF.diet,method="euclidean") #generate a distance matrix
noGDP3.diet.distG<-phyloseq::distance(dv.BO.gdp.noGDP3.clrG.diet,method="euclidean") #generate a distance matrix
noGDP3.diet.dist<-phyloseq::distance(dv.BO.gdp.noGDP3.clr.diet,method="euclidean") #generate a distance matrix

# Agglomerate and Transform at different taxonomic levels 
coi.P<-tax_glom(coi.dv.sub.noGDP3, taxrank = "Phylum", NArm = FALSE)
coi.C<-tax_glom(coi.dv.sub.noGDP3, taxrank = "Class", NArm = FALSE)
coi.O<-tax_glom(coi.dv.sub.noGDP3, taxrank = "Order", NArm = FALSE)
coi.F<-tax_glom(coi.dv.sub.noGDP3, taxrank = "Family", NArm = FALSE)
coi.G<-tax_glom(coi.dv.sub.noGDP3, taxrank = "Genus", NArm = FALSE)
coi.sp<-tax_glom(coi.dv.sub.noGDP3, taxrank = "Species", NArm = FALSE)


# CLR transformation 
coi.P.clr<-tran.clr.not.glom(coi.P) # 23 taxa and 27 samples
coi.C.clr<-tran.clr.not.glom(coi.C) # 57 taxa and 27 samples
coi.O.clr<-tran.clr.not.glom(coi.O) # 128 taxa and 27 samples 
coi.F.clr<-tran.clr.not.glom(coi.F) # 212 taxa and 27 samples
coi.G.clr<-tran.clr.not.glom(coi.G) # 266 taxa and 27 samples 
coi.S.clr<-tran.clr.not.glom(coi.sp) # 299 taxa and 27 samples

# Distances 
coi.P.dist<-phyloseq::distance(coi.P.clr,method="euclidean") #generate a distance matrix
coi.C.dist<-phyloseq::distance(coi.C.clr,method="euclidean") #generate a distance matrix
coi.O.dist<-phyloseq::distance(coi.O.clr,method="euclidean") #generate a distance matrix
coi.F.dist<-phyloseq::distance(coi.F.clr,method="euclidean") #generate a distance matrix
coi.G.dist<-phyloseq::distance(coi.G.clr,method="euclidean") #generate a distance matrix
coi.S.dist<-phyloseq::distance(coi.S.clr,method="euclidean") #generate a distance matrix


## Does the diet changes across the different GDP locations? ------------------------
# PERMDISP   
set.seed(1000)
permdisp(coi.P.clr,"GDP_point") # NS
permdisp(coi.C.clr,"GDP_point") # NS
permdisp(coi.O.clr,"GDP_point") # NS
permdisp(coi.F.clr,"GDP_point") # NS
permdisp(coi.G.clr,"GDP_point") # NS
permdisp(coi.S.clr,"GDP_point") # NS

# PERMANOVA
set.seed(1234)
res.permP<-permanova(coi.P.clr,"GDP_point") # Sign
res.permC<-permanova(coi.C.clr,"GDP_point") # Sign
res.permO<-permanova(coi.O.clr,"GDP_point") # Sign
res.permF<-permanova(coi.F.clr,"GDP_point") # Sign
res.permG<-permanova(coi.G.clr,"GDP_point") # Sign
res.permS<-permanova(coi.S.clr,"GDP_point") # Sign

# PAIRWISE 
meta.coi.sub<-subset(meta.coi.new,FishID %in% sample_names(coi.S.clr))

pairwise.adonis2(coi.P.dist~GDP_point, data = meta.coi.sub)
pairwise.adonis2(coi.C.dist~GDP_point, data = meta.coi.sub)
pairwise.adonis2(coi.O.dist~GDP_point, data = meta.coi.sub)
pairwise.adonis2(coi.F.dist~GDP_point, data = meta.coi.sub)
pairwise.adonis2(coi.G.dist~GDP_point, data = meta.coi.sub)
pairwise.adonis2(coi.S.dist~GDP_point, data = meta.coi.sub)



## Does the dissimilarities of the gut microbiota between fish correlate with their different diet? -----------
## Mantel test at different taxonomical levels 

# Phylum diet 
mantel(coi.P.dist,noGDP3.diet.distP,method = "spearman")# Not Sign
mantel(coi.P.dist,noGDP3.diet.distC,method = "spearman")# Not Sign
mantel(coi.P.dist,noGDP3.diet.distO,method = "spearman")# Not Sign
mantel(coi.P.dist,noGDP3.diet.distF,method = "spearman")# Not Sign
mantel(coi.P.dist,noGDP3.diet.distG,method = "spearman")# Not Sign
mantel(coi.P.dist,noGDP3.diet.dist,method = "spearman")# Not Sign

# Class diet
mantel(coi.C.dist,noGDP3.diet.distP,method = "spearman")# Not Sign
mantel(coi.C.dist,noGDP3.diet.distC,method = "spearman")# Not Sign
mantel(coi.C.dist,noGDP3.diet.distO,method = "spearman")# Not Sign
mantel(coi.C.dist,noGDP3.diet.distF,method = "spearman")# Not Sign
mantel(coi.C.dist,noGDP3.diet.distG,method = "spearman")# Not Sign
mantel(coi.C.dist,noGDP3.diet.dist,method = "spearman")# Not Sign

# Order diet
mantel(coi.O.dist,noGDP3.diet.distP,method = "spearman")# Not Sign
mantel(coi.O.dist,noGDP3.diet.distC,method = "spearman")# Not Sign
mantel(coi.O.dist,noGDP3.diet.distO,method = "spearman")# Not Sign
mantel(coi.O.dist,noGDP3.diet.distF,method = "spearman")# Not Sign
mantel(coi.O.dist,noGDP3.diet.distG,method = "spearman")# Not Sign
mantel(coi.O.dist,noGDP3.diet.dist,method = "spearman")# Not Sign

# Family diet
mantel(coi.F.dist,noGDP3.diet.distP,method = "spearman")# Not Sign
mantel(coi.F.dist,noGDP3.diet.distC,method = "spearman")# Not Sign
mantel(coi.F.dist,noGDP3.diet.distO,method = "spearman")# Not Sign
mantel(coi.F.dist,noGDP3.diet.distF,method = "spearman")# Not Sign
mantel(coi.F.dist,noGDP3.diet.distG,method = "spearman")# Not Sign
mantel(coi.F.dist,noGDP3.diet.dist,method = "spearman")# Not Sign

# Genus diet
mantel(coi.G.dist,noGDP3.diet.distP,method = "spearman")# Not Sign
mantel(coi.G.dist,noGDP3.diet.distC,method = "spearman")# Not Sign
mantel(coi.G.dist,noGDP3.diet.distO,method = "spearman")# Not Sign
mantel(coi.G.dist,noGDP3.diet.distF,method = "spearman")# Not Sign
mantel(coi.G.dist,noGDP3.diet.distG,method = "spearman")# Not Sign
mantel(coi.G.dist,noGDP3.diet.dist,method = "spearman")# Not Sign

# Species diet
mantel(coi.S.dist,noGDP3.diet.distP,method = "spearman")# Not Sign
mantel(coi.S.dist,noGDP3.diet.distC,method = "spearman")# Not Sign
mantel(coi.S.dist,noGDP3.diet.distO,method = "spearman")# Not Sign
mantel(coi.S.dist,noGDP3.diet.distF,method = "spearman")# Not Sign
mantel(coi.S.dist,noGDP3.diet.distG,method = "spearman")# Not Sign
mantel(coi.S.dist,noGDP3.diet.dist,method = "spearman") # Not Sign


# Is the structure of the micorbiota influenced by the prey diversity (Shannon's index') ----
# Calculate Alpha diversity of the diet 

alpha.dietP<-estimate_richness(coi.P,split= T, measures = c("Shannon")) %>%
  mutate(FishID = sample_names(coi.P)) # add sample IDs
colnames(alpha.dietP)[1]<-"Shannon.diet.Phylum"

alpha.dietC<-estimate_richness(coi.C,split= T, measures = c("Shannon")) %>%
  mutate(FishID = sample_names(coi.C)) # add sample IDs
colnames(alpha.dietC)[1]<-"Shannon.diet.Class"

alpha.dietO<-estimate_richness(coi.O,split= T, measures = c("Shannon")) %>%
  mutate(FishID = sample_names(coi.O)) # add sample IDs
colnames(alpha.dietO)[1]<-"Shannon.diet.Order"

alpha.dietF<-estimate_richness(coi.F,split= T, measures = c("Shannon")) %>%
  mutate(FishID = sample_names(coi.F)) # add sample IDs
colnames(alpha.dietF)[1]<-"Shannon.diet.Family"

alpha.dietG<-estimate_richness(coi.G,split= T, measures = c("Shannon")) %>%
  mutate(FishID = sample_names(coi.G)) # add sample IDs
colnames(alpha.dietG)[1]<-"Shannon.diet.Genus"

alpha.dietS<-estimate_richness(coi.sp,split= T, measures = c("Shannon")) %>%
  mutate(FishID = sample_names(coi.sp)) # add sample IDs
colnames(alpha.dietS)[1]<-"Shannon.diet.Species"

alpha.diet<-left_join(alpha.dietP,alpha.dietC,by="FishID")
alpha.diet<-left_join(alpha.diet,alpha.dietO,by="FishID")
alpha.diet<-left_join(alpha.diet,alpha.dietF,by="FishID")
alpha.diet<-left_join(alpha.diet,alpha.dietG,by="FishID")
alpha.diet<-left_join(alpha.diet,alpha.dietS,by="FishID")

meta.diet<-sample_data(dv.BO.gdp.noGDP3.diet)%>%data.frame()
meta.diet<-left_join(meta.diet,alpha.diet,by="FishID")


noGDP3.diet.distC
noGDP3.diet.distO
noGDP3.diet.distF
noGDP3.diet.distG
noGDP3.diet.dist

adonis2(noGDP3.diet.distP~meta.diet$Shannon.diet.Phylum)
adonis2(noGDP3.diet.distO~meta.diet$Shannon.diet.Phylum)
adonis2(noGDP3.diet.distC~meta.diet$Shannon.diet.Phylum)
adonis2(noGDP3.diet.distF~meta.diet$Shannon.diet.Phylum)
adonis2(noGDP3.diet.distG~meta.diet$Shannon.diet.Phylum)
adonis2(noGDP3.diet.dist~meta.diet$Shannon.diet.Phylum)

adonis2(noGDP3.diet.distP~meta.diet$Shannon.diet.Class)
adonis2(noGDP3.diet.distO~meta.diet$Shannon.diet.Class)
adonis2(noGDP3.diet.distC~meta.diet$Shannon.diet.Class)
adonis2(noGDP3.diet.distF~meta.diet$Shannon.diet.Class)
adonis2(noGDP3.diet.distG~meta.diet$Shannon.diet.Class)
adonis2(noGDP3.diet.dist~meta.diet$Shannon.diet.Class)

adonis2(noGDP3.diet.distP~meta.diet$Shannon.diet.Order)
adonis2(noGDP3.diet.distO~meta.diet$Shannon.diet.Order)
adonis2(noGDP3.diet.distC~meta.diet$Shannon.diet.Order)
adonis2(noGDP3.diet.distF~meta.diet$Shannon.diet.Order)
adonis2(noGDP3.diet.distG~meta.diet$Shannon.diet.Order)
adonis2(noGDP3.diet.dist~meta.diet$Shannon.diet.Order)

adonis2(noGDP3.diet.distP~meta.diet$Shannon.diet.Family)
adonis2(noGDP3.diet.distO~meta.diet$Shannon.diet.Family)
adonis2(noGDP3.diet.distC~meta.diet$Shannon.diet.Family)
adonis2(noGDP3.diet.distF~meta.diet$Shannon.diet.Family)
adonis2(noGDP3.diet.distG~meta.diet$Shannon.diet.Family)
adonis2(noGDP3.diet.dist~meta.diet$Shannon.diet.Family)

adonis2(noGDP3.diet.distP~meta.diet$Shannon.diet.Genus)
adonis2(noGDP3.diet.distO~meta.diet$Shannon.diet.Genus)
adonis2(noGDP3.diet.distC~meta.diet$Shannon.diet.Genus)
adonis2(noGDP3.diet.distF~meta.diet$Shannon.diet.Genus)
adonis2(noGDP3.diet.distG~meta.diet$Shannon.diet.Genus)
adonis2(noGDP3.diet.dist~meta.diet$Shannon.diet.Genus)

adonis2(noGDP3.diet.distP~meta.diet$Shannon.diet.Species)
adonis2(noGDP3.diet.distO~meta.diet$Shannon.diet.Species)
adonis2(noGDP3.diet.distC~meta.diet$Shannon.diet.Species)
adonis2(noGDP3.diet.distF~meta.diet$Shannon.diet.Species)
adonis2(noGDP3.diet.distG~meta.diet$Shannon.diet.Species)
adonis2(noGDP3.diet.dist~meta.diet$Shannon.diet.Species)

# Is there a correlation between the alpha diversity of the diet and that of the micorbitoa at all taxonomical levels?
# subset the alpha diversity tables for the samples analyzed only for the the diet
micro.gdp.alphaP.sub<-subset(micro.gdp.alphaP, FishID %in% meta.diet$FishID)
micro.gdp.alphaC.sub<-subset(micro.gdp.alphaC, FishID %in% meta.diet$FishID)
micro.gdp.alphaO.sub<-subset(micro.gdp.alphaO, FishID %in% meta.diet$FishID)
micro.gdp.alphaF.sub<-subset(micro.gdp.alphaF, FishID %in% meta.diet$FishID)
micro.gdp.alphaG.sub<-subset(micro.gdp.alphaG, FishID %in% meta.diet$FishID)
micro.gdp.alpha.sub<-subset(micro.gdp.alpha, FishID %in% meta.diet$FishID)

cor.test(micro.gdp.alphaP.sub$Shannon,meta.diet$Shannon.diet.Phylum,method = "spearman") # NS
cor.test(micro.gdp.alphaP.sub$Shannon,meta.diet$Shannon.diet.Class,method = "spearman") # NS
cor.test(micro.gdp.alphaP.sub$Shannon,meta.diet$Shannon.diet.Order,method = "spearman") # NS
cor.test(micro.gdp.alphaP.sub$Shannon,meta.diet$Shannon.diet.Family,method = "spearman") # NS
cor.test(micro.gdp.alphaP.sub$Shannon,meta.diet$Shannon.diet.Genus,method = "spearman") # NS
cor.test(micro.gdp.alphaP.sub$Shannon,meta.diet$Shannon.diet.Species,method = "spearman") # NS

cor.test(micro.gdp.alphaC.sub$Shannon,meta.diet$Shannon.diet.Phylum,method = "spearman") # NS
cor.test(micro.gdp.alphaC.sub$Shannon,meta.diet$Shannon.diet.Class,method = "spearman") # NS
cor.test(micro.gdp.alphaC.sub$Shannon,meta.diet$Shannon.diet.Order,method = "spearman") # NS
cor.test(micro.gdp.alphaC.sub$Shannon,meta.diet$Shannon.diet.Family,method = "spearman") # NS
cor.test(micro.gdp.alphaC.sub$Shannon,meta.diet$Shannon.diet.Genus,method = "spearman") # NS
cor.test(micro.gdp.alphaC.sub$Shannon,meta.diet$Shannon.diet.Species,method = "spearman") # NS

cor.test(micro.gdp.alphaO.sub$Shannon,meta.diet$Shannon.diet.Phylum,method = "spearman") # NS
cor.test(micro.gdp.alphaO.sub$Shannon,meta.diet$Shannon.diet.Class,method = "spearman") # NS
cor.test(micro.gdp.alphaO.sub$Shannon,meta.diet$Shannon.diet.Order,method = "spearman") # NS
cor.test(micro.gdp.alphaO.sub$Shannon,meta.diet$Shannon.diet.Family,method = "spearman") # NS
cor.test(micro.gdp.alphaO.sub$Shannon,meta.diet$Shannon.diet.Genus,method = "spearman") # NS
cor.test(micro.gdp.alphaO.sub$Shannon,meta.diet$Shannon.diet.Species,method = "spearman") # NS

cor.test(micro.gdp.alphaF.sub$Shannon,meta.diet$Shannon.diet.Phylum,method = "spearman") # NS
cor.test(micro.gdp.alphaF.sub$Shannon,meta.diet$Shannon.diet.Class,method = "spearman") # NS
cor.test(micro.gdp.alphaF.sub$Shannon,meta.diet$Shannon.diet.Order,method = "spearman") # NS
cor.test(micro.gdp.alphaF.sub$Shannon,meta.diet$Shannon.diet.Family,method = "spearman") # NS
cor.test(micro.gdp.alphaF.sub$Shannon,meta.diet$Shannon.diet.Genus,method = "spearman") # NS
cor.test(micro.gdp.alphaF.sub$Shannon,meta.diet$Shannon.diet.Species,method = "spearman") # NS

cor.test(micro.gdp.alphaG.sub$Shannon,meta.diet$Shannon.diet.Phylum,method = "spearman") # NS
cor.test(micro.gdp.alphaG.sub$Shannon,meta.diet$Shannon.diet.Class,method = "spearman") # NS
cor.test(micro.gdp.alphaG.sub$Shannon,meta.diet$Shannon.diet.Order,method = "spearman") # NS
cor.test(micro.gdp.alphaG.sub$Shannon,meta.diet$Shannon.diet.Family,method = "spearman") # NS
cor.test(micro.gdp.alphaG.sub$Shannon,meta.diet$Shannon.diet.Genus,method = "spearman") # NS
cor.test(micro.gdp.alphaG.sub$Shannon,meta.diet$Shannon.diet.Species,method = "spearman") # NS

cor.test(micro.gdp.alpha.sub$Shannon,meta.diet$Shannon.diet.Phylum,method = "spearman") # NS
cor.test(micro.gdp.alpha.sub$Shannon,meta.diet$Shannon.diet.Class,method = "spearman") # NS
cor.test(micro.gdp.alpha.sub$Shannon,meta.diet$Shannon.diet.Order,method = "spearman") # NS
cor.test(micro.gdp.alpha.sub$Shannon,meta.diet$Shannon.diet.Family,method = "spearman") # NS
cor.test(micro.gdp.alpha.sub$Shannon,meta.diet$Shannon.diet.Genus,method = "spearman") # NS
cor.test(micro.gdp.alpha.sub$Shannon,meta.diet$Shannon.diet.Species,method = "spearman") # NS
