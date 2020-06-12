library(tidyverse)
library(codyn)
library(rsq)
library(vegan)
source("Wd.R")


theme_set(theme_bw(12))

setwd("C:\\Users\\mavolio2\\Dropbox\\GEx_VirtualWorkshop_June2020\\")

dat<-read.csv("GEx_cleaned_11June2020.csv")

dat2<-dat%>%
  mutate(drop=ifelse(genus_species_use=="#N/A"|genus_species_use=="Dead unidentified"|genus_species_use=="Leaf.Litter"|genus_species_use=="cactus__dead_", 1, 0))%>%
  filter(drop!=1)

exage<-dat%>%
  select(site, exage)%>%
  unique()%>%
  filter(exage>10)

lsyear<-dat%>%
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage)%>%
  mutate(block = as.numeric(block))%>%
  mutate(site_block=paste(site, block, sep="##"))


#dropping probelmatic blocks and averaging over experiments where there are several blocks in a plot. 
#Averaging up to a block for Konza, Junner Koeland, and Jornada sites.

lsyear2<-lsyear%>%
  filter(site_block!="California_Sedgwick_Lisque##A"&
           site_block!="California_Sedgwick_Lisque##B"&
           site_block!="California_Sedgwick_Lisque##G"&
           site_block!="DesertLow##Alkali")%>%
           group_by(site, block, site_block, trt, exage, genus_species_use)%>%
  summarize(relcov=mean(relcov))

numreps<-lsyear%>%
  select(site, block)%>%
  unique()%>%
  group_by(site)%>%
  summarize(num=length(block))

lsyear2.2<-lsyear2%>%
  left_join(numreps)%>%
  filter(num>4)

gex_permanova<-data.frame()
sitelist<-unique(lsyear2.2$site)

for (i in sitelist){
  
  subset<-lsyear2.2%>%
    filter(site==i) %>%
    mutate(treatment=trt)
  
  wide<-subset%>%
    spread(genus_species_use, relcov, fill=0)
  
  #species<-wide[,7:ncol(wide)]
  #env<-wide[,1:6]
  out <- adonis(wide[,12:ncol(wide)]~trt, wide, strata = wide$block)

  
  # Assemble data frame
  perm_out <- data.frame(
    site = i,
    perm_Pvalue =  out$aov.tab$'Pr(>F)'[1],
    Tw_Pvalue = Tw_test_out$p.value,
    WdS_Pvalue = WdS_test_out$p.value
  )
  perm_out
  
  gex_permanova<-rbind(gex_permanova, perm_out)  
   
}


## look at NMDS to see gut check
msub<-lsyear2.2%>%
  filter(site=="AUS_Savernake")

mwide<-msub%>%
  spread(genus_species_use, relcov, fill=0)

species<-mwide[,7:ncol(mwide)]

env<-mwide[,1:6]


mds<-metaMDS(species, distance = "bray")

scores <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores
scores2<- cbind(env, scores) # binds the NMDS scores of year i to all years previously run

ggplot(scores2, aes(x=NMDS1, y=NMDS2, color=trt, shape=block))+
  geom_point(size=5)



###making comp diff shoowsh figure
comp<-read.csv("gex_multdiff_site_block.csv")

#getting CI for each site
gex_multdiff_CI<-comp%>%
  group_by(site) %>% 
  summarize_at(vars(composition_diff), funs(mean, sd, length)) %>% 
  mutate(se=sd/sqrt(length),
         ci=se*1.96)

#geting global mean comp diff
mean<-mean(gex_multdiff_CI$mean)
sd<-sd(gex_multdiff_CI$mean)
se<-sd/sqrt(252)
ci<-se*1.96
site<-"All Sites"

all<- data.frame(site, mean, sd, se, ci)

gex_multdiff_CI <- gex_multdiff_CI[order(gex_multdiff_CI$mean), ]  # sort
gex_multdiff_CI2<-gex_multdiff_CI%>%
  bind_rows(all)%>%
  mutate(colortrt=ifelse(site=="All Sites", 1,0))
gex_multdiff_CI2$site2 <- factor(gex_multdiff_CI2$site, levels = gex_multdiff_CI2$site)


ggplot(data=gex_multdiff_CI2, aes(x=site2, y=mean, color=as.factor(colortrt)))+
  geom_point()+
  scale_color_manual(values=c("black", "red"))+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci))+
  coord_flip()+
  xlab("Site")+
  ylab("Compositional Difference")+
  #scale_y_continuous(limits = c(0,1))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

###making comp diff sig figure based on permanovas

####flagging sig values
permsig <- gex_permanova %>%
  mutate(sig = ifelse(perm_Pvalue<.1, 1, 0))

toplot<-gex_multdiff_CI%>%
  right_join(permsig)

toplot$site2 <- factor(toplot$site, levels = toplot$site[order(toplot$mean)])

ggplot(data=toplot, aes(x=site2, y=mean, color=as.factor(sig)))+
  geom_point()+
  scale_color_manual(values=c("black", "red"), name="Sig. Permanova", labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci))+
  coord_flip()+
  xlab("Site")+
  ylab("Compositional Difference")+
  #scale_y_continuous(limits = c(0,1))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank())
