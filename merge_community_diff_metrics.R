################################################################################
##  merge_community_diff_metrics.R: Mergine files of community difference metrics.
##
##  Author: Kimberly Komatsu and Sally Koerner
##  Date created: April 25, 2019
################################################################################

###########Complete rewrite by Sally Koerner and Meghan Avolio ######################
###########Rewrite June 11, 2020
library(tidyverse)
library(codyn)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\data')
setwd("~/Dropbox/GEx_VirtualWorkshop_June2020")

#Calculate Simposons LRR
dat<-read.csv('GEx_cleaned_11June2020.csv')
dat2<-dat%>%
  mutate(drop=ifelse(genus_species_use=="#N/A"|genus_species_use=="Dead unidentified"|genus_species_use=="Leaf.Litter"|genus_species_use=="cactus__dead_", 1, 0))%>%
  filter(drop!=1)
lsyear<-dat2%>%
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage) %>% 
  mutate(site_block=paste(site, block, sep="##"))
  

#dropping probelmatic blocks and averaging over experiments where there are several blocks in a plot. 
#Averaging up to a block for Konza, Junner Koeland, and Jornada sites.

lsyear2<-lsyear%>%
  filter(site_block!="California_Sedgwick_Lisque##A"&
           site_block!="California_Sedgwick_Lisque##B"&
           site_block!="California_Sedgwick_Lisque##G"&
           site_block!="DesertLow##Alkali")%>%
  group_by(site, block, exage, year, trt, genus_species_use, site_block)%>%
  summarize(relcov=mean(relcov)) %>% 
  mutate(site_block_trt=paste(site_block, trt, sep="##"))

SimpD<-community_diversity(lsyear2, replicate.var="site_block_trt", abundance.var="relcov", metric="InverseSimpson") %>% 
  mutate(SimpD=1/InverseSimpson) %>% 
  separate(site_block_trt, into=c("site", "block", "trt"), sep="##") %>% 
  group_by(site, trt) %>% 
  summarise(SimpD=mean(SimpD)) %>% 
  spread(trt, SimpD) %>% 
  mutate(SimpD_rr=(log(G/U))) %>% 
  select(-G, -U)

exage<-lsyear2 %>% 
  group_by(site, exage) %>% 
  summarise(num=n()) %>% 
  select(-num)
  
RACdiff <- read.csv('gex_RACdiff_site_ave.csv') #RAC differences
compDiff <- read.csv('gex_multdiff_site_ave.csv') #composition difference
domIDdiff <- read.csv('Diff_BP_Dom_allyrs_11June2020.csv') %>%  #difference in identity of dominant spp
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage) %>% 
  summarise(bp_rr=mean(bp_rr), diff_sp=mean(diff_sp))

codomDiff <- read.csv('GEx_codominance_06112020.csv')%>% #codominance imported, need to calculate diff
  select(site, year, trt, num_codominants)%>%
  unique()%>%
  spread(key=trt, value=num_codominants)%>%
  mutate(num_codom_diff=log(G/U))%>%
  filter(!is.na(num_codom_diff)) %>% 
  ungroup()%>%
  group_by(site) %>% 
  mutate(lyear=max(year))%>%
  filter(year==lyear) %>% 
  select(site, year, num_codom_diff)
CmaxDiff <- read.csv('GEx_codominance_06112020.csv')%>% #Cmax imported, need to calculate diff
  select(site, year, trt, Cmax)%>%
  unique()%>%
  spread(key=trt, value=Cmax)%>%
  mutate(Cmax_diff=log(G/U))%>%
  filter(!is.na(Cmax_diff))%>%
  ungroup()%>%
  group_by(site) %>% 
  mutate(lyear=max(year))%>%
  filter(year==lyear) %>% 
  select(site, year, Cmax_diff) %>% 
  ungroup()

phylorealms <- read.csv('PhyloRealms.csv')%>%select(-X) #phylogenetic realms

blockNum <- dat%>%
  select(site, block)%>%
  unique()%>%
  group_by(site)%>%
  summarise(num_blocks=length(block))%>%
  ungroup()
photopath <- read.csv('percent_photosynthetic_pathway.csv') #photosynthetic pathways

###merge together
compDiffSite <- SimpD%>%
  left_join(compDiff)%>%
  left_join(RACdiff)%>%
  left_join(domIDdiff)%>%
  left_join(CmaxDiff)%>%
  left_join(codomDiff)%>%
  left_join(photopath)%>%
  left_join(phylorealms)%>%
  left_join(blockNum) %>% 
  left_join(exage)

###import data
All <- compDiffSite
Climate <- read.csv('climatePCAs.csv')
Meta <- read.csv('GEx-metadata-with-other-env-layers-v2.csv')
Herb<- read.csv('Meta_SEV2019_v2_with_body_size.csv')

All2<-All %>% 
  left_join(Meta)%>%
  left_join(Climate)%>%
  left_join(Herb) %>% 
  mutate(ALLC3=(C3)) %>% 
  mutate(ALLC4=(C4+C4.)) %>% 
  mutate(PhotoMix=abs(ALLC3-ALLC4)) 

write.csv(All2, 'community_difference_allmetrics_siteavg_12June2020c.csv', row.names=F)


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  test <- cor.test(x,y) 
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 1),
                   symbols = c("*", " "))
  
  
  text(0.5, 0.5, txt, cex = 2)
  text(0.8, 0.5, Signif, cex=5, col="red")
}

toplot<-All%>%
  filter(!is.na(bp_rr))

pairs(toplot[,c(3, 4, 5, 6, 7, 8, 2, 9, 11, 12)], pch = 21, font.labels=1, cex.labels=2,upper.panel=panel.cor)
