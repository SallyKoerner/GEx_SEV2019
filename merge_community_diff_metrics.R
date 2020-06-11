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
  
RACdiff <- read.csv('gex_RACdiff_site_ave.csv') #RAC differences
compDiff <- read.csv('gex_multdiff_site_ave.csv') #composition difference
domIDdiff <- read.csv('Diff_BP_Dom_allyrs_11June2020.csv') %>%  #difference in identity of dominant spp
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage) %>% 
  summarise(bp_rr=mean(bp_rr), diff_sp=mean(diff_sp))


### where we stopped on 11 June 2020


#####Kim's Old Code satrts here
###import data
domDiff <- read.csv('OnlineVersion_withSimpD_Apr2018.csv')%>%select(-X) #dominance difference, with a bunch of site level covariates

codomDiff <- read.csv('GEx_codominane_04242019.csv')%>% #codominance imported, need to calculate diff
  select(site, block, plot, trt, num_codominants)%>%
  unique()%>%
  spread(key=trt, value=num_codominants)%>%
  mutate(num_codom_diff=log(G/U))%>%
  filter(!is.na(num_codom_diff))%>%
  group_by(site, block)%>%
  summarise(num_codom_diff=mean(num_codom_diff))%>%
  ungroup()%>%
  group_by(site)%>%
  summarise(num_codom_diff=mean(num_codom_diff))%>%
  ungroup()
CmaxDiff <- read.csv('GEx_codominane_04242019.csv')%>% #Cmax imported, need to calculate diff
  select(site, block, plot, trt, Cmax)%>%
  unique()%>%
  spread(key=trt, value=Cmax)%>%
  mutate(Cmax_diff=log(G/U))%>%
  filter(!is.na(Cmax_diff))%>%
  group_by(site, block)%>%
  summarise(Cmax_diff=mean(Cmax_diff))%>%
  ungroup()%>%
  group_by(site)%>%
  summarise(Cmax_diff=mean(Cmax_diff))%>%
  ungroup()
photopath <- read.csv('percent_photosynthetic_pathway.csv') #photosynthetic pathways
phylorealms <- read.csv('PhyloRealms.csv')%>%select(-X) #phylogenetic realms
domRankDiff <- read.csv('sitelevel_domspecies_rankchange.csv')%>%select(-X) #rank change of dominant species
domPercentDiff <- read.csv('DomIdentityNumChangeLRR_bySite_SEV_April2019.csv') #percent abundance difference of dominant species
blockNum <- read.csv('GEx_cleaned_v3.csv')%>%
  select(site, block)%>%
  unique()%>%
  group_by(site)%>%
  summarise(num_blocks=length(block))%>%
  ungroup()
  

###merge together
compDiffSite <- domDiff%>%
  left_join(compDiff)%>%
  left_join(RACdiff)%>%
  left_join(domIDdiff)%>%
  left_join(CmaxDiff)%>%
  left_join(codomDiff)%>%
  left_join(photopath)%>%
  left_join(phylorealms)%>%
  left_join(domRankDiff)%>%
  left_join(domPercentDiff)%>%
  left_join(blockNum)%>%
  mutate(site_dom=(GDom+UDom)/2) #calculate average dominance between grazed and ungrazed areas

# write.csv(compDiffSite, 'community_difference_allmetrics_siteavg_04262019b.csv', row.names=F)


####Sallys Add in strats here --- reimport above write

###import data
All <- read.csv('community_difference_allmetrics_siteavg_04262019b.csv')
Climate <- read.csv('climatePCAs.csv')
Meta <- read.csv('GEx-metadata-with-other-env-layers-v2.csv')
Herb<- read.csv('Meta_SEV2019_v2_with_body_size.csv')

All2<-All %>% 
  left_join(Meta)%>%
  left_join(Climate)%>%
  left_join(Herb) %>% 
  mutate(ALLC3=(C3+C3.)) %>% 
  mutate(ALLC4=(C4+C4.)) %>% 
  mutate(PhotoMix=abs(ALLC3-ALLC4)) 

write.csv(All2, 'community_difference_allmetrics_siteavg_09June2020.csv', row.names=F)
