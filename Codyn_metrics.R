library(tidyverse)
library(codyn)
library(rsq)

setwd("C:\\Users\\mavolio2\\Dropbox\\GEx_VirtualWorkshop_June2020\\")

dat<-read.csv("GEx_cleaned_10June2020.csv")

lsyear<-dat%>%
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage)%>%
  mutate(site_block=paste(site, block, sep="##"))

#dropping probelmatic blocks and averaging over experiments where there are several blocks in a plot. 
#Averaging up to a block for Konza, Junner Koeland, and Jornada sites.

lsyear2<-lsyear%>%
  filter(site_block!="California_Sedgwick_Lisque##A"&
           site_block!="California_Sedgwick_Lisque##B"&
           site_block!="California_Sedgwick_Lisque##G"&
           site_block!="DesertLow##Alkali")%>%
  group_by(site, block, exage, year, trt, genus_species, site_block)%>%
  summarize(relcov=mean(relcov))

gex_multdiff<-data.frame()
siteblock<-unique(lsyear2$site_block)

for (i in 1:length(siteblock)){
  
  subset<-lsyear2%>%
    filter(site_block==siteblock[i]) %>%
    mutate(treatment=trt)
  
  out <- multivariate_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "trt", treatment.var = "treatment")
  
  out$site_block<-siteblock[i]
  
  gex_multdiff<-rbind(gex_multdiff, out)  
}

gex_multdiff_all<-gex_multdiff %>% 
  separate(site_block, into=c("site", "block"), sep="##")%>%
  select(site, block, composition_diff)

write.csv(gex_multdiff_all, "gex_multdiff_site_block.csv", row.names = F)

gex_multdiff_ave<-gex_multdiff_all%>%
  group_by(site) %>% 
  summarize(composition_diff=mean(composition_diff))

write.csv(gex_multdiff_ave, "gex_multdiff_site_ave.csv", row.names = F)


###doing RAC differences

gex_RACdiff<-data.frame()
siteblock<-unique(lsyear2$site_block)

for (i in 1:length(siteblock)){
  
  subset<-lsyear2%>%
    filter(site_block==siteblock[i]) %>%
    mutate(treatment=trt)
  
  out <- RAC_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "trt", treatment.var = "treatment")
  
  out$site_block<-siteblock[i]
  
  gex_RACdiff<-rbind(gex_RACdiff, out)  
}


gex_RACdiff_all<-gex_RACdiff %>% 
  separate(site_block, into=c("site", "block"), sep="##")%>%
  select(site, block, richness_diff, evenness_diff, rank_diff, species_diff)

write.csv(gex_RACdiff_all, "gex_RACdiff_site_block.csv", row.names = F)

gex_RACdiff_ave<-gex_RACdiff_all%>%
  group_by(site) %>% 
  summarize_at(vars(richness_diff, evenness_diff, rank_diff, species_diff), funs(mean))

write.csv(gex_RACdiff_ave, "gex_RACdiff_site_ave.csv", row.names = F)

###doing abundance differences
gex_abunddiff<-data.frame()
siteblock<-unique(lsyear2$site_block)

for (i in 1:length(siteblock)){
  
  subset<-lsyear2%>%
    filter(site_block==siteblock[i]) %>%
    mutate(treatment=trt)
  
  out <- abundance_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "trt", treatment.var = "treatment")
  
  out$site_block<-siteblock[i]
  
  gex_abunddiff<-rbind(gex_abunddiff, out)  
}

gex_abunddiff_all<-gex_abunddiff %>% 
  separate(site_block, into=c("site", "block"), sep="##")%>%
  select(site, block, genus_species, difference)

write.csv(gex_abunddiff_all, "gex_abund_diff_site_block.csv", row.names = F)

gex_abunddiff_ave<-gex_abunddiff_all%>%
  group_by(site, genus_species) %>% 
  summarize_at(vars(difference), funs(mean))

write.csv(gex_abunddiff_ave, "gex_abund_diff_site_ave.csv", row.names = F)

##not sure what this is for anymore.
###getting clean dataset
lsyear_export<-dat%>%
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage)%>%
  mutate(site_block=paste(site, block, sep="##")) %>%
  mutate(site_block_plot=paste(site_block, plot, sep="##"))

lsyear_export2<-lsyear_export%>%
  filter(site_block!="California_Sedgwick_Lisque##A"&
           site_block!="California_Sedgwick_Lisque##B"&
           site_block!="California_Sedgwick_Lisque##G"&
           site_block!="DesertLow##Alkali"&
           site_block_plot!="Jornada##east##30"&
           site_block_plot!="Jornada##west##18"&
           site_block_plot!="Jornada##west##19"&
           site_block_plot!="Jornada##west##20"&
           site_block_plot!="Jornada##west##21"&
           site_block_plot!="Jornada##west##22")%>%
  separate(site_block_plot, into=c("site","block","plot"), sep="##")%>%
  select(-X, -lyear, -mexage, -site_block)%>%
  select(site, block, plot, exage, year, trt, genus_species, relcov)
write.csv(lsyear_export2, "All_Cleaned_April2019_V2.csv", row.names=F)  


###looking at background variation in ungrazed and grazed comparisons
#first subset to sites with multiple plots
multreps<-lsyear%>%
  select(site, block)%>%
  unique()%>%
  group_by(site)%>%
  summarize(n=length(block))%>%
  filter(n>1)

lsyear_subset<-lsyear%>%
  right_join(multreps)

graze<-lsyear_subset%>%
  filter(site_block!="California_Sedgwick_Lisque##A"&
           site_block!="California_Sedgwick_Lisque##B"&
           site_block!="California_Sedgwick_Lisque##G"&
           site_block!="DesertLow##Alkali"&
           site!="Jornada"&
           site!="Junner Koeland"&
           site!="Konza")%>%
  filter(trt=="G")%>%
  separate(site_block, into=c('site', 'block'), sep="##")

comp_graze<-data.frame()
sites<-unique(graze$site)

for (i in 1:length(sites)){
  
  subset<-graze%>%
    filter(site==sites[i]) %>%
    mutate(block2=block)
  
  out <- multivariate_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "block", treatment.var = "block2")
  
  out$site<-sites[i]
  
  comp_graze<-rbind(comp_graze, out)  
}

ungraze<-lsyear_subset%>%
  filter(site_block!="California_Sedgwick_Lisque##A"&
           site_block!="California_Sedgwick_Lisque##B"&
           site_block!="California_Sedgwick_Lisque##G"&
           site_block!="DesertLow##Alkali"&
           site!="Jornada"&
           site!="Junner Koeland"&
           site!="Konza")%>%
  filter(trt=="U")%>%
  separate(site_block, into=c('site', 'block'), sep="##")

comp_ungraze<-data.frame()
sites<-unique(ungraze$site)

for (i in 1:length(sites)){
  
  subset<-ungraze%>%
    filter(site==sites[i]) %>%
    mutate(block2=block)
  
  out <- multivariate_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "block", treatment.var = "block2")
  
  out$site<-sites[i]
  
  comp_ungraze<-rbind(comp_ungraze, out)  
}

mean(comp_ungraze$composition_diff)#0.6005785
mean(comp_graze$composition_diff)#0.5654668