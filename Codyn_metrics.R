library(tidyverse)
library(codyn)

setwd("~/Dropbox/Dominance WG/")

dat<-read.csv("ALL_Cleaned_Apr2019.csv")

lsyear<-dat%>%
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage)%>%
  mutate(site_block=paste(site, block, sep="##"))

lsyear2<-lsyear%>%
  filter(site_block!="California_Sedgwick_Lisque_A"&
           site_block!="California_Sedgwick_Lisque_B"&
           site_block!="California_Sedgwick_Lisque_G"&
           site_block!="DesertLow_Alkali"&
           site!="Jornada"&
           site!="Junner Koeland"&
           site!="Konza")

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

##three sites, konza, jornada and junner koeland, had plots nested in block. need to use that informaiton here.

lsyear3<-lsyear%>%
  filter(site=="Jornada"|
         site=="Junner Koeland"|
         site=="Konza")%>%
  mutate(site_block_plot=paste(site_block, plot, sep="##"))%>%
  filter(site_block_plot!="Jornada_east_30"&
           site_block_plot!="Jornada_west_18"&
           site_block_plot!="Jornada_west_19"&
           site_block_plot!="Jornada_west_20"&
           site_block_plot!="Jornada_west_21"&
           site_block_plot!="Jornada_west_22")

gex_multdiff2<-data.frame()
siteblockplot<-unique(lsyear3$site_block_plot)

for (i in 1:length(siteblockplot)){
  
  subset<-lsyear3%>%
    filter(site_block_plot==siteblockplot[i]) %>%
    mutate(treatment=trt)
  
  out <- multivariate_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "trt", treatment.var = "treatment")
  
  un_site_block<-unique(subset$site_block)
  un_plot<-unique(subset$plot)
  
  out$site_block<-un_site_block
  out$plot<-un_plot
  
  gex_multdiff2<-rbind(gex_multdiff2, out)  
}

ave_gexmultdiff2<-gex_multdiff2%>%
  group_by(site_block, treatment, treatment2, trt_greater_disp, abs_dispersion_diff)%>%
  summarize(composition_diff=mean(composition_diff))

gex_multdiff_all<-gex_multdiff %>% 
  bind_rows(ave_gexmultdiff2)

gex_multdiff_ave<-gex_multdiff_all%>%
  separate(site_block, into=c("site", "block", sep="##"))
