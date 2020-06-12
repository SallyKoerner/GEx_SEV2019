####Code by M. Avolio June 2020 for virtual working group.
###this calculates log rr of berger parker differnces between grazed and ungrazed paired plots. 
##also determines whether there is a switch in the dom sp.


library(tidyverse)
library(codyn)
library(rsq)

theme_set(theme_bw(12))


###Meghan's Working Directory
setwd("C:\\Users\\mavolio2\\Dropbox\\GEx_VirtualWorkshop_June2020\\")
##Sally's Working Directory
setwd("~/Dropbox/GEx_VirtualWorkshop_June2020")


dat<-read.csv("GEx_cleaned_11June2020.csv")%>%
  mutate(site_block=paste(site, block, sep="##")) %>% 
  mutate(drop=ifelse(genus_species_use=="#N/A"|genus_species_use=="Dead unidentified"|genus_species_use=="Leaf.Litter"|genus_species_use=="cactus__dead_", 1, 0))%>%
  filter(drop!=1)

###all years of data
dat2<-dat%>%
  filter(site_block!="California_Sedgwick_Lisque##A"&
           site_block!="California_Sedgwick_Lisque##B"&
           site_block!="California_Sedgwick_Lisque##G"&
           site_block!="DesertLow##Alkali"&
           site!="Jornada"&
           site!="Junner Koeland"&
           site!="Konza")%>%
  select(-plot)


dat3<-dat%>%
  filter(site=="Jornada"|
           site=="Junner Koeland"|
           site=="Konza")%>%
  mutate(site_block_plot=paste(site_block, plot, sep="##"))%>%
  filter(site_block_plot!="Jornada##east##30"&
           site_block_plot!="Jornada##west##18"&
           site_block_plot!="Jornada##west##19"&
           site_block_plot!="Jornada##west##20"&
           site_block_plot!="Jornada##west##21"&
           site_block_plot!="Jornada##west##22")%>%
  group_by(site, year, exage, block, trt, genus_species_use, site_block)%>%
  summarise(relcov=mean(relcov))

dat4<-dat2%>%
  bind_rows(dat3)%>%
  group_by(site, block, year, trt)%>%
  mutate(dom=max(relcov))%>%
  filter(dom==relcov)

#we are dropping sites that have multiple dominant species in a block
keep<-dat4%>%
  group_by(site, block, year, trt)%>%
  summarize(n=length(relcov))%>%
  filter(n<2)

dat5<-dat4%>%
  right_join(keep)%>%
  select(-site_block, -dom, -n)

graze<-dat5%>%
  ungroup()%>%
  filter(trt=="G")%>%
  rename(graze_sp=genus_species_use,
         graze_bp=relcov)%>%
  select(-trt, -genus_species, - genus_species_clean)

ungraze<-dat5%>%
  ungroup()%>%
  filter(trt=="U")%>%
  rename(ungraze_sp=genus_species_use,
         ungraze_bp=relcov)%>%
  select(-trt, -genus_species, - genus_species_clean)

diffG_U<-ungraze%>%
  left_join(graze)%>%
  na.omit%>%
  mutate(bp_rr=log(graze_bp/ungraze_bp),
         diff_sp=ifelse(graze_sp==ungraze_sp, 0, 1))

##this is the file that beth and ben want
write.csv(diffG_U, "Diff_BP_Dom_allyrs.csv")


