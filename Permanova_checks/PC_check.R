library(tidyverse)
library(codyn)
library(rsq)
library(vegan)
library(prcomp)

setwd("/Users/avahoffman/Dropbox/Research/Grazing_consortium/2020/GEx_SEV2019")

dat<-read.csv("All_Cleaned_April2019_V2.csv")

lsyear<-dat%>%
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage)%>%
  mutate(block = as.numeric(block))%>%
  mutate(site_block=paste(site, block, sep="##"))

lsyear2<-
  lsyear%>%
  filter(site=="Argentina_ElPalmar")

wide<-
  lsyear2%>%
  spread(genus_species, relcov, fill=0)

species<-
  as.data.frame(wide[,12:ncol(wide)])

env<-
  as.data.frame(wide[,1:11])

pc1 <- 
  prcomp(t(species), scale. = T, center = T)$rotation[,1]

data.in <- 
  cbind(env,pc1)
