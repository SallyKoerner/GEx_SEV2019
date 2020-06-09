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