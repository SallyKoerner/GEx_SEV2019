### Calculate trait coverage per plot
###
### Author: Wilcox (kevin.wilcox@uwyo.edu)
### Created: April 25, 2019; updated: April 25, 2019

###
### Set up workspace
###
rm(list=ls())
library(tidyverse)
setwd("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\trait data\\")

relcov_data <- read.csv("GEx relcov_last year only_prob blocks removed.csv") %>%
  group_by(site, year, exage, plot, block, trt, clean_ejf) %>%
  summarize(relcov=mean(relcov, na.rm=F))

### Prepare trait data
trait_data_all <- read.csv("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\Merged trait database 5 traits Means2.csv", stringsAsFactors=F)
trait_data_try <- read.csv("try_data_traitsOnly.csv",  stringsAsFactors=F)
trait_data_all$

trait_data_try_short <- trait_data_try %>%
  filter(TraitID %in% c(3115, 14, 3106)) %>%
  mutate(TraitID = paste0("Trait",TraitID)) %>%
  dplyr::select(SpeciesName, TraitID, OrigValueStr) %>%
  mutate(trait_value = as.numeric(OrigValueStr)) %>%
  group_by(SpeciesName, TraitID) %>%
  summarise(mean_trait = mean(trait_value, na.rm=T)) %>%
  spread(key=TraitID, value=mean_trait)

n_trait_data <- trait_data_try %>%
  filter(TraitID==14) %>%
  mutate(unit = factor(UnitName))
levels(n_trait_data$unit)

### Combine GEx with trait means
relcov_traits <- relcov_data %>%
  left_join(trait_data_try_short, by=c("clean_ejf" = "SpeciesName"))

### Calculate trait_represented_cover
trait_cover_by_site <- relcov_traits %>%
  mutate(leafN_presence = ifelse(is.na(Trait14),0,1)) %>%
  mutate(height_presence = ifelse(is.na(Trait3106),0,1)) %>%
  mutate(SLA_presence = ifelse(is.na(Trait3115),0,1)) %>%
  mutate(leafN_relcov = leafN_presence*relcov) %>%
  mutate(height_relcov = height_presence*relcov) %>%
  mutate(SLA_relcov = SLA_presence*relcov) %>%
  group_by(site, year, exage, plot, block, trt) %>%
  summarize(leafN_site_abun = sum(leafN_relcov, na.rm=T),
            SLA_site_abun = sum(SLA_relcov, na.rm=T),
            height_site_abun = sum(height_relcov, na.rm=T)) %>%
  ungroup() %>%
  group_by(site) %>%
  summarize(leafN_site_abun = mean(leafN_site_abun, na.rm=T),
            SLA_site_abun = mean(SLA_site_abun, na.rm=T),
            height_site_abun = mean(height_site_abun, na.rm=T))

#ggplot(trait_cover_by_site, aes(x=1, y=trait_cover_by_plot$leafN_site_abun)

ggplot(trait_cover_by_site, aes(x=reorder(site,-leafN_site_abun), y=leafN_site_abun)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6))

ggplot(trait_cover_by_site, aes(x=reorder(site,-SLA_site_abun), y=SLA_site_abun)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6))

ggplot(trait_cover_by_site, aes(x=reorder(site,-height_site_abun), y=height_site_abun)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6))



z <- filter(trait_cover_by_site, leafN_site_abun>150)
zz <- filter(trait_cover_by_site, site=="Kruger_Satara")
