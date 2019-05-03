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

relcov_trait <- read.csv("GEx_All_Species_Trait.csv") %>%
  dplyr::select(site, block, plot, trt, clean_ejf, relcov, 
                Leaf_Area_Mean, Leaf_N_Mass_Mean, Lifespan_Mean, Height_Mean, SLA_Mean)

relcov_data <- read.csv("GEx relcov_last year only_prob blocks removed.csv") %>%
  group_by(site, year, exage, plot, block, trt, clean_ejf) %>%
  summarize(relcov=mean(relcov, na.rm=F))

### Prepare trait data
trait_data_all <- read.csv("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\Merged trait database 5 traits Means2.csv", stringsAsFactors=F) %>%
  dplyr::rename(SpeciesName = ï..SpeciesName)
trait_data_try <- read.csv("try_data_traitsOnly.csv",  stringsAsFactors=F)
sp_names <- data.frame(SpeciesName=unique(trait_data_all$SpeciesName), trait_dummy=1)
names(trait_data_all)

zz <- data.frame(spname=sort(unique(trait_data_try$AccSpeciesName)))

spname_key <- read.csv("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\trait data\\GEx_species_family_final.csv") %>%
  dplyr::select(clean_ejf) %>%
  mutate(key_dummy=1)

spname_merged <- sp_names %>%
  full_join(spname_key, by=c("SpeciesName" = "clean_ejf"))

write.csv(spname_merged, file="Merged aggregate species names and spname key.csv", row.names=F)

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
trait_cover_by_site <- relcov_trait %>%
  mutate(Leaf_N_presence = ifelse(is.na(Leaf_N_Mass_Mean),0,1)) %>%
  mutate(Height_presence = ifelse(is.na(Height_Mean),0,1)) %>%
  mutate(SLA_presence = ifelse(is.na(SLA_Mean),0,1)) %>%
  mutate(Leaf_Area_presence = ifelse(is.na(Leaf_Area_Mean),0,1)) %>%
  mutate(Lifespan_presence = ifelse(is.na(Lifespan_Mean),0,1)) %>%
  mutate(Leaf_N_relcov = Leaf_N_presence*relcov) %>%
  mutate(Height_relcov = Height_presence*relcov) %>%
  mutate(SLA_relcov = SLA_presence*relcov) %>%
  mutate(Leaf_Area_relcov = Leaf_Area_presence*relcov) %>%
  mutate(Lifespan_relcov = Lifespan_presence*relcov) %>%
  group_by(site, plot, block, trt) %>%
  summarize(Leaf_N_site_abun = sum(Leaf_N_relcov, na.rm=T),
            Leaf_Area_site_abun = sum(Leaf_Area_relcov, na.rm=T),
            Lifespan_site_abun = sum(Lifespan_relcov, na.rm=T),
            SLA_site_abun = sum(SLA_relcov, na.rm=T),
            Height_site_abun = sum(Height_relcov, na.rm=T)) %>%
  ungroup() %>%
  group_by(site) %>%
  summarize(Leaf_N_site_abun = mean(Leaf_N_site_abun, na.rm=T),
            SLA_site_abun = mean(SLA_site_abun, na.rm=T),
            Leaf_Area_site_abun = mean(Leaf_Area_site_abun, na.rm=T),
            Lifespan_site_abun = mean(Lifespan_site_abun, na.rm=T),
            Height_site_abun = mean(Height_site_abun, na.rm=T))

#ggplot(trait_cover_by_site, aes(x=1, y=trait_cover_by_plot$leafN_site_abun)

ggplot(trait_cover_by_site, aes(x=reorder(site,-Leaf_N_site_abun), y=Leaf_N_site_abun)) +
  geom_bar(stat="identity") +
  theme_classic() +
  geom_hline(yintercept=50, col="red", size=2) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6))
ggsave("figures\\leafN site abundance.png", width=10, height=4)

ggplot(trait_cover_by_site, aes(x=reorder(site,-SLA_site_abun), y=SLA_site_abun)) +
  geom_bar(stat="identity") +
  theme_classic() +
  geom_hline(yintercept=50, col="red", size=2) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6))
ggsave("figures\\SLA site abundance.png", width=10, height=4)

ggplot(trait_cover_by_site, aes(x=reorder(site,-Height_site_abun), y=Height_site_abun)) +
  geom_bar(stat="identity") +
  theme_classic() +
  geom_hline(yintercept=50, col="red", size=2) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6))
ggsave("figures\\Height site abundance.png", width=10, height=4)

ggplot(trait_cover_by_site, aes(x=reorder(site,-Leaf_Area_site_abun), y=Leaf_Area_site_abun)) +
  geom_bar(stat="identity") +
  theme_classic() +
  geom_hline(yintercept=50, col="red", size=2) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6))
ggsave("figures\\leafArea site abundance.png", width=10, height=4)

ggplot(trait_cover_by_site, aes(x=reorder(site,-Lifespan_site_abun), y=Lifespan_site_abun)) +
  geom_bar(stat="identity") +
  theme_classic() +
  geom_hline(yintercept=50, col="red", size=2) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6))
ggsave("figures\\lifespan site abundance.png", width=10, height=4)


###
###
### 





