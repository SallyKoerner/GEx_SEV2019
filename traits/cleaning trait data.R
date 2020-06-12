### Read trait data 
###
### Authors: kevin wilcox (kevin.wilcox@uwyo.edu)
### created: April 22, 2019;  last updated: April 22, 2019; 

file.choose()
### Set up workspace
rm(list=ls())
setwd("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\trait data\\")
library(tidyverse)
library(data.table)
library(readr)

### REad in data
try_data <- read.delim("try data_raw.txt")

try_data_clean1 <- try_data %>%
  filter(!is.na(TraitID))

write.csv(try_data_clean1, file="try_data_traitsOnly.csv", row.names=F)

trait_names <- unique(try_data_clean1$TraitName)

write.csv(trait_names, "trait names.csv", row.names=F)

trait_matrix <- try_data_clean1 %>%
  dplyr::select(AccSpeciesName,TraitID, OrigValueStr) %>%
  mutate(TraitID=paste0("trt",TraitID)) %>%
  group_by(AccSpeciesName, TraitID) %>%
  summarise(num_record = length(OrigValueStr)) %>%
  spread(key=TraitID, value=num_record)

trait_presence_counter <- trait_matrix %>%
  gather(key=TraitID, value=num_record, -AccSpeciesName) %>%
  group_by(TraitID) %>%
  summarise(num_present = length(which(!is.na(num_record))),
            tot_species = length(num_record),
            mean_num = mean(num_record, na.rm=T)) %>%
  ungroup() %>%
  mutate(pTrait_presence = num_present/tot_species) 

trait_presence_species <- trait_matrix %>%
  gather(key=TraitID, value=num_record, -AccSpeciesName) %>%
  group_by(AccSpeciesName) %>%
  summarise(num_present = length(which(!is.na(num_record))),
            tot_traits = length(num_record),
            mean_num = mean(num_record, na.rm=T)) %>%
  ungroup() %>%
  mutate(pSpecies_presence = num_present/tot_traits) 


ggplot(trait_presence_counter, aes(x=reorder(TraitID, -pTrait_presence), y=pTrait_presence)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=0))

ggplot(filter(trait_presence_counter, pTrait_presence>.18), aes(x=reorder(TraitID, -pTrait_presence), y=pTrait_presence)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=0))

ggplot(filter(trait_presence_species, pSpecies_presence >.3), aes(x=reorder(AccSpeciesName, -pSpecies_presence), y=pSpecies_presence)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=0))


###
### Subset dominant species only
###

dom_sp <- read.csv("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\analyses\\Getting trait data\\GEx_species_family.csv")

dom_trait_matrix <- trait_matrix %>%
  inner_join(dom_sp, by=c("AccSpeciesName" = "clean_ejf"))

dom_trait




###
### BIEN data
###

file.choose()
setwd("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\trait data\\")

trait_data <- read.csv("GEx_species_BEIN_trait.csv")
sp_names_gex <- read.csv("GEx_species_family.csv")

### Get trait list and write to file
trait_list_bien <- unique(trait_data$trait_name)
# write.csv(trait_list_bien, "BIEN_traitList.csv", row.names=F)

trait_matrix <- try_data_clean1 %>%
  dplyr::select(AccSpeciesName,TraitID, OrigValueStr) %>%
  mutate(TraitID=paste0("trt",TraitID)) %>%
  group_by(AccSpeciesName, TraitID) %>%
  summarise(num_record = length(OrigValueStr)) %>%
  spread(key=TraitID, value=num_record)

trait_presence_counter <- trait_matrix %>%
  gather(key=TraitID, value=num_record, -AccSpeciesName) %>%
  group_by(TraitID) %>%
  summarise(num_present = length(which(!is.na(num_record))),
            tot_species = length(num_record),
            mean_num = mean(num_record, na.rm=T)) %>%
  ungroup() %>%
  mutate(pTrait_presence = num_present/tot_species) 

trait_presence_species <- trait_matrix %>%
  gather(key=TraitID, value=num_record, -AccSpeciesName) %>%
  group_by(AccSpeciesName) %>%
  summarise(num_present = length(which(!is.na(num_record))),
            tot_traits = length(num_record),
            mean_num = mean(num_record, na.rm=T)) %>%
  ungroup() %>%
  mutate(pSpecies_presence = num_present/tot_traits) 




