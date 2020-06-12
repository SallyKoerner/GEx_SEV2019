### Read trait data from BIEN and TRY -- clean and merge
###
### Authors: kevin wilcox (kevin.wilcox@uwyo.edu)
### created: April 22, 2019;  last updated: Aug 13, 2019; 

file.choose()
### Set up workspace
rm(list=ls())
setwd("C:\\Users\\wilco\\Dropbox\\Cross_workstation_workspace\\Working groups\\GEx\\trait data\\11June2019Data\\")
library(tidyverse)
library(data.table)
library(readr)

### Read in data
try_data <- read.delim("6552_tryRaw_11June2019.txt")
bien_data <- read.csv ("..\\GExBIENtraits_Aug2019\\GExBIENtraits\\BIENtrait_GEx_July2019.csv")

trait_id_key <- read.csv("..\\GExBIENtraits_Aug2019\\GExBIENtraits\\BIEN_TRY_trait_key.csv")
spname_key <- read.csv("..\\GEx_species_family_final.csv")

try_data_clean1 <- try_data %>%
  filter(!is.na(TraitID))

write.csv(try_data_clean1, file="try_data_traitsOnly.csv", row.names=F)

## Remove potential erroneous data
try_data_clean2 <- try_data_clean1 %>%
  filter(ErrorRisk<4) %>%
  dplyr::select(AccSpeciesName,
                TraitName,
                TraitID,
                StdValue,
                UnitName,
                ObservationID) %>%
  rename( trait_name = TraitName,
          trait_value = StdValue,
          unit_name = UnitName) %>%
  mutate(trait_value = as.character(trait_value)) %>% ### for merging purposes
  left_join(dplyr::select(trait_id_key, TRY_trait_name, KW_trait_ID), by=c("trait_name"="TRY_trait_name")) %>%
  mutate(KW_trait_ID= coalesce(KW_trait_ID, TraitID)) ### fill KWTraitIDs with values from TRY traitID column

head(filter(try_data_clean2, trait_name == "Leaf length"))
tail(try_data_clean2)
missing <- filter(try_data_clean2, is.na(clean_ejf))

bien_data_clean1 <- bien_data %>%
  mutate(TraitID=NA) %>%
  dplyr::select(scrubbed_species_binomial,
                trait_name,
                TraitID,
                trait_value,
                unit,
                id) %>%
  rename(AccSpeciesName = scrubbed_species_binomial,
         unit_name = unit,
         ObservationID = id) %>%
  mutate(trait_value = replace(trait_value, trait_value=="*", NA)) %>%
  mutate(trait_value = replace(trait_value, trait_value=="<1", 0)) %>%
  mutate(trait_value = replace(trait_value, trait_value=="<5", 0)) %>%
  mutate(trait_value = replace(trait_value, trait_value==">100", 100)) %>%
  mutate(trait_value = replace(trait_value, trait_value==">10000", NA)) %>%
  mutate(trait_value = as.character(trait_value)) %>% ### for merging purposes%>%
  filter(!is.na(trait_name)) %>%
  left_join(dplyr::select(trait_id_key, BIEN_trait_name, KW_trait_ID), by=c("trait_name"="BIEN_trait_name"))

filter(bien_data_clean1, is.na(KW_trait_ID))
head(as.character(bien_data_clean1$trait_value))
head(bien_data_clean1$trait_value)
try_bien_data <- try_data_clean2 %>%
  bind_rows(bien_data_clean1)

species_names <- data.frame(x=1, sp_names=unique(try_bien_data$AccSpeciesName)) %>%
  arrange(sp_names)
  
trait_names <- unique(try_data_clean1$TraitName) 

bien_data_clean1$trait_value <- factor(bien_data_clean1$trait_value)
head(levels(bien_data_clean1$trait_value))

write.csv(trait_names, "trait names.csv", row.names=F)

trait_matrix <- try_data_clean2 %>%
  dplyr::select(AccSpeciesName,TraitID, StdValue) %>%
  mutate(TraitID=paste0("trt",TraitID)) %>%
  group_by(AccSpeciesName, TraitID) %>%
  summarise(num_record = length(StdValue)) %>%
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




