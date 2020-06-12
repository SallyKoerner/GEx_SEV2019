### Merge trait data with GEx
### Authors: Kevin Wilcox (kevin.wilcox@uwyo.edu)
###
### Created April 24, 2019; Last updated: April 24, 2019

### Set up workspace
rm(list=ls())
setwd("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\trait data\\")
library(tidyverse)

### Import data
relcov_full <- read.csv("..\\analyses\\cleaning_data\\ALL_Cleaned_Apr2019.csv")
trait_data <- read.csv("BIEN_traitFINAL.csv")
name_key <- read.csv("GEx_species_family_final.csv")

### compare trait names
spnames_from_trait_data <- data.frame(spname_trt = unique(trait_data$scrubbed_species_binomial))
name_merger <- spnames_from_trait_data %>%
  full_join(name_key, by=c("spname_trt" = "clean_ejf"))

write.csv(name_merger, file="matches between clean_ejf and spnames in trait data.csv", row.names=F)
  name_key %>%
  full_join(dplyr::select(trait_data, scrubbed_species_binomial, X))

### Merge TNRS names with names in GEx dataset
relcov_tnrs_names <- relcov_full %>%
  left_join(spname_key, by=) %>%
  dplyr::select(-genus_species)




### Select only last year of data and remove problem plots
relcov_last_year_only <- relcov_tnrs_names %>%
  group_by(site) %>%
  filter(year == max(year), exage == max(exage)) %>%
  filter(site !="California_Sedgwick_Lisque" | !block %in% c("A","B","G")) %>%
  filter(site !="DesertLow" | block != "Alkali") %>%
  filter(site !="Jornada" | block != "west" | !plot %in% 18:22) %>%
  filter(site !="Jornada" | block != "east" | plot != 30)

relcov_last_year_only <- relcov_full %>%
  group_by(site) %>%
  filter(year == max(year), exage == max(exage)) %>%
  filter(site !="California_Sedgwick_Lisque" | !block %in% c("A","B","G")) %>%
  filter(site !="DesertLow" | block != "Alkali") %>%
  filter(site !="Jornada" | block != "west" | !plot %in% 18:22) %>%
  filter(site !="Jornada" | block != "east" | plot != 30)



try_trait_data <- read.csv("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\trait data\\try_data_traitsOnly.csv",
                           stringsAsFactors=F)

try_lengths <- try_trait_data %>%
  mutate(trait_value_numeric=as.numeric(OrigValueStr)) %>%
  dplyr::select(AccSpeciesName, TraitName, OrigValueStr, trait_value_numeric) %>%
  group_by(AccSpeciesName, TraitName) %>%
  summarize(num_records = length(trait_value_numeric)) 

write.csv(try_lengths, file="number of records per species and trait_TRY.csv")


