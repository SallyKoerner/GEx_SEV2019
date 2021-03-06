### Merge species names to get tnrs names and select last year of data and greatest exage only
### Authors: Kevin Wilcox (kevin.wilcox@uwyo.edu)
###
### Created April 22, 2019; Last updated: April 24, 2019

### Set up workspace
rm(list=ls())
setwd("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\analyses\\cleaning_data\\")
library(tidyverse)

### Import data
relcov_full <- read.csv("ALL_Cleaned_Apr2019.csv")
spname_key <- read.csv()

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

