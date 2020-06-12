### Comparing 
### Authors: Kevin Wilcox (kevin.wilcox@uwyo.edu)
###
### Created April 24, 2019; Last updated: April 24, 2019

### Set up workspace
rm(list=ls())
setwd("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\trait data\\")
library(tidyverse)

### Import data
dom_trait_relcov <- read.csv("GEx_Dominant_Species_Trait_environmental.csv") %>%
  dplyr::select(site, block, plot, trt, clean_ejf, relcov, 
                Leaf_Area_Mean, Leaf_N_Mass_Mean, Lifespan_Mean, Height_Mean, SLA_Mean) %>%
  gather(key=trait_name, value=trait_value, -site:-relcov) %>%
  spread(key=trt, value=relcov) %>%
  mutate(U = replace(U, is.na(U), 0)) %>%
  mutate(G = replace(G, is.na(G), 0)) %>%
  mutate(GrResp_subtract = G-U,
         GrResp_pCh = abs(G-U)/U,
         GrResp_lnRR = log(G/U))

write.csv(dom_trait_relcov, file="GEx response ratios_cleanedUnits.csv", row.names=F)


trait_relcov_all <- read.csv("GEx_All_Sites_species_trait_environmental.csv")

site_chars <- trait_relcov_all %>%
  dplyr::select(site, CVpercBA:herb_map_2)

?read.csv
relcov_matrix <- trait_relcov_all %>%
  mutate(plot_id = paste(site, year, block, plot, trt, sep="_")) %>%
  dplyr::select(plot_id, clean_ejf, relcov) %>%
  group_by(plot_id, clean_ejf) %>%
  summarize(relcov=sum(relcov,na.rm=T)) %>%
  ungroup() %>%
  spread(key=clean_ejf, value=relcov)

relcov_matrix[c(3853,3860),]
trait_matrix <- trait_relcov_all %>%
  dplyr::select(clean_ejf,  Leaf_Area_Mean, Leaf_N_Mass_Mean, Lifespan_Mean, Height_Mean, SLA_Mean) %>%
  group_by(clean_ejf) %>%
  summarize(Leaf_Area_Mean = mean(Leaf_Area_Mean, na.rm=T),
            Leaf_N_Mass_Mean = mean(Leaf_N_Mass_Mean, na.rm=T),
            Height_Mean = mean(Height_Mean, na.rm=T),
            SLA_Mean = mean(SLA_Mean, na.rm=T))

write.csv(relcov_matrix, file="relcov matrix_forRobert.csv", row.names=F)
write.csv(trait_matrix, file="trait matrix_forRobert.csv", row.names=F)
### Clean up relcov data
relcov_last_year_only <- relcov_full %>%
  group_by(site) %>%
  filter(year == max(year), exage == max(exage)) %>%
  filter(site !="California_Sedgwick_Lisque" | !block %in% c("A","B","G")) %>%
  filter(site !="DesertLow" | block != "Alkali") %>%
  filter(site !="Jornada" | block != "west" | !plot %in% 18:22) %>%
  filter(site !="Jornada" | block != "east" | plot != 30)


###
### grazing impacts on community weighted means across sites
###

trait_comm_means <- read.csv("") 

trait_comm_responses %>% trait_comm_means

relcov_responses <- relcov_last_year_only %>%
  dplyr::select(-X, -X.1) %>%
  spread(key=trt, value=relcov) %>%
  mutate(U = replace(U, is.na(U), 0)) %>%
  mutate(G = replace(G, is.na(G), 0)) %>%
  mutate(GrResp_subtract = G-U,
         GrResp_pCh = abs(G-U)/U,
         GrResp_lnRR = log(G/U))
  

# trait_data_bien <- read.csv("BIEN_traitFINAL.csv") %>%
#   mutate(trait_name = replace(trait_name, trait_name=="leaf dry mass per leaf fresh mass","Leaf water content per leaf dry mass (not saturated)")) %>%
#   mutate(trait_name = replace(trait_name, trait_name=="leaf area per leaf dry mass","Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded")) %>%
#   mutate(trait_name = replace(trait_name, trait_name=="leaf nitrogen content per leaf dry mass","Leaf nitrogen (N) content per leaf dry mass")) %>%
#   mutate(trait_name = replace(trait_name, trait_name=="whole plant height","Plant height vegetative"))
  
trait_data_try <- read.csv("try_data_traitsOnly.csv",  stringsAsFactors=F)
#name_key <- read.csv("GEx_species_family_final.csv")

### Create common namdes for LDMC, SLA, Height, Leaf N across BIEN and TRY datasets
head(filter(trait_data_try, TraitID==3115))
trait_data_try_short <- trait_data_try %>%
  filter(TraitID %in% c(3115, 14, 3106)) %>%
  dplyr::select(SpeciesName, TraitName, TraitID, StdValue) %>%
  mutate(trait_value = as.numeric(StdValue)) %>%
  group_by(SpeciesName, TraitName, TraitID) %>%
  summarise(mean_trait = mean(trait_value, na.rm=T))

# LDMC = 3120, SLA (petiole excluded) 3115, lf N / dry mass = 14, plant height vegetative = 3106

filter(trait_data_bien, is.na(trait_name))

write.csv(relcov_last_year_only, file="GEx relcov_last year only_prob blocks removed.csv", row.names=F)


### 



### Combine trait data with GEx
relcov_traits <- relcov_responses %>%
  left_join(trait_data_try_short, by=c("clean_ejf" = "SpeciesName"))

###
levels(factor(relcov_traits$TraitName))
ggplot(filter(relcov_traits,TraitName=="Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded"),
       aes(x=mean_trait, y=GrResp_))

ggplot(dom_trait_relcov, aes(x=trait_value, y=GrResp_subtract)) +
  geom_point()+
  geom_smooth(method="lm") +
  facet_wrap(~trait_name, scales="free") +

ggplot(dom_trait_relcov, aes(x=trait_value, y=GrResp_lnRR)) +
  geom_point()+
  geom_smooth(method="lm") +
  facet_wrap(~trait_name, scales="free")

ggplot(relcov_traits, aes(x=mean_trait, y=GrResp_pCh)) +
  geom_point()+
  geom_smooth(method="lm") +
  facet_wrap(~TraitName, scales="free")
ggplot(relcov_traits, aes(x=mean_trait, y=GrResp_lnRR)) +
  geom_point()+
  geom_smooth(method="lm") +
  ylim(-4.5,4.5) +
  facet_wrap(~TraitName, scales="free")

ggplot(filter(dom_trait_relcov, trait_name=="Leaf_N_Mass_Mean"), aes(x=trait_value, y=GrResp_subtract)) +
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic() +
#  xlim(0,5)


boxplot(relcov_traits$mean_trait)

head(natdb-dump-GEx-2019-04-18)
natdb -> natdb-dump-GEx-2019-04-18
str(natdb-dump-GEx-2019-04-18)
write.csv()
file.choose()
natdb <- readRDS("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\trait data\\NATDB-20190425T172729Z-001\\NATDB\\natdb-dump-GEx-2019-04-18.RDS")
str(natdb)

head(natdb$categorical)
head(natdb$numeric)

write.csv(natdb$categorical, "natdb_data.csv", row.names=F)
write.csv(natdb$numeric, "natdb_data_2.csv", row.names=F)


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




try_trait_data <- read.csv("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\trait data\\try_data_traitsOnly.csv",
                           stringsAsFactors=F)

try_lengths <- try_trait_data %>%
  mutate(trait_value_numeric=as.numeric(OrigValueStr)) %>%
  dplyr::select(AccSpeciesName, TraitName, OrigValueStr, trait_value_numeric) %>%
  group_by(AccSpeciesName, TraitName) %>%
  summarize(num_records = length(trait_value_numeric)) 

write.csv(try_lengths, file="number of records per species and trait_TRY.csv")


