setwd("C://Users/Jason Donaldson/Documents/Workshops and conferences/GEx/")
library(dplyr)
##import plot level codominant species data and add region to datafile
GEx_codominant_species<-read.csv("GEx_codominane_04242019.csv")
GEx_metadata<-read.csv("Meta_SEV2019_v2.csv")

##Create clean species name file and put clean names into dominance file
GEx_raw_with_names<-read.csv("GEx_cleaned_v3.csv")
GEx_unique_species<-subset(GEx_raw_with_names,!duplicated(genus_species))
Clean_clean_ejfs<-GEx_unique_species[c(7,9)]

GEx_codominant_species$clean_ejf <- Clean_clean_ejfs$clean_ejf[match(GEx_codominant_species$genus_species, Clean_clean_ejfs$genus_species)]
                  
GEx_codominant_species$biogeographic.realm <- GEx_metadata$biogeographic.realm[match(GEx_codominant_species$site,GEx_metadata$site)]
GEx_codominant_species$location <- GEx_metadata$Location[match(GEx_codominant_species$site,GEx_metadata$site)]

##import trait data and subset to desired traits
#BIEN_Traits<-read.csv("BIEN_traitFINAL.csv")
#Desirable_traits<-read.csv("Trait list.csv")
#Filtered_BIEN<- BIEN_Traits[BIEN_Traits$trait_name %in% Desirable_traits$Trait, ]
#Unique_traits_by_species<-Filtered_BIEN%>%distinct(clean_ejf,trait_name,.keep_all = TRUE)
All_trait_data<-read.csv("Merged trait database 5 traits Means4.csv")
head(All_trait_data)

#Get all available trait means per species
Leaf_Area<-subset(All_trait_data,All_trait_data$TraitName == "Leaf Area")
Height<-subset(All_trait_data,All_trait_data$TraitName == "Height")
Leaf_N_Mass<-subset(All_trait_data,All_trait_data$TraitName == "Leaf N Mass")
Lifespan<-subset(All_trait_data,All_trait_data$TraitName == "Plant Lifespan")
SLA<-subset(All_trait_data,All_trait_data$TraitName == "SLA")

##link desired traits to dominant species names
#Leaf Area
GEx_codominant_species$Leaf_Area_Mean <- Leaf_Area$Mean.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Leaf_Area$SpeciesName)]
GEx_codominant_species$Leaf_Area_Std.Dev <- Leaf_Area$Std.Dev.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Leaf_Area$SpeciesName)]
GEx_codominant_species$Leaf_Area_Min <- Leaf_Area$Std.Dev.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Leaf_Area$SpeciesName)]
GEx_codominant_species$Leaf_Area_Max <- Leaf_Area$Max.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Leaf_Area$SpeciesName)]
GEx_codominant_species$Leaf_Area_n <- Leaf_Area$N.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Leaf_Area$SpeciesName)]
#Leaf_N_Mass
GEx_codominant_species$Leaf_N_Mass_Mean <- Leaf_N_Mass$Mean.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Leaf_N_Mass$SpeciesName)]
GEx_codominant_species$Leaf_N_Mass_Std.Dev <- Leaf_N_Mass$Std.Dev.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Leaf_N_Mass$SpeciesName)]
GEx_codominant_species$Leaf_N_Mass_Min <- Leaf_N_Mass$Std.Dev.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Leaf_N_Mass$SpeciesName)]
GEx_codominant_species$Leaf_N_Mass_Max <- Leaf_N_Mass$Max.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Leaf_N_Mass$SpeciesName)]
GEx_codominant_species$Leaf_N_Mass_n <- Leaf_N_Mass$N.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Leaf_N_Mass$SpeciesName)]
#Lifespan
GEx_codominant_species$Lifespan_Mean <- Lifespan$Mean.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Lifespan$SpeciesName)]
GEx_codominant_species$Lifespan_Std.Dev <- Lifespan$Std.Dev.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Lifespan$SpeciesName)]
GEx_codominant_species$Lifespan_Min <- Lifespan$Std.Dev.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Lifespan$SpeciesName)]
GEx_codominant_species$Lifespan_Max <- Lifespan$Max.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Lifespan$SpeciesName)]
GEx_codominant_species$Lifespan_n <- Lifespan$N.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Lifespan$SpeciesName)]
#Height
GEx_codominant_species$Height_Mean <- Height$Mean.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Height$SpeciesName)]
GEx_codominant_species$Height_Std.Dev <- Height$Std.Dev.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Height$SpeciesName)]
GEx_codominant_species$Height_Min <- Height$Std.Dev.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Height$SpeciesName)]
GEx_codominant_species$Height_Max <- Height$Max.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Height$SpeciesName)]
GEx_codominant_species$Height_n <- Height$N.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,Height$SpeciesName)]
#SLA
GEx_codominant_species$SLA_Mean <- SLA$Mean.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,SLA$SpeciesName)]
GEx_codominant_species$SLA_Std.Dev <- SLA$Std.Dev.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,SLA$SpeciesName)]
GEx_codominant_species$SLA_Min <- SLA$Std.Dev.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,SLA$SpeciesName)]
GEx_codominant_species$SLA_Max <- SLA$Max.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,SLA$SpeciesName)]
GEx_codominant_species$SLA_n <- SLA$N.Trait.Value.Cont.[match(GEx_codominant_species$clean_ejf,SLA$SpeciesName)]

##Create dominant only file
GEx_dominant_species<-subset(GEx_codominant_species,GEx_codominant_species$rank==1)

##Plot proportion of species covered for each trait by region
library(ggplot2)

gg_miss_var(GEx_codominant_species[c(16,17,22,27,32,37)], facet = location, show_pct = TRUE)


gg_miss_var(GEx_dominant_species[c(16,17,22,27,32,37)], facet = location, show_pct = TRUE)

#########################################
##For all species recorded in all plots##
#########################################

##import plot level codominant species data and add region to datafile
GEx_All_data<-read.csv("GEx_cleaned_v3.csv")
GEx_metadata<-read.csv("Meta_SEV2019_v2.csv")

GEx_All_data$biogeographic.realm <- GEx_metadata$biogeographic.realm[match(GEx_All_data$site,GEx_metadata$site)]
GEx_All_data$location <- GEx_metadata$Location[match(GEx_All_data$site,GEx_metadata$site)]

##import trait data and subset to desired traits
#Get all available trait means per species
Leaf_Area<-subset(All_trait_data,All_trait_data$TraitName == "Leaf Area")
Height<-subset(All_trait_data,All_trait_data$TraitName == "Height")
Leaf_N_Mass<-subset(All_trait_data,All_trait_data$TraitName == "Leaf N Mass")
Lifespan<-subset(All_trait_data,All_trait_data$TraitName == "Plant Lifespan")
SLA<-subset(All_trait_data,All_trait_data$TraitName == "SLA")

##link desired traits to dominant species names
#Leaf Area
GEx_All_data$Leaf_Area_Mean <- Leaf_Area$Mean.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Leaf_Area$SpeciesName)]
GEx_All_data$Leaf_Area_Std.Dev <- Leaf_Area$Std.Dev.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Leaf_Area$SpeciesName)]
GEx_All_data$Leaf_Area_Min <- Leaf_Area$Std.Dev.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Leaf_Area$SpeciesName)]
GEx_All_data$Leaf_Area_Max <- Leaf_Area$Max.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Leaf_Area$SpeciesName)]
GEx_All_data$Leaf_Area_n <- Leaf_Area$N.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Leaf_Area$SpeciesName)]
#Leaf_N_Mass
GEx_All_data$Leaf_N_Mass_Mean <- Leaf_N_Mass$Mean.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Leaf_N_Mass$SpeciesName)]
GEx_All_data$Leaf_N_Mass_Std.Dev <- Leaf_N_Mass$Std.Dev.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Leaf_N_Mass$SpeciesName)]
GEx_All_data$Leaf_N_Mass_Min <- Leaf_N_Mass$Std.Dev.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Leaf_N_Mass$SpeciesName)]
GEx_All_data$Leaf_N_Mass_Max <- Leaf_N_Mass$Max.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Leaf_N_Mass$SpeciesName)]
GEx_All_data$Leaf_N_Mass_n <- Leaf_N_Mass$N.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Leaf_N_Mass$SpeciesName)]
#Lifespan
GEx_All_data$Lifespan_Mean <- Lifespan$Mean.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Lifespan$SpeciesName)]
GEx_All_data$Lifespan_Std.Dev <- Lifespan$Std.Dev.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Lifespan$SpeciesName)]
GEx_All_data$Lifespan_Min <- Lifespan$Std.Dev.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Lifespan$SpeciesName)]
GEx_All_data$Lifespan_Max <- Lifespan$Max.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Lifespan$SpeciesName)]
GEx_All_data$Lifespan_n <- Lifespan$N.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Lifespan$SpeciesName)]
#Height
GEx_All_data$Height_Mean <- Height$Mean.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Height$SpeciesName)]
GEx_All_data$Height_Std.Dev <- Height$Std.Dev.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Height$SpeciesName)]
GEx_All_data$Height_Min <- Height$Std.Dev.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Height$SpeciesName)]
GEx_All_data$Height_Max <- Height$Max.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Height$SpeciesName)]
GEx_All_data$Height_n <- Height$N.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,Height$SpeciesName)]
#SLA
GEx_All_data$SLA_Mean <- SLA$Mean.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,SLA$SpeciesName)]
GEx_All_data$SLA_Std.Dev <- SLA$Std.Dev.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,SLA$SpeciesName)]
GEx_All_data$SLA_Min <- SLA$Std.Dev.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,SLA$SpeciesName)]
GEx_All_data$SLA_Max <- SLA$Max.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,SLA$SpeciesName)]
GEx_All_data$SLA_n <- SLA$N.Trait.Value.Cont.[match(GEx_All_data$clean_ejf,SLA$SpeciesName)]

##Plot proportion of species covered for each trait by region
library(ggplot2)

gg_miss_var(GEx_All_data[c(11,12,18,23,28,33)], facet = location, show_pct = TRUE)


#write.csv(GEx_codominant_species,"GEx_Codominant_Species_Trait.csv")
#write.csv(GEx_dominant_species,"GEx_Dominant_Species_Trait.csv")
#write.csv(GEx_All_data,"GEx_All_Species_Trait.csv")


#Add in Covariates by site
Covariate_data<-read.csv("GEx-metadata-with-other-env-layers-v2.csv")
head(Covariate_data)
Cov_dat_only<-Covariate_data[c(2,24:48)]

GEx_All_Sites_species_trait_environmental<-merge(GEx_All_data, Cov_dat_only, by.x = "site", by.y = "site")
GEx_Codominant_species_trait_environmental<-merge(GEx_codominant_species, Cov_dat_only, by.x = "site", by.y = "site")
GEx_Dominant_species_trait_environmental<-merge(GEx_dominant_species, Cov_dat_only, by.x = "site", by.y = "site")

write.csv(GEx_All_Sites_species_trait_environmental,"GEx_All_Sites_species_trait_environmental.csv")
write.csv(GEx_Codominant_species_trait_environmental,"GEx_Codominant_species_trait_environmental.csv")
write.csv(GEx_Dominant_species_trait_environmental,"GEx_Dominant_species_trait_environmental.csv")









