setwd("C://Users/Jason Donaldson/Documents/Workshops and conferences/GEx/")

##import plot level species data
speccomp<-read.csv("All_Cleaned_April2019_V2.csv")
##create list of unique species
uniquespec<-subset(speccomp,!duplicated(genus_species))
##import accepted species name list
Acc_spec_name<-read.csv("GEx_species_family.csv")


##get a dataframe with clean unique species names
uniquespec$scrubbed_species_binomial <- Acc_spec_name$clean_ejf[match(uniquespec$genus_species, Acc_spec_name$genus_species)]

##import BIEN species trait data
All_trait_dat<-read.csv("BIEN_traitFINAL.csv")
head(All_trait_dat)


Plot_Spec_Trait<-All_trait_dat[All_trait_dat$scrubbed_species_binomial %in% uniquespec$scrubbed_species_binomial, ,drop = F]

dat1[!dat1$id %in% dat2$id, , drop = FALSE]




