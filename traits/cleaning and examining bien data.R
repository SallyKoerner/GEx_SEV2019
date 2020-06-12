### Clean BIEN data and examine data coverage by traits and species
### 
### Authors: kevin wilcox (kevin.wilcox@uwyo.edu); Robert Griffin-Nolan (robertgn13@gmail.com)

file.choose()
setwd("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\trait data\\")

trait_data <- read.csv("GEx_species_BEIN_trait.csv")
sp_names_gex <- read.csv("GEx_species_family.csv")

unique(trait_data$trait_name)

head(filter(trait_data, trait_name=="maximum plant height"))

### Get trait list and write to file
trait_list_bien <- unique(trait_data$trait_name)
# write.csv(trait_list_bien, "BIEN_traitList.csv", row.names=F)

trait_matrix <- trait_data %>%
  dplyr::select(scrubbed_species_binomial,trait_name, trait_value) %>%
  group_by(scrubbed_species_binomial, trait_name) %>%
  summarise(num_record = length(trait_value)) %>%
  spread(key=trait_name, value=num_record)

trait_presence_counter <- trait_matrix %>%
  gather(key=trait_name, value=num_record, -scrubbed_species_binom) %>%
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


