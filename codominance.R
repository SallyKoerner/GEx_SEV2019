################################################################################
##  codominance.R: Function to calculate codominance of a community.
##
##  Author: Kimberly Komatsu
##  Date created: April 24, 2019
################################################################################

library(psych)
library(tidyverse)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\data')

###read in data
GEx <- read.csv('All_Cleaned_April2019_V2.csv')%>%
  rename(cover=relcov)%>%
  mutate(exp_unit=paste(site, block, plot, trt, year, sep='::'))

# #subset a test dataset of one plot to create a function
# test <- GEx%>%
#   filter(site=='Konza'&block=='A'&plot==1)


#############################################
#####calculate Cmax (codominance metric)#####

#calculate relative abundance
relCover <- GEx%>%
  group_by(exp_unit)%>%
  summarise(totcov=sum(cover))%>%
  ungroup()%>%
  right_join(GEx)%>%
  mutate(relcov=(cover/totcov)*100)%>%
  select(-cover, -totcov)

#generate rank of each species in each plot by relative cover, with rank 1 being most abundant
rankOrder <- relCover%>%
  group_by(exp_unit)%>%
  mutate(rank = as.numeric(order(order(relcov, decreasing=TRUE))))%>%
  ungroup()

###calculating harmonic means for all subsets of rank orders
#make a new dataframe with just the label
expUnit=GEx%>%
  select(exp_unit)%>%
  unique()

#makes an empty dataframe
harmonicMean=data.frame(row.names=1) 

#calculate harmonic means
for(i in 1:length(expUnit$exp_unit)) {
  
  #creates a dataset for each unique experimental unit
  subset <- rankOrder[rankOrder$exp_unit==as.character(expUnit$exp_unit[i]),]%>%
    select(exp_unit, genus_species, relcov, rank)
  
  for(j in 1:length(subset$rank)) {
    
    #creates a dataset for each series of ranks from 1 through end of the number of ranks
    subset2 <- subset[subset$rank<=j,]
    
    #calculate harmonic mean of values
    mean <- harmonic.mean(subset2$relcov)
    meanData <- data.frame(exp_unit=unique(subset2$exp_unit),
                           num_ranks=j, 
                           harmonic_mean=mean)
    
    harmonicMean=rbind(meanData, harmonicMean)
    
  }
  
}

differenceData <- harmonicMean%>%
  left_join(rankOrder)%>%
  filter(rank==num_ranks+1)%>% #only keep the next most abundant species after the number that went into the calculation of the harmonic mean
  mutate(difference=harmonic_mean-relcov) #calculates difference between harmonic mean and the relative cover of the next most abundant species
  
Cmax <- differenceData%>%
  group_by(exp_unit)%>%
  summarise(Cmax=max(difference))%>%
  ungroup()%>%
  left_join(differenceData)%>%
  filter(Cmax==difference)%>%
  rename(num_codominants=num_ranks)%>%
  select(exp_unit, Cmax, num_codominants)

codomSppList <- Cmax%>%
  left_join(rankOrder)%>%
  group_by(exp_unit)%>%
  filter(rank<=num_codominants)%>%
  ungroup()













