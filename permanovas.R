library(tidyverse)
library(codyn)
library(rsq)
library(vegan)

setwd("C:\\Users\\mavolio2\\Dropbox\\Dominance WG\\")

dat<-read.csv("All_Cleaned_April2019_V2.csv")

exage<-dat%>%
  select(site, exage)%>%
  unique()%>%
  filter(exage>10)

lsyear<-dat%>%
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage)%>%
  mutate(site_block=paste(site, block, sep="##"))

lsyear2<-lsyear%>%
  filter(site_block!="California_Sedgwick_Lisque##A"&
           site_block!="California_Sedgwick_Lisque##B"&
           site_block!="California_Sedgwick_Lisque##G"&
           site_block!="DesertLow##Alkali"&
           site!="Jornada"&
           site!="Junner Koeland"&
           site!="Konza")

numreps<-lsyear%>%
  select(site, block)%>%
  unique()%>%
  group_by(site)%>%
  summarize(num=length(block))

lsyear2.2<-lsyear2%>%
  left_join(numreps)%>%
  filter(num>4)

gex_permanova<-data.frame()
sitelist<-unique(lsyear2.2$site)


for (i in 1:length(sitelist)){
  
  subset<-lsyear2.2%>%
    filter(site==sitelist[i]) 
  
  wide<-subset%>%
    spread(genus_species, relcov, fill=0)
  
  species<-wide[,12:ncol(wide)]
  
  env<-wide[,1:11]
  out <- adonis(species~trt, data=env, method="bray", strata=env$block)
  #adonis(wide[,12:ncol(wide)]~trt, wide, strata = wide$block) this does the same thing
  
  perm_out <- data.frame(
    site = sitelist[i],
    perm_Pvalue =  out$aov.tab$'Pr(>F)'[1]
  )
  
   gex_permanova<-rbind(gex_permanova, perm_out)  
}

test2<-gex_permanova%>%
  left_join(numreps)

## look at NMDS to see gut check
msub<-lsyear2.2%>%
  filter(site=="Argentina_ElPalmar") 

mwide<-msub%>%
  spread(genus_species, relcov, fill=0)

species<-mwide[,12:ncol(mwide)]

env<-mwide[,1:11]


mds<-metaMDS(species, distance = "bray")

scores <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores 
scores2<- cbind(env, scores) # binds the NMDS scores of year i to all years previously run

ggplot(scores2, aes(x=NMDS1, y=NMDS2, color=trt, shape=block))+
  geom_point(size=5)

####ignore below
permanova_out_mod <- permanova_out_master %>%
  mutate(pval_flag_perm = ifelse(perm_Pvalue<.05, 1, 0)) %>%
  mutate(pval_flag_disp = ifelse(disp_Pvalue<.05, 1, 0)) %>%
  mutate(permanova="permanova") %>%
  group_by(site_project_comm, treatment) %>%
  summarise(tot_pval=sum(pval_flag)) %>%
  filter(tot_pval != 0)

##three sites, konza, jornada and junner koeland, had plots nested in block. need to use that informaiton here.

lsyear3<-lsyear%>%
  filter(site=="Jornada"|
           site=="Junner Koeland"|
           site=="Konza")%>%
  mutate(site_block_plot=paste(site_block, plot, sep="##"))%>%
  filter(site_block_plot!="Jornada##east##30"&
           site_block_plot!="Jornada##west##18"&
           site_block_plot!="Jornada##west##19"&
           site_block_plot!="Jornada##west##20"&
           site_block_plot!="Jornada##west##21"&
           site_block_plot!="Jornada##west##22")

gex_multdiff2<-data.frame()
siteblockplot<-unique(lsyear3$site_block_plot)
