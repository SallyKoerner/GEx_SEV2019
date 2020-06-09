library(tidyverse)
library(codyn)
library(rsq)
library(vegan)

setwd("/Users/avahoffman/Dropbox/Research/Grazing_consortium/2020/GEx_SEV2019")

dat<-read.csv("All_Cleaned_April2019_V2.csv")

lsyear<-dat%>%
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage)%>%
  mutate(block = as.numeric(block))%>%
  mutate(site_block=paste(site, block, sep="##"))

lsyear2<-
  lsyear%>%
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
  filter(num!=1)

gex_permanova<-data.frame()
sitelist<-unique(lsyear2.2$site)


for (i in sitelist){
  
  subset<-lsyear2.2%>%
    filter(site==i) %>%
    mutate(treatment=trt)
  
  wide<-
    subset%>%
    spread(genus_species, relcov, fill=0)
  
  species<-
    as.data.frame(wide[,12:ncol(wide)])
  
  env<-
    as.data.frame(wide[,1:11])
  
  # Perform permamova
  # Try with and without `strata =` option
  out <- 
    adonis(species~trt, data=env, method="bray")
  #adonis(wide[,12:ncol(wide)]~trt, wide, strata = wide$block)
  
  # Perform alternative W* test
  # First do distance matrix
  dm <- 
    as.matrix(dist(species))
  
  Tw_test_out <- 
    Tw2.test(dm, env$trt)
  WdS_test_out <- 
    WdS.test(dm, env$trt)
  
  # Assemble data frame
  perm_out <- data.frame(
    site = i,
    perm_Pvalue =  out$aov.tab$'Pr(>F)'[1],
    Tw_Pvalue = Tw_test_out$p.value,
    WdS_Pvalue = WdS_test_out$p.value
  )
  perm_out
  
  gex_permanova<-rbind(gex_permanova, perm_out)  
   
}

test2<-gex_permanova%>%
  left_join(numreps)

## look at NMDS to see gut check
test<-metaMDS(species, distance = "bray")

scores <- data.frame(scores(test, display="sites"))  # Extracts NMDS scores for year "i" #
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
