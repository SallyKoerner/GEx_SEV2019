library(tidyverse)
library(codyn)
library(vegan)

setwd("~/Dropbox/GEx_SEV/Data&Code")

dat<-read.csv("ALL_Cleaned_April2019_V2.csv")

lsyear<-dat%>%
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage)

R<-lsyear%>%
  group_by(site, year, exage, block, plot, trt)%>%
  summarise(rich=n())%>%
  spread(trt, rich)%>%
  rename(URich=U, GRich=G)
Rich<-R%>%
  mutate(Rich_LogRR=(log(GRich/URich)))%>%
  filter(Rich_LogRR!="NA")

E_q<-function(x){
  x1<-x[x!=0]
  if (length(x1)==1) {
    return(NA)
  }
  if (abs(max(x1) - min(x1)) < .Machine$double.eps^0.5) {##bad idea to test for zero, so this is basically doing the same thing testing for a very small number
    return(1)
  }
  r<-rank(x1, ties.method = "average")
  r_scale<-r/max(r)
  x_log<-log(x1)
  fit<-lm(r_scale~x_log)
  b<-fit$coefficients[[2]]
  2/pi*atan(b)
}


#function to calculate E1/D (inverse of Simpson's) from Smith and Wilson 1996
#' @S the number of species in the sample
#' @x the vector of abundances of each species
#' @N the total abundance
#' @p the vector of relative abundances of each species
E_simp<-function(x, S=length(x[x!=0]), N=sum(x[x!=0]), ps=x[x!=0]/N, p2=ps*ps ){
  D<-sum(p2)
  (1/D)/S
}

SimpD<-function(x, S=length(x[x!=0]), N=sum(x[x!=0]), ps=x[x!=0]/N, p2=ps*ps ){
  D<-sum(p2)
  D
}
#calculating gini coefficeint using the gini function in the reldist package
#' @x the vector of abundances of each species
#' this tive the inverse of other measures of evenness??
Gini<-function(x){
  x1<-x[x!=0]
  1-reldist::gini(x1)
}


newdom<-lsyear%>%
  group_by(site, year, exage, block, plot, trt)%>%
  summarise(S=n(), 
            div=diversity(relcov), 
            even=div/log(S),
            E_Q=E_q(relcov),
            Gini=Gini(relcov),
            SimpD=SimpD(relcov))%>%
  select (-div, -S, -E_Q, -Gini, -even)%>%
  spread(trt, SimpD)%>%
  rename(USimpD=U, GSimpD=G)
newdom2<-newdom%>%
  mutate(SimpD_LogRR=(log(GSimpD/USimpD)))%>%
  filter(SimpD_LogRR!="NA")

BPDom<-lsyear%>%
  group_by(site, year, exage, block, plot, trt)%>%
  summarise(dom=max(relcov))%>%
  spread(trt, dom)%>%
  rename(UDom=U, GDom=G)
BPDom2<-BPDom%>%
  mutate(BPDom_LogRR=(log(GDom/UDom)))%>%
  filter(BPDom_LogRR!="NA")

unique(Rich$site)
unique(newdom2$site)
RE<-merge(Rich, newdom2, by=c("site", "year", "exage", "block", "plot"))
RD<-merge(RE, BPDom2, by=c("site", "year", "exage", "block", "plot"))

RD2<-RD%>%
  group_by(site, year, exage, block)%>%
  summarise(URich=mean(URich), GRich=mean(GRich), UDom=mean(UDom), GDom=mean(GDom), 
            USimpD=mean(USimpD), GSimpD=mean(GSimpD), Rich_LogRR=mean(Rich_LogRR), 
            BPDom_LogRR=mean(BPDom_LogRR), SimpD_LogRR=mean(SimpD_LogRR))
unique(RD2$site)

##graphing richness differences and CI
gex_rich_CI<-RD2%>%
  group_by(site) %>% 
  summarize_at(vars(Rich_LogRR), funs(mean, sd, length)) %>% 
  mutate(se=sd/sqrt(length),
         ci=se*1.96)

mean<-mean(gex_rich_CI$mean)
sd<-sd(gex_rich_CI$mean)
se<-sd/sqrt(252)
ci<-se*1.96
site<-"All Sites"

all<- data.frame(site, mean, sd, se, ci)
gex_rich_CI2<-gex_rich_CI[order(gex_rich_CI$mean),]

gex_rich_CI2<-gex_rich_CI2%>%
  bind_rows(all)%>%
  mutate(colortrt=ifelse(site=="All Sites", 1,0))

gex_rich_CI2$site2<-factor(gex_rich_CI2$site, levels=gex_rich_CI2$site)

ggplot(data=gex_rich_CI2, aes(x=site2, y=mean, color=as.factor(colortrt)))+
  geom_point(stat="identity")+
  scale_color_manual(values=c("black", "red"))+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci))+
  coord_flip()+
  geom_hline(yintercept = 0)+
  xlab("Site")+
  ylab("Change in Richness (LRR)")+
  theme(axis.text.y=element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none")


##graphing BP Dom differences and CI
gex_dom_CI<-RD2%>%
  group_by(site) %>% 
  summarize_at(vars(BPDom_LogRR), funs(mean, sd, length)) %>% 
  mutate(se=sd/sqrt(length),
         ci=se*1.96)

mean<-mean(gex_dom_CI$mean)
sd<-sd(gex_dom_CI$mean)
se<-sd/sqrt(252)
ci<-se*1.96
site<-"All Sites"

all<- data.frame(site, mean, sd, se, ci)
gex_dom_CI2<-gex_dom_CI[order(gex_dom_CI$mean),]

gex_dom_CI2<-gex_dom_CI2%>%
  bind_rows(all)%>%
  mutate(colortrt=ifelse(site=="All Sites", 1,0))

gex_dom_CI2$site2<-factor(gex_dom_CI2$site, levels=gex_dom_CI2$site)

ggplot(data=gex_dom_CI2, aes(x=site2, y=mean, color=as.factor(colortrt)))+
  geom_point(stat="identity")+
  scale_color_manual(values=c("black", "red"))+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci))+
  coord_flip()+
  geom_hline(yintercept = 0)+
  xlab("Site")+
  ylab("Change in BP Dominance (LRR)")+
  theme(axis.text.y=element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none")



 
##### Numerical Dominant, identify, determine change
setwd("~/Dropbox/GEx_SEV/Data&Code")
dat<-read.csv("ALL_Cleaned_April2019_V2.csv")
test<-dat%>%
  select(site, block)%>%
  unique()

dat2<-dat%>%
  group_by(site, block, trt, genus_species)%>%
  summarise(relcov=mean(relcov))%>%
  ungroup()%>%
  group_by(site, block, trt)%>%
  filter(relcov==max(relcov))%>%
  ungroup()
#unique(dat2$site)

diff<-read.csv("gex_abund_diff_all.csv")
#unique(diff$site)

diff2<-left_join(dat2, diff, by=c("site", "block", "genus_species"))%>%
  mutate(site_block=paste(site, block, sep='::'))%>%
  mutate(exp_unit=paste(site, block, trt, sep='::'))

DomChangeNum<-diff2%>%
  filter(trt=="G")%>%
  select(site, block, genus_species)

#Dominant Identiy Change with % cover
DomChangeNum2<-left_join(DomChangeNum, dat2)%>%
  spread(trt, relcov)
  
  DomChangeNum2$U[is.na(DomChangeNum2$U)]<-0
  
DomChangeNum3<-DomChangeNum2 %>%
  mutate(U2=ifelse(U==0, U+1, U))%>%
  mutate(DomIdentityChange_LogRR=(log(G/U2)))
write.csv(DomChangeNum3, file="DomIdentityNumChangeLRR_byBlock_SEV_April2019.csv", row.names=FALSE)

DomChangeNum4<-DomChangeNum3 %>%
  group_by(site)%>%
  summarise(DomIdentityChange_LogRR=mean(DomIdentityChange_LogRR))
write.csv(DomChangeNum4, file="DomIdentityNumChangeLRR_bySite_SEV_April2019.csv", row.names=FALSE)


#generate rank of each species in each plot by relative cover, with rank 1 being most abundant
datrank<-dat%>%
  group_by(site, block, trt, genus_species)%>%
  summarise(relcov=mean(relcov))%>%
  ungroup()%>%
  mutate(exp_unit=paste(site, block, trt, sep='::'))
#unique(datrank$site)

rankOrder <- datrank%>%
  group_by(exp_unit)%>%
  mutate(ranks = rank(-relcov, ties.method = "min"))%>%
  ungroup()%>%
  select(site, block,trt, genus_species, ranks, relcov)
#unique(rankOrder$site)


# how does rank of dom in grazed change ---- could not Figure this out - Lauren and Tim made it work
gdom<-diff2%>%
  filter(trt=="G")%>%
  mutate(rankgrazed=1)%>%
  select(-trt, -relcov, -difference, -exp_unit, -site_block)
gdom2<- full_join(gdom, rankOrder, by=c("site", "block", "genus_species"))%>%
  filter(trt=="U")%>%
  filter(rankgrazed!="NA")

SpNum<-dat%>%
  filter(trt=="G")%>%
  group_by(site, year, exage, block)%>%
  summarise(numSp=n())%>%
  ungroup()

domrankshift<-merge(gdom2, SpNum, by=c("site", "block"), all=TRUE)%>%
  filter(rankgrazed!="NA")


#Create same/diff in dom species - cannot use simple spread because some sites have multiple co-dominants
overlapSpeciesDoubled <- data.frame(row.names=1)

diff3 <- diff2%>%
  mutate(site_block=paste(site, block, sep='::'))%>%
  select(-trt)%>%
  unique()

for(i in 1:length(diff3$site_block)) {
  subset <- diff2[diff2$site_block==diff3$site_block[i],]
  
  subsetG <- subset%>%
    filter(trt=='G')
  
  subsetU <- subset%>%
    filter(trt=='U')
  
  vector <- intersect(subsetG$genus_species,subsetU$genus_species)
  
  match <- data.frame(site=subsetG$site,
                      block=subsetG$block,
                      overlap=ifelse(is_empty(vector), 'diff', 'same'))
  
  overlapSpeciesDoubled=rbind(match, overlapSpeciesDoubled)

}

overlapSpecies <- unique(overlapSpeciesDoubled)
  
overlapNumber <- overlapSpecies%>%
  mutate(overlapNum=ifelse(overlap=="diff", 1, 0))%>%
    group_by(site)%>%
    summarise(num=mean(overlapNum))%>%
  ungroup()

write.csv(overlapSpecies, file="Dom_SameDiff_byBlock_SEV_April2019.csv", row.names=FALSE)
write.csv(overlapNumber, file="Dom_SameDiff_bySite_SEV_April2019.csv", row.names=FALSE)

