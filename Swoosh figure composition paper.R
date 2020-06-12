####created by M. Avolio June 2020 for virtual working group
###this code will create figure one for composition difference paper

library(tidyverse)
library(codyn)
library(rsq)
library(vegan)
library(gtable)
library(gridExtra)

theme_set(theme_bw(12))

setwd("C:\\Users\\mavolio2\\Dropbox\\GEx_VirtualWorkshop_June2020\\")

###making comp diff shoowsh figure
comp<-read.csv("gex_multdiff_site_block.csv")

#getting CI for each site
gex_multdiff_CI<-comp%>%
  group_by(site) %>% 
  summarize_at(vars(composition_diff), funs(mean, sd, length)) %>% 
  mutate(se=sd/sqrt(length),
         ci=se*1.96)

#geting global mean comp diff
mean<-mean(gex_multdiff_CI$mean)
sd<-sd(gex_multdiff_CI$mean)
se<-sd/sqrt(252)
ci<-se*1.96
site<-"All Sites"

all<- data.frame(site, mean, sd, se, ci)

gex_multdiff_CI <- gex_multdiff_CI[order(gex_multdiff_CI$mean), ]  # sort
gex_multdiff_CI2<-gex_multdiff_CI%>%
  bind_rows(all)%>%
  mutate(colortrt=ifelse(site=="All Sites", 1, ifelse(site=="AUS_Berry", 2, ifelse(site=="CPER", 3, ifelse(site=="Kruger_Satara", 4, 0)))))%>%
  mutate(size=ifelse(site=="All Sites"|site=="AUS_Berry"|site=="CPER"|site=="Kruger_Satara", 0.5, 0))

gex_multdiff_CI2$site2 <- factor(gex_multdiff_CI2$site, levels = gex_multdiff_CI2$site)

highlight<-gex_multdiff_CI2%>%
  filter(site=="AUS_Berry"|site=="CPER"|site=="Kruger_Satara")

shoosh<-
  ggplot(data=gex_multdiff_CI2, aes(x=site2, y=mean, color=as.factor(colortrt)))+
  geom_point(aes(size=as.factor(size)))+
  scale_size_manual(values=c(1, 3.5))+
  scale_color_manual(values=c("black", "green4", "dodgerblue", "red", "orange"))+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci))+
  geom_point(data=highlight, aes(x=site2, y= mean, color=as.factor(colortrt)), size=3.5)+
  coord_flip()+
  xlab("Site")+
  ylab("Compositional Difference")+
  #scale_y_continuous(limits = c(0,1))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

####makign the NMDS insets
dat<-read.csv("GEx_cleaned_11June2020.csv")

theme_set(theme_bw(10))

dat2<-dat%>%
  mutate(drop=ifelse(genus_species_use=="#N/A"|genus_species_use=="Dead unidentified"|genus_species_use=="Leaf.Litter"|genus_species_use=="cactus__dead_", 1, 0))%>%
  filter(drop!=1)

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

#dropping probelmatic blocks and averaging over experiments where there are several blocks in a plot. 
#Averaging up to a block for Konza, Junner Koeland, and Jornada sites.

lsyear2<-lsyear%>%
  filter(site_block!="California_Sedgwick_Lisque##A"&
           site_block!="California_Sedgwick_Lisque##B"&
           site_block!="California_Sedgwick_Lisque##G"&
           site_block!="DesertLow##Alkali")%>%
  group_by(site, block, site_block, trt, exage, genus_species_use)%>%
  summarize(relcov=mean(relcov))

site_list<-lsyear2%>%
  ungroup()%>%
  select(site)%>%
  unique()

## making NMDS figures
###High comp Diff example

hsub<-lsyear2%>%
  filter(site=="Junner Koeland")

hwide<-hsub%>%
  spread(genus_species_use, relcov, fill=0)

hspecies<-hwide[,6:ncol(hwide)]

henv<-hwide[,1:5]

hmds<-metaMDS(hspecies, distance = "bray")

hscores <- data.frame(scores(hmds, display="sites"))  # Extracts NMDS scores
hscores2<- cbind(henv, hscores) # binds the NMDS scores of year i to all years previously run

high<-ggplot(hscores2, aes(x=NMDS1, y=NMDS2, color=trt, fill=trt, shape=block))+
  geom_point(size=3)+
  scale_color_manual(name="Treatment", values=c("deepskyblue", "dodgerblue4"), labels=c("Grazed", 'Ungrazed'))+
  scale_fill_manual(name="Treatment", values=c("deepskyblue", "dodgerblue4"), labels=c("Grazed", 'Ungrazed'))+
  scale_shape_manual(name="Block Pair", values=c(21:25))+
  annotate("text", x=0.8, y = -0.5, label="0.082", size=3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(shape=F)


###Mid comp Diff example

msub<-lsyear2%>%
  filter(site=="California_RanchoMarino")

mwide<-msub%>%
  spread(genus_species_use, relcov, fill=0)

mspecies<-mwide[,6:ncol(mwide)]

menv<-mwide[,1:5]

mmds<-metaMDS(mspecies, distance = "bray")

mscores <- data.frame(scores(mmds, display="sites"))  # Extracts NMDS scores
mscores2<- cbind(menv, mscores) # binds the NMDS scores of year i to all years previously run

mid<-
ggplot(mscores2, aes(x=NMDS1, y=NMDS2, color=trt, fill=trt, shape=block))+
  geom_point(size=3)+
  scale_color_manual(name="Treatment", values=c("firebrick1", "firebrick4"), labels=c("Grazed", 'Ungrazed'))+
  scale_fill_manual(name="Treatment", values=c("firebrick1", "firebrick4"), labels=c("Grazed", 'Ungrazed'))+
  scale_shape_manual(name="Block Pair", values=c(21:25))+
  annotate("text", x=-0.4, y = -0.5, label="0.073", size=3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(shape = FALSE)

###Low comp Diff example

lsub<-lsyear2%>%
  filter(site=="Kruger_Satara")

lwide<-lsub%>%
  spread(genus_species_use, relcov, fill=0)

lspecies<-lwide[,6:ncol(wide)]

lenv<-lwide[,1:5]

lmds<-metaMDS(lspecies, distance = "bray")

lscores <- data.frame(scores(lmds, display="sites"))  # Extracts NMDS scores
lscores2<- cbind(lenv, lscores) # binds the NMDS scores of year i to all years previously run

low<-
  ggplot(lscores2, aes(x=NMDS1, y=NMDS2, color=trt, fill=trt, shape=block))+
  geom_point(size=3)+
  scale_color_manual(name="Treatment", values=c("orange", "darkorange3"), labels=c("Grazed", 'Ungrazed'))+
  scale_fill_manual(name="Treatment", values=c("orange", "darkorange3"), labels=c("Grazed", 'Ungrazed'))+
  scale_shape_manual(name="Block Pair", values=c(8,20,21:25))+
  annotate("text", x=0.6, y = -0.45, label="0.211", size=3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(shape = FALSE)

examples<-grid.arrange(high, mid, low, ncol=1)

figure1<-grid.arrange(shoosh, examples, ncol=2, widths=c(3,2))

ggsave("Figure1.jpg", figure1)
