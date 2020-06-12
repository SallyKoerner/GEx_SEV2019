library(tidyverse)
library(codyn)
library(rsq)

setwd("C:\\Users\\mavolio2\\Dropbox\\GEx_VirtualWorkshop_June2020\\")
comp_all<-read.csv("gex_multdiff_site_block.csv")
comp<-read.csv("gex_multdiff_site_ave.csv")
rac<-read.csv("gex_racdiff_site_ave.csv")
experiment<-read.csv("Meta_SEV2019.csv")%>%
  select(site, ANPP, PlotSize, domestic, dom, precip, Num_ex, grazing.pressure, Grazer_Domestic, herbivore.type, guild, sprich, X..herbivore.spp)
exage<-read.csv("All_Cleaned_April2019_V2.csv")%>%
  select(site, exage)%>%
  unique()

##graphing compositional differences and CI
gex_multdiff_CI<-comp_all%>%
  group_by(site) %>% 
  summarize_at(vars(composition_diff), funs(mean, sd, length)) %>% 
  mutate(se=sd/sqrt(length),
         ci=se*1.96)

mean<-mean(gex_multdiff_CI$mean)
sd<-sd(gex_multdiff_CI$mean)
se<-sd/sqrt(252)
ci<-se*1.96
site<-"All Sites"

all<- data.frame(site, mean, sd, se, ci)


gex_multdiff_CI <- gex_multdiff_CI[order(gex_multdiff_CI$mean), ]  # sort
gex_multdiff_CI2<-gex_multdiff_CI%>%
  bind_rows(all)%>%
  mutate(colortrt=ifelse(site=="All Sites", 1,0))
gex_multdiff_CI2$site2 <- factor(gex_multdiff_CI2$site, levels = gex_multdiff_CI2$site)

theme_set(theme_bw(12))
ggplot(data=gex_multdiff_CI2, aes(x=site2, y=mean, color=as.factor(colortrt)))+
  geom_point(stat="identity")+
  scale_color_manual(values=c("black", "red"))+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci))+
  coord_flip()+
  geom_hline(yintercept = 0.25)+
  xlab("Site")+
  ylab("Compositional Difference")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")


###doing a multiple regression, what is driving comp diff for racs?
comprac<-comp%>%
  left_join(rac)%>%
  mutate(abs_rich=abs(richness_diff),
         abs_even=abs(evenness_diff))%>%
  select(-richness_diff, -evenness_diff)

summary(m1<-lm(composition_diff~abs_rich+abs_even+rank_diff+species_diff, data=comprac))
rsq.partial(m1)

comprac2<-comprac %>% 
  na.omit()

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  test <- cor.test(x,y) 
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 1),
                   symbols = c("*", " "))
  
  
  text(0.5, 0.5, txt, cex = 2)
  text(0.8, 0.5, Signif, cex=5, col="red")
}
pairs(comprac2[,c(2, 5, 6, 3, 4)], pch = 21, labels=c("Compositional\nDiff","Abs. Richness\nDiff", "Abs. Evenness\nDiff","Rank\nDiff","Species\nDiff"), font.labels=1, cex.labels=2,upper.panel=panel.cor)

###doing multiple regression with experiment details
comp_herb<-comp%>%
  right_join(exage)

comp_herb2<-merge(comp_herb, experiment, by="site")%>%
  mutate(num_exc=as.numeric(as.character(Num_ex)),
        PlotSize=as.numeric(as.character(PlotSize)))

summary(m1<-lm(composition_diff~exage+num_exc+PlotSize, data=comp_herb2))
rsq.partial(m1)
comp_herb3<-comp_herb2%>%
  na.omit

pairs(comp_herb3[,c(2, 3, 5, 9)], pch = 21, labels=c("Comp\nDiff", "Exclosure\nAge","Plot Size", "Num\nExclosure"), font.labels=1, cex.labels=2,upper.panel=panel.cor)

###doing multiple regression with biolgoical details

comp_bio<-merge(comp, experiment, by="site")%>%
  mutate(precip=as.numeric(as.character(precip)),
         anpp=as.numeric(as.character(ANPP)))

summary(m1<-lm(composition_diff~precip+grazing.pressure+guild+X..herbivore.spp+dom+sprich, data=comp_bio))
rsq.partial(m1)

summary(m1<-lm(composition_diff~anpp+grazing.pressure+guild+X..herbivore.spp+dom+sprich, data=comp_bio))
rsq.partial(m1)

comp_bio2<-comp_bio%>%
  na.omit

pairs(comp_bio[,c(2, 15, 6,13, 7,9, 12,14)], pch = 21, labels=c("Comp\nDiff", "ANPP", "Dominance", "Richness", "Precip", "Graze\nPressure", "Guild", "Herb\nRich"), font.labels=1, cex.labels=2,upper.panel=panel.cor)

hist(comp_bio$grazing.pressure)
