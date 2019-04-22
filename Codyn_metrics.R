library(tidyverse)
library(codyn)
library(rsq)

setwd("~/Dropbox/Dominance WG/")

dat<-read.csv("ALL_Cleaned_Apr2019.csv")

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

gex_multdiff<-data.frame()
siteblock<-unique(lsyear2$site_block)

for (i in 1:length(siteblock)){
  
  subset<-lsyear2%>%
    filter(site_block==siteblock[i]) %>%
    mutate(treatment=trt)
  
  out <- multivariate_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "trt", treatment.var = "treatment")
  
  out$site_block<-siteblock[i]
  
  gex_multdiff<-rbind(gex_multdiff, out)  
}

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

for (i in 1:length(siteblockplot)){
  
  subset<-lsyear3%>%
    filter(site_block_plot==siteblockplot[i]) %>%
    mutate(treatment=trt)
  
  out <- multivariate_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "trt", treatment.var = "treatment")
  
  un_site_block<-unique(subset$site_block)
  un_plot<-unique(subset$plot)
  
  out$site_block<-un_site_block
  out$plot<-un_plot
  
  gex_multdiff2<-rbind(gex_multdiff2, out)  
}

ave_gexmultdiff2<-gex_multdiff2%>%
  group_by(site_block, treatment, treatment2, trt_greater_disp, abs_dispersion_diff)%>%
  summarize(composition_diff=mean(composition_diff))

gex_multdiff_all<-gex_multdiff %>% 
  bind_rows(ave_gexmultdiff2)%>%
  separate(site_block, into=c("site", "block"), sep="##")%>%
  select(site, block, composition_diff)

write.csv(gex_multdiff_all, "gex_multdiff_all.csv", row.names = F)

gex_multdiff_ave<-gex_multdiff_all%>%
  group_by(site) %>% 
  summarize(composition_diff=mean(composition_diff))

write.csv(gex_multdiff_ave, "gex_multdiff_ave.csv", row.names = F)

##graphing compositional differences and CI
gex_multdiff_CI<-gex_multdiff_all%>%
  separate(site_block, into=c("site", "block"), sep="##") %>% 
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

###doing RAC differences

gex_RACdiff<-data.frame()
siteblock<-unique(lsyear2$site_block)

for (i in 1:length(siteblock)){
  
  subset<-lsyear2%>%
    filter(site_block==siteblock[i]) %>%
    mutate(treatment=trt)
  
  out <- RAC_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "trt", treatment.var = "treatment")
  
  out$site_block<-siteblock[i]
  
  gex_RACdiff<-rbind(gex_RACdiff, out)  
}

gex_RACdiff2<-data.frame()
siteblockplot<-unique(lsyear3$site_block_plot)

for (i in 1:length(siteblockplot)){
  
  subset<-lsyear3%>%
    filter(site_block_plot==siteblockplot[i]) %>%
    mutate(treatment=trt)
  
  out <- RAC_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "trt", treatment.var = "treatment")
  
  un_site_block<-unique(subset$site_block)
  un_plot<-unique(subset$plot)
  
  out$site_block<-un_site_block
  out$plot<-un_plot
  
  gex_RACdiff2<-rbind(gex_RACdiff2, out)  
}

gexRACdiff2_ave<-gex_RACdiff2%>%
  group_by(site_block, trt, trt2, treatment, treatment2)%>%
  summarize_at(vars(richness_diff, evenness_diff, rank_diff, species_diff), funs(mean))

gex_RACdiff_all<-gex_RACdiff %>% 
  bind_rows(gexRACdiff2_ave)%>%
  separate(site_block, into=c("site", "block"), sep="##")%>%
  select(site, block, richness_diff, evenness_diff, rank_diff, species_diff)

write.csv(gex_RACdiff_all, "gex_RACdiff_all.csv", row.names = F)

gex_RACdiff_ave<-gex_RACdiff_all%>%
  group_by(site) %>% 
  summarize_at(vars(richness_diff, evenness_diff, rank_diff, species_diff), funs(mean))

write.csv(gex_RACdiff_ave, "gex_RACdiff_ave.csv", row.names = F)

###doing abundance differences

gex_abunddiff<-data.frame()
siteblock<-unique(lsyear2$site_block)

for (i in 1:length(siteblock)){
  
  subset<-lsyear2%>%
    filter(site_block==siteblock[i]) %>%
    mutate(treatment=trt)
  
  out <- abundance_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "trt", treatment.var = "treatment")
  
  out$site_block<-siteblock[i]
  
  gex_abunddiff<-rbind(gex_abunddiff, out)  
}

gex_abunddiff2<-data.frame()
siteblockplot<-unique(lsyear3$site_block_plot)

for (i in 1:length(siteblockplot)){
  
  subset<-lsyear3%>%
    filter(site_block_plot==siteblockplot[i]) %>%
    mutate(treatment=trt)
  
  out <- abundance_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "trt", treatment.var = "treatment")
  
  un_site_block<-unique(subset$site_block)
  un_plot<-unique(subset$plot)
  
  out$site_block<-un_site_block
  out$plot<-un_plot
  
  gex_abunddiff2<-rbind(gex_abunddiff2, out)  
}

gexabundiff2_ave<-gex_abunddiff2%>%
  group_by(site_block, trt, trt2, treatment, treatment2, genus_species)%>%
  summarize_at(vars(difference), funs(mean))

gex_abunddiff_all<-gex_abunddiff %>% 
  bind_rows(gexabundiff2_ave)%>%
  separate(site_block, into=c("site", "block"), sep="##")%>%
  select(site, block, genus_species, difference)

write.csv(gex_abunddiff_all, "gex_abund_diff_all.csv", row.names = F)

gex_abunddiff_ave<-gex_abunddiff_all%>%
  group_by(site, genus_species) %>% 
  summarize_at(vars(difference), funs(mean))

write.csv(gex_abunddiff_ave, "gex_abund_diff_ave.csv", row.names = F)

###doing a multiple regression, what is driving comp diff?
comprac<-gex_multdiff_ave%>%
  left_join(gex_RACdiff_ave)%>%
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


###getting clean dataset
lsyear_export<-dat%>%
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage)%>%
  mutate(site_block=paste(site, block, sep="##")) %>%
  mutate(site_block_plot=paste(site_block, plot, sep="##"))

lsyear_export2<-lsyear_export%>%
  filter(site_block!="California_Sedgwick_Lisque##A"&
           site_block!="California_Sedgwick_Lisque##B"&
           site_block!="California_Sedgwick_Lisque##G"&
           site_block!="DesertLow##Alkali"&
           site_block_plot!="Jornada##east##30"&
           site_block_plot!="Jornada##west##18"&
           site_block_plot!="Jornada##west##19"&
           site_block_plot!="Jornada##west##20"&
           site_block_plot!="Jornada##west##21"&
           site_block_plot!="Jornada##west##22")%>%
  separate(site_block_plot, into=c("site","block","plot"), sep="##")%>%
  select(-X, -lyear, -mexage, -site_block)%>%
  select(site, block, plot, exage, year, trt, genus_species, relcov)
write.csv(lsyear_export2, "All_Cleaned_April2019_V2.csv", row.names=F)  
