###code by M. Avolio June 2020 for virtual working group
##this code uses the berger-parker dominace file to create the swhoosh figure1 for dominance paper


library(tidyverse)
library(codyn)
library(rsq)
library(gridExtra)
theme_set(theme_bw(10))


###Meghan's Working Directory
#setwd("C:\\Users\\mavolio2\\Dropbox\\GEx_VirtualWorkshop_June2020\\")
setwd("/Users/skoerne/Dropbox/GEx/GEx_VirtualWorkshop_June2020")

diffG_U<-read.csv("Diff_BP_Dom_allyrs.11June2020_v2.csv")

###using the last year only for the figure
diffG_U.lyr<-diffG_U%>%
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage)

#how many see a change in dom sp?
sum(diffG_U.lyr$diff_sp)/561

diffG_U.lyr_mean<-diffG_U.lyr%>%
  group_by(site)%>%
  summarize(ave=mean(bp_rr))

diffG_U.lyr_mean$site2 <- factor(diffG_U.lyr_mean$site, levels = diffG_U.lyr_mean$site[order(diffG_U.lyr_mean$ave)])

toplot<-diffG_U.lyr%>%
  left_join(diffG_U.lyr_mean)%>%
  select(site, block, bp_rr, diff_sp, ave, site2)

#make plot of mean with mean a bigger dot
swoosh<-
ggplot(data=toplot, aes(x=site2, y=bp_rr, color=as.factor(diff_sp)))+
  geom_point()+
  scale_color_manual(name="Dominant\nSpecies Identity", values=c("blue", "red"), labels=c("No change", "Change"))+
  coord_flip()+
  geom_hline(yintercept = 0)+
  xlab("Site")+
  ylab("Difference in Dominance")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank())+
  annotate("text", x=200, y=-1.2, label="Herbivores decrease\ndominance")+
  annotate("text", x=60, y=1, label="Herbivores increase\ndominance")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##make a six bar chart for abs and raw, for no change, change, and global. This is the figure mendy wants
overall<-toplot%>%
  group_by(diff_sp)%>%
  summarize(abs=mean(abs(bp_rr)),
            global=mean(bp_rr),
            sd_abs=sd(abs(bp_rr)),
            sd_glob=sd(bp_rr),
            n=length(bp_rr))%>%
  mutate(se_abs=sd_abs/sqrt(n),
         se_glob=sd_glob/sqrt(n),
         ci_abs=se_abs*1.96,
         ci_glob=se_glob*1.96)%>%
  select(diff_sp, abs, global, ci_abs, ci_glob) %>% 
  ungroup()

overall2<-toplot%>%
  ungroup() %>% 
  summarize(abs=mean(abs(bp_rr)),
            global=mean(bp_rr),
            sd_abs=sd(abs(bp_rr)),
            sd_glob=sd(bp_rr),
            n=length(bp_rr))%>%
  mutate(se_abs=sd_abs/sqrt(n),
         se_glob=sd_glob/sqrt(n),
         ci_abs=se_abs*1.96,
         ci_glob=se_glob*1.96)%>%
  select(abs, global, ci_abs, ci_glob)%>%
  mutate(diff_sp=2)

toplot_all<-overall%>%
  bind_rows(overall2)%>%
  gather(type, mean, abs:global)%>%
  gather(ci_type, ci, ci_abs:ci_glob)%>%
  mutate(keep=ifelse(type=="abs"&ci_type=='ci_abs', 1, ifelse(type=="global"&ci_type=="ci_glob", 1, 0)))%>%
  filter(keep==1)

#replace with points? change black color.
means<-
ggplot(data=toplot_all, aes(x=type, y=mean, color=as.factor(diff_sp)))+
  geom_point(position = position_dodge(0.9), size=5)+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), position=position_dodge(0.9), width=0.5)+
  geom_hline(yintercept = 0)+
  xlab("")+
  ylab("Difference in Dominance")+
  scale_x_discrete(labels=c("Absolute Value", "Raw Value"))+
  scale_color_manual(name="Dominant\nSpecies Identity", values=c("blue", "red", "black"), labels=c("No Change", "Change", "All Data"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dominancefig1<-grid.arrange(swoosh, means, ncol=1)

ggsave("domfig1.jpeg",dominancefig1, width = 5.8, height=5.3, units="in")
