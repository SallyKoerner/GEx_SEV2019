library(tidyverse)
library(codyn)
library(rsq)

theme_set(theme_bw(12))

#do this with last year for the figure
#do for all years to export email to beth and ben, figure out if I am using the right folder
###Meghan's Working Directory
setwd("C:\\Users\\mavolio2\\Dropbox\\Dominance WG\\")
##Sally's Working Directory
setwd("~/Dropbox/GEx_VirtualWorkshop_June2020")


dat<-read.csv("GEx_cleaned_10June2020.csv")%>%
  mutate(site_block=paste(site, block, sep="##"))

###all years of data
dat2<-dat%>%
  filter(site_block!="California_Sedgwick_Lisque##A"&
           site_block!="California_Sedgwick_Lisque##B"&
           site_block!="California_Sedgwick_Lisque##G"&
           site_block!="DesertLow##Alkali"&
           site!="Jornada"&
           site!="Junner Koeland"&
           site!="Konza")%>%
  select(-plot)


dat3<-dat%>%
  filter(site=="Jornada"|
           site=="Junner Koeland"|
           site=="Konza")%>%
  mutate(site_block_plot=paste(site_block, plot, sep="##"))%>%
  filter(site_block_plot!="Jornada##east##30"&
           site_block_plot!="Jornada##west##18"&
           site_block_plot!="Jornada##west##19"&
           site_block_plot!="Jornada##west##20"&
           site_block_plot!="Jornada##west##21"&
           site_block_plot!="Jornada##west##22")%>%
  group_by(site, year, exage, block, trt, genus_species, site_block)%>%
  summarise(relcov=mean(relcov))

dat4<-dat2%>%
  bind_rows(dat3)%>%
  group_by(site, block, year, trt)%>%
  mutate(dom=max(relcov))%>%
  filter(dom==relcov)

#we are dropping sites that have multiple dominant species in a block
keep<-dat4%>%
  group_by(site, block, year, trt)%>%
  summarize(n=length(relcov))%>%
  filter(n<2)

dat5<-dat4%>%
  right_join(keep)%>%
  select(-site_block, -dom, -n)

graze<-dat5%>%
  ungroup()%>%
  filter(trt=="G")%>%
  rename(graze_sp=genus_species,
         graze_bp=relcov,
         graze_clean_ejf=clean_ejf)%>%
  select(-trt)

ungraze<-dat5%>%
  ungroup()%>%
  filter(trt=="U")%>%
  rename(ungraze_sp=genus_species,
         ungraze_bp=relcov, ungraze_clean_ejf=clean_ejf)%>%
  select(-trt)

diffG_U<-ungraze%>%
  left_join(graze)%>%
  na.omit%>%
  mutate(bp_rr=log(graze_bp/ungraze_bp),
         diff_sp=ifelse(graze_clean_ejf==ungraze_clean_ejf, 0, 1))

##this is the file that beth and ben want
write.csv(diffG_U, "Diff_BP_Dom_allyrs.csv")

###using the last year only for the figure
diffG_U.lyr<-diffG_U%>%
  group_by(site) %>% 
  mutate(lyear=max(year),
         mexage=max(exage))%>%
  filter(year==lyear, mexage==exage)

#how many see a change in dom sp?
sum(diffG_U.lyr$diff_sp)/573

diffG_U.lyr_mean<-diffG_U.lyr%>%
  group_by(site)%>%
  summarize(ave=mean(bp_rr))

diffG_U.lyr_mean$site2 <- factor(diffG_U.lyr_mean$site, levels = diffG_U.lyr_mean$site[order(diffG_U.lyr_mean$ave)])

toplot<-diffG_U.lyr%>%
  left_join(diffG_U.lyr_mean)%>%
  select(site, block, bp_rr, diff_sp, ave, site2)

#make plot of mean with mean a bigger dot

ggplot(data=toplot, aes(x=site2, y=bp_rr, color=as.factor(diff_sp)))+
  geom_point()+
  scale_color_manual(name="Dominant\nSpecies Identity", values=c("blue", "red"), labels=c("No change", "Change"))+
  coord_flip()+
  geom_hline(yintercept = 0)+
  xlab("Site")+
  ylab("Difference in Dominance")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank())+
  annotate("text", x=200, y=-1.25, label="Herbivores decrease\ndominance")+
  annotate("text", x=60, y=1.1, label="Herbivores increase\ndominance")

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

ggplot(data=toplot_all, aes(x=type, y=mean, fill=as.factor(diff_sp)))+
  geom_bar(stat = "identity", position = position_dodge(0.9))+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), position=position_dodge(0.9), width=0.5)+
  geom_hline(yintercept = 0)+
  xlab("")+
  ylab("Change in Dominance")+
  scale_x_discrete(labels=c("Absolute Value", "Raw Value"))+
  scale_fill_manual(name="Dominant\nSpecies Identity", values=c("blue", "red", "black"), labels=c("No Change", "Change", "All Data"))

#doing a different way, I don't like this as much.
ggplot(data=toplot_all, aes(x=diff_sp, y=mean, fill=type))+
  geom_bar(stat = "identity", position = position_dodge(0.9))+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), position=position_dodge(0.9), width=0.5)+
  geom_hline(yintercept = 0)+
  xlab("")+
  ylab("Change in Dominance")+
  #scale_x_discrete(labels=c("No Change", "Change", "All Data"))+ can't get this to work for some reason
  scale_fill_manual(name="Dominant\nSpecies Identity", values=c("blue", "red"), labels=c("Absolute Value", "Raw Value"))
