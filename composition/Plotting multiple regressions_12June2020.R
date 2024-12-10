### Plotting multiple regressions
###
### Created: June 12, 2020
### Last updated: June 12, 2020
### New updat: December 10, 2024

###
### Set up workspace
###
rm(list=ls())
#setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Working groups\\GEx\\community_paper\\")
setwd("/Users/skoerne/Dropbox/GEx/GEx_VirtualWorkshop_June2020")
library(tidyverse)
library(ggthemes)
library(ggpubr)

###
### Read in data, and clean! ( remove NMDS1 that is NULL and make grazing pressure a factor )
###

data_full <- read.csv("community_difference_allmetrics_siteavg_12June2020c.csv") %>%
  mutate(grazing.pressure = as.factor(grazing.pressure)) %>%
  drop_na(NMDS1) %>%
  mutate(grazing_pressure_name = ifelse(grazing.pressure==1, "Light", 
                                        ifelse(grazing.pressure==2, "Moderate",
                                               "Heavy"))) %>%
  mutate(grazing_pressure_name = factor(grazing_pressure_name, levels=c("Light","Moderate","Heavy")))

### Means and standard error by grazing pressure
grazing_pressure_summary <- data_full %>%
  group_by(grazing.pressure) %>%
  summarize(composition_diff_mean = mean(composition_diff),
            composition_diff_se = sd(composition_diff)/sqrt(length(composition_diff)))%>%
  mutate(grazing_pressure_name = ifelse(grazing.pressure==1, "Light", 
                                        ifelse(grazing.pressure==2, "Moderate",
                                               "Heavy"))) %>%
  mutate(grazing_pressure_name = factor(grazing_pressure_name, levels=c("Light","Moderate","Heavy")))


###
### Plot individual regressions and barplots
###

grazing_pressure_plot <- ggplot(grazing_pressure_summary, aes(x=grazing_pressure_name, y=composition_diff_mean,
                                     ymin=composition_diff_mean-composition_diff_se,
                                     ymax=composition_diff_mean+composition_diff_se)) +
  geom_col(fill="white", col="black") +
  geom_errorbar(width=0.1) +
  theme_few() +
  ylab("Composition Diff") +
  xlab("Grazing Pressure")

pca1_climate_plot <- ggplot(data_full, aes(x=Prin1.tropicsvscontinental, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Diff") +
  xlab("PCA1: Contin(-) vs Trop(+)")
    
pca2_climate_plot <- ggplot(data_full, aes(x=Prin2.hotvswet, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Diff") +
  xlab("PCA2: Wet (-) vs Hot (+)")

pca3_climate_plot <- ggplot(data_full, aes(x=Prin3.precipseasonality, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Diff") +
  xlab("PCA3: PPT seasonality")

nmds1_plot <- ggplot(data_full, aes(x=NMDS1, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Diff") +
  xlab("NMDS1 (Beth?)")

nmds2_plot <- ggplot(data_full, aes(x=NMDS2, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Diff") +
  xlab("NMDS2 (Beth?)")

dominance_plot <- ggplot(data_full, aes(x=dom, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Diff") +
  xlab("Dominance (Berger-Parker?)")

richness_plot <- ggplot(data_full, aes(x=sprich, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Diff") +
  xlab("Site-level Species Richness")

photomix_plot <- ggplot(data_full, aes(x=PhotoMix, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Diff") +
  xlab("C3-C4 Mix")

totalc4_plot <- ggplot(data_full, aes(x=ALLC4, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Diff") +
  xlab("Proportion C4?")

pca1_by_grazing_plot <- ggplot(data_full, aes(x=Prin1.tropicsvscontinental, y=composition_diff, col=grazing_pressure_name)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  theme_few() +
  ylab("Composition Diff") +
  xlab("PCA1: Contin(-) vs Trop(+)") +
  scale_color_manual(values=c("dodgerblue","darkorange2","firebrick4")) +
  theme(legend.position="none")

pca2_by_grazing_plot <- ggplot(data_full, aes(x=Prin2.hotvswet, y=composition_diff, col=grazing_pressure_name)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  theme_few() +
  ylab("Composition Diff") +
  xlab("PCA2: Wet (-) vs Hot (+)") +
  scale_color_manual(values=c("dodgerblue","darkorange2","firebrick4")) +
  theme(legend.position="none")

pca3_by_grazing_plot <- ggplot(data_full, aes(x=Prin3.precipseasonality, y=composition_diff, col=grazing_pressure_name)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  theme_few() +
  ylab("Composition Diff") +
  xlab("PCA3: PPT seasonality") +
  scale_color_manual(values=c("dodgerblue","darkorange2","firebrick4")) +
  theme(legend.position="none")

dominance_by_grazing_plot <- ggplot(data_full, aes(x=dom, y=composition_diff, col=grazing_pressure_name)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  theme_few() +
  ylab("Composition Diff") +
  xlab("Dominance (Berger-Parker(?))") +
  scale_color_manual(values=c("dodgerblue","darkorange2","firebrick4")) +
  theme(legend.position="none")

richness_by_grazing_plot <- ggplot(data_full, aes(x=sprich, y=composition_diff, col=grazing_pressure_name)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  theme_few() +
  ylab("Composition Diff") +
  xlab("Site-level Species Richness") +
  scale_color_manual(values=c("dodgerblue","darkorange2","firebrick4")) +
  theme(legend.position="none")

nmds1_by_grazing_plot <- ggplot(data_full, aes(x=NMDS1, y=composition_diff, col=grazing_pressure_name)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  theme_few() +
  ylab("Composition Diff") +
  xlab("NMDS1") +
  scale_color_manual(values=c("dodgerblue","darkorange2","firebrick4")) +
  theme(legend.position="none")

nmds2_by_grazing_plot <- ggplot(data_full, aes(x=NMDS2, y=composition_diff, col=grazing_pressure_name)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  theme_few() +
  ylab("Composition Diff") +
  xlab("NMDS2") +
  scale_color_manual(values=c("dodgerblue","darkorange2","firebrick4")) +
  theme(legend.position="none")

photomix_by_grazing_plot <- ggplot(data_full, aes(x=PhotoMix, y=composition_diff, col=grazing_pressure_name)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  theme_few() +
  ylab("Composition Diff") +
  xlab("C3-C4 Mix") +
  scale_color_manual(values=c("dodgerblue","darkorange2","firebrick4")) +
  theme(legend.position="none")

totalc4_by_grazing_plot <- ggplot(data_full, aes(x=ALLC4, y=composition_diff, col=grazing_pressure_name)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  theme_few() +
  ylab("Composition Diff") +
  xlab("Proportion C4") +
  scale_color_manual(values=c("dodgerblue","darkorange2","firebrick4")) +
  theme(legend.position="none")

###
### Stitch plots together into single panel
###

best_model_1_panels <- ggarrange(grazing_pressure_plot, 
                                 pca1_climate_plot,
                                 pca2_climate_plot,
                                 pca3_climate_plot,
                                 dominance_plot,
                                 richness_plot,
                                 nmds1_plot, 
                                 nmds2_plot,
                                 photomix_plot,
                    ncol = 3, nrow = 3)

best_model_2_panels <- ggarrange(
                                  pca1_by_grazing_plot,
                                  pca2_by_grazing_plot,
                                  pca3_by_grazing_plot,
                                  dominance_by_grazing_plot,
                                  richness_by_grazing_plot,
                                  nmds1_by_grazing_plot,
                                  nmds2_by_grazing_plot,
                                  photomix_by_grazing_plot,
                                  ncol = 3, nrow = 3)

pdf(file=paste0("Site factor bivariate plots_",Sys.Date(),"_1of2.pdf"), width=9, height=9, useDingbats = F)
print(best_model_1_panels)
dev.off()

png(file=paste0("Site factor bivariate plots_",Sys.Date(),"_1of2.png"), width=675, height=675, units = "px")
print(best_model_1_panels)
dev.off()

pdf(file=paste0("Site factor bivariate plots_",Sys.Date(),"_2of2.pdf"), width=9, height=9, useDingbats = F)
print(best_model_2_panels)
dev.off()

png(file=paste0("Site factor bivariate plots_",Sys.Date(),"_2of2.png"), width=675, height=675, units = "px")
print(best_model_2_panels)
dev.off()


#####Statistical Models
head(data_full)
#with PhyloRealms
summary(m1<-lm(composition_diff~grazing.pressure+Prin1.tropicsvscontinental+Prin2.hotvswet+Prin3.precipseasonality+NMDS1+NMDS2+dom + sprich + PhotoMix, data=data_full))
rsq.partial(m1)

#without PhyloRealms
summary(m1<-lm(composition_diff~grazing.pressure+Prin1.tropicsvscontinental+Prin2.hotvswet+Prin3.precipseasonality+dom + sprich + PhotoMix, data=data_full))
rsq.partial(m1)
 
library(rsq)
library(ggthemes)

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
pairs(data_full[,c(42, 81, 82, 83, 53, 54, 89)], pch = 21, labels=c("Grazing \nPressure","PCA1", "PCA2","PCA3","dom", "rich", "C3-C4 \nMix"), font.labels=1, cex.labels=2, upper.panel=panel.cor)
#export 1000x1000
###New Figure with new models (droping NMDS 1)
###
### Plot individual regressions and barplots
###

grazing_pressure_plot <- ggplot(grazing_pressure_summary, aes(x=grazing_pressure_name, y=composition_diff_mean,
                                                              ymin=composition_diff_mean-composition_diff_se,
                                                              ymax=composition_diff_mean+composition_diff_se)) +
  geom_col(fill="white", col="black") +
  geom_errorbar(width=0.1) +
  theme_few() +
  ylab("Composition Difference") +
  xlab("Grazing Pressure")

pca1_climate_plot <- ggplot(data_full, aes(x=Prin1.tropicsvscontinental, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Difference") +
  xlab("PCA1: Contin(-) vs Trop(+)")

pca2_climate_plot <- ggplot(data_full, aes(x=Prin2.hotvswet, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm", linetype="dotted") +
  theme_few() +
  ylab("Composition Difference") +
  xlab("PCA2: Wet (-) vs Hot (+)")

pca3_climate_plot <- ggplot(data_full, aes(x=Prin3.precipseasonality, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Difference") +
  xlab("PCA3: PPT seasonality")

dominance_plot <- ggplot(data_full, aes(x=dom, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm", linetype="dotted") +
  theme_few() +
  ylab("Composition Difference") +
  xlab("Site-level Dominance")

richness_plot <- ggplot(data_full, aes(x=sprich, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm", linetype="dotted") +
  theme_few() +
  ylab("Composition Difference") +
  xlab("Site-level Species Richness")

photomix_plot <- ggplot(data_full, aes(x=PhotoMix, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm", linetype="dotted") +
  theme_few() +
  ylab("Composition Difference") +
  xlab("C3-C4 Mix")

###
### Stitch plots together into single panel
###

full_model_1_panels <- ggarrange(grazing_pressure_plot, 
                                 pca1_climate_plot,
                                 pca2_climate_plot,
                                 pca3_climate_plot,
                                 dominance_plot,
                                 richness_plot,
                                 photomix_plot,
                                 ncol = 4, nrow = 2)


png(file=paste0("Site factor bivariate plots_",Sys.Date(),"_1of2.png"), width=750, height=475, units = "px")
print(full_model_1_panels)
dev.off()


##### DO same analysis and figure with mechanisms
head(data_full)
summary(m1<-lm(composition_diff~richness_diff+evenness_diff+rank_diff+species_diff+bp_rr+diff_sp, data=data_full))
rsq.partial(m1)

data_full_no0s<-data_full %>% 
  filter (bp_rr!="NA")
pairs(data_full_no0s[,c(4, 5, 6, 7, 8, 9)], pch = 21, labels=c("Richness\nDifference", "Evenness\nDifference","Rank\nDifference","Species\nDifference", "Dominance\nResponse Ratio", "Dominant \nIdentity Switch"), font.labels=1, cex.labels=2,upper.panel=panel.cor)
#export at 1100x1100

rich_diff <- ggplot(data_full, aes(x=richness_diff, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm", linetype="dotted") +
  theme_few() +
  ylab("Composition Difference") +
  xlab("Richness Difference")

even_diff <- ggplot(data_full, aes(x=evenness_diff, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm", linetype="dotted") +
  theme_few() +
  ylab("Composition Difference") +
  xlab("Evenness Difference")

rank_diff <- ggplot(data_full, aes(x=rank_diff, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Difference") +
  xlab("Rank Difference")

species_diff <- ggplot(data_full, aes(x=species_diff, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Difference") +
  xlab("Species Difference")

bp_diff <- ggplot(data_full, aes(x=bp_rr, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Difference") +
  xlab("Dominance Response Ratio")

diff_sp <- ggplot(data_full, aes(x=diff_sp, y=composition_diff)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_few() +
  ylab("Composition Difference") +
  xlab("Different Dominant Species")

### Stitch plots together into single panel
###

full_model_1_panels <- ggarrange(rich_diff, 
                                 even_diff,
                                 species_diff,
                                 rank_diff,
                                 bp_diff,
                                 diff_sp,
                                 ncol = 3, nrow = 2)


png(file=paste0("Mechanism bivariate plots_",Sys.Date(),"_1of2.png"), width=750, height=475, units = "px")
print(full_model_1_panels)
dev.off()

