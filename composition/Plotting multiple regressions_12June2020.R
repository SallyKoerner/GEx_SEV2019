### Plotting multiple regressions
###
### Created: June 12, 2020
### Last updated: June 12, 2020

###
### Set up workspace
###
rm(list=ls())
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Working groups\\GEx\\community_paper\\")
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





