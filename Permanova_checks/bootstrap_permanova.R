library(tidyverse)
library(codyn)
library(rsq)
library(vegan)
library(ggplot2)

setwd("/Users/avahoffman/Dropbox/Research/Grazing_consortium/2020/GEx_SEV2019")

dat <- read.csv("All_Cleaned_April2019_V2.csv")

lsyear <- dat %>%
  group_by(site) %>%
  mutate(lyear = max(year),
         mexage = max(exage)) %>%
  filter(year == lyear, mexage == exage) %>%
  mutate(site_block = paste(site, block, sep = "##"))

lsyear2 <-
  lsyear %>%
  filter(
    site_block != "California_Sedgwick_Lisque##A" &
      site_block != "California_Sedgwick_Lisque##B" &
      site_block != "California_Sedgwick_Lisque##G" &
      site_block != "DesertLow##Alkali" &
      site != "Jornada" &
      site != "Junner Koeland" &
      site != "Konza"
  )

numreps <- lsyear %>%
  select(site, block) %>%
  unique() %>%
  group_by(site) %>%
  summarize(num = length(block))

lsyear2.2 <- lsyear2 %>%
  left_join(numreps) %>%
  filter(num > 1)

run_grazing_bootstrap <-
  function(sample_threshold = 6) {
    
    gex_permanova <- data.frame()
    sitelist <- unique(lsyear2.2$site)
    
    for (i in sitelist) {
      # Subset dataset by site i
      subset <-
        lsyear2.2 %>%
        filter(site == i)
      
      # Convert long spp dataset to wide
      wide <-
        subset %>%
        spread(genus_species, relcov, fill = 0)
      
      # Species relative abundance matrix
      species <-
        as.data.frame(wide[, 12:ncol(wide)])
      
      # Treatment vector
      env <-
        as.data.frame(wide[, "trt"])
      
      # What is the site level replication?
      replication <- wide$num[1]
      
      # Calculate how many bootstrap samples needed
      if (replication < sample_threshold) {
        replication_needed <-
          sample_threshold - replication
      } else if (replication >= sample_threshold) {
        replication_needed <- 0
      }
      
      # Mean relative abundances of Grazing plots
      G_prob <-
        colMeans(as.data.frame((wide %>%
                                  filter(trt == "G"))[, 12:ncol(wide)]))
      
      # Mean relative abundances of Ungrazed plots
      U_prob <-
        colMeans(as.data.frame((wide %>%
                                  filter(trt == "U"))[, 12:ncol(wide)]))
      
      # Note: assumption of 100 individuals per plot
      # Bootstrap Grazed plots
      G_boot <-
        t(rmultinom(n = replication_needed,
                    size = 100,
                    prob = G_prob))
      # Bootstrap Ungrazed plots
      U_boot <-
        t(rmultinom(n = replication_needed,
                    size = 100,
                    prob = U_prob))
      
      species <- rbind(species, G_boot, U_boot)
      env <- rbind(env,
                   data.frame(trt = rep("G", replication_needed)),
                   data.frame(trt = rep("U", replication_needed)))
      
      # Perform permamova
      # Try with and without `strata =` option
      out <-
        adonis(species ~ trt, data = env, method = "bray")
      
      # Assemble data frame
      perm_out <- data.frame(
        site = i,
        perm_Pvalue =  out$aov.tab$'Pr(>F)'[1],
        real_n = replication,
        boot_n = replication_needed
      )
      
      # Append results from each site
      gex_permanova <- rbind(gex_permanova, perm_out)
    }
    
    return(gex_permanova)
  }

gex_permanova <- run_grazing_bootstrap() 


nmds_plot <-
  function() {
    ## Not run
    ## look at NMDS to see gut check
    test <- metaMDS(species, distance = "bray")
    
    scores <-
      data.frame(scores(test, display = "sites"))  # Extracts NMDS scores for year "i" #
    scores2 <-
      cbind(env, scores) # binds the NMDS scores of year i to all years previously run
    
    ggplot(scores2, aes(x = NMDS1, y = NMDS2, color = trt)) +
      geom_point(size = 5)
  }
