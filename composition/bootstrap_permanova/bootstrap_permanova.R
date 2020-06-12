# This script contains code for performing the bootstrap permanova analysis and output
###########################################################################################
###########################################################################################
# Load libraries
library(tidyverse)
library(codyn)
library(rsq)
library(vegan)
library(ggplot2)

## -- REPOSITORY DIR -- ##
setwd("/Users/avahoffman/Dropbox/Research/Grazing_consortium/2020/GEx_SEV2019")

source("composition/bootstrap_permanova/util.R")
source("composition/bootstrap_permanova/swooshplot.R")


process_comp_data <-
  function() {
    ## -- DATA CURRENTLY GITIGNORED -- SEE GOOGLE DRIVE FOR LATEST -- ##
    dat <- read.csv("GEx_cleaned_11June2020.csv")
    
    dat2 <-
      dat %>%
      mutate(
        drop = ifelse(
          genus_species_use == "#N/A" |
            genus_species_use == "Dead unidentified" |
            genus_species_use == "Leaf.Litter" |
            genus_species_use == "cactus__dead_",
          1,
          0
        )
      ) %>%
      filter(drop != 1)
    
    exage <-
      dat %>%
      select(site, exage) %>%
      unique() %>%
      filter(exage > 10)
    
    lsyear <-
      dat %>%
      group_by(site) %>%
      mutate(lyear = max(year),
             mexage = max(exage)) %>%
      filter(year == lyear, mexage == exage) %>%
      mutate(site_block = paste(site, block, sep = "##"))
    
    #dropping probelmatic blocks and averaging over experiments where there are several blocks in a plot.
    #Averaging up to a block for Konza, Junner Koeland, and Jornada sites.
    
    lsyear2 <-
      lsyear %>%
      filter(
        site_block != "California_Sedgwick_Lisque##A" &
          site_block != "California_Sedgwick_Lisque##B" &
          site_block != "California_Sedgwick_Lisque##G" &
          site_block != "DesertLow##Alkali"
      ) %>%
      group_by(site, block, site_block, trt, exage, genus_species_use) %>%
      summarize(relcov = mean(relcov))
    
    numreps <-
      lsyear %>%
      select(site, block) %>%
      unique() %>%
      group_by(site) %>%
      summarize(num = length(block))
    
    lsyear2.2 <-
      lsyear2 %>%
      left_join(numreps) %>%
      filter(num > 2)
    
    return(lsyear2.2)
  }


run_grazing_bootstrap <-
  function(lsyear2.2, sample_threshold = 6) {
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
        select(site,
               site_block,
               exage,
               block,
               trt,
               num,
               genus_species_use,
               relcov) %>%
        spread(genus_species_use, relcov, fill = 0)
      
      # Species relative abundance matrix
      species <-
        as.data.frame(wide[, 7:ncol(wide)])
      
      # Treatment vector
      env <-
        as.data.frame(wide[, c("trt", "block")])
      
      # Turn block into numeric var
      env$block <- as.integer(factor(as.factor(env$block)))
      
      # Add an indicator for original (not bootstrapped) samples
      env$type <- "original"
      
      # What is the site level replication?
      replication <- wide$num[1]
      
      # Calculate how many bootstrap samples needed
      if (replication < sample_threshold) {
        replication_needed <-
          sample_threshold - replication
        
        
        # Mean relative abundances of Grazing plots
        G_prob <-
          colMeans(as.data.frame((wide %>%
                                    filter(trt == "G"))[, 7:ncol(wide)]))
        
        # Mean relative abundances of Ungrazed plots
        U_prob <-
          colMeans(as.data.frame((wide %>%
                                    filter(trt == "U"))[, 7:ncol(wide)]))
        
        # Note: assumption of 100 individuals per plot
        # Bootstrap Grazed plots
        G_boot <-
          t(rmultinom(
            n = replication_needed,
            size = 100,
            prob = G_prob
          ))
        # Bootstrap Ungrazed plots
        U_boot <-
          t(rmultinom(
            n = replication_needed,
            size = 100,
            prob = U_prob
          ))
        
        # Bind bootstrapped communities to original
        species <- rbind(species, G_boot, U_boot)
        
        # Construct a data frame for the descriptions of the bootstrapped samples
        boot_env <-
          data.frame(trt = c(
            rep("G", replication_needed),
            rep("U", replication_needed)
          ),
          block = as.factor(c(rep(seq((max(env$block) + 1),
                                      (max(env$block) + replication_needed)
          ),
          2))),
          type = "boot")
        
        # Bind bootstrapped and original samples
        env <- rbind(env, boot_env)
        
      } else if (replication >= sample_threshold) {
        replication_needed <- 0
      }
      
      # Make plot for sense check
      nmds_plot(species, env, file_prefix = i)
      
      # Perform permamova
      # Try with and without `strata =` option
      out <-
        adonis(
          species ~ trt,
          data = env,
          method = "bray",
          strata = env$block
        )
      
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


gex_permanova <- run_grazing_bootstrap(process_comp_data())
write.csv(gex_permanova, file = "composition/bootstrap_permanova/out/permanova_bootstrap.csv")
swoosh_plot(gex_permanova)
