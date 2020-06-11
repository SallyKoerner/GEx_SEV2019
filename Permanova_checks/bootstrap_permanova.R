library(tidyverse)
library(codyn)
library(rsq)
library(vegan)
library(ggplot2)

setwd("/Users/avahoffman/Dropbox/Research/Grazing_consortium/2020/GEx_SEV2019")

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

nmds_plot <-
  function(species, env, file_prefix = "out") {
    ## look at NMDS to see gut check
    test <- metaMDS(species, distance = "bray")
    
    scores <-
      data.frame(scores(test, display = "sites"))  # Extracts NMDS scores for year "i" #
    scores2 <-
      cbind(env, scores) # binds the NMDS scores of year i to all years previously run
    
    ggplot(scores2,
           aes(
             x = NMDS1,
             y = NMDS2,
             color = trt,
             shape = type
           )) +
      geom_point(size = 5) +
      scale_shape_manual(values = c(1, 15))
    
    ggsave(
      paste("Permanova_checks/nmds/", file_prefix, ".pdf", sep = ""),
      height = 3,
      width = 4
    )
  }


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


gex_permanova <- run_grazing_bootstrap()
write.csv(gex_permanova, file = "Permanova_checks/permanova_bootstrap.csv")


swoosh_plot <-
  function(gex_permanova) {
    ###making comp diff shoowsh figure
    comp <- read.csv("gex_multdiff_site_block.csv")
    
    #getting CI for each site
    gex_multdiff_CI <- comp %>%
      group_by(site) %>%
      summarize_at(vars(composition_diff), funs(mean, sd, length)) %>%
      mutate(se = sd / sqrt(length),
             ci = se * 1.96)
    
    #geting global mean comp diff
    mean <- mean(gex_multdiff_CI$mean)
    sd <- sd(gex_multdiff_CI$mean)
    se <- sd / sqrt(252)
    ci <- se * 1.96
    site <- "All Sites"
    
    all <- data.frame(site, mean, sd, se, ci)
    
    gex_multdiff_CI <-
      gex_multdiff_CI[order(gex_multdiff_CI$mean),]  # sort
    gex_multdiff_CI2 <- gex_multdiff_CI %>%
      bind_rows(all) %>%
      mutate(colortrt = ifelse(site == "All Sites", 1, 0))
    gex_multdiff_CI2$site2 <-
      factor(gex_multdiff_CI2$site, levels = gex_multdiff_CI2$site)
    
    
    ggplot(data = gex_multdiff_CI2,
           aes(
             x = site2,
             y = mean,
             color = as.factor(colortrt)
           )) +
      geom_point() +
      scale_color_manual(values = c("black", "red")) +
      geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci)) +
      coord_flip() +
      xlab("Site") +
      ylab("Compositional Difference") +
      #scale_y_continuous(limits = c(0,1))+
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
      )
    
    ###making comp diff sig figure based on permanovas
    
    ####flagging sig values
    permsig <- gex_permanova %>%
      mutate(sig = ifelse(perm_Pvalue < .05, 1, 0))
    
    toplot <- gex_multdiff_CI %>%
      right_join(permsig)
    
    toplot$site2 <-
      factor(toplot$site, levels = toplot$site[order(toplot$mean)])
    
    ggplot(data = toplot, aes(
      x = site2,
      y = mean,
      color = as.factor(sig),
      shape = as.factor(boot_n)
    )) +
      theme_classic() +
      geom_point() +
      scale_color_manual(
        values = c("black", "red"),
        name = "Permanova p < 0.05",
        labels = c("No", "Yes")
      ) +
      scale_shape_manual(
        name = "Num. bootstrapped plots",
        values = c(1,2,3,4)
      ) +
      geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci)) +
      coord_flip() +
      xlab("Site") +
      ylab("Compositional Difference") +
      #scale_y_continuous(limits = c(0,1))+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      theme(legend.position = c(0.95,0.05),
            legend.justification = c(1,0))
    
    ggsave(file = "Permanova_checks/swooshplot_bootstrap.pdf",
           height = 6,
           width = 8)
  }
