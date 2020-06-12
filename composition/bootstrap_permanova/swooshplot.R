# This script contains code for making the "swooshplot" -- the magnitude of community
# composition change by site, highlighted by whether the bootstrapped permanova has p<0.05
###########################################################################################
###########################################################################################


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
      gex_multdiff_CI[order(gex_multdiff_CI$mean), ]  # sort
    gex_multdiff_CI2 <- gex_multdiff_CI %>%
      bind_rows(all) %>%
      mutate(colortrt = ifelse(site == "All Sites", 1, 0))
    gex_multdiff_CI2$site2 <-
      factor(gex_multdiff_CI2$site, levels = gex_multdiff_CI2$site)
    
    
    ####flagging sig values
    permsig <- gex_permanova %>%
      mutate(sig = ifelse(perm_Pvalue < .05, 1, 0),
             is_boot = ifelse(boot_n > 0, 1, 0))
    
    toplot <- gex_multdiff_CI %>%
      right_join(permsig)
    
    toplot$site2 <-
      factor(toplot$site, levels = toplot$site[order(toplot$mean)])
    
    ggplot(data = toplot,
           aes(
             x = site2,
             y = mean,
             color = as.factor(sig),
             shape = as.factor(is_boot)
           )) +
      theme_classic() +
      geom_point() +
      scale_color_manual(
        values = c("black", "red"),
        name = "Permanova p < 0.05",
        labels = c("No", "Yes")
      ) +
      scale_shape_manual(
        name = "Bootstrapped plots",
        values = c(16, 1),
        labels = c("No", "Yes")
      ) +
      geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci)) +
      coord_flip() +
      xlab("Site") +
      ylab("Compositional Difference") +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      theme(legend.position = c(0.95, 0.05),
            legend.justification = c(1, 0))
    
    ggsave(file = "composition/bootstrap_permanova/out/swooshplot_bootstrap.pdf",
           height = 6,
           width = 8)
  }
