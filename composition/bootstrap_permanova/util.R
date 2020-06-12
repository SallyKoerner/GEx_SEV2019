# This script contains utils for performing the bootstrap permanova analysis and output
###########################################################################################
###########################################################################################


nmds_plot <-
  function(species,
           env,
           file_prefix = "out") {
    # Generate NMDS plots for comparing grazed and ungrazed, with indicator
    # of bootstrapped or not
    #
    # species : dataframe, species matrix (plot X species)
    # env : dataframe, environmental vars per plot
    # file_prefix : should be the site name to use as the file name when saved
    
    # Perform NMDS
    test <- metaMDS(species, distance = "bray")
    
    # Gather scores
    scores <-
      data.frame(scores(test, display = "sites"))  # Extracts NMDS scores for year "i" #
    scores2 <-
      cbind(env, scores) # binds the NMDS scores of year i to all years previously run
    
    # Make plot
    ggplot(scores2,
           aes(
             x = NMDS1,
             y = NMDS2,
             color = trt,
             shape = type
           )) +
      geom_point(size = 5) +
      scale_shape_manual(values = c(1, 15))
    
    # Save image as site name
    ggsave(
      paste("composition/bootstrap_permanova/out/nmds/", file_prefix, ".pdf", sep = ""),
      height = 3,
      width = 4
    )
  }