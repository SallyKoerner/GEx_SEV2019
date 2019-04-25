################################################################################
##  modified_Codyn_metrics: Calculating community difference metrics, modified from codyn functions to give ln RR (making it all comparable to the dominance output).
##
##  Author: Kimberly Komatsu (modified from Meghan Avolio)
##  Date created: April 25, 2019
################################################################################

library(tidyverse)
library(codyn)
library(rsq)

#####IMPORTANT: must run utilities file first (see github), otherwise modifying codyn function won't work
##update SERSp to get lnRR----------------------
SERSp2 <- function (df, species.var, abundance.var, abundance.var2) 
{
  out <- c(species.var, "rank", "rank2", abundance.var, abundance.var2)
  out <- unique(df[!(names(df) %in% out)])
  df <- subset(df, df[[abundance.var]] != 0 | df[[abundance.var2]] != 
                 0)
  s_r1 <- S(df[[abundance.var]])
  e_r1 <- Evar(as.numeric(df[[abundance.var]]))
  s_r2 <- S(df[[abundance.var2]])
  e_r2 <- Evar(as.numeric(df[[abundance.var2]]))
  sdiff <- log(s_r2/s_r1)
  ediff <- log(e_r2/e_r1)
  spdiff <- df[df[[abundance.var]] == 0 | df[[abundance.var2]] == 
                 0, ]
  spdiffc <- nrow(spdiff)/nrow(df)
  rank_diff <- mean(abs(df[["rank"]] - df[["rank2"]]))/nrow(df)
  metrics <- data.frame(richness_diff = sdiff, evenness_diff = ediff, 
                        rank_diff = rank_diff, species_diff = spdiffc)
  return(cbind(out, metrics))
}

#update RAC_difference function to get lnRR
RAC_lnRR_difference <- function (df, time.var = NULL, species.var, abundance.var, replicate.var, 
                                 treatment.var = NULL, pool = FALSE, block.var = NULL, reference.treatment = NULL) 
{
  args <- as.list(match.call()[-1])
  df <- do.call(check_args, args, envir = parent.frame())
  if (pool) {
    rankdf <- pool_replicates(df, time.var, species.var, 
                              abundance.var, replicate.var, treatment.var)
  }
  else {
    by <- c(block.var, time.var)
    allsp <- split_apply_combine(df, by, FUN = fill_species, 
                                 species.var, abundance.var)
    by <- c(block.var, time.var, treatment.var, replicate.var)
    rankdf <- split_apply_combine(allsp, by, FUN = add_ranks, 
                                  abundance.var)
  }
  if (!is.null(block.var)) {
    cross.var <- treatment.var
  }
  else if (pool) {
    cross.var <- treatment.var
  }
  else {
    cross.var <- replicate.var
  }
  to_ordered = is.factor(rankdf[[cross.var]]) & !is.ordered(rankdf[[cross.var]])
  if (to_ordered) {
    class(rankdf[[cross.var]]) <- c("ordered", class(rankdf[[cross.var]]))
  }
  split_by <- c(block.var, time.var)
  merge_to <- !(names(rankdf) %in% split_by)
  cross.var2 <- paste(cross.var, 2, sep = "")
  if (is.null(reference.treatment)) {
    ranktog <- split_apply_combine(rankdf, split_by, FUN = function(x) {
      y <- x[merge_to]
      cross <- merge(x, y, by = species.var, suffixes = c("", 
                                                          "2"))
      idx <- cross[[cross.var]] < cross[[cross.var2]]
      cross[idx, ]
    })
  }
  else {
    ranktog <- split_apply_combine(rankdf, split_by, FUN = function(x) {
      y <- x[x[[treatment.var]] != reference.treatment, 
             merge_to]
      x <- x[x[[treatment.var]] == reference.treatment, 
             ]
      merge(x, y, by = species.var, suffixes = c("", "2"))
    })
  }
  if (to_ordered) {
    x <- class(ranktog[[cross.var]])
    class(ranktog[[cross.var]]) <- x[x != "ordered"]
    class(ranktog[[cross.var2]]) <- x[x != "ordered"]
  }
  idx <- is.na(ranktog[[abundance.var]])
  abundance.var2 <- paste(abundance.var, 2, sep = "")
  idx2 <- is.na(ranktog[[abundance.var2]])
  ranktog[idx, abundance.var] <- 0
  ranktog[idx2, abundance.var2] <- 0
  idx <- ranktog[[abundance.var]] != 0 | ranktog[[abundance.var2]] != 
    0
  ranktog <- ranktog[idx, ]
  split_by <- c(block.var, time.var, cross.var, cross.var2)
  output <- split_apply_combine(ranktog, split_by, FUN = SERSp2, 
                                species.var, abundance.var, abundance.var2)
  if (any(is.na(output$evenness_diff))) 
    warning(paste0("evenness_diff values contain NAs because there are plots", 
                   " with only one species"))
  output_order <- c(time.var, block.var, replicate.var, paste(replicate.var, 
                                                              2, sep = ""), treatment.var, paste(treatment.var, 2, 
                                                                                                 sep = ""), "richness_diff", "evenness_diff", "rank_diff", 
                    "species_diff")
  return(output[intersect(output_order, names(output))])
}



#kim's laptop--------------
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\data')


dat<-read.csv("All_Cleaned_April2019_V2.csv")

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


###doing RAC differences

gex_RACdiff<-data.frame()
siteblock<-unique(lsyear2$site_block)

for (i in 1:length(siteblock)){
  
  subset<-lsyear2%>%
    filter(site_block==siteblock[i]) %>%
    mutate(treatment=trt)
  
  out <- RAC_lnRR_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "trt", treatment.var = "treatment")
  
  out$site_block<-siteblock[i]
  
  gex_RACdiff<-rbind(gex_RACdiff, out)  
}

gex_RACdiff2<-data.frame()
siteblockplot<-unique(lsyear3$site_block_plot)

for (i in 1:length(siteblockplot)){
  
  subset<-lsyear3%>%
    filter(site_block_plot==siteblockplot[i]) %>%
    mutate(treatment=trt)
  
  out <- RAC_lnRR_difference(subset, species.var="genus_species", abundance.var = "relcov", replicate.var = "trt", treatment.var = "treatment")
  
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