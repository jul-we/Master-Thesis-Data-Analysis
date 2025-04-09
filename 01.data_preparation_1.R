####### Data preparation II ##############
# Only specific data is used in the data analysis. It's called foveated action goals,
# which means fixations on the agent. The data set is prepared and foveated action
# goals are filtered by cluster analysis.

library(dplyr)
library(tibble)
library(magrittr)
library(conflicted)
library(rio)
conflict_prefer("filter", "dplyr")
library(ggplot2)
library(QuClu)	# Quantile based clusteranalysis

# CLEAR WORKSPACE
rm(list=ls()); graphics.off(); cat("\014")

# WORKING DIRECTORY
try(docPath<-dirname(rstudioapi::documentPath()),T);try(setwd(docPath))


# Import ################################################################
ml3fix <- import("fixations_moonlander_iii.csv") %>% as_tibble()


## Prepare ########

# only fixations where obstacles are present
# only fixations with min 0.0125ms fixation duration
ml3fix %<>% filter(N_visible_obstacles > 0, fixation_duration >= 0.0125)


# level difficulty as factor
ml3fix$level_difficulty %<>%  factor(levels=c("medium", "hard"))


# add input noise as factor
ml3fix$input_noise_f <- factor(ml3fix$input_noise)
ml3fix$input_noise_fo <- factor(ml3fix$input_noise, ordered = TRUE)	# ordered factor


# select variables and omit NA's
ml3fix.clust <- ml3fix[c("distance_to_spaceship",
	     "N_visible_obstacles",
	     "input_noise",
	     "input_noise_f",
	     "input_noise_fo",
	     "ID",
	     "level",
	     "level_difficulty",
	     "fixation_duration")] %>% na.omit()


# Clusteranalysis {QuClu} ######################################################
set.seed(1)					# make it reproduceable
cla <- kquantiles(ml3fix.clust[1:3], k=2, method="CS")	# call cluster analysis
ml3fix.clust$cluster <- as.factor(cla$cl)		# add cluster to data



# Clusters Summarys
ml3fix.clust %>% filter(cluster==1) %>% pull(distance_to_spaceship) %>% summary()
ml3fix.clust %>% filter(cluster==2) %>% pull(distance_to_spaceship) %>% summary()


# Exclusion ###################################################################

# exclude cluster 1: fixations on agent
ml3fix.exlc <- ml3fix.clust %>% filter(cluster==2)
ml3fix.exlc$cluster %<>% droplevels()


# exclusion behind game display border
ml3fix.exlc %<>% filter(distance_to_spaceship <= 16.638)


# SAVE RDS #####################################################################
saveRDS(ml3fix.exlc, file="data_ml3fix_excl.rds")
