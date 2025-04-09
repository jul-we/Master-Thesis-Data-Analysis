####### Modelling with GAMLSS, truncated Gumbel Distribution #####
#* with "data_ml3fix_final.rds"
library(dplyr)
library(tibble)
library(magrittr)
library(conflicted)
conflict_prefer("filter", "dplyr")
library(gamlss)
library(gamlss.tr)

# CLEAR WORKSPACE
store <- c("store")
rm(list=ls()[-which(ls()==store)]); graphics.off()

# WORKING DIRECTORY
try(docPath<-dirname(rstudioapi::documentPath()),T);try(setwd(docPath))

# Data Import ##################################################################'
ml3fix <- readRDS("data_ml3fix_final.rds")

# Function #####################################################################
# with rerun, the model fitting process is optimized. Models are repeatedly fitted
# with increasing convergence critera.
rerun <- function(model, re.crit, cyc=30) {
	# model name
	model.name <- deparse(substitute(model))
	# init refit
	for (c in re.crit) {
		contr <- gamlss::gamlss.control(c.crit=c, n.cyc=cyc)
		i.contr <- gamlss::glim.control(cc=c)
		message("# CRIT = ", c)
		try(model <- refit(model))
		if (!model$converged) try(model <- gamlss::refit(model))
		message("# Converged:", model$converged, "\n")
	}
	# assign fitted model
	return(model)
}


# Set options for modelling #####################################################
# truncate Gumbel distribution
trunc.vals <- c(5.073984, 11.607166)			# important to use same values for truncation as in data preparation II
gen.trun(par = trunc.vals, type = "both", family=GU)


# gamlss Control options to optimize fitting
init.crit <- 0.9					# initial convergence criterium
re.crit <- c(0.05, 0.002)				# rerun convergence criterium
contr <- gamlss.control(c.crit=init.crit, n.cyc=30)
i.contr <- glim.control(cc=init.crit)
#############################################################################'

# Models ####
mods <- list()

#* input noise as categorical
# C1
mods$c1 <- gamlss(distance_to_spaceship ~ I(N_visible_obstacles^2)+input_noise_f +
	       	re(random=~N_visible_obstacles | ID, method = "REML"),
	       sigma.formula=~-1 + N_visible_obstacles + input_noise_f + level_difficulty + 
	       	re(random=~N_visible_obstacles | ID, method = "REML"),
	       family = GUtr, data = ml3fix, control=contr, i.control=i.contr)

mods$c1 %<>% rerun(re.crit=re.crit)


# C2 - ***
mods$c2 <- gamlss(distance_to_spaceship ~ I(N_visible_obstacles^2):input_noise_f +
	       	re(random=~N_visible_obstacles | ID, method = "REML"),
	       sigma.formula=~-1 + N_visible_obstacles + input_noise_f + level_difficulty + 
	       	re(random=~N_visible_obstacles | ID, method = "REML"),
	       family = GUtr, data = ml3fix, control=contr, i.control=i.contr)

mods$c2 %<>% rerun(re.crit=re.crit)


# C3
mods$c3 <- gamlss(distance_to_spaceship ~ I(N_visible_obstacles^2)*input_noise_f +
	       	re(random=~N_visible_obstacles | ID, method = "REML"),
	       sigma.formula=~-1 + N_visible_obstacles + input_noise_f + level_difficulty + 
	       	re(random=~N_visible_obstacles | ID, method = "REML"),
	       family = GUtr, data = ml3fix, control=contr, i.control=i.contr)

mods$c3 %<>% rerun(re.crit=re.crit)


#* input noise as metrix
# m1
mods$m1 <- gamlss(distance_to_spaceship ~ I(N_visible_obstacles^2) + input_noise +
	       	re(random=~N_visible_obstacles | ID, method = "REML"),
	       sigma.formula=~-1 + N_visible_obstacles + input_noise_f + level_difficulty + 
	       	re(random=~N_visible_obstacles | ID, method = "REML"),
	       family = GUtr, data = ml3fix, control=contr, i.control=i.contr)

mods$m1 %<>% rerun(re.crit=re.crit)

# m2
mods$m2 <- gamlss(distance_to_spaceship ~ I(N_visible_obstacles^2) + I(input_noise^2) +
	       	re(random=~N_visible_obstacles | ID, method = "REML"),
	       sigma.formula=~-1 + N_visible_obstacles + input_noise_f + level_difficulty + 
	       	re(random=~N_visible_obstacles | ID, method = "REML"),
	       family = GUtr, data = ml3fix, control=contr, i.control=i.contr)

mods$m2 %<>% rerun(re.crit=re.crit)


# Model comparisson
k = sqrt(log(nrow(ml3fix)))	# set k according to Stasinopoulus, S.379, 2.5 < AIC < 4
gaic.tab <- with(mods, GAIC.table(c1, c2, c3, m1, m2, k=k))

tab <- tibble(mod=rownames(gaic.tab), df=gaic.tab[,1], GAIC=gaic.tab[,2])
tab
arrange(tab, tab$GAIC) # sort by GAIC
lapply(mods, deviance) %>% unlist() %>% enframe() # deviance

LR.test(mods$c2, mods$c3)
LR.test(mods$m1, mods$c2)
#* LR Results
#  - c3 not better than c2 
#  - c2 better than m1 and m2

wp(mods$c2, ylim.all=0.2)
plot(mods$c2)


# -> C2 for bootstrapping

# Save Objects for Bootstrapping ##############################################
c2 <- mods$c2
save(c2,
     rerun,
     dGUtr, GUtr, pGUtr, qGUtr,
     init.crit,
     contr,
     i.contr,
     ml3fix,
     file="boot_input_data.RData")
