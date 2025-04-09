####### Data preparation II ######################################################
#* - exlusion of all outliers, infinite residuals & influencing values based on
#*   gumbel model residuals
#################################################################################'
#* input: foveated actions goals: only fixations on action goals -> "data_ml3fix_excl.rds"
#* output: final prepared data set -> "data_ml3fix_final.rds"

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


# Data Preparation ###########################################################
# Final Data Preparation for Gumbel based Models, incl. outlier analysis

# import from file
fix.pre <- readRDS("data_ml3fix_excl.rds")

# exclude outliers after 99th percentile
(excl.quantile <- quantile(fix.pre$distance_to_spaceship, .99))
fix.pre <- fix.pre[fix.pre$distance_to_spaceship <= excl.quantile, ]


# summary
summary(fix.pre$distance_to_spaceship); cat("rows:", nrow(fix.pre))
# plot(density(fix.pre$distance_to_spaceship))


# Modelling ####################################################################

# truncate GU distribution
gen.trun(par=range(fix.pre$distance_to_spaceship),
         type="both", family=GU)


# NULL Model for Outlier Analysis
mod.pre <- gamlss(distance_to_spaceship ~ 1, family = GUtr, data = fix.pre)


## Inspect pre-Fit ####
# plot(mod.pre)
# wp(resid=resid(mod.pre))


# Outlier Analysis & Exclusion #################################################
fix.pre$res <- resid(mod.pre)


# find infinte residuals & exclude
(inf <- which(is.infinite(fix.pre$res)))
if(length(inf)>0) fix.pre <- fix.pre[-inf,]

# linaear model with residuals
mod.resid <- lm(res ~ 1, data = fix.pre)

# find quantile outliers
(quant.out <- which(rstudent(mod.resid) <= -4 | rstudent(mod.resid) >= 4)) %>% unname()


# Influence Check by had & cooks distance
(infl <- car::influencePlot(mod.resid) %>% rownames() %>% as.numeric())
# car::infIndexPlot(mod.resid)

# Exclude influencial data
(out <- unique(c(quant.out, infl)))


# FINAL DATA
data_ml3fix_final <- fix.pre[-out,]

# SAVE DATA
saveRDS(data_ml3fix_final, file="data_ml3fix_final.rds")
