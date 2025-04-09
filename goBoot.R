####### Interrupted Bootsampling of model c2 ##############
# This script allows to interrupt the bootsamping process and to save results
# at any time because. The process can the be continued.
# it needs about 5 to 7 days.
library(dplyr)
library(tibble)
library(magrittr)
library(conflicted)
library(progress)
conflict_prefer("filter", "dplyr")
library(cli)
library(gamlss)

# CLEAR WORKSPACE
rm(list=ls()); graphics.off(); cat("\014")

# WORKING DIRECTORY
if (rstudioapi::isAvailable()) {
	try(docPath<-dirname(rstudioapi::documentPath()),T)
	try(setwd(docPath))
}


# README ######################################################################
#* RUN SCRIPT
#* 1. Call R in Terminal: R --vanilla
#* 2. run script with: source("goBoot")
#* 3. type 'sample' to start when asked
################################################################################'
#* - Object and file names are hard coded and should not be altered.


# DATA IMPORT #################################################################

# Define file names
input.data <- "boot_input_data.RData"		# file containing data to analyse. Other objects can also be part of that file
boot.data.file <- "boot_output_data.RData"	# file to store the resampled data 


# load existing data or create
data.obj <- list()

# Backup boot.data

if (file.exists(boot.data.file)) 
	cat("# Found ", boot.data.file, "!\n# Creating Backup now....", sep="")
	file.copy(from=boot.data.file, to=paste(boot.data.file,".bak",sep=""), overwrite = TRUE)


data.obj$input_data <- load(input.data) # loads objects and stores names
if (file.exists(boot.data.file)) data.obj$boot_data <- load(boot.data.file)[1] else {	# loads existing objects and stores it's names
	boot.data <- list()							# else:  creates object and saves name
	data.obj$boot_data <- deparse(substitute(boot.data))
}


# Sampling Function ############################################################
boot.sampling <- function(inp.data, out.data, K) {
	# inp.data - a data frame
	# out.data - a matrix as replicate return
	# K - number of bootstrap samples
	file <- file("boot.log", open = "a")
	sink(nullfile(), type = "output") 		# hide output
	
	# calculate K_pending samples
	if (is.null(dim(out.data$par))) K_achieved <- 0 else K_achieved <- dim(out.data$par)[2]
	K_pending <- K - K_achieved	# pending samples
	
	# Progressbar
	pbar <- progress_bar$new(total = K, clear = FALSE, width = 100,
			     format = "fitting Nr: :current | last fit: :what :spin [:bar] :percent finished in :eta")
	pbar$tick(K_achieved)
	
	if (K_pending==0) { 	# Abort if nothing to do
		sink()
		stop(paste("K_pending=0, Sampling Number", K, "reached!"), call.=FALSE)
	}
	
	out.data <- substitute(out.data) # store name bootstrap resamples object
	boot.t0 <<- proc.time()
	replicate(K_pending, {
		sink(file, append = TRUE, type = "message")	# to boot.log
		message(timestamp(quiet = TRUE))
		# INDICIES ##############################################'
		i <- sample(1:nrow(inp.data), replace = TRUE)
		boot.sample <<- inp.data[i, ] # gamlss search data argument in global env (bug)
		# CONTROLS ##############################################'
		init.crit <- 0.9
		contr <- gamlss.control(c.crit=init.crit, n.cyc=30)
		i.contr <- glim.control(cc=init.crit)
		# MODEL #################################################'
		try(boot.mod <- update(c2, data=boot.sample, start.from=c2, control=contr, i.control=i.contr), silent=FALSE)
		if (exists("boot.mod")) boot.mod <- rerun(boot.mod, re.crit=c(0.05, 0.005))
		#########################################################'
		# remove model if not converged
		try(if (!boot.mod$converged) {
			rm(boot.mod)
			message("# -> NOT converged & removed!")
		}, silent=TRUE)
		# Extract parameters ####################################'
		if(exists("boot.mod")) result <- boot.mod$mu.coefficients[-2]
		# Assign to out.data ####################################'
		if (exists("boot.mod")) {
			tmp.list <- eval(out.data) # copy boot.data list
			tmp.list$par <- cbind(tmp.list$par, result, deparse.level = 0) # cbind parameters to list
			# Calc Regression lines #################################'
			x <- seq(0, 20, 0.01)
			tmp.list$x <- x
			tmp.list$y0 <- cbind(tmp.list$y0, coef(boot.mod)[1] + coef(boot.mod)[3]*x^2)
			tmp.list$y0.5 <- cbind(tmp.list$y0.5, coef(boot.mod)[1] + coef(boot.mod)[4]*x^2)
			tmp.list$y1 <- cbind(tmp.list$y1, coef(boot.mod)[1] + coef(boot.mod)[5]*x^2)
			tmp.list$y1.5 <- cbind(tmp.list$y1.5, coef(boot.mod)[1] + coef(boot.mod)[6]*x^2)
			tmp.list$y2 <- cbind(tmp.list$y2, coef(boot.mod)[1] + coef(boot.mod)[7]*x^2)
			# assign tmp.list back to boot.data
			assign(deparse(out.data),
			       tmp.list,
			       envir=globalenv()
			       )
		}
		# Show Progressbar
		sink(type = "message"); pbar$tick(tokens = list(what = format(Sys.time(), "%H:%M:%S")))	# increase progress bar
	}
	)
	inter()
}



# rerun function ###############################################################
rm(rerun)
rerun <- function(model, re.crit, cyc=30) {
	# model name
	model.name <- deparse(substitute(model))
	# init refit
	for (c in re.crit) {
		# rerun options
		contr <- gamlss::gamlss.control(c.crit=c, n.cyc=cyc)
		i.contr <- gamlss::glim.control(cc=c)
		# 1. rerun
		message("# CRIT = ", c)
		model$converged <- FALSE # to catch singularity error
		try(model <- gamlss::refit(model), silent=FALSE)
		# POLICY: if not converged BREAK - never seccond chance
		message("# Converged with ", c, ": ", model$converged, "\n")
		if (!model$converged) break
	}
	
	# Final messages
	try(if (model$converged) message("# Finally Converged!\n") else message("# Finally NOT Converged!\n"))
	# return fitted model
	try(return(model))
}


# Interrupt Function ###########################################################
inter <- function(...) {
	# store bootsampling calculation time, create if not exist
	diff.t.last <- (proc.time()-boot.t0)[3]
	if (hasName(boot.data, "time")) boot.data$time %<>% + diff.t.last else boot.data$time <- diff.t.last
	for (s in sink.number()) sink()
	sink(type = "message")
	message("\n\n#### Interrupt by User! Wait until Data is saved...")
	# close(file)
	save(boot.data, file=boot.data.file)
	cli_h1(col_red("Data saved & Execution terminated"))
	# Show sampling Times
	show_last <- difftime(diff.t.last, 0, units="auto")
	cli_li(paste("Last Sampling Time:", show_last, attributes(show_last)[3]))
	show_total <- difftime(boot.data$time, 0, units="auto")
	cli_li(paste("Total Sampling Time:", show_total, attributes(show_total)[3]))
}


# SCRIPT LOOP ##################################################################
K = 9999

## Feedback ####
if (Sys.info()[1] == "Windows") shell("cls") else system2("clear")

cli_h1("Welcome to interrupted Bootstrapping")
cli_h2(col_cyan("Elemtents loaded"))
print(data.obj)

cli_h2(col_cyan("Boot Data Properties"))
cli_h3(col_br_blue("Planed Sampling Number K"))
cli_li(K)

if( is.null(boot.data$par)) cli_li(paste("boot.data$par", "is empty")) else {
	cli_h3(col_br_blue("Achived K-Samples"))
	cli_li(ncol(boot.data$par))
	cli_li(paste(ncol(boot.data$par)/K*100,"%", sep=""))
	cli_h3(col_br_blue("Parameter Names"))
	cli_li(rownames(boot.data$par))
	cli_h3(col_br_blue("Sum of sampling Time"))
	diff_t <- difftime(boot.data$time, 0, units="auto")
	cli_li(paste(diff_t, attributes(diff_t)[3]))
}
cli_rule()
input <- readline(col_yellow("-> Enter 'sample' to continue: "))
stopifnot(input=="sample")


# Script Call ####
# tryCatch(boot.sampling(studyData, boot.data, K), interrupt = inter)
tryCatch(boot.sampling(ml3fix, boot.data, K), interrupt = inter)

