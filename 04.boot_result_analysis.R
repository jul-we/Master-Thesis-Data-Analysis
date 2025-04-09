####### boot.data analysis ##############
# Bootstrapping was done in a separate process with k=9999 bootstrap resamples.
# A script was used to allow the process to be interrupted and resumed (goBoot.r).
# It took about one week computing time.
# This script performs data analysis on the bootstrap resamples.

library(dplyr)
library(tidyr)
library(tibble)
library(magrittr)
library(conflicted)
conflict_prefer("filter", "dplyr")
library(ggplot2)
library(gamlss)
library(patchwork)

# CLEAR WORKSPACE
rm(list=ls()); graphics.off(); cat("\014")

# WORKING DIRECTORY
try(docPath<-dirname(rstudioapi::documentPath()),T);try(setwd(docPath))

# MAIN ##########################################################################
load("boot_output_data.RData")
load("boot_input_data.RData")


# BOOTSTRAP #####################################################################

# bootstrap analysis
theta.mean <- apply(boot.data$parameters, 1, mean)
theta.ci <- apply(boot.data$parameters, 1, quantile, probs=c(0, 0.95)) # one-tailes test (uncorrected)

# Bonferroni holm: Order differences by size (differences from CI boundaries to test-values, because
# we have no p-values)
diff <- c(0, 0, theta.mean[2:5]) - theta.ci[2,]
names(diff) <- names(theta.mean)
diff
sort(diff[-1], T)

# calculate all BH Adjusment
theta.ci.bh <- apply(boot.data$parameters, 1, quantile,
	        probs=c(0, 1-0.05/5, 1-0.05/4, 1-0.05/3, 1-0.05/2, 1-0.05/1)) # one-tailed hypothesis, korrigiert nach Bonferroni

# select needed values, according to sorted sequence
theta.ci <- rbind(lower=theta.ci[1,], upper=theta.ci.bh[c(6, 8, 17, 22, 30, 33)])


# PLOT Effectcs ################################################################

## Plot Effects against 0
barplot(theta.mean[-1], ylim=c(-0.025, 0), main="One-tailed CI (95%) of input_noise\n(Bonferroni Adjustment)", cex.main=.9)
arrows(x0=seq(0.75, 5.75, 1.2),
       y0=theta.ci[1, -1],
       y1=theta.ci[2, -1], angle=90, code=3, lwd=2, col=c(4,1,1,1,4), length=.2)
abline(h=c(theta.ci[2, 6]), lty="dotted", col=4)


## Plot Effects according to succesive differences
adj.mean <- theta.mean[-1] - c(0, theta.mean[2:5])
ci.adj <- adj.mean - theta.mean[-1]
barplot(adj.mean, ylim=c(-0.018, 0.01), main="One-tailed CI (95%) of input_noise\n(Bonferroni Adjustment)", cex.main=.9)
arrows(x0=seq(0.75, 5.75, 1.2),
       y0=theta.ci[1, -1] + ci.adj,
       y1=theta.ci[2, -1] + ci.adj, angle=90, code=3, lwd=2, col=c(4,1,1,1,4), length=.2)


## + Effects with ggplot ####
comment <- c("int", "*1", "", "", "", "*2")
plot.dat <- tibble(Estimate=theta.mean, Parameters=names(theta.mean),
	         lower=theta.ci[1,], upper=theta.ci[2,], comment=comment)

# +++sequencial contrast correct version
library(ggtext)
plot.dat <- tibble(Estimate=c(0, adj.mean), Parameters=names(theta.mean),
	         lower=theta.ci[1,]+c(0,ci.adj), upper=theta.ci[2,]+c(0,ci.adj), comment=comment)

# pdf(file="results_parameters.pdf", width=7, height=5)
ggplot(plot.dat[-1,]) + geom_col(aes(x=Parameters, Estimate), show.legend=F, fill=4) +
	geom_hline(aes(yintercept=0), colour="slategray", size=0.8)+
	geom_errorbar(aes(x=Parameters, ymin=lower, ymax=upper, color=Parameters), linewidth=1.1, width=.5, 
		    show.legend = FALSE) +
	scale_color_manual(values=c("turquoise", rep("slategray",3), "turquoise")) +
	ylim(-0.018, .005) + geom_text(aes(label=comment, x=Parameters, y=-0.017), size=5) +
	labs(title="Estimated parameters", subtitle="(95% one-tailed CI, Bonferroni-Holm adjusted)",
	     caption="**Significant Effects**<br>
	     *1: 0 > input_noise_0<br>
	     *2: input_noise_1.5 > input_noise_2") +
	theme(plot.caption = element_markdown(lineheight = 1.3), legend.position = "bottom")

# dev.off()





# Regression Plot incl. Confidence Bands #######################################
y <- c("inp_noise_0", "inp_noise_0.5", "inp_noise_1", "inp_noise_1.5", "inp_noise_2")

reg.lines.df <- NULL
for (i in seq_along(y)) {
	reg.lines <- apply(X=boot.data[[y[i]]], MARGIN=1, FUN=quantile, probs=c(0.025, 0.5, 0.975))
	rownames(reg.lines)[2] <- 'distance'
	reg.lines.df <- bind_rows(reg.lines.df,
			      tibble(obstacles=boot.data$n_vis_obstacles,
			             as_tibble(t(reg.lines)), Level=y[i]))
}
# Rename confidence bands
attributes(reg.lines.df)$names[c(2,4)] <- c("lower", "upper")


# as ggplot
palette.colors(palette="ggplot2")
rib.alpha = .3

p1 <- ggplot(reg.lines.df) +
	geom_ribbon(aes(obstacles, ymin=lower, ymax=upper, group=Level, fill=Level),
		  alpha=rib.alpha, linewidth=.1, show.legend = FALSE) +
	scale_fill_manual(values=c("#F8766D", "NA" ,"NA", "NA", "NA")) + 
	geom_line(aes(obstacles, `distance`, col=Level), linewidth = 1) +
	theme(legend.position = "bottom")

p2 <- ggplot(reg.lines.df) +
	geom_ribbon(aes(obstacles, ymin=lower, ymax=upper, group=Level, fill=Level),
		  alpha=rib.alpha, linewidth=.1, show.legend = FALSE) +
	scale_fill_manual(values=c("NA", "#B79F00" ,"NA", "NA", "NA")) + 
	geom_line(aes(obstacles, `distance`, col=Level), linewidth = 1) +
	theme(legend.position = "bottom")

p3 <- ggplot(reg.lines.df) +
	geom_ribbon(aes(obstacles, ymin=lower, ymax=upper, group=Level, fill=Level),
		  alpha=rib.alpha, linewidth=.1, show.legend = FALSE) +
	scale_fill_manual(values=c("NA", "NA", "#00BA38", "NA", "NA")) + 
	geom_line(aes(obstacles, `distance`, col=Level), linewidth = 1) +
	theme(legend.position = "bottom")

p4 <- ggplot(reg.lines.df) +
	geom_ribbon(aes(obstacles, ymin=lower, ymax=upper, group=Level, fill=Level),
		  alpha=rib.alpha, linewidth=.1, show.legend = FALSE) +
	scale_fill_manual(values=c("NA", "NA", "NA", "#619CFF", "NA")) + 
	geom_line(aes(obstacles, `distance`, col=Level), linewidth = 1) +
	theme(legend.position = "bottom")

p5 <- ggplot(reg.lines.df) +
	geom_ribbon(aes(obstacles, ymin=lower, ymax=upper, group=Level, fill=Level),
		  alpha=rib.alpha, linewidth=.1, show.legend = FALSE) +
	scale_fill_manual(values=c("NA", "NA", "NA", "NA", "#F564E3")) + 
	geom_line(aes(obstacles, `distance`, col=Level), linewidth = 1) +
	theme(legend.position = "bottom")


mod.f <- latex2exp::TeX("$y_{i}=\\beta_{0}+input\\_noise*n\\_obstacles^2$ (input noise in 5 levels)")
# pdf(file="plots/results_model.pdf", width=8, height=5)
p1 + p2 + p3 + p4 + p5 +
	plot_annotation(title="Bootstrapped regression curves with confidence Bands", subtitle=mod.f) + 
	plot_layout(guides = "collect", ncol = 3) &
	theme(legend.position = "bottom")
# dev.off()
	
# Export to png is good with 5x7 inches


# + Plot of Bootstrap Distribution ##############################################
dat.par <- as_tibble(t(boot.data$parameters))
names(dat.par) <- c("intercept", "inp_noise_0", "inp_noise_0.5",
		"inp_noise_1", "inp_noise_1.5", "inp_noise_2")


# input noise
dat.par2 <- pivot_longer(dat.par[-1], 1:5)
names(dat.par2) <- c("Level", "value")
dat.par2$Level <- as.factor(dat.par2$Level)

library(ggridges)
library(ggtext)
# pdf(file="plots/results_bootstrap.pdf", width=8, height=5)
ggplot(dat.par2) +
	stat_density_ridges(aes(x=value, y=Level, fill=Level), geom="density_ridges",
			quantile_lines=TRUE, quantiles=c(0.5), alpha=.9, colour="darkslategray")+
	labs(title="Bootstrap distributions by level of input noise") +
	theme(legend.position = "bottom", plot.caption = element_markdown(lineheight = 1.3))
# dev.off()


