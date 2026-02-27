# Using BARTs to model flowering activity
# run locally
# last used/modified jby, 2026.02.26

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/jotr_bonus_blooms")

library("tidyverse")
library("embarcadero")
library("cowplot")

set.seed(19820604)

#-----------------------------------------------------------
# initial file loading

# set parameters as variables
taxnum <- 1595251 # Joshua trees
# Prunus ilicifolia = 57250

flow <- read.csv(paste("output/flowering_freq_climate_", taxnum, ".csv", sep="")) |> mutate(flr = prop_flr > 0) |> filter(year>=2016) # flowering/not flowering, biologically-informed candidate predictors

glimpse(flow) # 5,323; from 2016 on, 5070

hist(flow$prop_flr) # check the cutoff I've set for binary flowering

table(flow$year, flow$flr)
table(flow$year, flow$quarter)

flow.norm <- filter(flow, quarter <= 2)
flow.anom <- filter(flow, quarter > 2)

#-------------------------------------------------------------------------
# predictor selection

if(!dir.exists("output/models")) dir.create("output/models")

# predictors
xnames <- paste0(rep(c("tmax", "tmin", "ppt"),5),"_","Q",rep(0:4, each=3)) # weather data, all of it


# ALL YEAR, non-RI ------------------------------
# VARIMP variable importance across the whole predictor set .........
flow.varimp <- varimp.diag(y.data=as.numeric(flow[,"flr"]), x.data=flow[,xnames])

write_rds(flow.varimp, file=paste("output/models/bart.varimp.all.", taxnum, ".rds", sep="")) # save varimp() results
# flow.varimp <- read_rds(file=paste("output/models/bart.varimp.all.", taxnum, ".rds", sep=""))

# generate a better-organized varimp() figure
vi.fac <- flow.varimp$data |> group_by(variable) |> summarise(maxi = max(imp)) |> arrange(-maxi)

var_sel_compare <- flow.varimp$data |> mutate(trees = factor(trees, c(200, 150, 100, 50, 20, 10)), variable=factor(variable, vi.fac$variable))

levels(var_sel_compare$variable)

# colors for this: '#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#08589e'

predsel <- ggplot(data=filter(var_sel_compare, trees%in%c(10,20,50,100,200)), aes(x=variable, y=imp, color=trees, group=trees)) +
	geom_line(linewidth=0.5) + geom_point(size=2) +
	labs(x = "Predictor", y = "Importance (prop. splits)") +
	scale_color_manual(values=c('#ccebc5','#7bccc4','#4eb3d3','#2b8cbe','#08589e'), name="N trees") +
	theme_bw(base_size=12) + theme(axis.text.x=element_text(angle=75, hjust=1), legend.position="inside", legend.position.inside=c(0.9,0.725), legend.key.spacing.y=unit(1, "mm"), legend.key.size = unit(4, "mm"))


{cairo_pdf(paste0("output/figures/varimp_all_", taxnum, ".pdf"), width=5, height=4)

predsel

}
dev.off()

levels(var_sel_compare$variable)[1:7] 
# okay so that's "tmax_Q0" "ppt_Q1"  "tmax_Q1" "ppt_Q4"  "ppt_Q0"  "tmin_Q1" "ppt_Q3" 

# ALL YEAR, RI ----------------------------------
# VARIMP variable importance across the whole predictor set .........
flow.varimp <- varimp.diag(y.data=as.numeric(flow[,"flr"]), x.data=flow[,xnames], ri.data=flow[,"quarter"])

write_rds(flow.varimp, file=paste("output/models/ribart.varimp.all.", taxnum, ".rds", sep="")) # save varimp() results
# flow.varimp <- read_rds(file=paste("output/models/ribart.varimp.all.", taxnum, ".rds", sep=""))

# generate a better-organized varimp() figure
vi.fac <- flow.varimp$data |> group_by(variable) |> summarise(maxi = max(imp)) |> arrange(-maxi)

var_sel_compare <- flow.varimp$data |> mutate(trees = factor(trees, c(200, 150, 100, 50, 20, 10)), variable=factor(variable, vi.fac$variable))

levels(var_sel_compare$variable)

# colors for this: '#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#08589e'

predsel <- ggplot(data=filter(var_sel_compare, trees%in%c(10,20,50,100,200)), aes(x=variable, y=imp, color=trees, group=trees)) +
	geom_line(linewidth=0.5) + geom_point(size=2) +
	labs(x = "Predictor", y = "Importance (prop. splits)") +
	scale_color_manual(values=c('#ccebc5','#7bccc4','#4eb3d3','#2b8cbe','#08589e'), name="N trees") +
	theme_bw(base_size=12) + theme(axis.text.x=element_text(angle=75, hjust=1), legend.position="inside", legend.position.inside=c(0.85,0.725), legend.key.spacing.y=unit(1, "mm"), legend.key.size = unit(4, "mm"))


{cairo_pdf(paste0("output/figures/varimp_all_RI_", taxnum, ".pdf"), width=5, height=4)

predsel

}
dev.off()

levels(var_sel_compare$variable)[1:6]
# and that's "tmax_Q0" "ppt_Q4"  "ppt_Q1"  "ppt_Q0"  "tmax_Q1" "tmin_Q1"


#-------------------------------------------------------------------------
# MODEL FITTING 

# fill this in based on results of above section
topX <- levels(var_sel_compare$variable)[1:6] # adjust as needed, but use the set for RI


# vanilla BART for partials .....................
flr.mod.all <- bart(y.train=as.numeric(flow[,"flr"]), x.train=flow[,topX], keeptrees=TRUE, seed = 19861126)

invisible(flr.mod.all$fit$state)
write_rds(flr.mod.all, file=paste("output/models/bart.model.all.", taxnum, ".rds", sep="")) # save model
# flr.mod.all <- read_rds(paste("output/models/bart.model.all.", taxnum, ".rds", sep=""))

summary(flr.mod.all) # AUC reflects classification accuracy, how's that look? 0.8656455

mod_valid <- summary(flr.mod.all)$data %>% dplyr::select(fitted, observed) %>% mutate(type="Training data", classified=fitted>0.1540352)

rmse(mod_valid$classified, mod_valid$observed) # RMSE = 0.5094763

p <- partial(flr.mod.all, topX, trace=FALSE, smooth=5) # visualize partials; can't with RI?
varimp(flr.mod.all)

# Non-RI model to illustrate partials...

write_rds(p, file=paste("output/models/bart.model.all.partials.", taxnum, ".rds", sep=""))
# p <- read_rds(paste("output/models/bart.model.all.partials.", taxnum, ".rds", sep=""))

# Random-intercept, probably working model ......
flr.mod.all.ri <- rbart_vi(
	as.formula(paste(paste('flr', paste(topX, collapse=' + '), sep = ' ~ '), 'quarter', sep=' - ')),
	data = flow,
	group.by = flow[,'quarter'],
	n.chains = 1,
	k = 2,
	power = 2,
	base = 0.95,
	keepTrees = TRUE, 
	seed = 19861126)

invisible(flr.mod.all.ri$fit[[1]]$state)
write_rds(flr.mod.all.ri, file=paste("output/models/RIbart.model.all.", taxnum, ".rds", sep="")) # save model
# flr.mod.all.ri <- read_rds(paste("output/models/bart.model.all.", taxnum, ".rds", sep=""))

summary(flr.mod.all.ri) # AUC reflects classification accuracy, how's that look? 0.8553109

mod_valid <- summary(flr.mod.all.ri)$data %>% dplyr::select(fitted, observed) %>% mutate(type="Training data", classified=fitted>0.2928244)

rmse(mod_valid$classified, mod_valid$observed) # RMSE = 0.4901395

riplot <- plot.ri(flr.mod.all.ri, temporal=TRUE) + 
	geom_ribbon(fill="#cab2d6") + 
	geom_line(color="white") +
	scale_x_continuous(breaks=1:12) +
	labs(x="Month", y="RI effect") +
	theme_bw(base_size=12) +
	theme(axis.text.x=element_text(angle=0))


{cairo_pdf(paste("output/figures/RIbart_all_RI_effects_", taxnum, ".pdf", sep=""), width=4, height=3)

riplot

}
dev.off()


#-------------------------------------------------------------------------
# Partials and spartials in example years

# read back in, if necessary
flr.mod.all <- read_rds(paste("output/models/bart.model.all.", taxnum, ".rds", sep=""))

topX <- attr(flr.mod.all$fit$data@x, "term.labels")


# PARTIALS ------------------------------------------------
# OH and what we want here is normal-season vs all-data models overlaid
p <- read_rds(paste("output/models/bart.model.all.partials.", taxnum, ".rds", sep=""))

# reorganize raw data underlying the partials (gotta inspect to sort this out)
names(p) <- topX

partvals <- NULL # empty object to hold data

for(part in topX){
	partvals <- rbind(partvals, data.frame(predictor=part, p[[part]]$data))
}


glimpse(partvals)

partvals$predictor <- factor(partvals$predictor, topX)
levels(partvals$predictor) 

partvals$predtype <- "Precip"
partvals$predtype[grepl("tmax", partvals$predictor)] <- "Tmax"
partvals$predtype[grepl("tmin", partvals$predictor)] <- "Tmin"

table(partvals$predtype)



# generate a figure .............................

colors <- c("#b2df8a", "#fdbf6f", "#a6cee3")

partplots <- ggplot(partvals) + 
	geom_ribbon(aes(x=x, ymin=q05, ymax=q95, fill=predtype), alpha=1) + 
	geom_line(aes(x=x, y=med), color="white", alpha=1) + 
	facet_wrap("predictor", nrow=2, scale="free") + 
	labs(y="Marginal Pr(Flowers)", x="Predictor value") + 
	scale_fill_manual(values=colors, name="Predictor type") +
	theme_bw() + 
	theme(panel.spacing=unit(0.1,"in"), 
		legend.position="none", 
		legend.text=element_text(size=10),
		legend.title=element_text(size=11),
		legend.key.size=unit(3, "mm")
		)

partplots

{cairo_pdf(paste("output/figures/predictor_partials_compare_", taxnum, ".pdf", sep=""), width=6.5, height=4)

partplots

}
dev.off()

# predictor selection and partials together

{cairo_pdf(paste("output/figures/BART_predsel_partials_", taxnum, ".pdf", sep=""), width=9, height=4)

ggdraw() + draw_plot(predsel, 0, 0, 0.39, 1) + draw_plot(partplots, 0.39, 0, 0.6, 1) + draw_plot_label(label=c("A", "B"), x=c(0, 0.39), y=1)

}
dev.off()



