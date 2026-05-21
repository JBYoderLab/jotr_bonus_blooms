# Analyzing predicted historical flowering in Joshua tree
# next-level model consistency check, "LOO"
# Assumes local environment
# jby 2026.04.13

# starting up ------------------------------------------------------------

# rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/jotr_bonus_blooms")

library("tidyverse")
library("embarcadero")
library("cowplot")

library("sf")
library("prism")
library("terra")
library("geosphere")
library("barrks")

set.seed(18090212)

prism_set_dl_dir("../data/PRISM") # expect an error, oy

#-------------------------------------------------------------------------
# load up and organize historical data

# Jotr SDM range polygon
sdm <- read_sf("../data/Yucca/Jotr_SDM2023_range/Jotr_SDM2023_range_simple.shp")

# set parameters as variables
taxnum <- 1595251 # Joshua trees
# Prunus ilicifolia = 57250

flow <- read.csv(paste("output/flowering_freq_climate_", taxnum, ".csv", sep="")) |> mutate(flr = prop_flr > 0) |> filter(year>=2016) # flowering/not flowering, biologically-informed candidate predictors

glimpse(flow)

hist(flow$prop_flr) # check the cutoff I've set for binary flowering

table(flow$year, flow$flr)
table(flow$year, flow$quarter)


# raster files of predicted prFL
pred.files <- list.files("output/models/RIBART_all_year_predictions", pattern=".tiff", full=TRUE)
cf.files <- list.files("output/models/RIBART_counterfactual_predictions", pattern=".tiff", full=TRUE)


# useful bits and bobs
MojExt <- extent(-120, -112, 33, 39) # Mojave extent, maybe useful

# read the main model back in 
flr.mod.all.ri <- read_rds(paste("output/models/ribart.model.all.", taxnum, ".rds", sep=""))
topX <- attr(flr.mod.all.ri$fit[[1]]$data@x, "term.labels")

true.vector <- flr.mod.all.ri$fit[[1]]$data@y 
  
pred <- prediction(colMeans(pnorm(flr.mod.all.ri$yhat.train)), true.vector)
  
perf.tss <- performance(pred,"sens","spec")
tss.list <- (perf.tss@x.values[[1]] + perf.tss@y.values[[1]] - 1)
tss.df <- data.frame(alpha=perf.tss@alpha.values[[1]],tss=tss.list)
  
thresh <- min(tss.df$alpha[which(tss.df$tss==max(tss.df$tss))])
thresh


#-------------------------------------------------------------------------
# how well does the "all year" model predict anomalous events?

all.stack <- rast(pred.files)
names(all.stack) # confirm useful labeling

glimpse(filter(flow, quarter>2)) # here's our out of season records

flow.anom.preds <- flow |> filter(quarter>2) |> mutate(prFlr = NA, yrqu = paste0(year, "_Q", quarter))
glimpse(flow.anom.preds)

for(yq in unique(flow.anom.preds$yrqu)){

yqsub <- filter(flow.anom.preds, yrqu==yq)
flow.anom.preds$prFlr[flow.anom.preds$yrqu==yq] <- terra::extract(all.stack[[yq]], yqsub[,c("lon", "lat")])[,2]

}

flow.anom.preds$pred_flr <- flow.anom.preds$prFlr >= thresh

glimpse(flow.anom.preds) # check

write.table(flow.anom.preds, "output/RIBART_prediction_all_year_vs_anomalous_obs.csv", sep=",", col.names=TRUE, row.names=FALSE) # save


# okay and how well does observed flowering (flr) line up with predicted (pred_flr)?

table(flow.anom.preds$flr)
table(flow.anom.preds$pred_flr) # okay

t.test(prFlr~flr, data=flow.anom.preds) # THERE we go

table(flow.anom.preds$pred_flr, flow.anom.preds$flr) # okay?

performance(prediction(flow.anom.preds$prFlr, flow.anom.preds$flr), "auc")@y.values[[1]] # 0.94 # DANG


#-------------------------------------------------------------------------
# okay we have a model that predicts anomalous flowering!
# time to see what it says about early 20th c vs recent

pred.stack <- mask(rast(pred.files), st_transform(sdm, crs=crs(rast(pred.files[1]))), touches=TRUE)
pred.stack


# take a look at this ...
par(mfrow=c(2,2))

for(q in 1:4) plot(pred.stack[[q]]>thresh, main=paste0("1901-Q",q)) 

# So ... let's try some hashtag data science
pred.df <- as.data.frame(pred.stack) |> cbind(crds(pred.stack)) |> rename(lon=x, lat=y) |> pivot_longer(1:240, names_to="YQ", values_to="pred_PrFlr") |> mutate(year=as.numeric(gsub("(\\d+)_Q\\d", "\\1", YQ)), quarter=as.numeric(gsub("\\d+_Q(\\d)", "\\1", YQ)), season=if_else(quarter>=3, "Autumn", "Spring"), pred_flr = pred_PrFlr>thresh)

glimpse(pred.df)

# write out
write.table(pred.df, "output/RIBART_all_year_predicted_flowering_jtrange.csv", sep=",", col.names=TRUE, row.names=FALSE)
# pred.df <- read.csv("output/BART_all_year_predicted_flowering_jtrange.csv")

pred.flr.YQsumm <- pred.df |> group_by(year, quarter) |> summarize(tot=length(pred_flr), tot_flr=length(which(pred_flr))) |> mutate(freq_flr=tot_flr/tot) |> mutate(plotmonth=if_else(quarter==1, 2, if_else(quarter==2, 5, if_else(quarter==3, 8, 11))))

pred.flr.YQsumm # maybe useful, not quite what I want yet!

pred.flr.ESsumm <- pred.df |> group_by(year, quarter) |> 
			summarize(tot=length(pred_flr), tot_flr=length(which(pred_flr))) |> 
			mutate(freq_flr=tot_flr/tot) |>
			mutate(season=factor(if_else(quarter>=3, "Autumn", "Spring"), c("Spring", "Autumn")), era=if_else(year>1995, "1996-2025", "1901-1930"))

library("MASS")
library("lme4")
anova(glm(freq_flr~season*era, data=pred.flr.ESsumm, family=quasibinomial(link = "logit")))

ggplot(pred.flr.ESsumm, aes(y=freq_flr, x=era, fill=season)) + geom_boxplot()

{cairo_pdf("output/figures/pred_quarterly_flr_freq_by_era.pdf", width=4.5, height=4)

ggplot(pred.flr.ESsumm, aes(x=freq_flr, fill=era)) + geom_histogram(position="dodge") +
	facet_grid(season~era, scale="free_y") +
	scale_fill_manual(values=c("#1b7837", "#af8dc3"), guide="none") +
	labs(x="Proportion of range with predicted flowering", y="Quarters") +
	theme_bw(base_size=12) + 
	theme(panel.spacing=unit(5, "mm"))

}
dev.off()


wilcox.test(freq_flr~era, filter(pred.flr.ESsumm, season=="Autumn"), alternative="g") 
# early period greater, p = 0.0004759
pred.flr.ESsumm |> group_by(era, season) |> summarize(mdFlrFrq = median(freq_flr), mnFlrFrq = mean(freq_flr), NgtMin=length(which(freq_flr>0.05)))


{cairo_pdf("output/figures/pred_quarterly_flr_freq_by_year.pdf", width=6, height=4)

ggplot(pred.flr.summ, aes(x=ymd(paste(year,plotmonth, "14")), y=freq_flr, group=year, color=year<1995)) +
	geom_line(alpha=0.5, linewidth=0.75) +  
	scale_color_manual(values=c("#1b7837", "#af8dc3"), labels=c("1996-2025", "1901-1930"), name="Timeframe") +
	scale_x_continuous(breaks=c(3,6,9,12)) +
	labs(x = "Month", y = "Prop range with predicted flowering") + 
	theme_bw() + 
	theme(legend.position="inside", legend.position.inside=c(0.85, 0.75))

}
dev.off()

pred.flr.summ.summ <- pred.flr.summ |> mutate(norm.season=month<7, period=ifelse(year<1995, "early", "recent")) |> group_by(period, norm.season) |> summarize(mnPropFlr = mean(freq_flr, na.rm=TRUE), mdPropFlr = median(freq_flr, na.rm=TRUE), yrs_w_flr=length(which(freq_flr>0.001))/6)

pred.flr.summ.summ





