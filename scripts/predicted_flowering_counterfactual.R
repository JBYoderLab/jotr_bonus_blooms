# Analyzing predicted historical flowering in Joshua tree
# visualizing results of the counterfactual
# Assumes local environment
# jby 2026.05.12

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

flow <- read.csv(paste("output/flowering_freq_climate_", taxnum, ".csv", sep="")) %>% mutate(flr = prop_flr > 0) # flowering/not flowering, biologically-informed candidate predictors

glimpse(flow) # 5,323 since 2008

hist(flow$prop_flr) # check the cutoff I've set for binary flowering

table(flow$year, flow$flr)
table(flow$year, flow$month)



# useful bits and bobs
MojExt <- extent(-120, -112, 33, 39) # Mojave extent, maybe useful

# read the model back in
flr.mod.all.ri <- read_rds(paste("output/models/RIbart.model.all.", taxnum, ".rds", sep=""))
allX <- attr(flr.mod.all.ri$fit[[1]]$data@x, "term.labels")

summary(flr.mod.all.ri) #

true.vector <- flr.mod.all.ri$fit[[1]]$data@y 
  
pred <- prediction(colMeans(pnorm(flr.mod.all.ri$yhat.train)), true.vector)
  
perf.tss <- performance(pred,"sens","spec")
tss.list <- (perf.tss@x.values[[1]] + perf.tss@y.values[[1]] - 1)
tss.df <- data.frame(alpha=perf.tss@alpha.values[[1]],tss=tss.list)
  
thresh <- min(tss.df$alpha[which(tss.df$tss==max(tss.df$tss))])
thresh





# create some date pegs for the quarters
summout$month <- ifelse(summout$quarter==1, 2, ifelse(summout$quarter==2, 5, ifelse(summout$quarter==3, 8, 11)))
glimpse(summout)

shading <- data.frame(xmins=paste(2016:2025, "1 1"),xmaxes=paste(2016:2025, "6 30"))

{cairo_pdf("output/figures/quarterly_flowering_realvcf_vs_data.pdf", width=6.5, height=3)

ggplot() + 
	geom_rect(data=shading, aes(xmin=ymd(xmins), xmax=ymd(xmaxes), ymin=-Inf, ymax=Inf), fill="gray90") +
	geom_bar(data=filter(flr.raw.sums, year>=2016), aes(x=ymd(paste(year, month, 14)), y=PropFlr), fill="#66c2a4", stat="identity") +
	geom_bar(data=filter(flr.raw.sums, year %in% c(2018,2025), quarter==4), aes(x=ymd(paste(year, month, 14)), y=PropFlr), fill=NA, color="black", stat="identity", width=90) +
	geom_line(data=summout, aes(x=ymd(paste(year, month, 14)), y=FreqFlr, group=gcm, color=pred_type, alpha=pred_type), linewidth=0.75) +

	scale_color_manual(values=c("blue", "orange"), name="Predicted from", labels=c("Actual", "Counterfactuals")) + 
	scale_alpha_manual(values=c(0.75, 0.5), guide="none") + 

	labs(x="Date", y="Prop. range with flowering predicted\n or cells with flowering observed", title="Modeled quarterly flowering") +
	theme_bw() +
	theme(legend.position="inside", legend.position.inside=c(0.74,1.07), legend.text=element_text(lineheight=0.8), legend.key.height=unit(3, "mm"), legend.direction="horizontal")

}
dev.off()


#-------------------------------------------------------------------------
# how do the (posterior distribution of) 10 GCMs differ from actuals?

postActFiles <- list.files("output/models/RIBART_all_year_predictions/", pattern="posterior", full.names=TRUE)
postCfactFiles <- list.files("output/models/RIBART_counterfactual_predictions/", pattern="posterior", full.names=TRUE)

tsdm <- st_transform(sdm, crs=crs(rast(postActFiles[1])))

# set up data organization
propFlrOut <- data.frame(matrix(0,0,10))
names(propFlrOut) <-  c("pred_type", "mod_thresh", "year", "quarter", "gcm", "repl", "lon", "lat", "pred_PrFlr", "pred_Flr")

# LOOP over the individual stack files
for(quart in c(postActFiles, postCfactFiles)){

# quart=postActFiles[1]
# quart=postCfactFiles[1]

predFlr <- mask(rast(quart), tsdm, touches=TRUE)

# LOOP over replicates
for(repl in 1:300){

# extract and package info for the layer ...
repdf <- data.frame(
			pred_type = if(grepl("counterfactual", quart)) "counterfactual" else "actual", 
			mod_thresh = thresh, 
			year = as.numeric(gsub(".+_(\\d+)_Q\\d\\.tiff", "\\1", quart)),
			quarter = as.numeric(gsub(".+_\\d+_Q(\\d)\\.tiff", "\\1", quart)),
			gcm = if(grepl("counterfactual", quart)) as.numeric(gsub(".+_GCM(\\d+)_\\d+_Q.+\\.tiff", "\\1", quart)) else NA,
			repl = repl,
			nPredFlr = length(which(values(predFlr[[repl]]) > thresh)) ,
			propFlr = length(which(values(predFlr[[repl]]) > thresh))/length(which(!is.na(values(predFlr[[repl]]))))
			)

propFlrOut <- rbind(propFlrOut, repdf)

}

write.table(propFlrOut, "output/posterior_predictions_factuals_counterfactuals.csv", col.names=TRUE, row.names=FALSE, sep=",")

}

# propFlrOut <- read.csv("output/posterior_predictions_factuals_counterfactuals.csv")

#-------------------------------------------------------------------------
# quarterly precip and temp

tmax <- rast(pd_to_file(prism_archive_subset("tmax", "monthly", resolution="4km", years=2016:2025, mon=1:12)))
names(tmax) <- gsub(".+_(\\d+).*", "\\1", prism_archive_subset("tmax", "monthly", resolution="4km", years=2016:2025, mon=1:12))

tmax <- mask(tmax[[sort(names(tmax))]], tsdm, touches=TRUE)
tmaxQ <- tapp(tmax, rep(1:40, each=3), max)

tmin <- rast(pd_to_file(prism_archive_subset("tmin", "monthly", resolution="4km", years=2016:2025, mon=1:12)))
names(tmin) <- gsub(".+_(\\d+).*", "\\1", prism_archive_subset("tmin", "monthly", resolution="4km", years=2016:2025, mon=1:12))

tmin <- mask(tmin[[sort(names(tmin))]], tsdm, touches=TRUE)
tminQ <- tapp(tmin, rep(1:40, each=3), min)

ppt <- rast(pd_to_file(prism_archive_subset("ppt", "monthly", resolution="4km", years=2016:2025, mon=1:12)))
names(ppt) <- gsub(".+_(\\d+).*", "\\1", prism_archive_subset("ppt", "monthly", resolution="4km", years=2016:2025, mon=1:12))

ppt <- mask(ppt[[sort(names(ppt))]], tsdm, touches=TRUE)
pptQ <- tapp(ppt, rep(1:40, each=3), sum)


predsums <- data.frame(year = rep(2016:2025, each=12), month=1:12,
					pred = rep(c("tmax", "tmin", "ppt"), each=120),
					mn = sapply(c(tmax, tmin, ppt), function(x) mean(values(x), na.rm=TRUE)),
					md = sapply(c(tmax, tmin, ppt), function(x) median(values(x), na.rm=TRUE)),
					lo95 = sapply(c(tmax, tmin, ppt), function(x) quantile(values(x), 0.025, na.rm=TRUE)),
					up95 = sapply(c(tmax, tmin, ppt), function(x) quantile(values(x), 0.975, na.rm=TRUE))
					)

glimpse(predsums)

write.table(predsums, "output/JT_range_monthly_weather_summaries_2016_2025.csv", sep=",", col.names=TRUE, row.names=FALSE)
# predsums <- read.csv("output/JT_range_monthly_weather_summaries_2016_2025.csv")

predsums$month <- ifelse(sumModDiffs$quarter==1, 2, ifelse(sumModDiffs$quarter==2, 5, ifelse(sumModDiffs$quarter==3, 8, 11)))
glimpse(predsums)


# useful figure stuff ...
shading <- data.frame(xmins=paste(2016:2025, "1 1"),xmaxes=paste(2016:2025, "6 30"))
bonusshade <- data.frame(xmins=paste(c(2018,2025), "10 1"),xmaxes=paste(c(2018,2025), "12 31"))


MonTemps <- ggplot() + 
	geom_rect(data=shading, aes(xmin=ymd(xmins), xmax=ymd(xmaxes), ymin=-Inf, ymax=Inf), fill="gray90") +
	geom_rect(data=bonusshade, aes(xmin=ymd(xmins), xmax=ymd(xmaxes), ymin=-Inf, ymax=Inf), fill="slategray2") +

	geom_line(data=filter(predsums, pred=="tmax"), aes(x=ymd(paste(year, month, 14)), y=mn), color="#e31a1c") +
	geom_point(data=filter(predsums, pred=="tmax"), aes(x=ymd(paste(year, month, 14)), y=mn), shape=17, color="#e31a1c") +
	
	geom_line(data=filter(predsums, pred=="tmin"), aes(x=ymd(paste(year, month, 14)), y=mn), color="#ff7f00") +	
	geom_point(data=filter(predsums, pred=="tmin"), aes(x=ymd(paste(year, month, 14)), y=mn), shape=25, color="#ff7f00", fill="white") +

	labs(x="Date", title="A. Range-wide average monthly maximum and minimum temperature (°C)") +
	theme_bw() +
	theme(legend.position="none", axis.title=element_blank(), axis.text.x=element_blank())

MonTemps

MonPPTs <- ggplot() + 
	geom_rect(data=shading, aes(xmin=ymd(xmins), xmax=ymd(xmaxes), ymin=-Inf, ymax=Inf), fill="gray90") +
	geom_rect(data=bonusshade, aes(xmin=ymd(xmins), xmax=ymd(xmaxes), ymin=-Inf, ymax=Inf), fill="slategray2") +

	geom_line(data=filter(predsums, pred=="ppt"), aes(x=ymd(paste(year, month, 14)), y=mn), color="#1f78b4") +	
	geom_point(data=filter(predsums, pred=="ppt"), aes(x=ymd(paste(year, month, 14)), y=mn), shape=15, color="#1f78b4") +

	labs(x="Date", title="B. Range-wide average monthly total precipitation (mm)") +
	theme_bw() +
	theme(legend.position="none", axis.title=element_blank(), axis.text.x=element_blank())

MonPPTs


#-------------------------------------------------------------------------
# actual and modeled flowering

# summarize for pointranges
summods <- propFlrOut |> group_by(year, quarter, pred_type) |> summarize(mnPropFlr=mean(propFlr), mdPropFlr=median(propFlr), lo95=quantile(propFlr, 0.025), up95=quantile(propFlr, 0.975)) 

summods$month <- ifelse(summods$quarter==1, 2, ifelse(summods$quarter==2, 5, ifelse(summods$quarter==3, 8, 11)))

glimpse(summods)

# original data!
flr.raw.sums <- flow |> group_by(year, quarter) |> summarize(Ntot=length(flr), Nflr=length(which(flr)), PropFlr=Nflr/Ntot) |> mutate(month=if_else(quarter==1, 2, if_else(quarter==2, 5, if_else(quarter==3, 8, 11))))

glimpse(flr.raw.sums)

# useful figure stuff ...
shading <- data.frame(xmins=paste(2016:2025, "1 1"),xmaxes=paste(2016:2025, "6 30"))
bonusshade <- data.frame(xmins=paste(c(2018,2025), "10 1"),xmaxes=paste(c(2018,2025), "12 31"))

propFlrData <- ggplot() + 
	geom_rect(data=shading, aes(xmin=ymd(xmins), xmax=ymd(xmaxes), ymin=-Inf, ymax=Inf), fill="gray90") +
	geom_rect(data=bonusshade, aes(xmin=ymd(xmins), xmax=ymd(xmaxes), ymin=-Inf, ymax=Inf), fill="slategray2") +

	geom_bar(data=filter(flr.raw.sums, year>=2016), aes(x=ymd(paste(year, month, 14)), y=PropFlr), fill="#33a02c", stat="identity") +
#	geom_bar(data=filter(flr.raw.sums, year %in% c(2018,2025), quarter==4), aes(x=ymd(paste(year, month, 14)), y=PropFlr), fill=NA, color="black", stat="identity", width=90) +
	annotate("segment", x=ymd("2018-11-14"), y=0.57, xend=ymd("2018-11-14"), yend=0.37, arrow=arrow(length=unit(3, "mm"), angle=25), color="#ff7f00", linewidth=1) + 
	annotate("segment", x=ymd("2025-11-14"), y=0.43, xend=ymd("2025-11-14"), yend=0.23, arrow=arrow(length=unit(3, "mm"), angle=25), color="#ff7f00", linewidth=1) + 


	labs(x="Date", title="C. Prop. iNaturalist records with flowering observed, by quarter") +
	theme_bw() +
	theme(legend.position="none", axis.title=element_blank(), axis.text.x=element_blank())
	
propFlrData


propFlrPred <- ggplot() + 
	geom_rect(data=shading, aes(xmin=ymd(xmins), xmax=ymd(xmaxes), ymin=-Inf, ymax=Inf), fill="gray90") +
	geom_rect(data=bonusshade, aes(xmin=ymd(xmins), xmax=ymd(xmaxes), ymin=-Inf, ymax=Inf), fill="slategray2") +

	geom_linerange(data=filter(summods, pred_type=="actual"), aes(x=ymd(paste(year, month, 1)), ymin=lo95, ymax=up95), linewidth=0.75, color="#6a3d9a") +
	geom_point(data=filter(summods, pred_type=="actual"), aes(x=ymd(paste(year, month, 1)), y=mdPropFlr), color="#6a3d9a") +
	
	geom_linerange(data=filter(summods, pred_type=="counterfactual"), aes(x=ymd(paste(year, month, 28)), ymin=lo95, ymax=up95), linewidth=0.75, color="#ff7f00") +
	geom_point(data=filter(summods, pred_type=="counterfactual"), aes(x=ymd(paste(year, month, 28)), y=mdPropFlr), color="#ff7f00") +	
	
	labs(x="Date", title="D. Modeled quarterly actual and counterfactual prop. flowering (mean ± 95% CI)") +
	theme_bw() +
	theme(legend.position="none", axis.title=element_blank(), axis.text.x=element_blank())

propFlrPred



#-------------------------------------------------------------------------
# differences between 

actuals <- filter(propFlrOut, is.na(gcm)) |> dplyr::select(year, quarter, repl, propFlr) |> rename(ActPropFlr=propFlr)
glimpse(actuals)

moddiffs <- propFlrOut |> filter(!is.na(gcm)) |> left_join(actuals) |> mutate(cfdiff = ActPropFlr-propFlr)
glimpse(moddiffs)

sumModDiffs <- moddiffs |> group_by(year, quarter) |> summarize(mnCFdiff=mean(cfdiff), mdCFdiff=median(cfdiff), lo95=quantile(cfdiff, 0.025), up95=quantile(cfdiff, 0.975), sig=if_else((lo95<0 & up95<0) | (lo95>0 & up95>0), TRUE, FALSE) ) 
glimpse(sumModDiffs)

sumModDiffs$month <- ifelse(sumModDiffs$quarter==1, 2, ifelse(sumModDiffs$quarter==2, 5, ifelse(sumModDiffs$quarter==3, 8, 11)))
glimpse(sumModDiffs)

CFdiffs <- ggplot() + 
	geom_rect(data=shading, aes(xmin=ymd(xmins), xmax=ymd(xmaxes), ymin=-Inf, ymax=Inf), fill="gray90") +
	geom_rect(data=bonusshade, aes(xmin=ymd(xmins), xmax=ymd(xmaxes), ymin=-Inf, ymax=Inf), fill="slategray2") +

	geom_hline(yintercept=0, linetype=2, color="gray40") +

	geom_linerange(data=sumModDiffs, aes(x=ymd(paste(year, month, 14)), ymin=lo95, ymax=up95), color="#33a02c", linewidth=0.75) +

	geom_point(data=sumModDiffs, aes(x=ymd(paste(year, month, 14)), y=mdCFdiff, fill=sig), color="#33a02c", shape=21) +
	scale_fill_manual(values=c("white","#33a02c"), guide="none") +
	
	ylim(-1,1) +

	labs(x="Date", title="E. Factual-counterfactual difference in prop. flowering (mean ± 95% CI)") +
	theme_bw() +
	theme(legend.position="none", axis.title.y=element_blank())

CFdiffs


#-------------------------------------------------------------------------
# put it all together

{cairo_pdf("output/figures/Fig02_data_predictors_model.pdf", width=7.5, height=9)

MonTemps / MonPPTs / propFlrData / propFlrPred / CFdiffs

}
dev.off()


