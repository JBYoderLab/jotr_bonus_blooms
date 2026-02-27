# Analyzing predicted historical flowering in Joshua tree
# next-level model consistency check, "LOO"
# Assumes local environment
# jby 2026.02.24

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

flow <- read.csv(paste("output/flowering_freq_climate_", taxnum, ".csv", sep="")) %>% mutate(flr = prop_flr > 0) |> filter(year>=2016) # flowering/not flowering, biologically-informed candidate predictors

glimpse(flow) # 5,070 since 2016

hist(flow$prop_flr) # check the cutoff I've set for binary flowering

table(flow$year, flow$flr)
table(flow$year, flow$month)

flow.norm <- filter(flow, month <=6)
flow.anom <- filter(flow, month > 6)


# useful bits and bobs
MojExt <- extent(-120, -112, 33, 39) # Mojave extent, maybe useful

# read the model back in
flr.mod.all.ri <- read_rds(paste("output/models/RIbart.model.all.", taxnum, ".rds", sep=""))
allX <- attr(flr.mod.all.ri$fit[[1]]$data@x, "term.labels")

summary(flr.mod.all.ri) #

#-------------------------------------------------------------------------
# okay we have a model that predicts anomalous flowering!
# can a similar model do that if it doesn't know about one anomaly?

table(flow$year, flow$quarter)

# sequentially 
# (1) pull one year's months 7:12
# (2) train a new model
# (3) predict to observations from the left-out months 7:12

# First, train models and make predictions
for(yr in c(2018,2025)){

# yr = 2018

train <- filter(flow, !(year==yr & quarter>=3)) # train on all but fall of the LOO year
test <- filter(flow, year==yr, quarter>=3) # test on 

# testmod <- bart(y.train=as.numeric(train[,"flr"]), x.train=train[,allX], keeptrees=TRUE, seed=19820604)

testmod <- rbart_vi(
	as.formula(paste(paste('flr', paste(allX, collapse=' + '), sep = ' ~ '), 'quarter', sep=' - ')),
	data = train,
	group.by = train[,'quarter'],
	n.chains = 1,
	k = 2,
	power = 2,
	base = 0.95,
	keepTrees = TRUE, 
	seed = 19820604)

summary(testmod)

invisible(testmod$fit[[1]]$state)
write_rds(testmod, file=paste("output/models/RIbart.model.LOO", yr, ".", taxnum, ".rds", sep="")) # save model

for(predyr in 2016:2025){

	# predyr <- 2020; qu <- 2 # test condition
	
	yearpred <- NULL
	
	for(qu in 1:4){

	if(qu == 1){
		tmaxfiles <- c(prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr, mon=1:3),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr-1, mon=10:12),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr-1, mon=7:9),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr-1, mon=4:6),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr-1, mon=1:3) )
		
		tminfiles <- c(prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr, mon=1:3),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr-1, mon=10:12),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr-1, mon=7:9),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr-1, mon=4:6),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr-1, mon=1:3) )
			
		pptfiles <- c(prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr, mon=1:3),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr-1, mon=10:12),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr-1, mon=7:9),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr-1, mon=4:6),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr-1, mon=1:3) )
	}
	if(qu == 2){
		tmaxfiles <- c(prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr, mon=4:6),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr, mon=1:3),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr-1, mon=12:10),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr-1, mon=7:9),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr-1, mon=4:6) )
		
		tminfiles <- c(prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr, mon=4:6),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr, mon=1:3),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr-1, mon=10:12),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr-1, mon=7:9),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr-1, mon=4:6) )
			
		pptfiles <- c(prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr, mon=4:6),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr, mon=1:3),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr-1, mon=10:12),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr-1, mon=7:9),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr-1, mon=4:6) )
	}
	if(qu == 3){
		tmaxfiles <- c(prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr, mon=7:9),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr, mon=4:6),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr, mon=1:3),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr-1, mon=10:12),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr-1, mon=7:9) )
		
		tminfiles <- c(prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr, mon=7:9),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr, mon=4:6),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr, mon=1:3),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr-1, mon=10:12),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr-1, mon=7:9) )
			
		pptfiles <- c(prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr, mon=7:9),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr, mon=4:6),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr, mon=1:3),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr-1, mon=10:12),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr-1, mon=7:9) )
	}
	if(qu == 4){
		tmaxfiles <- c(prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr, mon=10:12),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr, mon=7:9),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr, mon=4:6),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr, mon=1:3),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=predyr-1, mon=10:12) )
		
		tminfiles <- c(prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr, mon=10:12),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr, mon=7:9),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr, mon=4:6),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr, mon=1:3),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=predyr-1, mon=10:12) )
			
		pptfiles <- c(prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr, mon=10:12),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr, mon=7:9),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr, mon=4:6),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr, mon=1:3),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=predyr-1, mon=10:12) )
	}	
	
	envs <- c(app(rast(pd_to_file(tmaxfiles[1:3])), "mean", na.rm=TRUE),
				app(rast(pd_to_file(tminfiles[1:3])), "mean", na.rm=TRUE),
				app(rast(pd_to_file(pptfiles[1:3])), "sum", na.rm=TRUE),
				app(rast(pd_to_file(tmaxfiles[4:6])), "mean", na.rm=TRUE),
				app(rast(pd_to_file(tminfiles[4:6])), "mean", na.rm=TRUE),
				app(rast(pd_to_file(pptfiles[4:6])), "sum", na.rm=TRUE),
				app(rast(pd_to_file(tmaxfiles[7:9])), "mean", na.rm=TRUE),
				app(rast(pd_to_file(tminfiles[7:9])), "mean", na.rm=TRUE),
				app(rast(pd_to_file(pptfiles[7:9])), "sum", na.rm=TRUE),
				app(rast(pd_to_file(tmaxfiles[10:12])), "mean", na.rm=TRUE),
				app(rast(pd_to_file(tminfiles[10:12])), "mean", na.rm=TRUE),
				app(rast(pd_to_file(pptfiles[10:12])), "sum", na.rm=TRUE),
				app(rast(pd_to_file(tmaxfiles[13:15])), "mean", na.rm=TRUE),
				app(rast(pd_to_file(tminfiles[13:15])), "mean", na.rm=TRUE),
				app(rast(pd_to_file(pptfiles[13:15])), "sum", na.rm=TRUE) )
	names(envs) <- paste0(rep(c("tmax", "tmin", "ppt"),5),"_","Q",rep(0:4, each=3))
		
	preds <- crop(envs, extent(-120, -112, 33, 39))

	# prediction with the RI predictor (year) removed
	pred.flr <- predict(testmod, stack(preds), splitby=20, ri.name="quarter", ri.data=NA, ri.pred=FALSE)
	
	# if(mo < 10) molab <- paste0(0,mo) else molab <- mo
	crs(pred.flr) <- crs(envs)
	
	yearpred <- c(yearpred, pred.flr)
	
	} # END loop over prediction quarters
	
	outyear <- rast(stack(yearpred))
	names(outyear) <- paste0(predyr,"_Q",1:4)
	
	writeRaster(outyear, paste0("output/models/LOO_test/RIBART_predicted_flowering_LOO",yr,"_",predyr,".tiff"), overwrite=TRUE)
	
} # END loop over prediction years

} # END loop over LOO years


#-------------------------------------------------------------------------
# okay, THEN ask what each LOO model predicts for out-of-season flowering in every year ...

# actually not possible to fit the whole pile of data in memory as such ... let's loop

LOOsumm <- data.frame(LOOyear=NULL, year=NULL, quarter=NULL, Ncells=NULL, NcellsFlr=NULL, FreqFlr=NULL)

# LOOP over LOO years
for(looyr in c(2018,2025)){

# looyr <- 2020 # trial run

# load files
loofiles <- list.files("output/models/LOO_test", pattern=paste0("LOO",looyr,".*.tiff"), full.names=TRUE)
loostack <- mask(rast(loofiles), st_transform(sdm, crs=crs(rast(loofiles[1]))), touches=TRUE)

# everything I want to do here requires getting the model-specific cutoff, oy
loomod <- read_rds(paste("output/models/RIbart.model.LOO", looyr, ".", taxnum, ".rds", sep="")) # saved model

true.vector <- loomod$fit[[1]]$data@y 
  
pred <- prediction(colMeans(pnorm(loomod$yhat.train)), true.vector)
  
perf.tss <- performance(pred,"sens","spec")
tss.list <- (perf.tss@x.values[[1]] + perf.tss@y.values[[1]] - 1)
tss.df <- data.frame(alpha=perf.tss@alpha.values[[1]],tss=tss.list)
  
thresh <- min(tss.df$alpha[which(tss.df$tss==max(tss.df$tss))])

# now ...
loolong <- as.data.frame(loostack) |> cbind(crds(loostack)) |> pivot_longer(1:40, names_to="yrqu", values_to="pred_PrFlr") |> mutate(LOO_mod_thresh=thresh, pred_Flr=pred_PrFlr>thresh, LOOyear=looyr, year=as.numeric(gsub("(\\d+)_Q\\d", "\\1", yrqu)), quarter=as.numeric(gsub("\\d+_Q(\\d)", "\\1", yrqu))) |> rename(lon=x, lat=y) |> dplyr::select(LOOyear, LOO_mod_thresh, year, quarter, lon, lat, pred_PrFlr, pred_Flr)

# stash that
write.table(loolong, paste0("output/models/RIbart_LOO", looyr, "_predictions.csv"), sep=",", col.names=TRUE, row.names=FALSE)

# summarize that
loosumout <- loolong |> group_by(LOOyear, year, quarter) |> summarize(Ncells=length(pred_Flr), NcellsFlr=length(which(pred_Flr)), FreqFlr=NcellsFlr/Ncells)

LOOsumm <- rbind(LOOsumm, loosumout)

# stash that
write.table(LOOsumm, paste0("output/models/RIbart_LOO_predictions_summary.csv"), sep=",", col.names=TRUE, row.names=FALSE)

} # END loop over years


#-------------------------------------------------------------------------
# and now, a figure

inat_pheno_data <- read.csv(paste0("data/inat_phenology_data_", taxnum, ".csv"), h=TRUE)

flr.raw.sums <- flow |> group_by(year, quarter) |> summarize(Ntot=length(flr), Nflr=length(which(flr)), PropFlr=Nflr/Ntot) |> mutate(month=if_else(quarter==1, 2, if_else(quarter==2, 5, if_else(quarter==3, 8, 11))))

glimpse(flr.raw.sums)

# create some date pegs for the quarters
LOOsumm$month <- ifelse(LOOsumm$quarter==1, 2, ifelse(LOOsumm$quarter==2, 5, ifelse(LOOsumm$quarter==3, 8, 11)))
glimpse(LOOsumm)

shading <- data.frame(xmins=paste(2016:2025, "1 1"),xmaxes=paste(2016:2025, "6 30"))

{cairo_pdf("output/figures/LOO_quarterly_flowering_vs_data.pdf", width=6.5, height=3)

ggplot() + 
	geom_rect(data=shading, aes(xmin=ymd(xmins), xmax=ymd(xmaxes), ymin=-Inf, ymax=Inf), fill="gray90") +
	geom_bar(data=flr.raw.sums, aes(x=ymd(paste(year, month, 14)), y=PropFlr), fill="#66c2a4", stat="identity") +
	geom_line(data=LOOsumm, aes(x=ymd(paste(year, month, 14)), y=FreqFlr, group=LOOyear, color=factor(LOOyear)), alpha=0.75, linewidth=0.75) +
	scale_color_manual(values=c("blue", "orange"), name="LOO year") + 
	scale_fill_manual(values=c("#66c2a4"), name="Prop records\nflowering") + 
	labs(x="Date", y="Prop. range with flowering predicted\n or cells with flowering observed", title="Modeled monthly flowering") +
	theme_bw() +
	theme(legend.position="inside", legend.position.inside=c(0.8,1.1), legend.text=element_text(lineheight=0.8), legend.key.height=unit(3, "mm"), legend.direction="horizontal")

}
dev.off()


