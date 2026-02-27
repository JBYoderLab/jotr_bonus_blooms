# Using trained BART models to predict historical Joshua tree flowering
# last used/modified jby, 2026.02.26

# rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/jotr_bonus_blooms")

library("tidyverse")
library("embarcadero")
library("cowplot")

library("prism")
library("terra")
library("geosphere")
library("barrks")

set.seed(18090212)

prism_set_dl_dir("../data/PRISM") # expect an error, oy

#-----------------------------------------------------------
# initial file loading

# set parameters as variables
taxnum <- 1595251 # Joshua trees
# Prunus ilicifolia = 57250

flow <- read.csv(paste("output/flowering_freq_climate_", taxnum, ".csv", sep="")) |> mutate(flr = prop_flr > 0) |> filter(year>=2016)# flowering/not flowering, biologically-informed candidate predictors

glimpse(flow) # 5070 from 2016 on

hist(flow$prop_flr) # check the cutoff I've set for binary flowering

table(flow$year, flow$flr)
table(flow$year, flow$quarter)

flow.norm <- filter(flow, quarter <= 2)
flow.anom <- filter(flow, quarter > 2)


#-----------------------------------------------------------
# predictions from PRISM archives: "all-year" model, current and historical

# make a place to stash it all
if(!dir.exists("output/models/RIBART_all_year_predictions")) dir.create("output/models/RIBART_all_year_predictions")

# read the model back in, if necessary
flr.mod.all <- read_rds(paste("output/models/RIbart.model.all.", taxnum, ".rds", sep=""))
topX <- attr(flr.mod.all$fit[[1]]$data@x, "term.labels")

# LOOP over "early" years: three decades, 1901-1930
for(yr in c(1901:1930, 1996:2025)){

# yr <- 2020; qu <- 2 # test condition
	
	yearpred <- NULL

	for(qu in 1:4){
	
	# hmm goddammit, I have to figure out how to account for backing into a second year
	if(qu == 1){
		tmaxfiles <- c(prism_archive_subset("tmax", "monthly", resolution="4km", years=yr, mon=1:3),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr-1, mon=10:12),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr-1, mon=7:9),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr-1, mon=4:6),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr-1, mon=1:3) )
		
		tminfiles <- c(prism_archive_subset("tmin", "monthly", resolution="4km", years=yr, mon=1:3),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr-1, mon=10:12),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr-1, mon=7:9),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr-1, mon=4:6),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr-1, mon=1:3) )
			
		pptfiles <- c(prism_archive_subset("ppt", "monthly", resolution="4km", years=yr, mon=1:3),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr-1, mon=10:12),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr-1, mon=7:9),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr-1, mon=4:6),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr-1, mon=1:3) )
	}
	if(qu == 2){
		tmaxfiles <- c(prism_archive_subset("tmax", "monthly", resolution="4km", years=yr, mon=4:6),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr, mon=1:3),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr-1, mon=12:10),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr-1, mon=7:9),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr-1, mon=4:6) )
		
		tminfiles <- c(prism_archive_subset("tmin", "monthly", resolution="4km", years=yr, mon=4:6),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr, mon=1:3),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr-1, mon=10:12),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr-1, mon=7:9),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr-1, mon=4:6) )
			
		pptfiles <- c(prism_archive_subset("ppt", "monthly", resolution="4km", years=yr, mon=4:6),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr, mon=1:3),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr-1, mon=10:12),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr-1, mon=7:9),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr-1, mon=4:6) )
	}
	if(qu == 3){
		tmaxfiles <- c(prism_archive_subset("tmax", "monthly", resolution="4km", years=yr, mon=7:9),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr, mon=4:6),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr, mon=1:3),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr-1, mon=10:12),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr-1, mon=7:9) )
		
		tminfiles <- c(prism_archive_subset("tmin", "monthly", resolution="4km", years=yr, mon=7:9),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr, mon=4:6),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr, mon=1:3),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr-1, mon=10:12),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr-1, mon=7:9) )
			
		pptfiles <- c(prism_archive_subset("ppt", "monthly", resolution="4km", years=yr, mon=7:9),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr, mon=4:6),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr, mon=1:3),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr-1, mon=10:12),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr-1, mon=7:9) )
	}
	if(qu == 4){
		tmaxfiles <- c(prism_archive_subset("tmax", "monthly", resolution="4km", years=yr, mon=10:12),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr, mon=7:9),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr, mon=4:6),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr, mon=1:3),
			prism_archive_subset("tmax", "monthly", resolution="4km", years=yr-1, mon=10:12) )
		
		tminfiles <- c(prism_archive_subset("tmin", "monthly", resolution="4km", years=yr, mon=10:12),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr, mon=7:9),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr, mon=4:6),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr, mon=1:3),
			prism_archive_subset("tmin", "monthly", resolution="4km", years=yr-1, mon=10:12) )
			
		pptfiles <- c(prism_archive_subset("ppt", "monthly", resolution="4km", years=yr, mon=10:12),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr, mon=7:9),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr, mon=4:6),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr, mon=1:3),
			prism_archive_subset("ppt", "monthly", resolution="4km", years=yr-1, mon=10:12) )
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
	pred.flr <- predict(flr.mod.all, stack(preds), splitby=20, ri.name="quarter", ri.data=NA, ri.pred=FALSE)

	pred.flr # BOOM

	crs(pred.flr) <- crs(envs)
	
	yearpred <- c(yearpred, pred.flr)
		
	} # end loop over QUARTERS

outyear <- rast(stack(yearpred))
names(outyear) <- paste0(yr,"_Q",1:4)

writeRaster(outyear, paste0("output/models/RIBART_all_year_predictions/RIBART_predicted_flowering_",yr,".tiff"), overwrite=TRUE)

} # end loop over YEARS









