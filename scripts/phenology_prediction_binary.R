# Using trained BART models to predict historical Joshua tree flowering
# last used/modified jby, 2026.05.11

# rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/jotr_bonus_blooms")

library("tidyverse")
library("embarcadero")
#library("dbarts")
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
# predictions from PRISM archives: all years 1900-present

# make a place to stash it all
if(!dir.exists("output/models/RIBART_historic_predictions")) dir.create("output/models/RIBART_historic_predictions")

# read the model back in, if necessary
flr.mod.all <- read_rds(paste("output/models/RIbart.model.all.", taxnum, ".rds", sep=""))
topX <- attr(flr.mod.all$fit[[1]]$data@x, "term.labels")

# LOOP over "early" years: three decades, 1901-1930
for(yr in c(1900:2025)){

# yr <- 2016; qu <- 2 # test condition

pred.yr <- NULL
	
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

	preds <- crop(envs, ext(-120, -112, 33, 39))
	
	pred.flr <- rast(predict(flr.mod.all, stack(preds), ri.pred=FALSE, ri.data=qu, ri.name="quarter"))
	crs(pred.flr) <- crs(preds[[1]])
	
	# plot(pred.flr)
	
	pred.yr <- c(pred.yr, pred.flr)
		
	} # end loop over QUARTERS
	
	pred.yr <- rast(pred.yr)
	names(pred.yr) <- paste0(yr, "_Q", 1:4)
	
	writeRaster(pred.yr, paste0("output/models/RIBART_historic_predictions/RIBART_predicted_flowering_",yr,".tiff"), overwrite=TRUE)

} # end loop over YEARS

#-----------------------------------------------------------
# predictions from PRISM archives: posteriors of "all-year" model, current and historical

# make a place to stash it all
if(!dir.exists("output/models/RIBART_all_year_predictions")) dir.create("output/models/RIBART_all_year_predictions")

# read the model back in, if necessary
flr.mod.all <- read_rds(paste("output/models/RIbart.model.all.", taxnum, ".rds", sep=""))
topX <- attr(flr.mod.all$fit[[1]]$data@x, "term.labels")

# LOOP over years
for(yr in c(2016:2025)){

# yr <- 2016; qu <- 2 # test condition
	
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

	preds <- crop(envs, ext(-120, -112, 33, 39))
	preds.df <- as.data.frame(preds)
	
	# prediction with the RI predictor (year) removed
	pred.flr.mat <- pnorm(dbarts:::predict.rbart(flr.mod.all, as.matrix(preds.df[,topX]), group.by=rep(1,nrow(preds.df)), type="bart")) # pnorm gets us back to 0,1 bounds
		
	pred.flr.df <- cbind(as.data.frame(preds[[1]], xy=TRUE)[,c("x","y")], t(pred.flr.mat))
	names(pred.flr.df)[3:302] <- paste0("tree", c(paste0("00",1:9), paste0("0",10:99), 100:300))
	
	# convert it into a SpatRaster
	pred.rast.post <- rast(pred.flr.df, type="xyz", crs=crs(preds[[1]]))
	
	# plot(pred.rast.post[[150]])
	
	writeRaster(pred.rast.post, paste0("output/models/RIBART_all_year_predictions/RIBART_posterior_predicted_flowering_",yr,"_Q",qu,".tiff"), overwrite=TRUE)

		
	} # end loop over QUARTERS

} # end loop over YEARS

#-----------------------------------------------------------
# predictions for counterfactuals


# make a place to stash it all
if(!dir.exists("output/models/RIBART_counterfactual_predictions")) dir.create("output/models/RIBART_counterfactual_predictions")

# read the model back in, if necessary
flr.mod.all <- read_rds(paste("output/models/RIbart.model.all.", taxnum, ".rds", sep=""))
topX <- attr(flr.mod.all$fit[[1]]$data@x, "term.labels")

stanrast <- rast(pd_to_file(prism_archive_subset("tmax", "monthly", resolution="4km", years=2020, mon=1)))

diffs <- rast("output/cmip6_temperature_precipitation_delta_conus_4km_quarterly_2000-2019.nc") |> project(crs(stanrast)) |> resample(stanrast)
diffs[[1]]
plot(diffs[[1]])
names(diffs) # temp 

temp_diff_q1 <- diffs[[0:9*4+1]]
temp_diff_q2 <- diffs[[0:9*4+2]]
temp_diff_q3 <- diffs[[0:9*4+3]]
temp_diff_q4 <- diffs[[0:9*4+4]]

ppt_diff_q1 <- diffs[[0:9*4+41]]
ppt_diff_q2 <- diffs[[0:9*4+42]]
ppt_diff_q3 <- diffs[[0:9*4+43]]
ppt_diff_q4 <- diffs[[0:9*4+44]]

# oh yeah and let's generate a figure of these for the SI?
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library("sf")

# map elements
states <- ne_states(country="united states of america", returnclass="sf")
countries <- ne_countries(scale=10, continent="north america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")

temp_deltas <- crop(c(app(temp_diff_q1,"mean"), app(temp_diff_q2,"mean"), app(temp_diff_q3,"mean"), app(temp_diff_q4,"mean")), extent(-120, -112, 33, 39))
names(temp_deltas) <- paste0("Q",1:4)

temp_deltas_ln <- cbind(as.data.frame(temp_deltas), crds(temp_deltas)) |> pivot_longer(1:4, names_to="quarter", values_to="delta") |> rename(lon=x, lat=y) 
glimpse(temp_deltas_ln)

{cairo_pdf("output/figures/SI_temp_deltas.pdf", width=6.5, height=3)

ggplot() + 

geom_tile(data=temp_deltas_ln, aes(x=lon, y=lat, fill=delta)) + 

#geom_sf(data=coast, color="slategray2", linewidth=3) + 
geom_sf(data=countries, fill=NA, color="white") + 
geom_sf(data=states, fill=NA, color="white") + 
#geom_sf(data=filter(states, name=="California"), fill="cornsilk3", color="antiquewhite4") + 

#geom_sf(data=sdm, fill='antiquewhite3', color=NA, linewidth=0.3, linetype=2) + 

#geom_text(data=filter(statelabs, state!="UT"), aes(label=state, x=lon, y=lat), color="white", size=20, alpha=0.75) +


#geom_sf(data=states, fill=NA, color="antiquewhite3") + 
		
facet_wrap("quarter", nrow=1) +

scale_fill_gradient(low="#4575b4", high="#d73027", name="Temp delta (°C)") + 
labs(x="Longitude", y="Latitude") + 
		
coord_sf(xlim = c(-119, -112.75), ylim = c(33.25, 38.25), expand = TRUE) +

labs(x="Longitude", y="Latitude") +
	
theme_minimal(base_size=12) + 

theme(legend.position="bottom",
	  axis.text=element_blank(), 
	  axis.title=element_blank(),
	  legend.key.height=unit(5, "mm"),
	  legend.key.width=unit(10, "mm"),
	  plot.margin=unit(c(0.01,0.01,0.01,0.01), "inches"), 
	  panel.background=element_rect(fill="slategray3", color=NA), 
	  panel.grid=element_blank()
	  )

}
dev.off()


ppt_deltas <- crop(c(app(ppt_diff_q1,"mean"), app(ppt_diff_q2,"mean"), app(ppt_diff_q3,"mean"), app(ppt_diff_q4,"mean")), extent(-120, -112, 33, 39))
names(ppt_deltas) <- paste0("Q",1:4)

ppt_deltas_ln <- cbind(as.data.frame(ppt_deltas), crds(ppt_deltas)) |> pivot_longer(1:4, names_to="quarter", values_to="delta") |> rename(lon=x, lat=y) 
glimpse(ppt_deltas_ln)

{cairo_pdf("output/figures/SI_ppt_deltas.pdf", width=6.5, height=3)

ggplot() + 

geom_tile(data=ppt_deltas_ln, aes(x=lon, y=lat, fill=delta)) + 

#geom_sf(data=coast, color="slategray2", linewidth=3) + 
geom_sf(data=countries, fill=NA, color="white") + 
geom_sf(data=states, fill=NA, color="white") + 
#geom_sf(data=filter(states, name=="California"), fill="cornsilk3", color="antiquewhite4") + 

#geom_sf(data=sdm, fill='antiquewhite3', color=NA, linewidth=0.3, linetype=2) + 

#geom_text(data=filter(statelabs, state!="UT"), aes(label=state, x=lon, y=lat), color="white", size=20, alpha=0.75) +


#geom_sf(data=states, fill=NA, color="antiquewhite3") + 
		
facet_wrap("quarter", nrow=1) +

scale_fill_gradient2(low="#8c510a", mid="white", high="#01665e", name="Precip delta (mm)") + 
labs(x="Longitude", y="Latitude") + 
		
coord_sf(xlim = c(-119, -112.75), ylim = c(33.25, 38.25), expand = TRUE) +

labs(x="Longitude", y="Latitude") +
	
theme_minimal(base_size=12) + 

theme(legend.position="bottom",
	  axis.text=element_blank(), 
	  axis.title=element_blank(),
	  legend.key.height=unit(5, "mm"),
	  legend.key.width=unit(10, "mm"),
	  plot.margin=unit(c(0.01,0.01,0.01,0.01), "inches"), 
	  panel.background=element_rect(fill="slategray3", color=NA), 
	  panel.grid=element_blank()
	  )

}
dev.off()






# LOOP over in-sample years: 2016-2025
for(yr in c(2016:2017)){

# yr <- 2020; qu <- 2 # test condition
	
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
	
	names(envs) <- paste0(rep(c("tmax", "tmin", "ppt"),5),"_Q",rep(0:4, each=3))
	
	# loop over the individual GCS differentials
	for(gcm in 1:10){
	# now counterfactual envs
		if(qu == 1){
		cfenvs <- c(envs[["tmax_Q0"]] - temp_diff_q1[[gcm]], 
				envs[["tmin_Q0"]] - temp_diff_q1[[gcm]],
				envs[["ppt_Q0"]] - ppt_diff_q1[[gcm]],
				envs[["tmax_Q1"]] - temp_diff_q2[[gcm]], 
				envs[["tmin_Q1"]] - temp_diff_q2[[gcm]],
				envs[["ppt_Q1"]] - ppt_diff_q2[[gcm]],
				envs[["tmax_Q2"]] - temp_diff_q3[[gcm]], 
				envs[["tmin_Q2"]] - temp_diff_q3[[gcm]],
				envs[["ppt_Q2"]] - ppt_diff_q3[[gcm]],
				envs[["tmax_Q3"]] - temp_diff_q4[[gcm]], 
				envs[["tmin_Q3"]] - temp_diff_q4[[gcm]],
				envs[["ppt_Q3"]] - ppt_diff_q4[[gcm]],
				envs[["tmax_Q4"]] - temp_diff_q1[[gcm]], 
				envs[["tmin_Q4"]] - temp_diff_q1[[gcm]],
				envs[["ppt_Q4"]] - ppt_diff_q1[[gcm]]
				)	
		}
		if(qu == 2){
		cfenvs <- c(envs[["tmax_Q0"]] - temp_diff_q2[[gcm]], 
				envs[["tmin_Q0"]] - temp_diff_q2[[gcm]],
				envs[["ppt_Q0"]] - ppt_diff_q2[[gcm]],
				envs[["tmax_Q1"]] - temp_diff_q3[[gcm]], 
				envs[["tmin_Q1"]] - temp_diff_q3[[gcm]],
				envs[["ppt_Q1"]] - ppt_diff_q3[[gcm]],
				envs[["tmax_Q2"]] - temp_diff_q4[[gcm]], 
				envs[["tmin_Q2"]] - temp_diff_q4[[gcm]],
				envs[["ppt_Q2"]] - ppt_diff_q4[[gcm]],
				envs[["tmax_Q3"]] - temp_diff_q1[[gcm]], 
				envs[["tmin_Q3"]] - temp_diff_q1[[gcm]],
				envs[["ppt_Q3"]] - ppt_diff_q1[[gcm]],
				envs[["tmax_Q4"]] - temp_diff_q2[[gcm]], 
				envs[["tmin_Q4"]] - temp_diff_q2[[gcm]],
				envs[["ppt_Q4"]] - ppt_diff_q2[[gcm]]
				)	
		}
		if(qu == 3){
		cfenvs <- c(envs[["tmax_Q0"]] - temp_diff_q3[[gcm]], 
				envs[["tmin_Q0"]] - temp_diff_q3[[gcm]],
				envs[["ppt_Q0"]] - ppt_diff_q3[[gcm]],
				envs[["tmax_Q1"]] - temp_diff_q4[[gcm]], 
				envs[["tmin_Q1"]] - temp_diff_q4[[gcm]],
				envs[["ppt_Q1"]] - ppt_diff_q4[[gcm]],
				envs[["tmax_Q2"]] - temp_diff_q1[[gcm]], 
				envs[["tmin_Q2"]] - temp_diff_q1[[gcm]],
				envs[["ppt_Q2"]] - ppt_diff_q1[[gcm]],
				envs[["tmax_Q3"]] - temp_diff_q2[[gcm]], 
				envs[["tmin_Q3"]] - temp_diff_q2[[gcm]],
				envs[["ppt_Q3"]] - ppt_diff_q2[[gcm]],
				envs[["tmax_Q4"]] - temp_diff_q3[[gcm]], 
				envs[["tmin_Q4"]] - temp_diff_q3[[gcm]],
				envs[["ppt_Q4"]] - ppt_diff_q3[[gcm]]
				)	
		}
		if(qu == 4){
		cfenvs <- c(envs[["tmax_Q0"]] - temp_diff_q4[[gcm]], 
				envs[["tmin_Q0"]] - temp_diff_q4[[gcm]],
				envs[["ppt_Q0"]] - ppt_diff_q4[[gcm]],
				envs[["tmax_Q1"]] - temp_diff_q1[[gcm]], 
				envs[["tmin_Q1"]] - temp_diff_q1[[gcm]],
				envs[["ppt_Q1"]] - ppt_diff_q1[[gcm]],
				envs[["tmax_Q2"]] - temp_diff_q2[[gcm]], 
				envs[["tmin_Q2"]] - temp_diff_q2[[gcm]],
				envs[["ppt_Q2"]] - ppt_diff_q2[[gcm]],
				envs[["tmax_Q3"]] - temp_diff_q3[[gcm]], 
				envs[["tmin_Q3"]] - temp_diff_q3[[gcm]],
				envs[["ppt_Q3"]] - ppt_diff_q3[[gcm]],
				envs[["tmax_Q4"]] - temp_diff_q4[[gcm]], 
				envs[["tmin_Q4"]] - temp_diff_q4[[gcm]],
				envs[["ppt_Q4"]] - ppt_diff_q4[[gcm]]
				)	
		}	

		preds <- crop(cfenvs, extent(-120, -112, 33, 39))
		
		preds.df <- as.data.frame(preds)
	
		# prediction with the RI predictor (year) removed
		pred.flr.mat <- pnorm(dbarts:::predict.rbart(flr.mod.all, as.matrix(preds.df[,topX]), group.by=rep(NA,nrow(preds.df)), type='bart'))
		pred.flr.df <- cbind(as.data.frame(preds[[1]], xy=TRUE)[,c("x","y")], t(pred.flr.mat))
		names(pred.flr.df)[3:302] <- paste0("tree", c(paste0("00",1:9), paste0("0",10:99), 100:300))
	
		# convert it into a SpatRaster
		pred.rast.post <- rast(pred.flr.df, type="xyz", crs=crs(preds[[1]]))
	
		writeRaster(pred.rast.post, paste0("output/models/RIBART_counterfactual_predictions/RIBART_posterior_predicted_flowering_GCM", gcm, "_", yr, "_Q", qu, ".tiff"), overwrite=TRUE)
		
		} # end loop over GCMs
		
	} # end loop over QUARTERS
	
} # end loop over YEARS


#-------------------------------------------------------------------------
# figure illustrating quarterly predictions

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


# need to figure out 

Qpreds <- rast(paste0("output/models/RIBART_all_year_predictions/RIBART_predicted_flowering_", 2016:2025, ".tiff"))

tsdm <- st_transform(sdm, crs=crs(Qpreds))

Qpreds.mask <- mask(Qpreds, tsdm)

Qpreds.mask.ln <- cbind(as.data.frame(Qpreds.mask), crds(Qpreds.mask)) |> pivot_longer(1:40, names_to="yrqu", values_to="pred_prFlr") |> rename(lon=x, lat=y) |> mutate(year=gsub("(\\d+)_Q\\d", "\\1", yrqu), quarter=gsub("\\d+_(Q\\d)", "\\1", yrqu))
glimpse(Qpreds.mask.ln)

# h=530, w=705

{cairo_pdf(paste("output/figures/predicted_distribution_quarterly_", taxnum, ".pdf", sep=""), width=14, height=5.7)

ggplot() + 

#geom_sf(data=coast, color="slategray2", linewidth=3) + 
geom_sf(data=countries, fill="antiquewhite3", color="antiquewhite4") + 
geom_sf(data=states, fill="antiquewhite2", color="antiquewhite3") + 
#geom_sf(data=filter(states, name=="California"), fill="cornsilk3", color="antiquewhite4") + 

#geom_sf(data=sdm, fill='antiquewhite3', color=NA, linewidth=0.3, linetype=2) + 

#geom_text(data=filter(statelabs, state!="UT"), aes(label=state, x=lon, y=lat), color="white", size=20, alpha=0.75) +

geom_tile(data=Qpreds.mask.ln, aes(x=lon, y=lat, fill=pred_prFlr>thresh)) + 

#geom_sf(data=states, fill=NA, color="antiquewhite3") + 
		
facet_grid(quarter~year, switch="y") +

scale_fill_manual(values=c("#7fbc41", "#762a83"), name="Flowering observed") + 
labs(x="Longitude", y="Latitude") + 
		
coord_sf(xlim = c(-119, -112.75), ylim = c(33.25, 38.25), expand = TRUE) +

labs(x="Longitude", y="Latitude") +
	
theme_minimal(base_size=20) + 

theme(legend.position="none",
	  axis.text=element_blank(), 
	  axis.title=element_blank(), 
	  plot.margin=unit(c(0.01,0.01,0.01,0.01), "inches"), 
	  panel.background=element_rect(fill="slategray3", color=NA), 
	  panel.grid=element_blank()
	  )

}
dev.off()

#-------------------------------------------------------------------------
# figure comparing factual and counterfactual JUST for anomalies

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


# need to figure out 

Qpreds <- rast(paste0("output/models/RIBART_all_year_predictions/RIBART_predicted_flowering_", 2016:2025, ".tiff"))

tsdm <- st_transform(sdm, crs=crs(Qpreds))

Qpreds.mask <- mask(Qpreds, tsdm)

Qpreds.mask.ln <- cbind(as.data.frame(Qpreds.mask), crds(Qpreds.mask)) |> pivot_longer(1:40, names_to="yrqu", values_to="pred_prFlr") |> rename(lon=x, lat=y) |> mutate(year=gsub("(\\d+)_Q\\d", "\\1", yrqu), quarter=gsub("\\d+_(Q\\d)", "\\1", yrqu))
glimpse(Qpreds.mask.ln)

# h=530, w=705

{cairo_pdf(paste("output/figures/predicted_distribution_quarterly_", taxnum, ".pdf", sep=""), width=14, height=5.7)

ggplot() + 

#geom_sf(data=coast, color="slategray2", linewidth=3) + 
geom_sf(data=countries, fill="antiquewhite3", color="antiquewhite4") + 
geom_sf(data=states, fill="antiquewhite2", color="antiquewhite3") + 
#geom_sf(data=filter(states, name=="California"), fill="cornsilk3", color="antiquewhite4") + 

#geom_sf(data=sdm, fill='antiquewhite3', color=NA, linewidth=0.3, linetype=2) + 

#geom_text(data=filter(statelabs, state!="UT"), aes(label=state, x=lon, y=lat), color="white", size=20, alpha=0.75) +

geom_tile(data=Qpreds.mask.ln, aes(x=lon, y=lat, fill=pred_prFlr>thresh)) + 

#geom_sf(data=states, fill=NA, color="antiquewhite3") + 
		
facet_grid(quarter~year, switch="y") +

scale_fill_manual(values=c("#7fbc41", "#762a83"), name="Flowering observed") + 
labs(x="Longitude", y="Latitude") + 
		
coord_sf(xlim = c(-119, -112.75), ylim = c(33.25, 38.25), expand = TRUE) +

labs(x="Longitude", y="Latitude") +
	
theme_minimal(base_size=20) + 

theme(legend.position="none",
	  axis.text=element_blank(), 
	  axis.title=element_blank(), 
	  plot.margin=unit(c(0.01,0.01,0.01,0.01), "inches"), 
	  panel.background=element_rect(fill="slategray3", color=NA), 
	  panel.grid=element_blank()
	  )

}
dev.off()

