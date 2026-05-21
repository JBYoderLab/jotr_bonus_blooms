# working with phenology-annotated iNat observations
# Assumes local environment 
# jby 2026.05.05

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/jotr_bonus_blooms")

library("tidyverse")
library("terra")
library("sf")
library("prism")

prism_set_dl_dir("../data/PRISM") # expect an error, oy

# set parameters as variables
taxnum <- 1595251 # Joshua trees

#-------------------------------------------------------------------------
# read in iNat observations compiled using inat_phenology_obs.R

inat <- read.csv(paste0("data/inat_phenology_data_", taxnum, ".csv"), h=TRUE) %>% mutate(observed_on = ymd(observed_on), quarter = if_else(month%in%1:3, 1, if_else(month%in%4:6, 2, if_else(month%in%7:9, 3, 4))))

glimpse(inat) # how many raw observations?
table(inat$phenology) # by phenophase
table(inat$phenology, inat$year) # by phenophase and year
table(inat$phenology, inat$quarter) # by phenophase and quarter
table(inat$year, inat$quarter)

#-------------------------------------------------------------------------
# organize iNat observations for extraction of summarized PRISM data

prism_temp_rast <- rast(paste("../data/PRISM/quarterlies/tmax_cropped_2010Q1.bil", sep="")) # raster grid base

# data structure setup
flowering <- data.frame(matrix(0,0,6))
names(flowering) <- c("lon","lat", "year", "quarter", "prop_flr", "n_obs")

# then LOOP over years in the raw data ....
for(yr in unique(inat$year)){
	for(qu in sort(unique(filter(inat, year==yr)$quarter))){
	
	# yr <- 2020; qu <- 3
	
	Qdat <- dplyr::filter(inat, year==yr, quarter==qu)
	
	obs <- rasterize(as.matrix(Qdat[,c("longitude","latitude")]), prism_temp_rast, fun=length, background=NA)
	
	if(length(which(Qdat$phenology%in%c("Flower Budding", "Flowering")))>0){
		flr <- rasterize(as.matrix(dplyr::filter(Qdat, phenology%in%c("Flower Budding", "Flowering"))[,c("longitude","latitude")]), prism_temp_rast, fun=length, background=0)
		
		flrfrq <- flr/obs
		
		}else{
		
		flrfrq <- 0/obs
		
		}
		
	outp <- data.frame(lon=crds(flrfrq)[,"x"], lat=crds(flrfrq)[,"y"], year=yr, quarter=qu, prop_flr=as.data.frame(flrfrq)$values, n_obs=as.data.frame(obs)$values)
	
	# put it all together
	flowering <- rbind(flowering, outp)
	
	cat("Done with year", yr, "quarter", qu, "\n")
	} # END LOOP over quarters
} # END LOOP over years

head(flowering)
glimpse(flowering) # okay okay okay!
table(flowering$year) # smiling serenely
table(flowering$quarter)
hist(flowering$prop_flr)
hist(flowering$n_obs)

write.table(flowering, paste("output/flowering_freq_rasterized_", taxnum,".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE)
# flowering <- read.csv(paste("output/flowering_freq_rasterized_", taxnum,".csv", sep=""))

#-------------------------------------------------------------------------
# attach PRISM data to flowering/not flowering observations

# we're going to start with PRISM data for the four months prior to the month of observation
# this follows Brenskelle, more or less (they used 120 days)

prism_temp_rast <- rast(paste("../data/PRISM/quarterlies/tmax_cropped_2010Q1.bil", sep="")) # raster grid base

# new data structure
flr.clim <- data.frame(matrix(0,0,ncol(flowering)+15))
names(flr.clim) <- c(colnames(flowering), paste0(rep(c("tmax", "tmin", "ppt"),5),"_","Q",rep(0:4, each=3)))

# LOOP over years and months ...
for(yr in sort(unique(flowering$year))){
	for(qu in sort(unique(filter(flowering, year==yr)$quarter))){
		
	# yr <- 2020; qu <- 3 # test condition

	# pull subset of flowering observations for year
	flsub <- filter(flowering, year==yr, quarter==qu)
	
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

	flsubenvs <- cbind(flsub, terra::extract(envs, flsub[,c("lon","lat")])[,-1])
	
	flr.clim <- rbind(flr.clim, flsubenvs) 
	
	write.table(flr.clim, paste0("output/flowering_freq_climate_", taxnum, ".csv"), sep=",", col.names=TRUE, row.names=FALSE)
	
	cat("Done with year", yr, "quarter", qu, "\n")

	} # END LOOP over quarters
} # END LOOP over years

glimpse(flr.clim)

table(flr.clim$year, flr.clim$quarter)


# I do wonder about additional predictors, specifically day length ...

library("chillR") # would give this ...

flr.clim.daylen <- flr.clim |> mutate(tmean=(tmax+tmin)/2, Jdate = as.POSIXlt(paste(year,month,01,sep="-"))$yday)
flr.clim.daylen$daylen <- mapply(function(L,D) daylength(L, D)$Daylength, flr.clim.daylen$lat, flr.clim.daylen$Jdate)

glimpse(flr.clim.daylen) # boom

table(flr.clim$year, flr.clim$month)

write.table(flr.clim.daylen, paste0("output/flowering_freq_climate_", taxnum, ".csv"), sep=",", col.names=TRUE, row.names=FALSE)

# and that's generated a data file we can feed into embarcadero ... in the next script!

#-------------------------------------------------------------------------
# map gridded records

flr.clim <- read.csv(paste0("output/flowering_freq_climate_", taxnum, ".csv"))
glimpse(flr.clim)

flr.clim.summed <- flr.clim |> filter(year>=2016) |> group_by(lat, lon) |> summarize(tot_obs = sum(n_obs))
glimpse(flr.clim.summed)

library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library("sf")

# map elements
states <- ne_states(country="united states of america", returnclass="sf")
countries <- ne_countries(scale=10, continent="north america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")

sdm <- read_sf("../data/Yucca/Jotr_SDM2023_range/jotr_SDM2023_range_simple.shp")

statelabs <- data.frame(state=c("CA","NV","UT","AZ"), lon=c(-116.75,-115.75,-113.5,-113.45), lat=c(35.25,37.65,37.5,35.5))


{cairo_pdf(paste("output/figures/record_distribution_map_", taxnum, ".pdf", sep=""), width=6, height=6)

ggplot() + 

geom_sf(data=coast, color="slategray2", linewidth=3) + 
geom_sf(data=countries, fill="antiquewhite3", color="antiquewhite4") + 
geom_sf(data=states, fill="antiquewhite2", color="antiquewhite4") + 
#geom_sf(data=filter(states, name=="California"), fill="cornsilk3", color="antiquewhite4") + 

geom_sf(data=sdm, fill='antiquewhite3', color=NA, linewidth=0.3, linetype=2) + 

geom_text(data=filter(statelabs, state!="UT"), aes(label=state, x=lon, y=lat), color="white", size=20, alpha=0.75) +

geom_tile(data=flr.clim.summed, aes(x=lon, y=lat, fill=log10(tot_obs))) + 

geom_sf(data=states, fill=NA, color="antiquewhite4") + 
		
scale_fill_gradient(low="#a1d99b", high="#00441b", name="Records per cell", breaks=c(0,1,2), labels=c(1, 10, 100)) + 
labs(x="Longitude", y="Latitude") + 
		
coord_sf(xlim = c(-119, -112.75), ylim = c(33.25, 38.25), expand = TRUE) +

labs(x="Longitude", y="Latitude") +

annotation_scale(location = "bl", width_hint = 0.3) + 
annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.15, "in"), pad_y = unit(0.25, "in"), style = north_arrow_fancy_orienteering, height=unit(1, "in"), width=unit(0.65, "in")) +
	
theme_minimal(base_size=14) + 

theme(legend.position="inside", 
	  legend.position.inside=c(0.7,0.07), 
	  legend.key.width=unit(0.3, "inches"), 
	  legend.key.height=unit(0.1, "in"), 
	  legend.direction="horizontal", 
	  legend.box.spacing=unit(0.01,"inches"), 
	  legend.box="horizontal", 
	  legend.text=element_text(size=12), 
	  legend.title=element_text(size=14, margin=margin(0,5,0,2, unit='mm')), 
	  legend.background=element_rect(fill="white", color=NA),
	  axis.text=element_blank(), 
	  axis.title=element_blank(), 
	  plot.margin=unit(c(0.01,0.01,0.01,0.01), "inches"), 
	  panel.background=element_rect(fill="slategray3", color="black"), 
	  panel.grid=element_blank()
	  )

}
dev.off()




