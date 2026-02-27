# Scraping phenology-annotated iNat observations
# Assumes local environment 
# jby 2026.02.23

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/jotr_bonus_blooms")

library("tidyverse")
library("sf")
library("terra")
library("rinat")


#-------------------------------------------------------------------------
# Pull down iNat observations of Joshua tree with specific phenology code

# term id: 12 for Plant Phenology then term_id_value: 13 =Flowering, 14 =Fruiting, 15 =Flower Budding

# trial run, to make sure it works as expected
taxnum <- 1595251


test <- get_inat_obs(quality="research", place_id=53170, taxon_id=taxnum, annotation=c(12,13), year=2021, maxresults=1e4) 

glimpse(test)

# years to read in for the loop
years <- 2008:2025
verbose <- TRUE # read out progress
max_loc_uncertainty <- 2000
write.out <- TRUE

# set up dataframe
inat_pheno_data <- data.frame(matrix(0,0,8))
names(inat_pheno_data) <- c("scientific_name", "latitude", "longitude", "url", "image_url", "observed_on", "phenology", "year")

# LOOP over years, downloading by phenophase
for(y in years){

    bud.y <- try(rinat::get_inat_obs(quality="research", taxon_id=taxnum, annotation=c(12, 15), year=y, maxresults=1e4))
    Sys.sleep(5) # throttling under the API limit, maybe?
    flo.y <- try(rinat::get_inat_obs(quality="research", taxon_id=taxnum, annotation=c(12, 13), year=y, maxresults=1e4))
    Sys.sleep(5)
    fru.y <- try(rinat::get_inat_obs(quality="research", taxon_id=taxnum, annotation=c(12, 14), year=y, maxresults=1e4))
    Sys.sleep(5)
    non.y <- try(rinat::get_inat_obs(quality="research", taxon_id=taxnum, annotation=c(12, 21), year=y, maxresults=1e4))
    Sys.sleep(5)


    if(class(bud.y)=="data.frame") bud.o <- bud.y |> dplyr::filter(captive_cultivated=="false", positional_accuracy < max_loc_uncertainty) |> dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) |> dplyr::mutate(phenology="Flower Budding", year=gsub("(\\d{4})-.+","\\1", observed_on)) else bud.o <- NULL

    if(class(flo.y)=="data.frame") flo.o <- flo.y |> dplyr::filter(captive_cultivated=="false", positional_accuracy < max_loc_uncertainty) |> dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) |> dplyr::mutate(phenology="Flowering", year=gsub("(\\d{4})-.+","\\1", observed_on)) else flo.o <- NULL

    if(class(fru.y)=="data.frame") fru.o <- fru.y |> dplyr::filter(captive_cultivated=="false", positional_accuracy < max_loc_uncertainty) |> dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) |> dplyr::mutate(phenology="Fruiting", year=gsub("(\\d{4})-.+","\\1", observed_on)) else fru.o <- NULL

    if(class(non.y)=="data.frame") non.o <- non.y |> dplyr::filter(captive_cultivated=="false", positional_accuracy < max_loc_uncertainty) |> dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) |> dplyr::mutate(phenology="No Evidence of Flowering", year=gsub("(\\d{4})-.+","\\1", observed_on)) else non.o <- NULL

    inat_pheno_data <- rbind(inat_pheno_data, bud.o, flo.o, fru.o, non.o)

    if(write.out){
      # insert a check here for a valid directory and filename?
      if(!file.exists("data")) dir.create("data") # make sure there's a folder to write to!
      utils::write.table(inat_pheno_data, paste("data/inat_phenology_data_", taxnum, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
    }

    # provide some indication of progress
    if(verbose) cat("\nDownloaded", nrow(rbind(bud.o,flo.o,fru.o,non.o)), "records from", y, "\n\n")

}

# expect error messages if searches return zero obs with a given phenology status; this may not be a problem, but see what the final data table looks like
glimpse(inat_pheno_data) # 14,189 records on 2026.01.24
filter(inat_pheno_data, year>=2008) %>% glimpse()
table(inat_pheno_data$year, inat_pheno_data$phenology)


#-------------------------------------------------------------------------
# More adjustments

inat_pheno_data$month <- as.numeric(gsub(".+-(\\d{2})-.+","\\1", inat_pheno_data$observed_on))

glimpse(inat_pheno_data)

write.table(inat_pheno_data, paste0("data/inat_phenology_data_", taxnum, ".csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)


#-------------------------------------------------------------------------
# visualize, if you like
inat_pheno_data <- read.csv(paste0("data/inat_phenology_data_", taxnum, ".csv"), h=TRUE)

flr.raw.tots <- inat_pheno_data |> group_by(year, month) |> summarize(Ntot=length(phenology))
flr.raw.flow <- inat_pheno_data |> filter(phenology%in%c("Flowers Budding", "Flowering")) |> group_by(year, month) |> summarize(Nflr=length(phenology))

flr.raw.sums <- flr.raw.tots |> left_join(flr.raw.flow)
flr.raw.sums$Nflr[is.na(flr.raw.sums$Nflr)] <- 0
flr.raw.sums$PropFlr <- flr.raw.sums$Nflr/flr.raw.sums$Ntot

glimpse(flr.raw.sums)


{cairo_pdf("output/figures/iNat_obs_raw.pdf", width=4.5, height=9)

ggplot(flr.raw.sums, aes(x=month, y=PropFlr)) + 
	geom_bar(stat="identity") +
	facet_grid(rows="year") +
	scale_x_continuous(breaks=1:12) + scale_y_continuous(breaks=c(0,0.5,1)) +
	labs(y="Proportion of records with flowers or buds", x="Month of observation") +
	theme_bw()
	
}
dev.off()

flr.sums.ln <- flr.raw.tots |> left_join(flr.raw.flow) |> mutate(Nnon = Ntot-Nflr) |> pivot_longer(all_of(c("Ntot", "Nflr", "Nnon")), names_to="count_of", values_to="count") |> mutate(count_of=factor(count_of, c("Ntot", "Nnon", "Nflr"))) |> ungroup()

shading <- data.frame(xmins=ymd(paste(2008:2025, "1 1")),xmaxes=ymd(paste(2008:2025, "6 30")) )

{cairo_pdf("output/figures/iNat_obs_raw.pdf", width=6.5, height=3)

ggplot() + 
	geom_rect(data=shading, aes(xmin=xmins, xmax=xmaxes, ymin=-Inf, ymax=Inf), fill="gray90") +
	geom_bar(data=filter(flr.sums.ln, count_of!="Ntot"), aes(x=ymd(paste(year, month, 1)), y=count, fill=count_of), position="stack", stat="identity") +
	annotate("text", label="in season", angle=90, color="white", x=ymd("2008 03 15"), y=200) +
	scale_fill_manual(values=c("#a6cee3", "#33a02c"), labels=c("No flowers, or fruits", "Flowers or buds"), name=NULL) +
	labs(y="Records", x="Date observed") +
	xlim(ymd("2008 01 01"),ymd("2025 12 31")) +	
	theme_bw() +
	theme(legend.position="inside", legend.position.inside=c(0.175,0.775), legend.key.size=unit(3, "mm"), legend.background=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
	
}
dev.off()

