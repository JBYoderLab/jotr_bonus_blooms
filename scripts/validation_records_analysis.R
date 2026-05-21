# Comparing predicted historical anomalous flowering to independent validation records
# Assumes local environment
# jby 2026.05.11

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/jotr_bonus_blooms")

library("tidyverse")

library("terra")
library("sp")
library("sf")
library("hexbin")
library("embarcadero")
library("cowplot")

source("../../Active_projects/shared/Rscripts/base.R") # my special mix of personal functions
source("../../Active_projects/shared/Rscripts/base_graphics.R") # my special mix of personal functions

set.seed(18090212)

taxnum <- 1595251


#-------------------------------------------------------------------------
# Load up and inspect NPN records

# US NPN data
npn <- read.csv("data/validation/NPN_datasheet_1778265717004/status_intensity_observation_data.csv")

glimpse(npn)
length(unique(npn$Site_ID)) # 17 unique sites
length(unique(paste(npn$Longitude, npn$Latitude))) # confirm 17 locations at resolution of this data

table(year(ymd(npn$Observation_Date)), npn$Site_ID) # records for 16 years, 2011-2026
sum(table(year(ymd(npn$Observation_Date)), npn$Site_ID)>0) # works out to 115 location-year combos

table(npn$Phenophase_Category) # ONLY flowers or fruits recorded

table(npn$Intensity_Category_ID)
# 23 = fruits present
# 24 = ripe fruits present
# 25 = recent fruit drop
# 31 = open flowers (peak)
# 35 = flowers or flower heads present
# 48 = flowers and flower buds 
# 50 = open flowers percentage (individual)
# 56 = fruits present
# 58 = ripe fruit percentage
# 59 = recent fruit or seed drop

# How are these organized?
filter(npn, month(ymd(Observation_Date))>=6, grepl("flower", Phenophase_Description)) %>% dplyr::select(Observation_Date, Phenophase_Description) # Holy shit I wasn't even trying to find that

# where are these located?
sdm.pres <- read_sf("../data/Yucca/Jotr_SDM2023_range/Jotr_SDM2023_range_simple.shp")
ggplot() + geom_sf(data=sdm.pres) + geom_point(data=npn, aes(x=Longitude, y=Latitude))

filter(npn, Longitude < -122) # HUH
filter(npn, Longitude > -113) %>% group_by(year(ymd(Observation_Date))) %>% summarize(n=length(Observation_Date))
filter(npn, Longitude > -113, year(ymd(Observation_Date))==2018) %>% group_by(Phenophase_Category) %>% summarize(n=length(Observation_Date))

# okay, let's try this ... natural range records only
npn.vobs <- npn %>% filter(Longitude > -120, Longitude < -114) |> 
				   mutate(year = year(ymd(Observation_Date)),
				   		  month = month(ymd(Observation_Date)),
				   		  quarter = if_else(month < 4, 1, if_else(month < 7, 2, if_else(month < 10, 3, 4)))
				   		) |> 
				   dplyr::select(Longitude, Latitude, Site_ID, Observation_Date, year, month, quarter, Phenophase_Description, Phenophase_Status) |>
				   rename(lon=Longitude, lat=Latitude, location=Site_ID) 

glimpse(npn.vobs)

npn.val <- npn.vobs |> filter(grepl("flower", Phenophase_Description), Phenophase_Status>=0) |> 
				group_by(year, quarter, location, lat, lon) |> 
				summarize(nobs = sum(Phenophase_Status)) |>  # Phenophase_Status = 1 if present, 0 if absent
				mutate(flr = nobs>0, obs_by="npn", location=as.character(location)) |> 
				dplyr::select(lat, lon, location, year, obs_by, flr)

glimpse(npn.val)

table(npn.val$year, npn.val$quarter) # Oookay, this is solid, I think
table(npn.val$flr)

#-------------------------------------------------------------------------
# Herbarium records

# Cal Consortium of Herbaria
cch <- read.csv("data/validation/SymbOutput_2024-02-11_233150_DwC-A/occurrences_cleaned.csv")

glimpse(cch)
table(cch$year) # LMAO back to 1876??
hist(cch$coordinateUncertaintyInMeters) # hoo boy

cch.val <- cch %>% filter(year>=1900, coordinateUncertaintyInMeters<=4000) |> 
 	mutate(quarter = if_else(month < 4, 1, if_else(month < 7, 2, if_else(month < 10, 3, 4)))) |> 
	dplyr::select(year, quarter, location, lat, lon, obs_by, flr)

glimpse(cch.val) # okay nice
table(cch.val$year, cch.val$quarter)
table(cch.val$flr)

# TORCH (and compatible)
torch <- read.csv("data/validation/SymbOutput_2024-02-17_183940_DwC-A/occurrences_cleaned.csv")

glimpse(torch)
table(torch$year) 
hist(torch$coordinateUncertaintyInMeters) # hoo boy

torch.val <- torch %>% filter(year>=1900, coordinateUncertaintyInMeters<=4000) |> 
 	mutate(quarter = if_else(month < 4, 1, if_else(month < 7, 2, if_else(month < 10, 3, 4)))) |> 
	dplyr::select(year, quarter, location, lat, lon, obs_by, flr)

glimpse(torch.val) # okay nice
table(torch.val$year, torch.val$quarter)
table(torch.val$flr)


#-------------------------------------------------------------------------
# Merge validation sources

all.val <- rbind(npn.val, cch.val, torch.val) 

glimpse(all.val)
table(all.val$obs_by, all.val$year)
table(all.val$flr, all.val$year)
table(all.val$flr, all.val$obs_by)

write.table(all.val, "data/validation_cleaned.csv", sep=",", col.names=TRUE, row.names=FALSE)
# all.val <- read.csv("data/validation_cleaned.csv")

#-------------------------------------------------------------------------
# Link validation records to predictions

# read the model back in
flr.mod.all.ri <- read_rds(paste("output/models/RIbart.model.all.", taxnum, ".rds", sep=""))
allX <- attr(flr.mod.all.ri$fit[[1]]$data@x, "term.labels")

summary(flr.mod.all.ri) #

true.vector <- flr.mod.all.ri$fit[[1]]$data@y 
  
pred <- prediction(colMeans(pnorm(flr.mod.all.ri$yhat.train)), true.vector)
  
perf.tss <- performance(pred,"sens","spec")
tss.list <- (perf.tss@x.values[[1]] + perf.tss@y.values[[1]] - 1)
tss.df <- data.frame(alpha=perf.tss@alpha.values[[1]],tss=tss.list)
  
cutoff <- min(tss.df$alpha[which(tss.df$tss==max(tss.df$tss))])
cutoff # confirm this works

glimpse(all.val) # confirm this is in memory

valid <- NULL # initialize the data object, in a lazy way

valyears <- sort(unique(filter(all.val, year>=1900, year<2026)$year))

# LOOP over years and quarters ...
for(yr in valyears){

# yr <- 1999

pred <- rast(paste0("output/models/RIBART_historic_predictions/RIBART_predicted_flowering_", yr, ".tiff"))


for(q in unique(filter(all.val, year==yr)$quarter)){

yqsub <- filter(all.val, year==yr, quarter==q)

yqsub$pred_PrFlr <- terra::extract(pred[[q]], yqsub[,c("lon", "lat")])[,2]
yqsub$pred_Flr <- yqsub$pred_PrFlr >= cutoff

valid <- rbind(valid, yqsub)

} # end loop over quarters

} # end loop over years

glimpse(valid) 
table(valid$obs_by, useNA="ifany")
head(valid) # cool
tail(valid)

write.table(valid, "output/validation_cleaned_predictions.csv", sep=",", col.names=TRUE, row.names=FALSE)

# valid <- read.csv("output/validation_cleaned_predictions.csv")

# full data set, all quarters -----
performance(prediction(valid$pred_PrFlr, valid$flr), "auc")@y.values[[1]] # 0.81 overall

t.test(pred_PrFlr~flr, data=valid) # p < 2.2e-16 okay!

# for .. how many? ... Q3 and Q4 observations ---
valid |> filter(quarter>2) # quite a lot for out of flowering season
valid |> filter(quarter>2, flr) # huh okay THREE, but the *model* says these aren't flowering

performance(prediction(filter(valid, quarter>2)$pred_PrFlr, filter(valid, quarter>2)$flr), "auc")@y.values[[1]] # 0.59

t.test(pred_PrFlr~flr, data=filter(valid, quarter>2)) # eek?

# for Q4 only?
valid |> filter(quarter==4) # and even just for Q4
valid |> filter(quarter==4, flr) # and even just for Q4

performance(prediction(filter(valid, quarter==4)$pred_PrFlr, filter(valid, quarter==4)$flr), "auc")@y.values[[1]] # 0.44

t.test(pred_PrFlr~flr, data=filter(valid, quarter==4)) # p = 0.002 but reverse-oriented, LOL

# Let's write these out for notes and have a look ...
valid |> filter(quarter>2, flr) # okay cleaned some up with closer inspection, didn't get everything



#-------------------------------------------------------------------------
# newspaper records

news <- read.csv("data/news_accounts/news_reports.csv")
glimpse(news)
table(news$year)

# For a first attempt, let's just look at range-wide average ...

jotr.flr <- read.csv("output/historic_flowering_reconst_jotr.csv")
glimpse(jotr.flr)

# What we really need, though, is to use what geographic information we have
jotr.histStack <- raster::stack("output/BART/jotr_BART_predicted_flowering_1900-2023.grd")

# relevant polygons
JTNP <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_parks_and_protected_lands_scale_rank") %>% filter(unit_code=="JOTR") # the national park, when that's specified
Lanc <- read_sf(dsn = "../data/spatial/CA_Counties/", lay= "CA_Counties_TIGER2016") %>% filter(NAME=="Los Angeles") # LA County as a proxy for "Lancaster" or "Cajon" or similar
Kern <- read_sf(dsn = "../data/spatial/CA_Counties/", lay= "CA_Counties_TIGER2016") %>% filter(NAME=="Kern") # and Kern
Riverside <- read_sf(dsn = "../data/spatial/CA_Counties/", lay= "CA_Counties_TIGER2016") %>% filter(NAME=="Riverside") # and Kern

# then mask
JTNP.histFlr <- mask(jotr.histStack, st_transform(JTNP, crs=4269), touches=TRUE) # go generous
Lanc.histFlr <- mask(jotr.histStack, st_transform(Lanc, crs=4269), touches=TRUE) # go generous
Kern.histFlr <- mask(jotr.histStack, st_transform(Kern, crs=4269), touches=TRUE) # go generous
Riverside.histFlr <- mask(jotr.histStack, st_transform(Kern, crs=4269), touches=TRUE) # go generous

# then this gets complex
news.prFLs <- news %>% mutate(prFL.loc=NA, prFL.all=NA)

for(i in 1:nrow(news.prFLs)){

if(news.prFLs$location[i]=="JTNP") loc.yr <- JTNP.histFlr[[paste("prFL", news.prFLs$year[i], sep=".")]]
if(news.prFLs$location[i]=="Lancaster") loc.yr <- Lanc.histFlr[[paste("prFL", news.prFLs$year[i], sep=".")]]
if(news.prFLs$location[i]=="Kern County") loc.yr <- Kern.histFlr[[paste("prFL", news.prFLs$year[i], sep=".")]]
if(news.prFLs$location[i]=="Riverside County") loc.yr <- Kern.histFlr[[paste("prFL", news.prFLs$year[i], sep=".")]]

news.prFLs$prFL.all[i] <- cellStats(jotr.histStack[[paste("prFL", news.prFLs$year[i], sep=".")]], "mean")
news.prFLs$prFL.loc[i] <- cellStats(loc.yr, "mean")

}

glimpse(news.prFLs)

# (guh, the quotes bork this every way I try it)
write.table(news.prFLs, "output/news_reports_prFL.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# news.prFLs <- read.csv("output/news_reports_prFL.tsv", sep="\t")

t.test(prFL.loc~flr, data=news.prFLs) # p = 0.0004, YES
t.test(prFL.all~flr, data=news.prFLs) # p = 0.003


news.prFLs %>% group_by(flr) %>% summarize(mnPrLoc = mean(prFL.loc), sePrLoc = sd(prFL.loc)/sqrt(length(prFL.loc)), mnPrAll = mean(prFL.all), sePrAll = sd(prFL.all)/sqrt(length(prFL.all)))



ggplot() + geom_histogram(data=prFL.sum, aes(x=mn.prFL), fill="gray60") + geom_vline(data=news.prFLs, aes(xintercept=prFL.loc, color=flr)) + geom_vline(xintercept=0.26, linetype=2)


# visualize news reports and other validation sources ................
valid_srcs <- valid |> dplyr::select(obs_by, jotr_prFlr, flr) |> rename(fitted=jotr_prFlr, observed=flr) |> mutate(observed=as.numeric(observed), type="Validation records", classified=fitted>0.2635843)
glimpse(valid_srcs)

accu <- rbind(obs_valid, exp_valid) |> mutate(type=factor(type, c("Training data", "Validation records")), classified=fitted>0.2635843)

ggplot(valid_srcs, aes(y=as.factor(observed), x=fitted)) + geom_vline(xintercept=0.26, linewidth=0.5, color="gray60") + geom_jitter(height=0.15, aes(fill=classified, color=classified, shape=classified, size=classified, alpha=classified)) + geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + 

facet_wrap("obs_by", ncol=1) + 

annotate(geom="text", x=0.27, y=1.45, label="Classification cutoff", hjust=0, color="gray40", size=3) +

scale_y_discrete(labels=c("False", "True")) + labs(y="Flowers observed", x="Predicted Pr(Flowers)") +

scale_fill_manual(values=c("#ffffcc", NA)) +
scale_color_manual(values=c("gray60", "#253494")) +
scale_shape_manual(values=c(21, 19)) + 
scale_size_manual(values=c(1, 1.2)) + 
scale_alpha_manual(values=c(0.75, 0.5)) + 

theme_minimal(base_size=10) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))


{cairo_pdf("output/figures/FigXX_validation_sources.pdf", width=8, height=6)

ggdraw() + draw_plot(mod_varimp, 0, 0.53, 0.45, 0.44) + draw_plot(mod_accuracy, 0, 0, 0.45, 0.55) + draw_plot(partplot, 0.45, 0, 0.55, 1) + draw_plot_label(label=c("A", "B", "C"), x=c(0, 0, 0.45), y=c(1,0.55,1))

}
dev.off()



#-------------------------------------------------------------------------
# visualize accuracy with original data and validations 

jotr.mod <- read_rds("output/BART/bart.model.Jotr.rds")
names(jotr.mod)

obs_valid <- summary(jotr.mod)$data |> dplyr::select(fitted, observed) |> mutate(type="Training data", classified=fitted>0.2635843)
glimpse(obs_valid) 


# model accuracy figure .........................
mod_accuracy <- ggplot(filter(obs_valid, type=="Training data"), aes(y=as.factor(observed), x=fitted)) + 
	geom_vline(xintercept=0.26, linewidth=0.5, color="gray60") + geom_jitter(height=0.15, aes(fill=as.factor(observed), color=as.factor(observed), shape=as.factor(observed), size=as.factor(observed), alpha=as.factor(observed))) + 
	geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + 

	# facet_wrap("type", nrow=2) + 

	annotate(geom="text", x=0.27, y=1.45, label="Classification cutoff", hjust=0, color="gray40", size=4) +

	scale_y_discrete(labels=c("False", "True")) + labs(y="Flowers observed", x="Predicted Pr(flowers)") +

	scale_fill_manual(values=c("#ffffcc", NA)) +
	scale_color_manual(values=c("gray60", "#253494")) +
	scale_shape_manual(values=c(21, 19)) + 
	scale_size_manual(values=c(1, 1.2)) + 
	scale_alpha_manual(values=c(0.75, 0.5)) + 

	theme_minimal(base_size=12) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))

mod_accuracy

# and then also the predictor selection figure ............
jotr.varimp <- read_rds("output/BART/bart.varimp.Jotr.rds")

jotr.varimp$data <- jotr.varimp$data |> mutate(trees = factor(trees, c(10,20,50,100,150,200)))

levels(jotr.varimp$data$names) <- c("Delta[Y1-2]*PPT", "Delta[Y0-1]*PPT", "Max*VPD[Y0]", "Delta[Y0-1]*Min*VPD", "Min*Temp[Y0]", "Delta[Y0-1]*Max*Temp", "PPT[Y0]", "PPT[Y1]", "Min*VPD[Y0]", "Delta[Y0-1]*Max*VPD", "Max*Temp[Y0]", "PPT[Y2]", "Delta[Y0-1]*Min*Temp")

jotr.varimp$labels$group <- "Trees"
jotr.varimp$labels$colour <- "Trees"

label_parse <- function(breaks){ parse(text=breaks) } # need this, for reasons

mod_varimp <- ggplot(jotr.varimp$data, aes(x=names, y=imp, color=trees, group=trees)) + 
	geom_line(linewidth=0.75) + geom_point(size=1.5, color="gray30") + 
	scale_color_manual(values=c('#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#0c2c84'), name="Trees") + 
	
	labs(y="Relative contribution") + scale_x_discrete(label=label_parse) + 
	
	theme_minimal() + theme(legend.position=c(0.8, 0.75), axis.text=element_text(size=9), axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank(), legend.text=element_text(size=8), legend.title=element_text(size=8), legend.key.size=unit(0.15, "in"))  # okay nice

mod_varimp

# and finally the partial effects ...............
p <- read_rds("output/BART/bart.model.Jotr.partials.rds")

partvals <- rbind(
				data.frame(predictor="Delta[Y1-2]*PPT~(mm)", p[[1]]$data),
				data.frame(predictor="Delta[Y0-1]*PPT~(mm)", p[[2]]$data),
				data.frame(predictor="Max*VPD[Y0]~(hPa)", p[[3]]$data),
				data.frame(predictor="Delta[Y0-1]*Min*VPD~(hPa)", p[[4]]$data),
				data.frame(predictor="Min*Temp[Y0]~(degree*C)", p[[5]]$data),
				data.frame(predictor="Delta[Y0-1]*Max*Temp~(degree*C)", p[[6]]$data)
				) |> mutate(predictor=factor(predictor, c("Delta[Y1-2]*PPT~(mm)", "Delta[Y0-1]*PPT~(mm)", "Max*VPD[Y0]~(hPa)", "Delta[Y0-1]*Min*VPD~(hPa)", "Min*Temp[Y0]~(degree*C)", "Delta[Y0-1]*Max*Temp~(degree*C)")))

partplot <- ggplot(partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, labeller="label_parsed", scale="free") + labs(y="Marginal Pr(Flowers)") + theme_minimal() + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

partplot

# put it all together
{cairo_pdf("output/figures/Fig02_predictors_accuracy_partials.pdf", width=8, height=6)

ggdraw() + draw_plot(mod_varimp, 0, 0.4, 0.45, 0.6) + draw_plot(mod_accuracy, 0, 0, 0.45, 0.4) + draw_plot(partplot, 0.45, 0, 0.55, 1) + draw_plot_label(label=c("A", "C", "B"), x=c(0, 0, 0.45), y=c(1,0.42,1))

}
dev.off()


#-------------------------------------------------------------------------
# NEW validation viz

form_valid <- valid |> dplyr::select(jotr_prFlr, flr) |> rename(fitted=jotr_prFlr, observed=flr) |> mutate(observed=as.numeric(observed)) |> mutate(classified=fitted>0.2635843)
glimpse(form_valid)

valid$src <- valid$obs_by
valid$src[valid$obs_by%in%c("CCH","TORCH")] <- "Herbaria"
valid$src[valid$obs_by%in%c("CIS","Ray Yeager")] <- "Notes"


# accuracy scatterplot for formal records
form_acc <- ggplot(valid, aes(y=flr, x=jotr_prFlr)) + 
	geom_vline(xintercept=0.26, linewidth=0.5, color="gray60") + geom_jitter(height=0.2, aes(color=src), alpha=0.75, size=1) + 
	geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + 

	annotate(geom="text", x=0.27, y=1.45, label="Classification cutoff", hjust=0, color="gray40", size=3.5) +

	scale_y_discrete(labels=c("False", "True")) + labs(y="Flowering recorded", x="Predicted Pr(flowers)", title="Formal records") +

	scale_color_manual(values=c('#b2df8a','#33a02c','#a6cee3','#1f78b4')[4:1]) +

	theme_minimal(base_size=12) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))

form_acc

# accuracy scatterplot for informal records
news_valid <- news.prFLs %>% dplyr::select(-qt) %>% rename(fitted=prFL.loc, observed=flr) %>% mutate(classified=fitted>0.2635843)
glimpse(news_valid)

news_acc <- ggplot(news_valid, aes(y=as.factor(observed), x=fitted)) + 
	geom_vline(xintercept=0.26, linewidth=0.5, color="gray60") + geom_jitter(height=0.15, color='#fb9a99', size=2) + 
	geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + 
 
 	annotate(geom="text", x=0.27, y=1.45, label="Classification cutoff", hjust=0, color="gray40", size=3.5) +

	scale_y_discrete(labels=c("Poor", "Intense")) + labs(y="Flowering described", x="Mean Pr(flowers) for locale", title="Newspaper accounts") +

	theme_minimal(base_size=12) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))
	
news_acc

# and now a timeline
table(valid$obs_by)
table(valid$year)

val.ln <- table(valid$year, valid$src) %>% as.data.frame() %>% rename(year=Var1, src=Var2, records=Freq) %>% mutate(year=as.numeric(as.character(year)), src=as.character(src))
glimpse(val.ln)

val_counts <- rbind(val.ln, data.frame(year=as.numeric(names(table(news_valid$year))), records=c(table(news_valid$year)), src="Newspapers") ) %>% mutate(src=factor(src, c("Herbaria", "Notes", "NPN", "StClair&Hoines", "Newspapers")))

glimpse(val_counts)

# bar plot over the years .............
val_obs <- ggplot(val_counts, aes(x=year, y=records, fill=src)) + 
	geom_bar(stat="identity", position="stack") + 
	
	labs(x="Year of observation", y="Records") +
	
	scale_fill_manual(values=c('#fb9a99','#b2df8a','#33a02c','#a6cee3','#1f78b4')[5:1], labels=c("Herbarium collections", "Field notes", "USA-NPN data", "St. Clair & Hoines (2018)", "Newspapers"), name="Validation record source") + 

	theme_minimal(base_size=14) + theme(legend.position="inside", legend.position.inside=c(0.3, 0.6), legend.key.size=unit(0.15,"in"), plot.margin=unit(c(0.1,0.1,0.05,0.1),"in"), axis.text.x = element_text(angle=18), legend.background=element_rect(fill="white", color=NA))

val_obs

# and piece it together
{cairo_pdf(file="output/figures/Fig04_validation.pdf", width=6.5, height=4.5)

ggdraw() + draw_plot(val_obs, 0, 0.45, 1, 0.55) + draw_plot(form_acc, 0.05, 0, 0.45, 0.45) + draw_plot(news_acc, 0.55, 0, 0.45, 0.45) + draw_plot_label(c("A", "B", "C"), x=c(0, 0, 0.5), y=c(1, 0.45, 0.45))

}
dev.off()

# and piece it together, alt format

{cairo_pdf(file="output/figures/tall_validation.pdf", width=6.5, height=6)

ggdraw() + draw_plot(val_obs, 0, 0.45, 1, 0.55) + draw_plot(form_acc, 0.05, 0, 0.45, 0.45) + draw_plot(news_acc, 0.55, 0, 0.45, 0.45) + draw_plot_label(c("A", "B", "C"), x=c(0, 0, 0.5), y=c(1, 0.45, 0.45))

}
dev.off()



# figure to do the above but with the two species separately -------------

jotr.preds <- jotr.preds <- c("pptY1Y2", "pptY0Y1", "vpdmaxY0", "vpdminY0Y1", "tminY0", "tmaxY0Y1")

# first YUJA 
yuja.mod <- read_rds("output/BART/bart.model.yuja.rds")
names(yuja.mod)

yuja_obs_valid <- summary(yuja.mod)$data |> dplyr::select(fitted, observed) |> mutate(type="Training data", classified=fitted>0.1589413)
glimpse(yuja_obs_valid) 


# model accuracy figure
yuja_accuracy <- ggplot(yuja_obs_valid, aes(y=as.factor(observed), x=fitted)) + 
	geom_vline(xintercept=0.16, linewidth=0.5, color="gray60") + 
	geom_jitter(height=0.15, aes(fill=as.factor(observed), color=as.factor(observed), shape=as.factor(observed), size=as.factor(observed), alpha=as.factor(observed))) + 
	geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + 
	
	annotate(geom="text", x=0.17, y=1.45, label="Classification cutoff", hjust=0, color="gray40", size=3) +

	scale_y_discrete(labels=c("False", "True")) + labs(y="Flowers observed", x="Predicted Pr(Flowers)") +

	scale_fill_manual(values=c("#ffffcc", NA)) +
	scale_color_manual(values=c("gray60", "#253494")) +
	scale_shape_manual(values=c(21, 19)) + 
	scale_size_manual(values=c(1, 1.2)) + 
	scale_alpha_manual(values=c(0.75, 0.5)) + 

	labs(title="Model accuracy, YUJA") + theme_minimal(base_size=10) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))

yuja_accuracy

yuja_p <- partial(yuja.mod, jotr.preds, trace=FALSE, smooth=5)

yuja_partvals <- rbind(
				data.frame(predictor="Delta[Y1-2]*PPT~(mm)", yuja_p[[1]]$data),
				data.frame(predictor="Delta[Y0-1]*PPT~(mm)", yuja_p[[2]]$data),
				data.frame(predictor="Max*VPD[Y0]~(hPa)", yuja_p[[3]]$data),
				data.frame(predictor="Delta[Y0-1]*Min*VPD~(hPa)", yuja_p[[4]]$data),
				data.frame(predictor="Min*Temp[Y0]~(degree*C)", yuja_p[[5]]$data),
				data.frame(predictor="Delta[Y0-1]*Max*Temp~(degree*C)", yuja_p[[6]]$data)
				) |> mutate(predictor=factor(predictor, c("Delta[Y1-2]*PPT~(mm)", "Delta[Y0-1]*PPT~(mm)", "Max*VPD[Y0]~(hPa)", "Delta[Y0-1]*Min*VPD~(hPa)", "Min*Temp[Y0]~(degree*C)", "Delta[Y0-1]*Max*Temp~(degree*C)")))

yuja_partplot <- ggplot(yuja_partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, labeller="label_parsed", scale="free") + labs(y="Marginal Pr(Flowers)", title="Partial effects, YUJA") + theme_minimal() + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

yuja_partplot


# now YUBR
yubr.mod <- read_rds("output/BART/bart.model.yubr.rds")
names(yubr.mod)

yubr_obs_valid <- summary(yubr.mod)$data |> dplyr::select(fitted, observed) |> mutate(type="Training data", classified=fitted>0.2776204)
glimpse(yubr_obs_valid) 


# model accuracy figure
yubr_accuracy <- ggplot(yubr_obs_valid, aes(y=as.factor(observed), x=fitted)) + 
	geom_vline(xintercept=0.28, linewidth=0.5, color="gray60") + 
	geom_jitter(height=0.15, aes(fill=as.factor(observed), color=as.factor(observed), shape=as.factor(observed), size=as.factor(observed), alpha=as.factor(observed))) + 
	geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + 
	
	annotate(geom="text", x=0.29, y=1.45, label="Classification cutoff", hjust=0, color="gray40", size=3) +

	scale_y_discrete(labels=c("False", "True")) + labs(y="Flowers observed", x="Predicted Pr(Flowers)") +

	scale_fill_manual(values=c("#ffffcc", NA)) +
	scale_color_manual(values=c("gray60", "#253494")) +
	scale_shape_manual(values=c(21, 19)) + 
	scale_size_manual(values=c(1, 1.2)) + 
	scale_alpha_manual(values=c(0.75, 0.5)) + 

	labs(title="Model accuracy, YUBR") + theme_minimal(base_size=10) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))

yubr_accuracy

yubr_p <- partial(yubr.mod, jotr.preds, trace=FALSE, smooth=5)

yubr_partvals <- rbind(
				data.frame(predictor="Delta[Y1-2]*PPT~(mm)", yubr_p[[1]]$data),
				data.frame(predictor="Delta[Y0-1]*PPT~(mm)", yubr_p[[2]]$data),
				data.frame(predictor="Max*VPD[Y0]~(hPa)", yubr_p[[3]]$data),
				data.frame(predictor="Delta[Y0-1]*Min*VPD~(hPa)", yubr_p[[4]]$data),
				data.frame(predictor="Min*Temp[Y0]~(degree*C)", yubr_p[[5]]$data),
				data.frame(predictor="Delta[Y0-1]*Max*Temp~(degree*C)", yubr_p[[6]]$data)
				) |> mutate(predictor=factor(predictor, c("Delta[Y1-2]*PPT~(mm)", "Delta[Y0-1]*PPT~(mm)", "Max*VPD[Y0]~(hPa)", "Delta[Y0-1]*Min*VPD~(hPa)", "Min*Temp[Y0]~(degree*C)", "Delta[Y0-1]*Max*Temp~(degree*C)")))

yubr_partplot <- ggplot(yubr_partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, labeller="label_parsed", scale="free") + labs(y="Marginal Pr(Flowers)", title="Partial effects, YUBR") + theme_minimal() + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

yubr_partplot


{cairo_pdf("output/figures/SIFig_subspecies_accuracy_partials.pdf", width=6.5, height=6.5)

ggdraw() + draw_plot(yubr_partplot, 0, 0.3, 0.5, 0.7) + draw_plot(yuja_partplot, 0.5, 0.3, 0.5, 0.7) + draw_plot(yubr_accuracy, 0, 0, 0.5, 0.3) + draw_plot(yuja_accuracy, 0.5, 0, 0.5, 0.3) + draw_plot_label(label=c("A", "B", "C", "D"), x=c(0, 0, 0.5, 0.5), y=c(1, 0.3, 1, 0.3))

}
dev.off()

