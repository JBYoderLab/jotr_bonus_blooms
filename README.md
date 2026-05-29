README: Joshua tree "bonus blooms"
==================================

This repository contains code and data in support of the article

> Yoder JB, CJ Carlson, and CW Callahan. 2026. Limited role of climate change in out-of-season Joshua tree flowering. bioRxiv. [doi.org/10.1101/2026.05.27.727294](https://doi.org/10.1101/2026.05.27.727294)

From the article's abstract

> Last winter, community scientists recorded anomalous winter flowering by Joshua trees (*Yucca brevifolia* and
*Y. jaegeriana*), and some speculated the bloom was caused by climate change. We trained machine learning models that reliably identify weather triggers for seasonal Joshua tree flowering, then applied frontier methods from climate science to simulate flowering in a counterfactual world without human-caused climate change. Surprisingly, we found winter blooms in 2018–9 and 2025–6 were driven by high winter rainfall, not rising temperatures—and therefore, are probably the result of natural weather variability, not climate change.


Directory structure and contents
--------------------------------

These are the major contents and subdirectories of the project folder. 

- `scripts`: R code for downloading and managing records from [iNaturalist](https://www.inaturalist.org), linking them to weather data from the [PRISM](https://prism.oregonstate.edu) database, BART model training and analysis
- `data`: mostly cleaned, organized iNaturalist records and associated weather predictors; larger weather and deltas to simulate counterfactual weather data are not uploaded here
	- `validation`: data from [USA-NPN](https://www.usanpn.org/data/observational) and digital herbarium records ([CCH](https://www.cch2.org/portal/index.php) and [TORCH](https://portal.torcherbaria.org/portal/index.php)), used for validation of model predictions
- `output`: anything generated in the course of analysis
	- `figures`: figures specifically


Key dependencies
----------------

- [r-project.org](http://www.r-project.org) --- for all analyses
- [`rinat`](https://cran.r-project.org/web/packages/rinat/index.html) --- for downloading observations from the iNaturalist API
- [`prism`](https://cran.r-project.org/web/packages/prism/index.html) --- for downloading and managing PRISM data layers
- [`embarcadero`](https://github.com/cjcarlson/embarcadero) --- for training and interpreting BART models


External data
-------------

- `Jotr_SDM2023_range` --- folder with `jotr_SDM2023_range.shp` and `jotr_SDM2023_range_simple.shp`, shapefiles (and accessory files) that provide range map polygons for Joshua tree (eastern and western together) derived from the [Esque *et al.* (2023)](https://doi.org/10.3389/fevo.2023.1266892) data set.
- `PRISM` --- folder of raster-formatted monthly weather data from 1900 to the most recent full year, from the [PRISM](https://prism.oregonstate.edu) database
