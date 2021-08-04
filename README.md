
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MMASN: Multi-scale Morphometric Analysis of Stream Networks

<!-- badges: start -->

<!-- badges: end -->

The goal of MMASN is to facilitate the measurement of morphometric
(i.e. geometric) characteristics and spatial composition and
configuration of habitat patches within stream networks

## Installation

The development version can be installed from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dkopp3/MMASN")
```

## Example

In this example I’ll the demonstrate the functionality of MMASN by 1)
Downloading Data from the National Hydrography Dataset Plus Version 2
and ancillary data sources, and 2) measuring Watershed Scale, COMID
Scale and Reach Scale Morphometry. Then using these values, I’ll use a
clustering analysis to identify habitat patches with in stream networks
baised on their morphometry. Finally, I’ll quantify the spatial pattern
(Composition & Configuration) of these habitat pataches with landscape
metrics.

``` r
library(MMASN)
```

``` r
#Other packages for the example
library(archive) #7zip software workaround 
library(foreign) #reads .dbf files
library(sf) #workhorse
library(lwgeom) #line stuff 
library(raster) #external data
library(cluster) #patch identification 

library(tidyverse) #landscape metrics #tidyverse masks tidygraph, load first 
library(igraph) #landscape metrics
library(tidygraph) #landscape metrics

library(ggplot2)#lazy plots :) 
library(gridExtra)
```

# Data Preparation

## Download data

``` r
#use NHDPlusV2 data from VPU 17 
#https://www.epa.gov/waterdata/nhdplus-pacific-northwest-data-vector-processing-unit-17

#Flowlines - hydropgraphy 
download.file("https://s3.amazonaws.com/edap-nhdplus/NHDPlusV21/Data/NHDPlusPN/NHDPlusV21_PN_17_NHDSnapshot_08.7z",
               "Data/NHDPLUSV2/NHDPlusV21_PN_17_NHDSnapshot_08.7z", mode = "wb")
archive_extract(archive("Data/NHDPLUSV2/NHDPlusV21_PN_17_NHDSnapshot_08.7z"), dir = "Data/NHDPLUSV2")

#Value added attribute table 
#PlusFlow table - navigation
download.file("https://s3.amazonaws.com/edap-nhdplus/NHDPlusV21/Data/NHDPlusPN/NHDPlusV21_PN_17_NHDPlusAttributes_10.7z",
              "Data/NHDPLUSV2/NHDPlusV21_PN_17_NHDPlusAttributes_10.7z", mode = "wb")
archive_extract(archive("Data/NHDPLUSV2/NHDPlusV21_PN_17_NHDPlusAttributes_10.7z"), dir = "Data/NHDPLUSV2")

#catchment files were unavailable at NHDPlusV2 page
#download.file(mode = "wb")
#archive_extract(archive(), dir = "Data/NHDPLUSV2")

#download StreamCat (filename = NULL, downloads all files)
#https://www.epa.gov/national-aquatic-resource-surveys/streamcat-metrics-and-definitions
StreamCat_download(path = "Data/StreamCat", vpu = "17", filename = "NRSA_PredictedBioCondition")
StreamCat_download(path = "Data/StreamCat", vpu = "17", filename = "Lithology")
StreamCat_download(path = "Data/StreamCat", vpu = "17", filename = "NLCD2006RipBuf100_Region")

#download NRSA site information 
download.file("https://www.epa.gov/sites/production/files/2021-04/nrsa_1819_site_information_-_data.csv",
              "Data/nrsa_1819_site_information_-_data.csv", mode = "wb")

#download Ancillary data - Accumulated Degree days  >32degF
url <- paste0("http://geoserver.usanpn.org/geoserver/gdd/wcs?service=WCS&",
              "version=2.0.1&",
              "request=GetCoverage&CoverageId=gdd:30yr_avg_agdd&",
              "subset=http://www.opengis.net/def/axis/OGC/0/elevation(", 365,")&",
              "format=image/geotiff")
download.file(url, destfile = paste0("Data/AGDD_30yrAVG_", 365) , method = "libcurl", mode = "wb")
```

## Read NHDPlusV2 Files

``` r
flow = foreign::read.dbf("Data/NHDPLUSV2/NHDPlusPN/NHDPlus17/NHDPlusAttributes/PlusFlow.dbf")
vaa = foreign::read.dbf("Data/NHDPLUSV2/NHDPlusPN/NHDPlus17/NHDPlusAttributes/PlusFlowlineVAA.dbf")
elevslope = foreign::read.dbf("Data/NHDPLUSV2/NHDPlusPN/NHDPlus17/NHDPlusAttributes/elevslope.dbf")
NHDFlowline = sf::read_sf("Data/NHDPLUSV2/NHDPlusPN/NHDPlus17/NHDSnapshot/Hydrography/NHDFlowline.shp")
NHDFlowline <- sf::st_zm(NHDFlowline, what = "ZM", drop = T)
NHDCatchments <- sf::read_sf("Data/NHDPLUSV2/NHDPlusPN/NHDPlus17/NHDPlusCatchment/Catchment.shp")
```

## Identify network outlets

``` r
outlets <- read.csv("Data/nrsa_1819_site_information_-_data.csv")
outlets <- unique(subset(outlets, select = c("UNIQUE_ID","COMID","LON_DD83","LAT_DD83")))
outlets <- sf::st_as_sf(outlets, coords = c("LON_DD83","LAT_DD83"), crs = 4269)

vpu_poly <- sf::read_sf("Data/VPU_NAD83.shp")
vpu_ls <- sf::st_contains(st_transform(vpu_poly, crs = 5070), st_transform(outlets, crs = 5070))
names(vpu_ls) <- vpu_poly$VPUID

#select outlets in vpu 17
outlets <- outlets[vpu_ls[["17"]],]

ggplot()+
geom_sf(data = vpu_poly, color = "blue") +
  geom_sf(data = outlets, color = "black")
```

\#find NHDPlus COMID’s near sampling point

``` r
#returns all nhdplus COMID within the maxdist of a buffer
comids <- find_comid(pts = outlets, NHDFlowline, maxdist = 200)

#select closest
comids <- do.call(rbind, lapply(comids, function(x) x[which.min(x$dist_m), ]))

outlets$COMID[comids$index] <- comids$COMID
head(comids)
```

\#Delineate network upstream of a comid

``` r
upstr_comid_1 <- net_delin(comid = 24423157, flow = flow, vaa = vaa)
upstr_comid_2 <- net_delin(comid = 947100116, flow = flow, vaa = vaa)
```

\#quick plot of Stream Networks

``` r
#select outlets
sites_1 <- outlets[outlets$COMID == 24423157, ] 
sites_2 <- outlets[outlets$COMID == 947100116, ] 

#extract flowlines from NHDPlus
network_1 <- NHDFlowline[NHDFlowline$COMID%in%upstr_comid_1, ]
network_2 <- NHDFlowline[NHDFlowline$COMID%in%upstr_comid_2, ]

#plot
p1 <- ggplot()+
  geom_sf(data = network_1, color = "blue") +
  geom_sf(data = sites_1, color = "black")

p2<-ggplot()+
  geom_sf(data = network_2, color = "blue")+
  geom_sf(data = sites_2, color = "black")
  
#install.packages("gridExtra")
gridExtra::grid.arrange(p1, p2, nrow = 1)
```

# Watershed Scale Morphometry

Watershed scale morphometry refers to the geometric characteristics of
the watershed (i.e. basin shape) or the stream network itself
(i.e. arrangements of the stream channels).

Basin Shape can be measured with a variety of indicied. Compactness
coefficient (CmpC) is defined as the ratio of the watershed perimeter to
the circumference of equivalent circular area. Elongation ratio (ElnR)
is defined as the ratio of diameter of a circle of the same area as the
watershed to the maximum watershed length. The numerical value varies
from 0 (in highly elongated shape) to 1 (in circular shape). Circulatory
ratio (CrcR) is defined as the ratio of watershed area to the area of
the circle having the same perimeter as the watershed perimeter. The
numeric value may vary in between 0 (in line) and 1 (in a circle). see
<http://www.jnkvv.org/PDF/04042020192203Geomorphology%20of%20Watershed.pdf>

``` r
#calculate basin shape metrics
cats_1 <- cat_shp(network = network_1, NHDCatchments)
cats_2 <- cat_shp(network = network_2, NHDCatchments)
cats_1
```

the shape of the stream network involves measuring the arrangement of
the channels. Bifurcation Ratio (Rb) is defined as the ratio of number
of streams of a particular order (Nu) to number of streams of next
higher order (Nu+i). Length Ratio (Rl) is defined as the ratio of mean
stream length (Lu) of a particular stream order to mean stream length of
the next lower order (Lu-1). Area Ratio: It is defined as the ratio of
mean catchment area (Lu) of a particular stream order to mean catchment
area of the next lower order (Lu-1). For a stream network I estimate
these values from geometric relationship between the ratio and a given
order.

``` r
#calcualte stream network characteristics
netshp_1 <- net_shp(network = network_1, vaa = vaa)
netshp_2 <- net_shp(network = network_2, vaa = vaa)
netshp_1
```

\#view data

``` r

WS_SCALE <- data.frame(rbind(cbind(st_set_geometry(cats_1, NULL), netshp_1),
                  cbind(st_set_geometry(cats_2, NULL), netshp_2)), row.names = NULL)
head(WS_SCALE)
```

\#quick plot of watersheds

``` r
p1 <- ggplot()+
  geom_sf(data = cats_1, color = "green")+
  geom_sf(data = network_1, color = "blue") +
  geom_sf(data = sites_1, color = "black")

p2 <- ggplot()+
  geom_sf(data = cats_2, color = "green")+
  geom_sf(data = network_2, color = "blue")+
  geom_sf(data = sites_2, color = "black")

#install.packages("gridExtra")
gridExtra::grid.arrange(p1, p2, nrow = 1)
```
