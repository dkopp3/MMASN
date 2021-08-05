
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MMASN: Multi-scale Morphometric Analysis of Stream Networks

<!-- badges: start -->

<!-- badges: end -->

The goal of MMASN is to facilitate morphometric analysis of stream
networks and provide landscape metrics that quantify spatial pattern of
habitat patches in dendritic networks. These functions rely on data
provide by the [National Hydrography Dataset Plus Version 2;
NHDPlusv2](https://www.epa.gov/waterdata/get-nhdplus-national-hydrography-dataset-plus-data).

# Installation

This Package is still under development but can be installed from
[GitHub](https://github.com/):

``` r
# install.packages("devtools")
devtools::install_github("dkopp3/MMASN")
```

# Example

Here, I download data from the NHDPlusV2 and optional data sources. I
then use MMASN to measure watershed-, NHDPlus COMID- and reach- scale
morphometry and demonstrate how these values can be used to identify
habitat patches in stream networks with hierarchical clustering.
Finally, I quantify the spatial pattern of these habitat patches using
landscape metrics developed for dendritic stream networks.

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

NHDPlusv2 (Required, digital stream network and catchments)

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
```

StreamCat (optional, COMID Scale values)

``` r
#download StreamCat (filename = NULL, downloads all files)
#https://www.epa.gov/national-aquatic-resource-surveys/streamcat-metrics-and-definitions
StreamCat_download(path = "Data/StreamCat", vpu = "17", filename = "NRSA_PredictedBioCondition")
StreamCat_download(path = "Data/StreamCat", vpu = "17", filename = "Lithology")
StreamCat_download(path = "Data/StreamCat", vpu = "17", filename = "NLCD2006RipBuf100_Region")
```

Average accumulated growing degree days (Optional, Reach scale dataset
values)

``` r
#download Ancillary data - Accumulated Degree days >32degF
#https://www.usanpn.org/data/agdd_maps
url <- paste0("http://geoserver.usanpn.org/geoserver/gdd/wcs?service=WCS&",
              "version=2.0.1&",
              "request=GetCoverage&CoverageId=gdd:30yr_avg_agdd&",
              "subset=http://www.opengis.net/def/axis/OGC/0/elevation(", 365,")&",
              "format=image/geotiff")
download.file(url, destfile = paste0("Data/AGDD_30yrAVG_", 365) , method = "libcurl", mode = "wb")
```

Stream network outlets
(Required)

``` r
#download Sampling locations from the National Rivers and Streams Assessment
download.file("https://www.epa.gov/sites/production/files/2021-04/nrsa_1819_site_information_-_data.csv",
              "Data/nrsa_1819_site_information_-_data.csv", mode = "wb")
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

## Select stream network outlets

``` r
#read in file 
outlets <- read.csv("Data/nrsa_1819_site_information_-_data.csv")

#select unique values 
outlets <- unique(subset(outlets, 
                         select = c("UNIQUE_ID", "COMID",
                                    "LON_DD83","LAT_DD83")))
#create sf object
outlets <- sf::st_as_sf(outlets, 
                        coords = c("LON_DD83","LAT_DD83"), 
                        crs = 4269)
```

``` r
#read in NHDPlusV2 vector processing unit (vpu)
vpu_poly <- sf::read_sf("Data/VPU_NAD83.shp")

#associate outlets with vpu polugons
vpu_ls <- sf::st_contains(st_transform(vpu_poly, crs = 5070), 
                          st_transform(outlets, crs = 5070))
names(vpu_ls) <- vpu_poly$VPUID

#select outlets in vpu 17
outlets <- outlets[vpu_ls[["17"]],]
```

``` r
#quick plot
ggplot()+
geom_sf(data = vpu_poly, color = "blue") +
  geom_sf(data = outlets, color = "black")
```

## Find Nearest NHDPlus COMID

In many instances the sampling point may not completely align with the
NHDPlusV2 flowline. The find\_comid() function returns a list of length
outlets, with each element containing the COMID’s within a given buffer
distance.

``` r
#returns all nhdplus COMID within the maxdist buffer
comids <- find_comid(pts = outlets, NHDFlowline, maxdist = 200)

#select closest COMID for this ecample
comids <- do.call(rbind, lapply(comids, function(x) x[which.min(x$dist_m), ]))

#assign COMID value
outlets$COMID[comids$index] <- comids$COMID
head(comids)
```

## Delineate Upstream Network

net\_delin extracts all the COMID’s upstream of a given COMID. In this
example we delineate the stream network upstream of two comids

``` r
upstr_comid_1 <- net_delin(comid = 24423157, flow = flow, vaa = vaa)
upstr_comid_2 <- net_delin(comid = 947100116, flow = flow, vaa = vaa)
```

``` r
#quick plot of Stream Networks

#select outlets
sites_1 <- outlets[outlets$COMID == 24423157, ] 
sites_2 <- outlets[outlets$COMID == 947100116, ] 

#extract flowlines from NHDPlus via simple matching
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

Watershed scale morphometry refers to indices that quantify the shape of
the watershed (i.e. land area surrounding the entire stream network) or
the stream network itself.

## Basin Shape

Basin Shape includes: Compactness coefficient (CmpC), defined as the
ratio of the watershed perimeter to the circumference of equivalent
circular area.

Elongation ratio (ElnR), defined as the ratio of diameter of a circle of
the same area as the watershed to the maximum watershed length. The
numerical value varies from 0 (in highly elongated shape) to 1 (in
circular shape).

Circulatory ratio (CrcR), defined as the ratio of watershed area to the
area of the circle having the same perimeter as the watershed perimeter.
The numeric value may vary in between 0 (in line) and 1 (in a circle).
[see for more
information](http://www.jnkvv.org/PDF/04042020192203Geomorphology%20of%20Watershed.pdf)

``` r
#calculate basin shape metrics
cats_1 <- cat_shp(network = network_1, NHDCatchments)
cats_2 <- cat_shp(network = network_2, NHDCatchments)
cats_1
```

## Stream Network Shape

Stream network shape includes: Bifurcation Ratio (Rb) is defined as the
ratio of number of streams of a particular order (Nu) to number of
streams of next higher order (Nu+i).

Length Ratio (Rl) is defined as the ratio of mean stream length (Lu) of
a particular stream order to mean stream length of the next lower order
(Lu-1).

Area Ratio is defined as the ratio of mean catchment area (Lu) of a
particular stream order to mean catchment area of the next lower order
(Lu-1).

These values are estimated from the geometric relationship between the
ratio and a given order.

``` r
#calcualte stream network characteristics for two stream networks
netshp_1 <- net_shp(network = network_1, vaa = vaa)
netshp_2 <- net_shp(network = network_2, vaa = vaa)
netshp_1
```

``` r
# Quick plot of watersheds
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

``` r
# Combine watershed scale morphometry into dataframe for later use 
WS_SCALE <- data.frame(rbind(
  cbind(st_set_geometry(cats_1, NULL), netshp_1),
  cbind(st_set_geometry(cats_2, NULL), netshp_2)), 
  row.names = NULL)

head(WS_SCALE)
```

# COMID Scale Metrics

COMID scale metrics refer to attributes of individual stream channels

## Channel Shape

channel shape includes the mean elevation and slope of the COMID
provided by NHDPlusV2 and the sinuosity. Sinuosity is measured as the
length of the stream channel divided by the straight line distance
between the COMID inlet/outlet.

``` r
planform <- chn_geom(network = network_1, elevslope = elevslope)

#drop the NHDPlusV2 smoothed values
planform <- subset(planform, select = -c(MAXELEVSMO, MINELEVSMO))
head(planform)
```

## Confluence Attributes

Confluence attributes include: Tributary angle defined as the angle that
the tributaty and mainstem meet

Area ratio is the catchment area of the tributary divided by the
catchment area of the mainstem

Confluence area is the catchment area upstream of the confluence

Confluence order is the stream order directly downstream of the
confluence

Confluence class is the concatenation of the tributary and mainstem
stream orders

``` r
confl_attrs <- net_confl(network = network_1, vaa = vaa, flow = flow)
head(confl_attrs)
```

``` r
#quick plot of confluence
ggplot()+
  geom_sf(data = cats_1, color = "green")+
  geom_sf(data = network_1, color = "blue") +
  geom_sf(data = network_1[network_1$COMID %in%confl_attrs$COMID,], color = "black")+
  geom_sf(data = network_1[network_1$COMID %in% sites_1$COMID,], color = "purple", lwd = 1.25)+
  geom_sf(data = sites_1, color = "black") 

confl_attrs[confl_attrs$COMID==sites_1$COMID, ]
```

## Adding StreamCat Data

The get StreamCat function extracts values from a given StreamCat data
table for all COMIDs within a network

``` r
#get lithology values for each COMID in the network
lithology <- get_streamcat(comid = network_1$COMID, 
                           strcat_path = "Data/StreamCat",
                           csvfile = "Lithology", vpu = "17")

#clean output: select columns with values
#many lithology values are 0 for all COMIDs
lithology <- subset(lithology, 
                    select = c("COMID", grep("Cat", names(lithology), value = T)))
lithology <- lithology[,apply(lithology, 2, sum, na.rm = T)>0]
```

``` r
#get predicted biological condition
biocond <- get_streamcat(comid = network_1$COMID, 
                         strcat_path = "Data/StreamCat",
                         csvfile = "NRSA_PredictedBioCondition", 
                         vpu = "17")
biocond <- biocond[,c("COMID", "prG_BMMI")]
```

``` r
#clean COMID scale data
COMID_SCALE <- Reduce(function(x,y) merge(x, y, by = "COMID", all.x = T), 
                      list(planform, confl_attrs, lithology, biocond))

#some reaches are not directly downstream of confluence 
x<-COMID_SCALE[,names(confl_attrs)]
#suppress invalid factor level on tributary class 
suppressWarnings(x[is.na(x)]<-0)
#update levels 
levels(x$Confl_class)<-c(levels(x$Confl_class), "0.0")
x$Confl_class[is.na(x$Confl_class)] <- "0.0"
COMID_SCALE[,names(confl_attrs)]<-x
rm(x)
head(COMID_SCALE)
```

# Reach scale metrics

In some instances users may want to use data at a finer resolution than
COMID. We define reaches as divisions of NHPLUS Comids and provide a
function that splits the COMIDs within a stream network into equal
parts, identifies the point and creates lateral transects. These sf
objects (reaches, midpoints, and transects) can be used to extract new
data at finer resolutions.

## Creating Stream Network Reaches

``` r
#split NHDPlusV2 COMIDs into 2 reaches (i.e. n=2)
net_rch <- net_sample(network = network_1, vaa = vaa, n = 2, what = "reaches", keep.divergent = F)

#create point at midpoint of reaches
net_pts <- net_sample(network_1, vaa, n = 2, what = "midpoints", keep.divergent = F)

#create lateral transects at each midpoint
net_trn <- net_sample(network_1, vaa, n = 2, what = "transects", 
                      NHDCatchments = NHDCatchments, keep.divergent = F)
```

``` r
#quick plot of watersheds
p1 <- ggplot()+
  geom_sf(data = cats_1, color = "green")+
  geom_sf(data = net_rch, color = "blue") +
  geom_sf(data = net_pts, color = "black")

p2 <- ggplot()+
  geom_sf(data = cats_1, color = "green")+
  geom_sf(data = net_rch, color = "blue") +
  geom_sf(data = net_trn, color = "black")

#install.packages("gridExtra")
gridExtra::grid.arrange(p1, p2, nrow = 1)
```

## extract data at midpoints and along transects

``` r
#Use annual growing degree days from the National Phenology Network
r <- raster("Data/AGDD_30yrAVG_365")
agdd = raster::extract(r,
                       st_coordinates(
                         st_transform(net_pts, crs = crs(r))))

agdd <- data.frame(COMID = net_pts$COMID, 
                   PointID = net_pts$PointID,
                   agdd)

#Valley width is length of lateral transect
valley_width<-data.frame(COMID = net_trn$COMID,
                         PointID = net_trn$PointID, 
                         vln = units::set_units(st_length(net_trn),NULL))
```

``` r
#merge data
REACH_SCALE <- merge(agdd, valley_width, by = c("COMID", "PointID"), all.x = T)
head(REACH_SCALE)
```

# Idenfitying “Patches” within a Stream Network

## combine Watershed, COMID, and Reach data

``` r
#combine multiscale data
StrNet_data <- cbind(WS_SCALE[1,], 
                     merge(COMID_SCALE, 
                           REACH_SCALE, 
                           by = "COMID", 
                           all.x = T))

#select variables for clustering
StrNet_vars <- names(StrNet_data)[!names(StrNet_data) %in% 
                                    c("COMID", "PointID", 
                                      "prG_BMMI", "Class", 
                                      "CatPctFull","Rl_rsq",
                                      "Ra_rsq","Rb_rsq")]
#Use complete cases
StrNet_data <- StrNet_data[complete.cases(StrNet_data[,StrNet_vars]),]

head(StrNet_data[,c("COMID", "PointID", StrNet_vars)])
```

## Hierarchical Clustering Analysis

``` r
#calculate dissimiilarity 
d <- daisy(StrNet_data[,StrNet_vars], metric = "gower")

#Clustering using wards method   
cl <- agnes(d, method = "ward")
```

``` r
#view dendrogram
plot(cl, which.plots = 2)

#define 5 groups of reaches
StrNet_data$grps <- cutree(cl, 5)
```

``` r
#plot Cluster Results
#merge results with the sampled reaches 
#retain all.x = T because some reaches were missing 
#data remained unclassified
net <- merge(net_rch, 
             StrNet_data[,c("PointID", "COMID", "prG_BMMI", "grps")], 
             by = c("COMID", "PointID"), all.x = T)

#update sites that could not be included in classification because
#missing data to 0
net$grps[is.na(net$grps)] <- 0

p1 <- ggplot()+
  geom_sf(data = net, color = net$grps + 1)

gridExtra::grid.arrange(p1, nrow = 1)
```

# Quantifying Spatial Pattern

Spatial pattern can refer to composition of configuration of habitat
patches within the steam network. Below, each landscape metric is
calculated for each patch type (i.e. class scale).

## Composition

``` r
#merge adjacent habitat patches of the same type into continuous segment
net_HGP <- merge_lines(lines = st_transform(net, 5070), 
                       groupName = "grps")

#create example composition metrics network
#claculate the number and total and mean length,
#of a given patch type
HGP_cmp <- do.call(rbind,
                      lapply(split(net_HGP, net_HGP$grps), 
                             function(x) 
                               data.frame(patch = unique(x$grps), 
                                          num = nrow(x), 
                                          tlen = sum(st_length(x)),
                                          meanlen = mean(st_length(x)))))
HGP_cmp
```

## Configuration

``` r
#function sf_to_tidygraph adapted from #https://www.r-spatial.org/r/2019/09/26/spatial-networks.html 
#Enables measurement of watercourse distance between 
#habitat patches within the stream network
graph <- sf_to_tidygraph(st_transform(net, 5070), directed = F)
#here we calculate pairwise watercourse distances 
#separating patches of type 2
HGP_cnf <- patch_dist(graph, net_HGP, groupName = "grps", patch = "2")

#claculate a modified dendritic connectivity index baised on watercourse 
#distance between the patches 
DCI_dist(plen_i = HGP_cnf$patch_i_len, 
         plen_j = HGP_cnf$patch_j_len, 
         dist_ij = HGP_cnf$d, mu = 12000)
```
