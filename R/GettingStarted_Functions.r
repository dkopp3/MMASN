#' StreamCat Region Download
#'
#' Downloads zip files for VPU from StreamCat webpage
#' \url{https://www.epa.gov/national-aquatic-resource-surveys/streamcat}
#'
#' @param path Directory for downloaded files
#' @param vpu NHDPlusV2 Vector Processing Unit
#' @param filename name of streamcat file to download
#'
#' @return StreamCat data files are downloaded to "StreamCat" directory
#'
#'
#' @export

StreamCat_download <- function (path = getwd(),
                                vpu = "17",
                                filename = "NRSA_PredictedBioCondition"){

  if(!is.character(vpu)){
    stop("vpu must be character")
  }

  #list all file names on StreamCat ftp directory
  ftpPath <- "ftp://newftp.epa.gov/EPADataCommons/ORD/NHDPlusLandscapeAttributes/StreamCat/HydroRegions/"
  url = ftpPath
  filenames = RCurl::getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  filenames <- strsplit(filenames, "\r\n")
  filenames = unlist(filenames)

  if(!is.null(filename)){
    filenames <- grep(paste0(filename), filenames, value = T)
    files.download <- grep(paste0("Region", vpu), filenames, value = T)
  } else {
    #paste file extensions
    files.download <- grep(paste0("Region", vpu), filenames, value = T)
  }

  print(paste("downloading", length(files.download),
              "StreamCat files for Region", vpu, sep = " "))

  for (i in files.download){
    utils::download.file(paste(ftpPath, i, sep= ""),
                         paste(path, i, sep="/"),
                         mode = "wb")
    utils::unzip(paste(path, i, sep = "/"), exdir = paste0(path,"/Region",vpu))
  }
}




#' @title Find NHDPlusV2 COMID's
#'
#' Identifies NHDPlusV2 COMID's within a search radius
#'
#' @param pts Sampling points
#' @param NHDFlowline Flowline coverage from NHDPlusV2
#' @param maxdist search radius around points (m)
#'
#' @return list of length points with elements of COMID
#' within search radius
#'
#' @export
#'

find_comid <- function(pts, NHDFlowline, maxdist = 200){

  #prepare sample and flowlines
  #####
  sample_points <- sf::st_transform(pts, crs = 5070)
  NHDFlowline <- sf::st_transform(NHDFlowline, crs = (5070))

  geom <- sf::st_geometry(NHDFlowline)
  sf::st_geometry(NHDFlowline) <- NULL
  names(NHDFlowline) <- toupper(names(NHDFlowline))
  NHDFlowline <- sf::st_sf(NHDFlowline, geom = geom)
  ###################

  #find flowlines close to sample point
  #####
  rad <- sf::st_buffer(sample_points, dist = maxdist)
  z <- sf::st_intersects(rad, NHDFlowline)

  NoFlolines <- unlist(lapply(z, function(x) length(x)==0))

  if(any(NoFlolines)){
    warning("flowlines not found for some points, check index")
  }

  COMIDs <- sapply(which(!NoFlolines), simplify = F,function(x)
    data.frame(index = x,
               COMID = NHDFlowline$COMID[z[[x]]],
               dist_m = t(sf::st_distance(sample_points[x, ], NHDFlowline[z[[x]], ]))))

  return(COMIDs)
}
