#' Network Delineation
#'
#' Identifies flowlines upstream from a COMID
#'
#' @param comid COMID at network outlet
#' @param flow NHDPlusV2 flow table, PlusFlow.dbf
#' @param vaa NHDPlusV2 value added attributes table, PlusFlowlineVAA.dbf"

#' @return vector of upstream COMIDs
#'
#' @export

net_delin <- function(comid, flow = flow, vaa = vaa){

  names(flow) <- toupper(names(flow))
  names(vaa) <- toupper(names(vaa))

  fcomid <- comid
  VPUOUT <- vaa[vaa$COMID %in% fcomid, c("VPUOUT")]
  net <- data.frame(fcomid, VPUOUT)

  # delineate upstream until the full net is delineated or 5 wbs are found
  while (any(flow$TOCOMID %in% fcomid) & #the fcomid has somthing upstream
         all(!is.na(net$fcomid)) & #na might occur if missmatch
         all(net$VPUOUT == 0))# if network leaves the vpu stop the loop
  {
    #upstream comid's
    fcomid <- flow[flow[, "TOCOMID"] %in% fcomid, "FROMCOMID"] #upstream comids
    fcomid <- unique(fcomid[fcomid != 0])

    #vaa of upstream COMID
    VPUOUT <- vaa[vaa$COMID %in% fcomid, c("VPUOUT")]

    if(length(VPUOUT) > 0){
      #Identify whether lake network travels outside vpu
      if(any(VPUOUT == 1)){VPUOUT = 1} else {VPUOUT = 0}
    }

    temp <- data.frame(VPUOUT, fcomid)
    net <- rbind(temp[!temp$fcomid %in% net$fcomid, ], net)

    if (nrow(net) > 5e+05) {
      stop(print("There is a problem"))
    }
  }

  if(any(net$VPUOUT != 0)){
    warning("Stream Network Left VPU")
  }

  net <- net$fcomid

  return(net)
}



#' Basin Shape Metrics
#'
#' Calculates metrics the describe the shape of a networks waterhsed.
#'
#' @param network a single network extracted from NHDPlusV2
#' @param NHDCatchments NHDPlusV2 catchment coverage, Catchment.shp
#'
#' @return sf object of watershed attributed with
#'  width (\code{width}), maximum length (\code{mxln}),
#'  area (\code{area}), parimeter (\code{prmt}),
#'  compactedness coefficient (\code{CmpC}),
#'  circularoty ratio (\code{CrcR}), and elongation
#'  ratio (\code{ElnR}). Where applicable, units = m
#'
#' @details Compactness coefficient (CmpC) is defined as the ratio of the
#' watershed perimeter to the circumference of equivalent circular area.
#'
#' Elongation ratio (ElnR) is defined as the ratio of diameter of a circle of the same area as
#' the watershed to the maximum watershed length. The numerical value varies
#' from 0 (in highly elongated shape) to 1 (in circular shape).
#'
#' Circulatory ratio (CrcR) is defined as the ratio of watershed area to the area of the circle
#' having the same perimeter as the watershed perimeter. The numeric value
#' may vary in between 0 (in line) and 1 (in a circle).
#'
#' @export

cat_shp <- function(network, NHDCatchments){

  cats <- NHDCatchments[NHDCatchments$FEATUREID%in%network$COMID,]
  #project to albers

  cats <- sf::st_transform(cats, 5070)
  cats <- sf::st_cast(sf::st_union(cats), "POLYGON")

  box <- sf::st_bbox(cats)

  vert <- sf::st_linestring(as.matrix(rbind(c(box$xmin,box$ymin), c(box$xmin,box$ymax))))
  vert <- sf::st_length(vert)

  horiz <- sf::st_linestring(as.matrix(rbind(c(box$xmin,box$ymin),c(box$xmax,box$ymin))))
  horiz <- sf::st_length(horiz)

  width <- min(horiz, vert)
  mxln <- max(horiz, vert)
  area <- sf::st_area(cats)
  prmt <- lwgeom::st_perimeter(cats)

  #Perimeter of watershed/Perimeter of circle of watershed area
  #Circle Perimeter (P) = √(4πA)
  #circle area (A) = P^2 / 4π
  #circle diameter D = 2r; r=√(A/π)

  #Compactness coefficient is defined as the
  #ratio of the watershed perimeter to the
  #circumference of equivalent circular area.
  CmpC <- prmt/sqrt(4* base::pi *area)

  #Circulatory ratio is defined as the
  #ratio of watershed area to the area of
  #the circle having the same perimeter as the watershed perimeter
  CrcR <- area/(prmt^2 / 4 * base::pi)

  #Elongation ratio is defined as the
  #ratio of diameter of a circle of the same area as the watershed
  #to the maximum watershed length.
  #The numerical value varies from 0 (in highly elongated shape) to
  #1 (in circular shape). These values can be
  ElnR <- (2*sqrt(area / base::pi))/mxln

  cats <- sf::st_as_sf(data.frame(cbind(width = width/1000, mxln = mxln/1000,
                                        area = area/1e6, prmt = prmt/1000,
                                        CmpC, CrcR, ElnR), geometry = cats),
                       crs = 5070)

  #transform back into network crs
  cats <- sf::st_transform(cats, crs = sf::st_crs(network))

  return(cats)
}



#' Network Shape Metrics
#'
#' Calculates metrics the describe the shape of a network
#'
#' @param network a single network extracted from NHDPlusV2
#' @param vaa NHDPlusV2 value added attributes table, PlusFlowlineVAA.dbf"
#'
#' @return data.frame with network stream order (\code{Ord}),
#' Mainstem Length (\code{Ord}), total length (\code{TL}), Drainage
#' Density (\code{Dd}), Horton Laws: bifurcation ratio (\code{Rb}),
#' length Ratio (\code{Rl}) and area Ratio (\code{Ord})
#'
#' @details Horton ratios estimated as semilog relationships
#'
#' @export

net_shp<-function (network, vaa){

  net <- sf::st_set_geometry(network, NULL)
  names(net) <- toupper(names(net))
  names(vaa) <- toupper(names(vaa))


  #remove diveregences
  vaa <- vaa[vaa$STREAMORDE == vaa$STREAMCALC & vaa$COMID %in% net$COMID,]

  #split each path - more accurately represents stream reaches in calcuations
  #z<-dat[[1]][,c("COMID", "STREAMORDE", "LEVELPATHI")]
  #plot(st_geometry(network))
  #plot(st_geometry(network[network$COMID%in%z$COMID,]), col = "red", add = T)

  dat <- split(vaa, as.character(vaa$LEVELPATHI))

  #calculate Mainstem Length (while we are here)
  MSL <- max(unlist(lapply(dat, function(x) sum(x$LENGTHKM))))

  dat <- do.call(rbind, lapply(dat, function(z)
    do.call(rbind, lapply(split(z, z$STREAMORDE), function(x)
      data.frame(STREAMORDE = unique(x$STREAMORDE),
                 TOTDASQKM = max(x$TOTDASQKM),
                 LENGTHKM = sum(x$LENGTHKM))))))


  #summarise by stream order
  dat <- split(dat, as.character(dat$STREAMORDE))
  dat <- do.call(rbind, lapply(dat, function(x)
    data.frame(STREAMORDE = unique(x$STREAMORDE),
               N = nrow(x),logN = log(nrow(x)) ,
               A = mean(x$TOTDASQKM), logA = log(mean(x$TOTDASQKM)),
               L = mean(x$LENGTHKM), logL = log(mean(x$LENGTHKM)))))


  #y = log(c(60,13,9,4,1))
  #x = c(1:5)
  #1/exp(coef(lm(y~x))[2])

  # estimate log(Rb) as slope
  #plot(dat$STREAMORDE,dat$logN)
  #abline(lm(dat$logN~dat$STREAMORDE))
  #fit semilog model
  lmRb <- stats::lm(logN ~ STREAMORDE, data = dat)
  # slove for Rb w/base e
  Rb <- 1/exp(stats::coef(lmRb)["STREAMORDE"])
  Rb_rsq <- summary(lmRb)$r.squared

  #length
  #plot(dat$STREAMORDE, dat$logL)
  #abline(lm(dat$logL~dat$STREAMORDE))
  lmRl <- stats::lm(logL ~ STREAMORDE, data = dat)
  # slove for Rb w/base e
  Rl <- exp(stats::coef(lmRl)["STREAMORDE"])
  Rl_rsq <- summary(lmRl)$r.squared

  #area
  #plot(dat$STREAMORDE, dat$logA)
  #abline(lm(dat$logA~dat$STREAMORDE))
  lmRa <- stats::lm(logA ~ STREAMORDE, data = dat)
  # slove for Rb w/base e
  Ra <- exp(stats::coef(lmRa)["STREAMORDE"])
  Ra_rsq <- summary(lmRa)$r.squared

  TL <- sum(vaa$LENGTHKM)
  Dd <- sum(vaa$LENGTHKM)/max(vaa$TOTDASQKM)
  Ord <- max(vaa$STREAMORDE)

  out <- data.frame(Ord, MSL, TL, Dd, Rb,Rb_rsq, Rl, Rl_rsq, Ra, Ra_rsq, row.names = NULL)

  return(out)
}

