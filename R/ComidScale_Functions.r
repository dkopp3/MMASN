
#' Channel Shape
#'
#' Calculates metrics that describe shape of the COMIDs wihtin a network
#'
#' @param network a single network extracted from NHDPlusV2
#' @param elevslope NHDPlusV2 elevation slope table, elevslope.dbf"
#'
#' @return data.frame of NHDPlusV2 attributes: slope (\code{SLOPE}), smoothed
#' maximum \code{MAXELEVSMO} and minimum elevation (\code{MINELEVSMO}) and
#' COMID length (\code{LENGTHKM}), mean elevation (\code{elev}) and channel sinuosity
#' (\code{sinuosity}).
#'
#' @details
#' Elevation is the mean of the smoothed max and min elevation
#'
#' Sinuosity is COMID length divided by straight line length
#'
#'
#' @export

chn_geom <- function (network, elevslope){
  slope <- elevslope
  names(slope) <- toupper(names(slope))
  net <- sf::st_transform(network, crs = 5070)


  slopeelev <- slope[slope$COMID%in%network$COMID, c("COMID", "SLOPE", "MAXELEVSMO", "MINELEVSMO")]
  slopeelev$elev <- apply(slopeelev[,c("MAXELEVSMO", "MINELEVSMO")],1,mean)

  # to ensure the input is linestring geometry -
  # vpu 17 read in a multiline and caused error
  net <- suppressWarnings(sf::st_cast(net, "LINESTRING"))
  tot.len <- sf::st_length(net)

  xy <- sf::st_coordinates(net)
  xy <- data.frame(xy)
  xy <- split(xy, xy$L1)
  #create line string from 1st and lask point of each comid
  xy <- lapply(xy, function(x) sf::st_linestring(as.matrix(rbind(x[1, c("X", "Y")], x[nrow(x), c("X", "Y")]))))

  #name list elements
  names(xy) <- net$COMID

  #check plots
  #plot(st_geometry(net[net$COMID==23963061,]))
  #plot(xy[["23963061"]], add = T)

  #calculate channel sinuosity
  sinuosity <- do.call(rbind, sapply(net$COMID, simplify = F, function(x)
    data.frame(COMID = x,
               LENGTHKM = net$LENGTHKM[net$COMID==x],
               sinuosity = sf::st_length(net[net$COMID==x,]) / sf::st_length(xy[[as.character(x)]]))))
  sinuosity$sinuosity <- units::set_units(sinuosity$sinuosity, NULL)
  out <- merge(slopeelev, sinuosity, by = "COMID")

  return(out)
}



#' Confluence Attributes
#'
#' Calculates metrics related to confluences in a network
#'
#' @param network a single network extracted from NHDPlusV2
#' @param vaa NHDPlusV2 value added attributes table, PlusFlowlineVAA.dbf"
#' @param flow NHDPlusV2 flow table, PlusFlow.dbf
#'
#' @return data.frame of confluence angle (\code{Confl_angle}),
#' area ratio (\code{AreaRatio}), total upstream area \code{Confl_AREA},
#' order of the tributary (\code{Confl_order}), and tributary classification
#' \code{Confl_class}
#'
#' @details Lambert Conformal Cone (crs = 102004) used to
#' preserve angles Seybold et al. (2017) is no longer available
#'
#' in some cases tributary angle is NAN because average tributary direction is near verticle.
#' the function cannot determine starts/end of line needed to calculate direction position
#' vector direction. this could be could fix this by changing the maximum verticies
#'
#' Complex confluences >2 tributaries are incorrect.
#'
#' Confluence angle is the angle two tributarys meet. Calculated as the intersection between
#' two slopes orthoganal regression slopes reflecting each tributary's average direction
#'
#' Area ratio is the ratio of the tributary catchment area to the mainstem catchment area. Mainstem
#' is identified as tributary with the largest catchment area
#'
#' confluence area is the total area upstream of the confluence
#'
#' confluence order is the stream order directly downstream of the confluence
#'
#' confluence class is the stream order of each tributary in the confluence
#'
#' @export

net_confl <- function (network, vaa, flow){

  names(vaa) <- toupper(names(vaa))
  names(flow) <- toupper(names(flow))

  NHD <- network
  #NHD <- sf::st_cast(NHD,"LINESTRING")
  #Lambert Conformal Cone - preserves angles as per seybold et al 2017
  #NHD <- sf::st_transform(NHD, crs = 102009)

  NHD <- merge(NHD, vaa[, c("COMID", "STREAMORDE", "STREAMCALC", "STREAMLEVE", "TOTDASQKM")], by = "COMID")

  #remove divergent flow paths
  NHD <- NHD[NHD$STREAMORDE == NHD$STREAMCALC, ]

  #identify confluences (>1 flowline upstream)
  confluences <- lapply(split(NHD, as.character(NHD$COMID)),
                        function(x) flow[flow$TOCOMID == x$COMID, "FROMCOMID"])

  #ind <- unlist(lapply(confluences, function(x) length(x) > 1))

  confluences <- lapply(confluences, function(x) NHD[NHD$COMID%in%x,])
  ind <- unlist(lapply(confluences, function(x) nrow(x))>1)
  confluences <- confluences[ind]

  angles <- do.call(rbind,lapply(confluences, function(q){
    coords <- list(q[1,], q[2,])
    coords <- lapply(coords, function(z) sf::st_coordinates(z))
    pos <- lapply(coords,function(x) pos_vec(x, diagnostics = F))
    angle_deg <- vec_angle(pos[[1]], pos[[2]])
    return(angle_deg)
  }))

  #angles <- do.call(rbind,lapply(confluences, function(q) confluence_angle(q[1,],q[2,])))
  colnames(angles) <- "Confl_angle"

  #summarize vaa table valuse from upstream confluences
  attrssum <- do.call(rbind,lapply(confluences, function(t)
    data.frame(AreaRatio = min(t$TOTDASQKM) / max(t$TOTDASQKM),
               Confl_AREA = sum(t$TOTDASQKM),
               Confl_ord = ifelse(anyDuplicated(t$STREAMORDE),
                                  unique(t$STREAMORDE) + 1,
                                  max(t$STREAMORDE)),
               Confl_class = paste(t$STREAMORDE[order(t$STREAMORDE)],
                                   collapse = "."))))

  out <- merge(angles, attrssum, by = 0)
  names(out)[1] <- "COMID"

  return(out)
}


#' Average line direction
#'
#' Calculates average direction of a line used for concluence angle
#'
#' @param p an points of an sf LINESTRING
#' @param maxvert maximum number of line verticies used to fit regression
#' @param diagnostics logical used for plotting
#'
#' @return position vector of a line unless diagnostics = T
#'
#' @details this function uses orthogonal regression to find the
#' average direction of a series of points and returns
#' the position vector.
#'
#' max vert is the maximum number of verticies to fit average line
#'
#' in some instances the function will return NAN when
#' line is near verticle. Could fix this by changing the
#' maximum verticies
#'
#' @example
#' #two lines
#' x <- confluences[[54]][1,]
#' y <- confluences[[54]][2,]
#'
#' #coordinates for each line
#' coords <- list(x, y)
#' coords <- lapply(coords, function(z) sf::st_coordinates(z))
#'
#' #apply each function over
#' pos <- lapply(coords,function(x) pos_vec(x, diagnostics = T))
#'
#' #plot diagnostics
#' plot(st_coordinates(confluences[[54]]))
#'
#' abline(pos[[1]]$y_int1, pos[[1]]$beta1)
#' abline(pos[[2]]$y_int1, pos[[2]]$beta1)
#'
#' points(pos[[1]]$start, pch = 19, col = "green")
#' points(pos[[1]]$end, pch = 19, col = "red")
#'
#' points(pos[[2]]$start, pch = 19, col = "green")
#' points(pos[[2]]$end, pch = 19, col = "red")
#' @export


pos_vec <- function(p , maxvert = nrow(p), diagnostics = T){
  #p <- st_coordinates(sfline)
  #plot(st_geometry(st_as_sf(data.frame(p), coords=c("X","Y"), crs = st_crs(sfline))), add = T)
  #maxvert = nrow(p)-2
  #p <- coords[[2]]

  #only considering the 1st 5 vertices from confluence
  if(nrow(p)>maxvert){
    p <- p[(nrow(p)-maxvert):nrow(p),]
  }

  x1 <- p[, c(1)]
  y1 <- p[, c(2)]

  #plot(x1,y1, ylim = c((min(y1) - 0.002), (max(y1)+0.003)),
  #     xlim = c((min(x1)-0.002), (max(x1)+0.003)))
  #average tributary direction orthogonal regression (total) on x y coords
  #this is the average direction of the flowline, orthogonal regression
  #minimizes perpendicular error
  v <- stats::prcomp(cbind(x1, y1))$rotation
  beta1 <- v[2,1] / v[1,1]

  #given the average direction of the trib, find intercept.
  #force lines to pass through confluence
  #the last row is start in coords subset is start (confluence)

  #y_int1 <- y1[length(y1)] - beta1 * x1[length(x1)]
  #new_xy1_1 <- cbind(x1[length(x1)], y1[length(y1)])

  y_int1 <- mean(y1) - beta1 * mean(x1)
  #abline(y_int1, beta1)
  new_xy1_1 <- data.frame(X = x1, Y = beta1 * (x1) + y_int1)

  closest <- as.matrix(stats::dist(rbind(data.frame(X = x1[length(x1)],
                                                    Y = y1[length(y1)]),
                                         new_xy1_1)))[1,]
  new_xy1_1 <- new_xy1_1[which.min(closest[-1]),]
  #points(new_xy1_1, col = "green")
  #the first row is the terminous
  #find the predicted point that is closest to the
  #most upstream end of the flowline. needed to define
  #the correct direction of the vector because some
  new_xy1_2 <- data.frame(X = x1, Y = beta1 * (x1) + y_int1)

  closest <- as.matrix(stats::dist(rbind(data.frame(X = x1[1],
                                                    Y = y1[1]),
                                         new_xy1_2)))[1,]
  new_xy1_2 <- new_xy1_2[which.min(closest[-1]),]
  #points(new_xy1_2, col = "red")
  #position vector is terminus minus start
  #makes vector start at orgin (0,0)
  pos_vec_xy1 <- new_xy1_2 - new_xy1_1

  if(diagnostics){
    return(list(pos_vec_xy1 = pos_vec_xy1,
                start = new_xy1_1,
                end = new_xy1_2,
                y_int1 = y_int1,
                beta1 = beta1))
  } else {

    return(pos_vec_xy1)

  }
}




#' Intersection Angle
#'
#' @param pos1 position vector of line 1
#' @param pos2 position vector of line 2
#' @param deg logical return angle in degrees
#'
#' @return intersection angle of two positon vectors
#'
#' @example
#' pos1 <- c(2,2)
#' pos2 <- c(0,3)
#'
#' dat<-rbind(pos1, c(0,0), pos2, c(0,0))
#'
#' vec_angle(pos1, pos2)
#'
#' plot(dat, xlim = c(-3,3), ylim = c(-3,3))
#'
#' lines(dat)
#'
#' @export


vec_angle <- function(pos1, pos2, deg =T){
  #dot product divided by magnitude
  alpha <- acos(sum(pos1 * pos2) / (sqrt(sum(pos1^2)) * sqrt(sum(pos2^2))))

  if(deg){
    alpha<-(alpha  * 180) / (base::pi)#in degs
  }
  return(alpha)
}



#' Network StreamCat Data
#'
#' Extracts StreamCat variables for Network COMID
#'
#' @param comid NHDPlusV2 COMID
#' @param strcat_path directory containing StreamCat data
#' @param csvfile Name of StreamCat File
#' @param vpu NHDPlusV2 Vector Processing Unit
#'
#' @return \code{data.frame} of streamCat variables
#'
#' @export
#'

get_streamcat <- function (comid, strcat_path = "Data/StreamCat",
                           csvfile = "Lithology", vpu = "01"){

  Region_dir <- grep(vpu, list.dirs(strcat_path), value = T)
  f_name <- grep(csvfile, list.files(Region_dir, full.names = T), value = T)

  if(length(f_name) != 1){
    stop("multiple StreamCat files, give more specific name")
  }

  strcat <- utils::read.csv(file = f_name)

  out <- strcat[strcat$COMID %in% comid, ]

  return(out)
}
