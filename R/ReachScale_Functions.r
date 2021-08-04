
#' Create Network Reaches, Points and Transects
#'
#' Splits network COMIDs into reaches, identifies midpoints
#' and creates lateral transects
#'
#' @param network a single network extracted from NHDPlusV2
#' @param vaa NHDPlusV2 value added attributes table, PlusFlowlineVAA.dbf"
#' @param n number of reaches per network COMID
#' @param what one of c("reaches", "midpoints", "transects")
#' @param NHDCatchments NHDPlusV2 catchment coverage, Catchment.shp, required if what = "transects"
#' @param keep.divergent logical keep divergent flow paths
#'
#' @return sf object of either lines or points. Transects are perpendicular to the reach
#'
#' @export


net_sample <- function(network, vaa, n = 5, what = c("reaches", "midpoints", "transects"),
                       NHDCatchments = NULL, keep.divergent =  F){

  names(vaa)<-toupper(names(vaa))
  #need to supply catchment
  if (what == "transects" & is.null(NHDCatchments)){
    stop("must supply NHDCatchments for transects")
  }

  #drop divergent paths
  if(!keep.divergent){
    mainpath <- vaa[vaa$STREAMORDE == vaa$STREAMCALC, "COMID"]
    network <- network[network$COMID %in% mainpath, ]
  }

  splits <- seq(0, 1, by = 1/n)
  start <- splits[-length(splits)]
  end <- splits[-1]

  lines <- sf::st_cast(network[["geometry"]], "LINESTRING")

  #lines<-network[1,]
  lines <- sf::st_transform(lines, 5070)

  #splits flowlines into equal parts (line segments)
  lines <- lapply(lines, function(z)
    sf::st_sfc(sapply(1:n, simplify = F, function(x)
      lwgeom::st_linesubstring(z, start[x], end[x], 0.001)),
      crs = 5070))

  if(what == "reaches"){
    #name output
    names(lines) <- network$COMID
    lines <- sapply(names(lines), simplify = F, function(x)
      sf::st_as_sf(data.frame(PointID = 1:length(lines[[x]]), COMID = x, lines[[x]])))
    lines <- do.call(rbind,lines)
    lines <- sf::st_transform(lines, crs = sf::st_crs(network))
    return(lines)
  }

  #find midpoints of each line
  coords <- lapply(lines, function(g)
    sf::st_sfc(lapply(g, function(x)
      sf::st_point(get_mids(sf::st_coordinates(x)[,c(1,2)]))),
      crs = 5070))

  if (what == "midpoints"){
    #name output
    names(coords) <- network$COMID
    coords <- sapply(names(coords), simplify = F, function(x)
      sf::st_as_sf(data.frame(PointID = 1:length(coords[[x]]), COMID = x, coords[[x]])))
    coords <- do.call(rbind,coords)
    coords <- sf::st_transform(coords, crs = sf::st_crs(network))
    return(coords)
  }

  #need to supply catchment
  if (what == "transects"){

    cats <- NHDCatchments[NHDCatchments$FEATUREID %in% network$COMID, ]
    cats <- sf::st_transform(cats, crs = 5070)

    #some times the flowlines do not have catchments
    network2 <- network[which(network$COMID%in%cats$FEATUREID),]
    lines <- lines[which(network$COMID%in%cats$FEATUREID)]
    coords <- coords[which(network$COMID%in%cats$FEATUREID)]


    q <- sapply(1:length(network2$COMID), simplify = F, function(p){
      #p<-667
      #for(p in 1:length(network2$COMID)){
      #print(p)}
      cat <- cats[cats$FEATUREID==network2$COMID[p],]
      lns <- lines[[p]]
      mid_coord <- coords[[p]]

      #uses bounding box of catchment to determine transect length
      box <- sf::st_bbox(cat)

      vert <- sf::st_linestring(as.matrix(rbind(c(box$xmin,box$ymin), c(box$xmin,box$ymax))))
      vert <- sf::st_length(vert)

      horiz <- sf::st_linestring(as.matrix(rbind(c(box$xmin,box$ymin),c(box$xmax,box$ymin))))
      horiz <- sf::st_length(horiz)

      mxln <- max(horiz, vert)

      z <- sapply(1:n, simplify = F, function (s){
        #s<-4
        #coordinates for straingt line connecting start & end points of reach
        tcoords <- sf::st_coordinates(lns[s])
        pos_vec <- tcoords[1, c("X","Y")] - tcoords[nrow(tcoords), c("X","Y")]


        #perpendicular lines have dot product sum(pos_vec*opp_pos) = 0
        #transpose and make one negative
        opp_pos <- c(X = pos_vec["Y"], Y = -1*pos_vec["X"])

        #divied by magnitude such that its a unit vector...
        #unitvectors have magnutude == 1; dist(rbind(c(0,0),unitvector))
        #sqrt((opp_pos[1]^2) + (opp_pos[2]^2)) is distance formula (magnitude)
        #starting form 0,0
        unitvector <- opp_pos/sqrt((opp_pos[1]^2) + (opp_pos[2]^2))

        mpt <- mid_coord[s]
        mpt <- sf::st_coordinates(mpt)

        #calculate transect on either sidt of midpoint
        zneg <- c(mpt - (mxln * as.numeric(opp_pos)))
        zpos <- c(mpt + (mxln * as.numeric(opp_pos)))
        transect <- sf::st_sfc(sf::st_linestring(rbind(zpos,zneg)), crs = 5070)
        transect <- sf::st_cast(sf::st_intersection(transect, cat), "LINESTRING")

        #sometimes the transect can intersect catchment multipl times
        #choose the line that intersects the straight channel
        #that its peprpendicular to

        if(length(transect)>1){
          #stop()
          tcoords <- sf::st_linestring(rbind(tcoords[1,c("X","Y")],tcoords[nrow(tcoords),c("X","Y")]))
          transect <- transect[sf::st_intersects(transect, tcoords, sparse = F),]
          #plot(transect, add = T, col = "green")
        }

        return (sf::st_as_sf(transect))
      })

      if(any(unlist(lapply(z, function(x) nrow(x)))>0)){
        z <- sf::st_as_sf(data.frame(COMID = network2$COMID[p],
                                     PointID = which(unlist(lapply(z, function(x) nrow(x)))>0),
                                     sf::st_geometry(do.call(rbind, z))))
      }
      #return(z)

    })

    q <- do.call(rbind, q)
    q <- sf::st_transform(q, crs = sf::st_crs(network))

    return(q)
  }
}


#' INTERNAL FUNCTION identify midpoints
#'
#' identify mid point from a set of points
#' used in net sample points
#'
#' @param coords coordinates to identify midpoint
#'
#' @details adapted from View(maptools:::getMidpoints)
#'
#' @return  midpoint returned from set of points

get_mids <- function (coords) {
  dist <- sqrt((diff(coords[, 1])^2 + (diff(coords[, 2]))^2))
  dist_mid <- sum(dist)/2
  dist_cum <- c(0, cumsum(dist))
  end_index <- which(dist_cum > dist_mid)[1]
  start_index <- end_index - 1
  start <- coords[start_index, ]
  end <- coords[end_index, ]
  dist_remaining <- dist_mid - dist_cum[start_index]
  mid <- start + (end - start) * (dist_remaining/dist[start_index])
  return(mid)
}




#flowlines <- net_p1[[1]]

#' @title merge adjacent line segements of the same class
#'
#' @description function takes a series of flowlines and merges adjacent line segments
#' of the same class
#'
#' @param lines an sf object of LINESTRING
#' @param groupName field containing classification
#'
#' @return sf object with adjacent lines of the same class merged together
#'
#' @export


merge_lines <- function(lines, groupName = "grps"){


  q<-lapply(split(lines, sf::st_set_geometry(lines[,groupName], NULL)), function(x){

    #q<-split(lines, st_set_geometry(lines[,groupName], NULL)
    #x<-q[[3]]

    #intersecting lines with themselves
    net_int <- suppressWarnings(sf::st_intersection(x))

    #merge the list elements that contain the same index
    inlist <- net_int$origins

    outlist <- list()
    count <- 1

    #inlist elements are iterariively removed
    while(length(inlist) > 0){
      j <- inlist[[1]]

      #iteratively match across list elements
      while(any(do.call(rbind, lapply(inlist, function (x) any(x %in% j))))){
        #index of matching elements
        ind <- do.call(rbind, lapply(inlist, function (x) any(x %in% j)))
        #new values to combine
        j <- unique(append(j, unique(unlist(inlist[ind])), length(j)))
        #update to remove elements that need to be combined
        inlist <- inlist[!ind]
      }

      #update new list with values
      outlist[[count]] <- j
      count <- count + 1
    }


    if (length(setdiff(1:nrow(x), sort(unlist(outlist))))!=0|
        any(duplicated(sort(unlist(outlist))))){

      stop("something is wrong")
    }

    nmatrix <- lapply(outlist, function(z) sf::st_line_merge(sf::st_combine(x[z,])))

    nmatrix <- lapply(nmatrix, function(w)
      sf::st_sf(data.frame(groupName = unique(sf::st_set_geometry(x[,groupName],NULL))),
                geometry = w))
    nmatrix <- do.call(rbind, nmatrix)
    nmatrix$patchid <- 1:nrow(nmatrix)

    return(nmatrix)
  })



  return(do.call(rbind,q))
}


#' @title sf to tidygraph
#'
#' @description creates tidy graph from sf object
#'
#' @param x sf network
#' @param directed logical
#'
#' @return super useful graph
#'
#' @details adapted from https://www.r-spatial.org/r/2019/09/26/spatial-networks.html embedded within class metrics
#'
#' @importFrom magrittr %>%
#'
#' @export
#'

sf_to_tidygraph <- function(x, directed = TRUE) {
  #x<-splitln
  L1<-edgeID<-.<-xy<-start_end<-nodeID<-geometry<-NULL
  edges <- x %>%
    dplyr::mutate(edgeID = c(1:dplyr::n()))

  nodes <- edges %>%
    sf::st_coordinates() %>%
    dplyr::as_tibble() %>%
    dplyr::rename(edgeID = L1) %>%
    dplyr::group_by(edgeID) %>%
    dplyr::slice(c(1, dplyr::n())) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(start_end = rep(c('start', 'end'), times = dplyr::n()/2)) %>%
    dplyr::mutate(xy = paste(.$X, .$Y)) %>%
    dplyr::mutate(nodeID = tidygraph::group_indices(., factor(xy, levels = unique(xy)))) %>%
    dplyr::select(-xy)

  source_nodes <- nodes %>%
    dplyr::filter(start_end == 'start') %>%
    dplyr::pull(nodeID)

  target_nodes <- nodes %>%
    dplyr::filter(start_end == 'end') %>%
    dplyr::pull(nodeID)

  edges = edges %>%
    dplyr::mutate(from = source_nodes, to = target_nodes)

  nodes <- nodes %>%
    tidygraph::distinct(nodeID, .keep_all = TRUE) %>%
    dplyr::select(-c(edgeID, start_end)) %>%
    sf::st_as_sf(coords = c('X', 'Y')) %>%
    sf::st_set_crs(sf::st_crs(edges))

  graph <- tidygraph::tbl_graph(nodes = nodes, edges = dplyr::as_tibble(edges), directed = directed)

  graph %>%
    tidygraph::activate(edges) %>%
    dplyr::mutate(length = sf::st_length(geometry))

}


#' @title average distance between patches
#'
#' @description embedded within class metrics
#'
#' @param graph result from sf_todygraph
#' @param net_HGP a classified stream network
#' @param groupName name of field containing classification result
#' @param patch calss value
#'
#' @return long pairwise patch distances
#'
#' @importFrom magrittr %>%
#' @export


patch_dist <- function(graph, net_HGP, groupName, patch = "2"){
  edges<-geometry<-NULL
  plen <- split(net_HGP, sf::st_set_geometry(net_HGP[,groupName], NULL))
  plen <- plen[[patch]]

  nodes <- graph %>%
    tidygraph::activate(nodes)%>%
    dplyr::as_tibble() %>%
    sf::st_as_sf()

  #nodes where patch intersects a matrix
  intnodes <- suppressWarnings(sf::st_intersection(plen, nodes))

  #calculate length of edges and pairwise distances to all nodes
  graph <- graph %>% tidygraph::activate(edges) %>% dplyr::mutate(length = sf::st_length(geometry))
  dist_tbl <- igraph::distances(graph = graph, weights = graph %>% tidygraph::activate(edges) %>% dplyr::pull(length))

  i <- split(intnodes$nodeID,intnodes$patchid)
  j <- split(intnodes$nodeID,intnodes$patchid)

  out<-lapply(i, function(p)
    lapply(j, function(q) {
      m<-dist_tbl[p, q]
      #closest vertex
      d <- min(m)
      indx <- which(m == d, arr.ind = TRUE)
      data.frame(to_node = p[indx[1,1]],
                 from_node = q[indx[1,2]], d)
    }))

  out <- do.call(rbind, lapply(out, function(z) do.call(rbind, z)))
  pnam <- stats::setNames(data.frame(do.call(rbind, strsplit(row.names(out), "\\."))), c("patch_i", "patch_j"))
  out<-cbind(pnam, out)
  out<-merge(data.frame(patch_j = plen$patchid, patch_j_len = sf::st_length(plen)), out, by = "patch_j")
  out<-merge(data.frame(patch_i = plen$patchid, patch_i_len = sf::st_length(plen)), out, by = "patch_i")

  return(out)
}


#' @title Distance-Baised Dendritic Connectivity Index
#'
#' @description calculates modified dentritic connectivity
#' index baised on watercourse distance
#'
#' @param plen_i length of patch i
#' @param plen_j length of patch j
#' @param dist_ij distance separating patch i and patch j
#' @param mu mean dispersal distance
#'
#' @return long pairwise patch distances
#'
#' @details index modified from cote et al 2009
#'
#' median parent offspring dispersal distance is
#' 12km = 12,000m (Comte and Olden 2017)
#'
#' @example
#' mu = 1000; mean(rexp(100000, rate = 1/mu)) ~ 1,000
#'
#' @export



DCI_dist <- function(plen_i, plen_j, dist_ij, mu = 12000){
  #DCI Calculation

  tlen <- sum(unique(plen_i))
  li <- (plen_i / tlen)
  lj <- (plen_j / tlen)

  #probability is dist_tblce based on exponential distribution
  #i.e. the probability of making it past a given dist_tblce
  p <- 1 - stats::pexp(dist_ij, rate = 1/mu)

  #DCI calculation
  DCI <- sum(p * li * lj) * 100

  return(units::set_units(DCI,NULL))

}
