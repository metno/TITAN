#+ conversion of coordinates between CRSs
spatconv <- function( argv, 
                      data, 
                      extent ) {
#============================================================================
  if (argv$spatconv) {
    cat( "conversion of spatial coordinates\n")
    cat( paste(" from CRS:",argv$proj4_input_obsfiles,"\n"))
    cat( paste("   to CRS:",argv$proj4_where_dqc_is_done,"\n"))
    coord <- SpatialPoints( cbind( data$lon, data$lat),
                            proj4string = CRS( argv$proj4_input_obsfiles))
    coord.new <- spTransform( coord, CRS( argv$proj4_where_dqc_is_done))
    xy.new <- coordinates( coord.new)
  #  x  <- round( xy.new[,1],0)
  #  y  <- round( xy.new[,2],0)
    x  <- xy.new[,1]
    y  <- xy.new[,2]
    xp <- expand.grid( c( extent[1], extent[2]),
                       c( extent[3], extent[4]))
    coord <- SpatialPoints( xp,
                            proj4string = CRS(argv$proj4_input_obsfiles))
    coord.new <- spTransform( coord, CRS( argv$proj4_where_dqc_is_done))
    # define the extent for the SCT grid
    e <- extent(coord.new)
    xl <- e[1:2]
    yl <- e[3:4]
  } else {
    x  <- data$lon
    y  <- data$lat
    xl <- c( extent[1], extent[2])
    yl <- c( extent[3], extent[4])
    e  <- extent( c( xl, yl))
  }
  if (argv$debug) 
    save( argv, data, extent, x, y, xl, yl, e,
          file.path(argv$debug.dir,"input_data.RData")) 
  cat("+---------------------------------+\n")
  #
  return( list( x=x, y=y, xl=xl, yl=yl, e=e))
}
