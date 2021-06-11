#+ check on the metadata
metadata_check_r <- function( argv, 
                              data, 
                              z,
                              extent, 
                              dqcflag) {
#==============================================================================
# input arguments:
# argv. list with the input arguments (command line and/or config file)
# data. input data structure
# z. elevation
# extent.
# dqcflag. data quality control flags.
#
# returned values:
# dqcflag. updated data quality control flags.
#==============================================================================
# NOTE: keep-listed stations could be flagged here

  ix <- which( is.na( dqcflag) | dqcflag == argv$code.keep)

  if ( length( ix) > 0) {
    meta <- is.na( data$lat[ix]) | 
            is.na( data$lon[ix]) |
            is.na( z[ix]) | 
            z[ix] < argv$zmin | 
            z[ix] > argv$zmax |
            is.na( data$value[ix]) 
    if ( argv$dqc_inbox_only) 
      meta <- meta | ( data$lat[ix] < extent[3] | 
                       data$lat[ix] > extent[4] |
                       data$lon[ix] < extent[1] | 
                       data$lon[ix] > extent[2] )
    if ( any( meta)) dqcflag[ix[which(meta)]] <- argv$code.nometa
  } else {
    print("no valid observations left, no metadata check")
  }

  if (argv$verbose) {
    flagaux<-dqcflag==argv$code.nometa & !is.na(dqcflag)
    print("test for no metdata, statistics over the whole dataset")
    print(paste("# observations lacking metadata and/or NAs=",
          length(which(flagaux))))
    print(paste("  # NAs                 =",
          length(which(flagaux & is.na(data$value))))) # coincides with all the NAs
    if (argv$dqc_inbox_only) {
      print(paste("  # lon-lat missing (*) =",
            length(which(flagaux & (is.na(data$lat) | 
                                    is.na(data$lon) | 
                                    data$lat < extent[3] | 
                                    data$lat > extent[4] |
                                    data$lon < extent[1] | 
                                    data$lon > extent[2] ) ))))
    } else {
      print(paste("  # lon-lat missing =",
            length(which(flagaux & (is.na(data$lat) | 
                                    is.na(data$lon) ))))) 
    }
    print(paste("  # z missing           =",
          length(which(flagaux & is.na(z)))))
    print(paste("  # z out of range      =",
          length(which(flagaux & !is.na(z) & 
                       (z<argv$zmin | z>argv$zmax) ))))
    if (argv$dqc_inbox_only) 
      print("(*) or outside the specified box")
    rm(flagaux)
    print("+---------------------------------+")
  }

  #
  return(dqcflag)

}
