#+ Isolation test
isolation_test <- function( argv,
                            data, 
                            dqcflag){
#==============================================================================

  cat( paste0( "isolation test (code=", argv$code.iso, ")\n"))

  ndata <- length(data$lat)

  # number of observation providers
  M <- nfin

  # set doit vector
  if ( length( argv$doit.iso) != M) 
    argv$doit.iso <- rep( argv$doit.iso[1], M)
  doit   <-vector( length=ndata, mode="numeric")
  doit[] <- NA
  for (f in 1:M) 
    doit[data$prid==argv$prid[f]] <- argv$doit.iso[f]

  #
  ix <- which( is.na(dqcflag) & doit!=0)
  if ( length(ix) > 0) {
    flag <- isolation_check( points = Points( data$lat[ix], data$lon[ix], rep( 0, length(ix))),
                             argv$n.iso,
                             argv$dr.iso,
                             Inf) 

    # suspect if: 
    sus <- which( flag == 1 &
                  is.na(dqcflag[ix]) &
                  doit[ix] == 1 )

    # set dqcflag
    if ( length(sus) > 0) dqcflag[ix[sus]] <- argv$code.iso

  } else {
    cat( "no valid observations left, no check for isolated observations\n")
  }

  cat( paste( "# isolated observations=", length( which( dqcflag == argv$code.iso)), "\n"))
  cat( "+---------------------------------+\n")

  #
  return(dqcflag)
}
