#+ Final decision
final_decision <- function( data, 
                            dqcflag){
#==============================================================================

  nfin <- length( argv$input.files)

  # if an observation makes it this far, then it is good
  dqcflag[is.na(dqcflag)] <- 0

  # print summary
  nsus       <- length( which( dqcflag != 0))
  nok        <- length( which( dqcflag == 0))
  nsus_notna <- length( which( dqcflag != 0 & !is.na( data$value)))
  nok_notna  <- length( which( dqcflag == 0 & !is.na( data$value)))
  nnotna     <- length( which( !is.na( data$value)))
  nna        <- length( which(  is.na( data$value)))

  cat( "summary:\n")
  cat( " #  NAs, number of observations with no value (NAs)\n")
  cat( " #  sus, number of suspicious observations or no-metadata\n")
  cat( " # good, number of good observations\n")
  cat( " NOTE for sus and good, the statistics consider only observations not NAs\n")
  cat( "summary:\n")
  cat( paste0( " #  NAs= ", nna,"\n"))
  cat( paste0( " #  sus= ", nsus_notna," [",round(100*nsus_notna/nnotna,0),"%]\n"))
  cat( paste0( " # good= ", nok," [", round(100*nok_notna/nnotna,0), "%]\n"))

  # number of observation providers
  M <- nfin

  if ( M > 1) {
    for (f in 1:M) {

      faux <- data$prid == argv$prid[f]

      nsus       <- length( which( dqcflag != 0 & faux))
      nok        <- length( which( dqcflag == 0 & faux))
      nsus_notna <- length( which( dqcflag != 0 & !is.na( data$value) & faux))
      nok_notna  <- length( which( dqcflag == 0 & !is.na( data$value) & faux))
      nnotna     <- length( which( !is.na( data$value) & faux))
      nna        <- length( which( is.na(data$value) & faux))

      cat(  paste( "--> summary provider", argv$prid[f], "\n"))
      cat( paste0( "  #  NAs= ", nna, "\n"))
      cat( paste0( "  #  sus= ", nsus_notna, " [", round(100*nsus_notna/nnotna,0), "%]\n"))
      cat( paste0( "  # good= ", nok, " [", round(100*nok_notna/nnotna,0), "%]\n"))

    }
    cat("+---------------------------------+\n")
  }

  return(dqcflag)

}

