#+ plausibility test
plausibility_test <- function( argv, data, dqcflag){
#..............................................................................
# use range_check from titanlib
#==============================================================================

  if ( length( ix <- which( ( is.na(dqcflag) | dqcflag==argv$code.keep) & 
                              range_check( as.numeric( data$value), 
                                           as.numeric(  argv$vmin),
                                           as.numeric(  argv$vmax)) == 1)) > 0) 
    dqcflag[ix] <- argv$code.p

  # verbose
  cat(paste("plausibility test (",argv$code.p,")\n"))
  cat(paste(" min/max thresholds =",argv$vmin,argv$vmax,"\n"))
  cat(paste(" # <min=", length( which( dqcflag==argv$code.p &
                                       data$value<argv$vmin)),"\n"))
  cat(paste(" # >max=", length( which( dqcflag==argv$code.p &
                                       data$value>argv$vmax)),"\n"))
  cat(paste(" # suspect observations=",length(which(dqcflag==argv$code.p &
                                                    !is.na(dqcflag))),"\n"))
  cat("+---------------------------------+\n")

  #
  return(dqcflag)
}
