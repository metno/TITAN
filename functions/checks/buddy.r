#+ check against first-guess (deterministic)
buddy <- function( argv, 
                   data,
                   z, 
                   dqcflag) {
#------------------------------------------------------------------------------
# test deviations between observations and spatial trend in a region
# See the wiki 
#==============================================================================

  cat( paste0( "buddy-check (code=", argv$code.buddy, ")\n"))

  nfin  <- length( argv$input.files)

  ndata <- length( data$lat)

  # number of observation providers
  M <- nfin

  # number of tests
  if ( ( N <- length( argv$inner_radius.buddy)) == 0) return( FALSE)

  nsus <- vector( mode="numeric", length=N)

  debug <- FALSE

  background_values <- -999
  background_uncertainties <- -999

  if ( is.null( argv$transf.buddy)) argv$transf.buddy <- F

  if ( length( argv$doit.buddy) != M) 
    argv$doit.buddy <- rep( argv$doit.buddy[1], M)
  if ( length( argv$prio.buddy) != M) 
    argv$prio.buddy <- rep( argv$prio.buddy[1], M)

  if ( length( argv$tpos.buddy) != (M*N)) 
    argv$tpos.buddy <- rep( argv$tpos.buddy[1], N*M)
  if ( length( argv$tneg.buddy) != (M*N)) 
    argv$tneg.buddy <- rep( argv$tneg.buddy[1], N*M)

  if ( length( argv$background_elab_type.buddy) != N) 
    argv$background_elab_type.buddy <- rep( argv$background_elab_type.buddy[1], N)
  if ( length( argv$num_min_outer.buddy) != N) 
    argv$num_min_outer.buddy <- rep( argv$num_min_outer.buddy[1], N)
  if ( length( argv$num_max_outer.buddy) != N) 
    argv$num_max_outer.buddy <- rep( argv$num_max_outer.buddy[1], N)
  if ( length( argv$outer_radius.buddy) != N) 
    argv$outer_radius.buddy <- rep( argv$outer_radius.buddy[1], N)
  if ( length( argv$num_min_prof.buddy) != N) 
    argv$num_min_prof.buddy <- rep( argv$num_min_prof.buddy[1], N)
  if ( length( argv$min_elev_diff.buddy) != N) 
    argv$min_elev_diff.buddy <- rep( argv$min_elev_diff.buddy[1], N)

  #............................................................................
  # set doit/prio vectors
  doit <- vector( length=ndata, mode="numeric"); doit[]<-NA
  prio <- vector( length=ndata, mode="numeric"); prio[]<-NA
  for (f in 1:M) {
    ix <- which( data$prid == argv$prid[f])
    if ( length(ix) == 0) next
    doit[ix] <- argv$doit.buddy[f]
    prio[ix] <- argv$prio.buddy[f]
  }
  prio_unique <- sort( unique( prio, na.rm=T), decreasing=T)
  rm( ix)

  #............................................................................
  # prepare vectors of valid and admissible values
  if (argv$variable == "T") {
    values_mina <- data$value - 20
    values_maxa <- data$value + 20
    values_minv <- data$value - 1
    values_maxv <- data$value + 1
  } else if (argv$variable == "RR") {
    values_mina <- pmin( pmax( data$value - 10, 0), 
                         pmax( data$value - 0.5 * data$value, 0))
    values_maxa <- pmax( data$value + 10, data$value + 0.5 * data$value)
    values_minv <- pmin( pmax( data$value - 1, 0), 
                         pmax( data$value - 0.1 * data$value, 0))
    values_maxv <- pmax( data$value + 1, data$value + 0.1 * data$value)
  }

  #............................................................................
  # data transformation
  if (argv$transf.buddy) {
    values_mina <- boxcox( x=values_mina, lambda=argv$boxcox.lambda)
    values_maxa <- boxcox( x=values_maxa, lambda=argv$boxcox.lambda)
    values_minv <- boxcox( x=values_minv, lambda=argv$boxcox.lambda)
    values_maxv <- boxcox( x=values_maxv, lambda=argv$boxcox.lambda)
    data$value  <- boxcox( x=data$value,  lambda=argv$boxcox.lambda)
  }

  #............................................................................
  # test
  for (i in 1:argv$i.buddy) {

    nsus[]<-0

    for (j in 1:N) {

      tpos <- vector( length=ndata, mode="numeric"); tpos[]<-NA
      tneg <- vector( length=ndata, mode="numeric"); tneg[]<-NA

      for (f in 1:M) {
        ix <- which( data$prid == argv$prid[f])
        if ( length(ix) == 0) next
        tpos[ix] <- argv$tpos.buddy[(j-1)*M+f]
        tneg[ix] <- argv$tneg.buddy[(j-1)*M+f]
      }
      prio_unique <- sort( unique( prio, na.rm=T), decreasing=T)
      rm( ix)

      for (k in 1:length(prio_unique)) {

        # use only (probably) good observations with doit!=0
        flag_aux<-( ( is.na( dqcflag) | dqcflag == argv$code.keep) &
                    !is.na( data$lon) & !is.na( data$lat) & 
                    !is.na( data$value) &
                    doit != 0 & !is.na( prio))
        if (argv$variable == "T") flag_aux <- flag_aux & !is.na(z)
        ix <- which( flag_aux); rm( flag_aux)
        t0a <- Sys.time()


        obsToCheck_n<-length(ix)
        if (obsToCheck_n>0) {
          # define global 1D vector used in statSpat (1D for fast access)
          obsToCheck_lon  <- data$lon[ix]
          obsToCheck_lat  <- data$lat[ix]
          obsToCheck_z    <- z[ix]
          obsToCheck_val  <- data$value[ix]
          obsToCheck_mina <- values_mina[ix]
          obsToCheck_maxa <- values_maxa[ix]
          obsToCheck_minv <- values_minv[ix]
          obsToCheck_maxv <- values_maxv[ix]
          obsToCheck_tpos <- tpos[ix]
          obsToCheck_tneg <- tneg[ix]
          obsToCheck_chk  <- rep( 0, obsToCheck_n)
          # check only those observations with priorities geq than this
          obsToCheck_chk[prio[ix]>=prio_unique[k]] <- 1


          res <- fgt( points = Points( obsToCheck_lat, 
                                       obsToCheck_lon, 
                                       obsToCheck_z),
                      obsToCheck_val,
                      obsToCheck_chk,
                      background_values,
                      background_uncertainties,
                      argv$background_elab_type.buddy[j],
                      argv$num_min_outer.buddy[j],
                      argv$num_max_outer.buddy[j],
                      argv$inner_radius.buddy[j],
                      argv$outer_radius.buddy[j],
                      100,
                      argv$num_min_prof.buddy[j],
                      argv$min_elev_diff.buddy[j],
                      obsToCheck_mina,
                      obsToCheck_maxa,
                      obsToCheck_minv,
                      obsToCheck_maxv,
                      obsToCheck_tpos,
                      obsToCheck_tneg,
                      debug)
 
          flag <- res[[1]]

          # suspect if: 
          sus<-which( flag[1:obsToCheck_n] == 1 &
                      is.na(dqcflag[ix]) &
                      doit[ix]==1 )

          # set dqcflag
          if (length(sus)>0) dqcflag[ix[sus]] <- argv$code.buddy

        } else {
          cat( "no valid observations left, no buddy check\n")
        }
        nsus[j] <- ifelse( exists("sus"), length(sus), 0)
        t1a  <- Sys.time()
        str  <- " (#TOT "
        str1 <- ""
        for (f in 1:M) {
          if (f>1) str<-paste0(str,"; ")
          str <- paste0( str, "prid", argv$prid[f], "=", 
                  length( which( dqcflag==argv$code.buddy & 
                                 data$prid==argv$prid[f])))
          str1 <- paste0( str1, "::prid ", argv$prid[f]," ",
                                "prio=", argv$prio.buddy[f],",",
                                "tpos=", argv$tpos.buddy[(j-1)*M+f],",",
                                "tneg=", argv$tneg.buddy[(j-1)*M+f])
        }
        str<-paste0(str,")")
        cat( paste0( "iteration=", i,
                     "/test=", j,
                     "/prio>=", prio_unique[k],
                     "/dqc param:",
                     "inner_rad=", argv$inner_radius.buddy[j],",",
                     "outer_rad=", argv$outer_radius.buddy[j],",",
                     str1,
                     "/time ",round(t1a-t0a,1),attr(t1a-t0a,"unit"),"\n"))
        cat( paste0(nsus[j]," new suspect observations",str,"\n"))
        rm(str)
      } # end for k
    } # end for j
    if ( sum(nsus) <= argv$break.buddy) break
  }  # end for i
  cat("+---------------------------------+\n")
  #
  return(dqcflag)
}

