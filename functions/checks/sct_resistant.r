#+ spatial consistency test resistant to outlier in the data
sct_resistant_r <- function( argv, 
                             data,
                             z, 
                             dqcflag) {
#------------------------------------------------------------------------------
# test deviations between observations and spatial trend in a region
# See the wiki 
#==============================================================================

  cat( paste0( "sct (code=", argv$code.sct, ")\n"))

  ndata <- length( data$lat)

  # number of observation providers
  M  <- length( argv$input.files)

  # number of tests
  if ( ( N <- length( argv$inner_radius.sct)) == 0) return( FALSE)

  nsus <- vector( mode="numeric", length=N)

  debug <- FALSE

  background_values <- -999
  background_uncertainties <- -999

  if ( is.null( argv$transf.sct)) argv$transf.sct <- F

  if ( length( argv$doit.sct) != M) 
    argv$doit.sct <- rep( argv$doit.sct[1], M)
  if ( length( argv$prio.sct) != M) 
    argv$prio.sct <- rep( argv$prio.sct[1], M)

  if ( length( argv$tpos.sct) != (M*N)) 
    argv$tpos.sct <- rep( argv$tpos.sct[1], N*M)
  if ( length( argv$tneg.sct) != (M*N)) 
    argv$tneg.sct <- rep( argv$tneg.sct[1], N*M)
  if ( length( argv$eps2.sct) != (M*N)) 
    argv$eps2.sct <- rep( argv$eps2.sct[1], N*M)

  if ( length( argv$background_elab_type.sct) != N) 
    argv$background_elab_type.sct <- rep( argv$background_elab_type.sct[1], N)
  if ( length( argv$num_min_outer.sct) != N) 
    argv$num_min_outer.sct <- rep( argv$num_min_outer.sct[1], N)
  if ( length( argv$num_max_outer.sct) != N) 
    argv$num_max_outer.sct <- rep( argv$num_max_outer.sct[1], N)
  if ( length( argv$outer_radius.sct) != N) 
    argv$outer_radius.sct <- rep( argv$outer_radius.sct[1], N)
  if ( length( argv$num_min_prof.sct) != N) 
    argv$num_min_prof.sct <- rep( argv$num_min_prof.sct[1], N)
  if ( length( argv$min_elev_diff.sct) != N) 
    argv$min_elev_diff.sct <- rep( argv$min_elev_diff.sct[1], N)
  if ( length( argv$min_horizontal_scale.sct) != N) 
    argv$min_horizontal_scale.sct <- rep( argv$min_horizontal_scale.sct[1], N)
  if ( length( argv$max_horizontal_scale.sct) != N) 
    argv$max_horizontal_scale.sct <- rep( argv$max_horizontal_scale.sct[1], N)
  if ( length( argv$kth_closest_obs_horizontal_scale.sct) != N) 
    argv$kth_closest_obs_horizontal_scale.sct <- rep( argv$kth_closest_obs_horizontal_scale.sct[1], N)
  if ( length( argv$vertical_scale.sct) != N) 
    argv$vertical_scale.sct <- rep( argv$vertical_scale.sct[1], N)

  #............................................................................
  # set doit/prio vectors
  doit <- vector( length=ndata, mode="numeric"); doit[]<-NA
  prio <- vector( length=ndata, mode="numeric"); prio[]<-NA
  for (f in 1:M) {
    ix <- which( data$prid == argv$prid[f])
    if ( length(ix) == 0) next
    doit[ix] <- argv$doit.sct[f]
    prio[ix] <- argv$prio.sct[f]
  }
  prio_unique <- sort( unique( prio, na.rm=T), decreasing=T)
  rm( ix)

  #............................................................................
  # prepare vectors of valid and admissible values
  if (argv$variable == "T") {
    values_mina <- data$value - argv$a_delta.sct
    values_maxa <- data$value + argv$a_delta.sct
    values_minv <- data$value - argv$v_delta.sct
    values_maxv <- data$value + argv$v_delta.sct
  } else if (argv$variable == "RR") {
    values_mina <- pmin( pmax( data$value - argv$a_delta.sct, 0), 
                         pmax( data$value - argv$a_fact.sct * data$value, 0))
    values_maxa <- pmax( data$value + argv$a_delta.sct,
                         data$value + argv$a_fact.sct * data$value)
    values_minv <- pmin( pmax( data$value - argv$v_delta.sct, 0), 
                         pmax( data$value - argv$v_fact.sct * data$value, 0))
    values_maxv <- pmax( data$value + argv$v_delta.sct,
                         data$value + argv$v_fact.sct * data$value)
  }

  #............................................................................
  # data transformation
  if (argv$transf.sct) {
    values_mina <- boxcox( x=values_mina, lambda=argv$boxcox.lambda)
    values_maxa <- boxcox( x=values_maxa, lambda=argv$boxcox.lambda)
    values_minv <- boxcox( x=values_minv, lambda=argv$boxcox.lambda)
    values_maxv <- boxcox( x=values_maxv, lambda=argv$boxcox.lambda)
    data$value  <- boxcox( x=data$value,  lambda=argv$boxcox.lambda)
  }

  #............................................................................
  # test
  for (i in 1:argv$i.sct) {

    nsus[]<-0

    for (j in 1:N) {

      tpos <- vector( length=ndata, mode="numeric"); tpos[]<-NA
      tneg <- vector( length=ndata, mode="numeric"); tneg[]<-NA
      eps2 <- vector( length=ndata, mode="numeric"); eps2[]<-NA

      for (f in 1:M) {
        ix <- which( data$prid == argv$prid[f])
        if ( length(ix) == 0) next
        tpos[ix] <- argv$tpos.sct[(j-1)*M+f]
        tneg[ix] <- argv$tneg.sct[(j-1)*M+f]
        eps2[ix] <- argv$eps2.sct[(j-1)*M+f]
      }
      prio_unique <- sort( unique( prio, na.rm=T), decreasing=T)
      rm( ix)

      for (k in 1:length(prio_unique)) {

        # use only (probably) good observations with doit!=0
        flag_aux<-( ( is.na( dqcflag) | dqcflag == argv$code.keep) &
                    !is.na( data$lat) & !is.na( data$lon) &
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
          obsToCheck_eps2 <- eps2[ix]
          obsToCheck_chk  <- rep( 0, obsToCheck_n)
          # check only those observations with priorities geq than this
          obsToCheck_chk[prio[ix]>=prio_unique[k]] <- 1

          cat("\n")
          res <- sct_resistant(
                      points = Points( obsToCheck_lat, 
                                       obsToCheck_lon, 
                                       obsToCheck_z),
                      obsToCheck_val,
                      obsToCheck_chk,
                      background_values,
                      argv$background_elab_type.sct[j],
                      argv$num_min_outer.sct[j],
                      argv$num_max_outer.sct[j],
                      argv$inner_radius.sct[j],
                      argv$outer_radius.sct[j],
                      100,
                      argv$num_min_prof.sct[j],
                      argv$min_elev_diff.sct[j],
                      argv$min_horizontal_scale.sct[j],
                      argv$max_horizontal_scale.sct[j],
                      argv$kth_closest_obs_horizontal_scale.sct[j],
                      argv$vertical_scale.sct[j],
                      obsToCheck_mina,
                      obsToCheck_maxa,
                      obsToCheck_minv,
                      obsToCheck_maxv,
                      obsToCheck_eps2,
                      obsToCheck_tpos,
                      obsToCheck_tneg,
                      debug,
                      argv$basic.sct)
 
          flag <- res[[1]]

          # suspect if: 
          sus<-which( flag[1:obsToCheck_n] == 1 &
                      is.na(dqcflag[ix]) &
                      doit[ix]==1 )

          # set dqcflag
          if (length(sus)>0) dqcflag[ix[sus]] <- argv$code.sct

        } else {
          cat( "no valid observations left, no sct check\n")
        }
        nsus[j] <- ifelse( exists("sus"), length(sus), 0)
        t1a  <- Sys.time()
        str  <- " (#TOT "
        str1 <- ""
        for (f in 1:M) {
          if (f>1) str<-paste0(str,"; ")
          str <- paste0( str, "prid", argv$prid[f], "=", 
                  length( which( dqcflag==argv$code.sct & 
                                 data$prid==argv$prid[f])))
          str1 <- paste0( str1, "::prid ", argv$prid[f]," ",
                                "prio=", argv$prio.sct[f],",",
                                "eps2=", argv$eps2.sct[(j-1)*M+f],",",
                                "tpos=", argv$tpos.sct[(j-1)*M+f],",",
                                "tneg=", argv$tneg.sct[(j-1)*M+f])
        }
        str<-paste0(str,")")
        cat( paste0( "++++>> SCT iteration=", i,
                     "/test=", j,
                     "/prio>=", prio_unique[k],
                     "/dqc param:",
                     "inner_rad=", argv$inner_radius.sct[j],",",
                     "outer_rad=", argv$outer_radius.sct[j],",",
                     "nmin_outer=", argv$num_min_outer.sct[j],",",
                     "nmax_outer=", argv$num_max_outer.sct[j],",",
                     "backg_type=", argv$background_elab_type.sct[j],",",
                     "nmin_prof=", argv$num_min_prof.sct[j],",",
                     "min_elev_diff=", argv$min_elev_diff.sct[j],",",
                     "min_hscale=", argv$min_horizontal_scale.sct[j],",",
                     "max_hscale=", argv$max_horizontal_scale.sct[j],",",
                     "kth_dh=", argv$kth_closest_obs_horizontal_scale.sct[j],",",
                     "dz=", argv$vertical_scale.sct[j],",",
                     str1,
                     "/time ",round(t1a-t0a,1),attr(t1a-t0a,"unit"),"\n"))
        cat( paste0(nsus[j]," new suspect observations",str,"\n"))
        rm(str)
      } # end for k
    } # end for j
    if ( sum(nsus) <= argv$break.sct) break
  }  # end for i
  cat("+---------------------------------+\n")
  #
  return(dqcflag)
}

