#+ spatial consistency test for dichotomous (yes/no) variables with background fields
sct_fg_dual_r <- function( argv, 
                           data,
                           x,
                           y,
                           z, 
                           dqcflag) {
#------------------------------------------------------------------------------
#==============================================================================

  ndata       <- length( data$lat)

  nfin        <- length( argv$input.files)

  cat( paste0( "sct_fg_dual with external first guess (code=", argv$code.sct_fg_dual, ")\n"))

  # number of observation providers
  M <- nfin

  # number of tests
  if ( ( N <- length( argv$fglab.sct_fg_dual)) == 0) return( FALSE)

  # number of background fields
  if ( ( B <- length( fg_env$fg)) == 0) return( FALSE)

  debug <- FALSE

  if ( length( argv$doit.sct_fg_dual) != M) 
    argv$doit.sct_fg_dual <- rep( argv$doit.sct_fg_dual[1], M)
  if ( length( argv$prio.sct_fg_dual) != M) 
    argv$prio.sct_fg_dual <- rep( argv$prio.sct_fg_dual[1], M)

  if ( length( argv$event_thresholds.sct_fg_dual) != (M*N)) 
    argv$event_thresholds.sct_fg_dual <- rep( argv$event_thresholds.sct_fg_dual[1], N*M)
  if ( length( argv$conditions.sct_fg_dual) != N) 
    argv$conditions.sct_fg_dual <- rep( argv$conditions.sct_fg_dual[1], N)

  if ( length( argv$fg_event_thresholds.sct_fg_dual) != N) 
    argv$fg_event_thresholds.sct_fg_dual <- rep( argv$fg_event_thresholds.sct_fg_dual[1], N)

  if ( length( argv$thr_relinfo.sct_fg_dual) != (M*N)) 
    argv$thr_relinfo.sct_fg_dual <- rep( argv$thr_relinfo.sct_fg_dual[1], N*M)

  if ( length( argv$num_min_outer.sct_fg_dual) != N)
    argv$num_min_outer.sct_fg_dual <- rep( argv$num_min_outer.sct_fg_dual[1], N)
  if ( length( argv$num_max_outer.sct_fg_dual) != N)
    argv$num_max_outer.sct_fg_dual <- rep( argv$num_max_outer.sct_fg_dual[1], N)
  if ( length( argv$inner_radius.sct_fg_dual) != N) 
    argv$inner_radius.sct_fg_dual <- rep( argv$inner_radius.sct_fg_dual[1], N)
  if ( length( argv$outer_radius.sct_fg_dual) != N) 
    argv$outer_radius.sct_fg_dual <- rep( argv$outer_radius.sct_fg_dual[1], N)
  if ( length( argv$min_horizontal_scale.sct_fg_dual) != N) 
    argv$min_horizontal_scale.sct_fg_dual <- rep( argv$min_horizontal_scale.sct_fg_dual[1], N)
  if ( length( argv$max_horizontal_scale.sct_fg_dual) != N) 
    argv$max_horizontal_scale.sct_fg_dual <- rep( argv$max_horizontal_scale.sct_fg_dual[1], N)
  if ( length( argv$kth_closest_obs_horizontal_scale.sct_fg_dual) != N) 
    argv$kth_closest_obs_horizontal_scale.sct_fg_dual <- rep( argv$kth_closest_obs_horizontal_scale.sct_fg_dual[1], N)
  if ( length( argv$vertical_scale.sct_fg_dual) != N) 
    argv$vertical_scale.sct_fg_dual <- rep( argv$vertical_scale.sct_fg_dual[1], N)

  #............................................................................
  # set doit/prio vectors
  doit <- vector( length=ndata, mode="numeric"); doit[]<-NA
  prio <- vector( length=ndata, mode="numeric"); prio[]<-NA
  for (f in 1:M) {
    ix <- which( data$prid == argv$prid[f])
    if ( length(ix) == 0) next
    doit[ix] <- argv$doit.sct_fg_dual[f]
    prio[ix] <- argv$prio.sct_fg_dual[f]
  }
  prio_unique <- sort( unique( prio, na.rm=T), decreasing=T)
  rm( ix)

  #............................................................................
  # prepare fg

  # loop over background fields 
  for (ifg in 1:B) {

    if (!is.null( fg_env$fg[[ifg]]$r_main)) {

      # prepare for the accumulation of flags
      dqcflag_acc <- rep( 0, ndata)

      if ( ( Bi <- length( ixfglab <- which( argv$fglab.sct_fg_dual == ifg))) == 0 ) next

      nens <- nlayers(fg_env$fg[[ifg]]$r_main)

      nsus <- vector( mode="numeric", length=Bi)

      # loop over ensembles 
      for (ens in 1:nens) {
 
        dqcflag_f   <- dqcflag

        if ( class( fg_env$fg[[ifg]]$r_main) == "RasterLayer") {
          rfg <- fg_env$fg[[ifg]]$r_main
        } else {
          rfg <- raster( fg_env$fg[[ifg]]$r_main, ens)
        }

        if ( !is.null( fg_env$fg[[ifg]]$r_aux1)) rfgdem <- fg_env$fg[[ifg]]$r_aux1

        fg_x    <- integer(0)
        fg_y    <- integer(0)
        fg_lat  <- integer(0)
        fg_lon  <- integer(0)
        fg_z    <- integer(0)
        fg_val  <- integer(0)
        t0a <- Sys.time()
        if ( is.null(rfg)) boom( "ERROR in sct_fg_dual")
        dfg <- getValues( rfg)
        if ( !is.null( rfgdem)) { dfgdem <- getValues( rfgdem) } else
                                { dfgdem <- rep( 0, length( dfg)) }

        # get coordinates into CRS 
        ixx <- which( !is.na(dfg) & !is.na(dfgdem))
        if ( length( ixx) > 0) {
          fgxy         <- as.data.frame( xyFromCell( rfg, ixx))
          names( fgxy) <- c( "x", "y")
          coordinates( fgxy) <- c( "x", "y")
          proj4string( fgxy) <- CRS( fg_env$fg[[ifg]]$main.proj4)
          fgxy_transf <- as.data.frame( spTransform( fgxy, crs(argv$proj4_where_dqc_is_done)))
          fg_x        <- fgxy_transf[,1]
          fg_y        <- fgxy_transf[,2]
          fgll_transf <- as.data.frame( spTransform( fgxy, crs("+proj=longlat +datum=WGS84")))
          fg_lon      <- fgll_transf[,1]
          fg_lat      <- fgll_transf[,2]
          fg_z        <- dfgdem[ixx]
          fg_val      <- dfg[ixx]
          rm( fgxy, fgxy_transf, fgll_transf)
        }

        # test
        flagok_pre <- !is.na( data$lat) & !is.na( data$lon) & !is.na( data$value) &
                      doit != 0 & !is.na( prio)

        for (i in 1:argv$i.sct_fg_dual) {

          nsus[]<-0

          for (jj in 1:Bi) {

            j <- ixfglab[jj]

            t     <- vector( length=ndata, mode="numeric"); t[]<-NA
            r_eve <- vector( length=ndata, mode="numeric"); r_eve[]<-NA

            for (f in 1:M) {
              ix <- which( data$prid == argv$prid[f])
              if ( length(ix) == 0) next
              t[ix]     <- argv$thr_relinfo.sct_fg_dual[(j-1)*M+f]
              r_eve[ix] <- argv$event_thresholds.sct_fg_dual[(j-1)*M+f]
            }
            rm( ix)

            first_k <- T
            for (k in 1:length(prio_unique)) {
              # use only (probably) good observations with doit!=0
              flag_aux<-( (is.na(dqcflag_f) | dqcflag_f==argv$code.keep) &
                          flagok_pre)

              if (argv$variable == "T") flag_aux <- flag_aux & !is.na(z) 

              ix <- which( flag_aux); rm( flag_aux)
              t0a <- Sys.time()

              obsToCheck_n   <- length( ix)
              obsToCheck_fgn <- length( fg_lon)

              if (obsToCheck_n>0) {
                # define global 1D vector used in statSpat (1D for fast access)
                obsToCheck_lon  <- c( data$lon[ix], fg_lon)
                obsToCheck_lat  <- c( data$lat[ix], fg_lat)
                obsToCheck_z    <- c(        z[ix], fg_z)
                obsToCheck_val  <- c( data$value[ix], fg_val)
                obsToCheck_t    <- c( t[ix], rep( 0, obsToCheck_fgn))
                obsToCheck_r    <- c( r_eve[ix],
                                      rep( argv$fg_event_thresholds.sct_fg_dual[j],
                                           obsToCheck_fgn))
                obsToCheck_chk  <- rep( 0, obsToCheck_n)
                # check only those observations with priorities geq than this
                obsToCheck_chk[prio[ix]>=prio_unique[k]] <- 1
                # do not check the fg
                obsToCheck_chk  <- c( obsToCheck_chk, rep(0,obsToCheck_fgn))

                cat("\n") 

                res <- sct_dual(
                            points = Points( obsToCheck_lat, 
                                             obsToCheck_lon, 
                                             obsToCheck_z),
                            obsToCheck_val,
                            obsToCheck_chk,
                            obsToCheck_r,
                            argv$conditions.sct_fg_dual[j],
                            argv$num_min_outer.sct_fg_dual[j],
                            argv$num_max_outer.sct_fg_dual[j],
                            argv$inner_radius.sct_fg_dual[j],
                            argv$outer_radius.sct_fg_dual[j],
                            100,
                            argv$min_horizontal_scale.sct_fg_dual[j],
                            argv$max_horizontal_scale.sct_fg_dual[j],
                            argv$kth_closest_obs_horizontal_scale.sct_fg_dual[j],
                            argv$vertical_scale.sct_fg_dual[j],
                            obsToCheck_t,
                            debug)
       
                flag <- res

                # suspect if: 
                sus<-which( flag[1:obsToCheck_n] == 1 &
                            is.na(dqcflag_f[ix]) &
                            doit[ix]==1 )

                # set dqcflag
                if (length(sus)>0) dqcflag_f[ix[sus]] <- argv$code.sct_fg_dual

              } else {
                cat( "no valid observations left, no sct_fg_dual\n")
              }
              nsus[jj] <- ifelse( exists("sus"), length(sus), 0)
              t1a  <- Sys.time()
              str  <- " (#TOT "
              str1 <- ""
              for (f in 1:M) {
                if (f>1) str<-paste0(str,"; ")
                str <- paste0( str, "prid", argv$prid[f], "=", 
                        length( which( dqcflag_f==argv$code.sct_fg_dual & 
                                       data$prid==argv$prid[f])))
                str1 <- paste0( str1, "::prid ", argv$prid[f]," ",
                                      "prio=", argv$prio.sct_fg_dual[f],",",
                                      "t=", argv$thr_relinfo.sct_fg_dual[(j-1)*M+f],",",
                                      "r=", argv$event_thresholds.sct_fg_dual[(j-1)*M+f])
              }
              str<-paste0(str,")")
              cat( paste0( "++++>> SCT_fg_dual-TMP ifg=",ifg,"of",B,
                           "/ens=", ens, "of", nens,
                           "/iteration=", i,
                           "/test=", j,
                           "/prio>=", prio_unique[k],
                           "/dqc param:",
                           "inner_rad=", argv$inner_radius.sct_fg_dual[j],",",
                           "outer_rad=", argv$outer_radius.sct_fg_dual[j],",",
                           "nmin_outer=", argv$num_min_outer.sct_fg_dual[j],",",
                           "nmax_outer=", argv$num_max_outer.sct_fg_dual[j],",",
                           "min_hscale=", argv$min_horizontal_scale.sct_fg_dual[j],",",
                           "max_hscale=", argv$max_horizontal_scale.sct_fg_dual[j],",",
                           "kth_dh=", argv$kth_closest_obs_horizontal_scale.sct_fg_dual[j],",",
                           "dz=", argv$vertical_scale.sct_fg_dual[j],",",
                           "cond=", argv$conditions.sct_fg_dual[j],",",
                           "fg_r=", argv$fg_event_thresholds.sct_fg_dual[j],",",
                           str1,
                           "/time ",round(t1a-t0a,1),attr(t1a-t0a,"unit"),"\n"))
              cat( paste0( nsus[jj]," new suspect observations", str, "\n"))
              rm(str)

            } # end for k - priorities

          } # end for j - number of test based on the j-th fg field

          if ( sum(nsus) <= argv$break.sct_fg_dual) { cat("++++>> BREAK <<++++\n"); break}

        }   # end for i - iterations i.sct_fg_dual

        # accumulate dqcflag
        dqcflag_acc[which(dqcflag_f==argv$code.sct_fg_dual)] <- dqcflag_acc[which(dqcflag_f==argv$code.sct_fg_dual)] + 1

      } # end for ens

      # set dqcflag
      sus  <- which( dqcflag_acc >= (nens/2))
      nsus <- length( sus)
      if ( nsus > 0) dqcflag[sus] <- argv$code.sct_fg_dual
      cat( paste0( "\n++++>> SCT_fg_dual-DEF ifg=", ifg,", total number of suspect observations=", nsus, "\n"))

    } # end if

  } # end for first-guesses

  cat("+---------------------------------+\n")

  #
  return(dqcflag)

}

