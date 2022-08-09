#+ spatial consistency test, using background fields
sct_fg_resistant <- function( argv, 
                              data,
                              x,
                              y,
                              z, 
                              dqcflag) {
#------------------------------------------------------------------------------
#==============================================================================

  require(RANN)

  cat( paste0( "sct with external first guess (code=", argv$code.sct_fg, ")\n"))

  nfin  <- length( argv$input.files)

  ndata <- length( data$lat)

  # number of observation providers
  M <- nfin

  # number of tests
  if ( ( N <- length( argv$fglab.sct_fg)) == 0) return( FALSE)

  # number of background fields
  if ( ( B <- length( fg_env$fg)) == 0) return( FALSE)

  debug <- FALSE
  background_elab_type <- "External"
  num_min_prof  <- -999
  min_elev_diff <- 1

  if ( is.null( argv$transf.sct_fg)) argv$transf.sct_fg <- F

  if ( length( argv$doit.sct_fg) != M) 
    argv$doit.sct_fg <- rep( argv$doit.sct_fg[1], M)
  if ( length( argv$prio.sct_fg) != M) 
    argv$prio.sct_fg <- rep( argv$prio.sct_fg[1], M)

  if ( length( argv$tpos.sct_fg) != (M*N)) 
    argv$tpos.sct_fg <- rep( argv$tpos.sct_fg[1], N*M)
  if ( length( argv$tneg.sct_fg) != (M*N))
    argv$tneg.sct_fg <- rep( argv$tneg.sct_fg[1], N*M)
  if ( length( argv$eps2.sct_fg) != (M*N)) 
    argv$eps2.sct_fg <- rep( argv$eps2.sct_fg[1], N*M)

  if ( length( argv$num_min_outer.sct_fg) != N)
    argv$num_min_outer.sct_fg <- rep( argv$num_min_outer.sct_fg[1], N)
  if ( length( argv$num_max_outer.sct_fg) != N)
    argv$num_max_outer.sct_fg <- rep( argv$num_max_outer.sct_fg[1], N)
  if ( length( argv$num_max_aggn.sct_fg) != N) 
    argv$num_max_aggn.sct_fg <- rep( argv$num_max_aggn.sct_fg[1], N)
  if ( length( argv$aggn_radius.sct_fg) != N) 
    argv$aggn_radius.sct_fg <- rep( argv$aggn_radius.sct_fg[1], N)
  if ( length( argv$outer_radius.sct_fg) != N) 
    argv$outer_radius.sct_fg <- rep( argv$outer_radius.sct_fg[1], N)
  if ( length( argv$min_horizontal_scale.sct_fg) != N) 
    argv$min_horizontal_scale.sct_fg <- rep( argv$min_horizontal_scale.sct_fg[1], N)
  if ( length( argv$max_horizontal_scale.sct_fg) != N) 
    argv$max_horizontal_scale.sct_fg <- rep( argv$max_horizontal_scale.sct_fg[1], N)
  if ( length( argv$kth_closest_obs_horizontal_scale.sct_fg) != N) 
    argv$kth_closest_obs_horizontal_scale.sct_fg <- rep( argv$kth_closest_obs_horizontal_scale.sct_fg[1], N)
  if ( length( argv$vertical_scale.sct_fg) != N) 
    argv$vertical_scale.sct_fg <- rep( argv$vertical_scale.sct_fg[1], N)

  #............................................................................
  # set doit/prio vectors
  doit <- vector( length=ndata, mode="numeric"); doit[]<-NA
  prio <- vector( length=ndata, mode="numeric"); prio[]<-NA
  for (f in 1:M) {
    ix <- which( data$prid == argv$prid[f])
    if ( length(ix) == 0) next
    doit[ix] <- argv$doit.sct_fg[f]
    prio[ix] <- argv$prio.sct_fg[f]
  }
  prio_unique <- sort( unique( prio, na.rm=T), decreasing=T)
  rm( ix)

  #............................................................................
  # prepare vectors of valid and admissible values
  if (argv$variable == "T") {
    values_mina <- data$value - argv$a_delta.sct_fg
    values_maxa <- data$value + argv$a_delta.sct_fg
    values_minv <- data$value - argv$v_delta.sct_fg
    values_maxv <- data$value + argv$v_delta.sct_fg
  } else if (argv$variable == "RR") {
    values_mina <- pmin( pmax( data$value - argv$a_delta.sct_fg, 0), 
                         pmax( data$value - argv$a_fact.sct_fg * data$value, 0))
    values_maxa <- pmax( data$value + argv$a_delta.sct_fg,
                         data$value + argv$a_fact.sct_fg * data$value)
    values_minv <- pmin( pmax( data$value - argv$v_delta.sct_fg, 0), 
                         pmax( data$value - argv$v_fact.sct_fg * data$value, 0))
    values_maxv <- pmax( data$value + argv$v_delta.sct_fg,
                         data$value + argv$v_fact.sct_fg * data$value)
  }

  # data transformation
  if (argv$transf.sct_fg) {
    values_mina <- boxcox( x=values_mina, lambda=argv$boxcox.lambda)
    values_maxa <- boxcox( x=values_maxa, lambda=argv$boxcox.lambda)
    values_minv <- boxcox( x=values_minv, lambda=argv$boxcox.lambda)
    values_maxv <- boxcox( x=values_maxv, lambda=argv$boxcox.lambda)
    data$value <- boxcox( x=data$value, lambda=argv$boxcox.lambda)
  }

  #............................................................................
  # prepare fg

  # loop over background fields 
  for (ifg in 1:B) {

    if (!is.null( fg_env$fg[[ifg]]$r_main)) {

      # prepare for the accumulation of flags
      dqcflag_acc <- rep( 0, ndata)

      if ( ( Bi <- length( ixfglab <- which( argv$fglab.sct_fg == ifg))) == 0 ) next

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
        if ( is.null(rfg)) boom( "ERROR in SCT_fg")
        dfg <- getValues( rfg)
        if ( !is.null( rfgdem)) { dfgdem <- getValues( rfgdem) } else
                                { dfgdem <- rep( 0, length( dfg)) }

        # data transformation
        if (argv$transf.sct_fg) dfg <- boxcox( x=dfg, lambda=argv$boxcox.lambda)

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

        # functions
        spatagg <- function(i) { mean( fg_val[nnix[i,]]) }
        demspatagg <- function(i) { mean( fg_z[nnix[i,]]) }

        # test
        flagok_pre <- !is.na( x) & !is.na( y) & !is.na( data$value) &
                      doit != 0 & !is.na( prio)

        for (i in 1:argv$i.sct_fg) {

          nsus[]<-0

          for (jj in 1:Bi) {

            j <- ixfglab[jj]

            tpos <- vector( length=ndata, mode="numeric"); tpos[]<-NA
            tneg <- vector( length=ndata, mode="numeric"); tneg[]<-NA
            eps2 <- vector( length=ndata, mode="numeric"); eps2[]<-NA
            for (f in 1:M) {
              ix <- which( data$prid == argv$prid[f])
              if ( length(ix) == 0) next
              tpos[ix] <- argv$tpos.sct_fg[(j-1)*M+f]
              tneg[ix] <- argv$tneg.sct_fg[(j-1)*M+f]
              eps2[ix] <- argv$eps2.sct_fg[(j-1)*M+f]
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

              obsToCheck_n<-length(ix)

              if (obsToCheck_n>0) {
                # define global 1D vector used in statSpat (1D for fast access)
                obsToCheck_lon <- data$lon[ix]
                obsToCheck_lat <- data$lat[ix]
                obsToCheck_z   <- z[ix]
                obsToCheck_val <- data$value[ix]
                obsToCheck_mina <- values_mina[ix]
                obsToCheck_maxa <- values_maxa[ix]
                obsToCheck_minv <- values_minv[ix]
                obsToCheck_maxv <- values_maxv[ix]
                obsToCheck_tpos <- tpos[ix]
                obsToCheck_tneg <- tneg[ix]
                obsToCheck_eps2 <- eps2[ix]
                obsToCheck_chk <- rep( 0, obsToCheck_n)
                # check only those observations with priorities geq than this
                obsToCheck_chk[prio[ix]>=prio_unique[k]] <- 1

                if (first_k) {
                  nn2 <- nn2( cbind( fg_x, fg_y), 
                              query = cbind( x[ix], y[ix]), 
                              k = argv$num_max_aggn.sct_fg[j], 
                              searchtype = "radius", radius = argv$aggn_radius.sct_fg[j])
                  nnix <- nn2[[1]]
                  if (!is.na(argv$cores)) {
                    aux <- t( mcmapply( spatagg,
                                        1:length(ix),
                                        mc.cores = argv$cores,
                                        SIMPLIFY = T))
                  # no-multicores
                  } else {
                    aux <- t( mapply( spatagg,
                                      1:length(ix),
                                      SIMPLIFY = T))
                  }
                  fgstat_at_opoint   <- vector( mode="numeric", length=ndata)
                  fgstat_at_opoint[] <- NA
                  fgstat_at_opoint[ix] <- aux

                  if ( !is.null( rfgdem)) {
                    if (!is.na(argv$cores)) {
                      aux <- t( mcmapply( demspatagg,
                                          1:length(ix),
                                          mc.cores = argv$cores,
                                          SIMPLIFY = T))
                    # no-multicores
                    } else {
                      aux <- t( mapply( demspatagg,
                                        1:length(ix),
                                        SIMPLIFY = T))
                    }
                    demstat_at_opoint <- vector( mode="numeric", length=ndata)
                    demstat_at_opoint[] <- NA
                    demstat_at_opoint[ix] <- aux 
                    rm( aux)
                  }

                  first_k <- F
                }

                background_values <- fgstat_at_opoint[ix]

                if ( !is.null( rfgdem)) 
                  background_values <- background_values + 
                   argv$gamma.standard * ( z[ix] - demstat_at_opoint[ix])

                cat("\n")  
                res <- sct_resistant(
                            points = Points( obsToCheck_lat, 
                                             obsToCheck_lon, 
                                             obsToCheck_z),
                            obsToCheck_val,
                            obsToCheck_chk,
                            background_values,
                            background_elab_type,
                            argv$num_min_outer.sct_fg[j],
                            argv$num_max_outer.sct_fg[j],
                            argv$inner_radius.sct_fg[j],
                            argv$outer_radius.sct_fg[j],
                            100,
                            num_min_prof,
                            min_elev_diff,
                            argv$min_horizontal_scale.sct_fg[j],
                            argv$max_horizontal_scale.sct_fg[j],
                            argv$kth_closest_obs_horizontal_scale.sct_fg[j],
                            argv$vertical_scale.sct_fg[j],
                            obsToCheck_mina,
                            obsToCheck_maxa,
                            obsToCheck_minv,
                            obsToCheck_maxv,
                            obsToCheck_eps2,
                            obsToCheck_tpos,
                            obsToCheck_tneg,
                            debug,
                            argv$basic.sct_fg)
       
                flag <- res[[1]]

                # suspect if: 
                sus<-which( flag[1:obsToCheck_n] == 1 &
                            is.na(dqcflag_f[ix]) &
                            doit[ix]==1 )

                # set dqcflag
                if (length(sus)>0) dqcflag_f[ix[sus]] <- argv$code.sct_fg

              } else {
                cat( "no valid observations left, no SCT_fg\n")
              }
              nsus[jj] <- ifelse( exists("sus"), length(sus), 0)
              t1a  <- Sys.time()
              str  <- " (#TOT "
              str1 <- ""
              for (f in 1:M) {
                if (f>1) str<-paste0(str,"; ")
                str <- paste0( str, "prid", argv$prid[f], "=", 
                        length( which( dqcflag_f==argv$code.sct_fg & 
                                       data$prid==argv$prid[f])))
                str1 <- paste0( str1, "::prid ", argv$prid[f]," ",
                                      "prio=", argv$prio.sct_fg[f],",",
                                      "eps2=", argv$eps2.sct_fg[(j-1)*M+f],",",
                                      "tpos=", argv$tpos.sct_fg[(j-1)*M+f],",",
                                      "tneg=", argv$tneg.sct_fg[(j-1)*M+f])
              }
              str<-paste0(str,")")
              cat( paste0( "++++>> SCT_fg-TMP ifg=",ifg,"of",B,
                           "/ens=", ens, "of", nens,
                           "/iteration=", i,
                           "/test=", j,
                           "/prio>=", prio_unique[k],
                           "/dqc param:",
                           "inner_rad=", argv$inner_radius.sct_fg[j],",",
                           "outer_rad=", argv$outer_radius.sct_fg[j],",",
                           "nmin_outer=", argv$num_min_outer.sct_fg[j],",",
                           "nmax_outer=", argv$num_max_outer.sct_fg[j],",",
                           "min_hscale=", argv$min_horizontal_scale.sct_fg[j],",",
                           "max_hscale=", argv$max_horizontal_scale.sct_fg[j],",",
                           "kth_dh=", argv$kth_closest_obs_horizontal_scale.sct_fg[j],",",
                           "dz=", argv$vertical_scale.sct_fg[j],",",
                           "nmax_agg=", argv$num_max_aggn.sct_fg[j],",",
                           "agg_radius=", argv$aggn_radius.sct_fg[j],",",
                           str1,
                           "/time ",round(t1a-t0a,1),attr(t1a-t0a,"unit"),"\n"))
              cat( paste0( nsus[jj]," new suspect observations", str, "\n"))
              rm(str)
            } # end for k
          } # end for j
          if ( sum(nsus) <= argv$break.sct_fg) break
        }   # end for i

        # accumulate dqcflag
        dqcflag_acc[which(dqcflag_f==argv$code.sct_fg)] <- dqcflag_acc[which(dqcflag_f==argv$code.sct_fg)] + 1

      } # end for ens

      # set dqcflag
      sus  <- which( dqcflag_acc >= (nens/2))
      nsus <- length( sus)
      if ( nsus > 0) dqcflag[sus] <- argv$code.sct_fg
      cat( paste0( "\n++++>> SCT_fg-DEF ifg=", ifg,", total number of suspect observations=",nsus,"\n"))

    } # end if

  } # end for first-guesses

  cat("+---------------------------------+\n")

  #
  return(dqcflag)

}

