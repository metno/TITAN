#+ check against first-guess fields
fgt_r <- function( argv, 
                   data,
                   x,
                   y,
                   z, 
                   dqcflag) {
#------------------------------------------------------------------------------
# aggregation neighbourhoods (aggn) 
#                      k = argv$num_max_aggn.fgt[i], 
#                      searchtype = "radius", radius = argv$aggn_radius.fgt[i])
# test circle (circle_radius)
#
# aggregate the background over aggregation neighbourhoods
# test deviations between observations and aggregated backgrounds over test circles
#==============================================================================

  require(RANN)

  cat( paste0( "fgt (code=", argv$code.fgt, ")\n"))

  nfin  <- length(argv$input.files)

  ndata <- length(data$lat)

  # number of observation providers
  M <- nfin

  # number of tests
  if ( ( N <- length( argv$fglab.fgt)) == 0) return( FALSE)

  # number of background fields
  if ( ( B <- length( fg_env$fg)) == 0) return( FALSE)


  debug <- FALSE
  background_elab_type <- "External"
  num_min_prof  <- -999
  min_elev_diff <- 1

  if ( is.null( argv$transf.fgt)) argv$transf.fgt <- F

  if ( length( argv$doit.fgt) != M) 
    argv$doit.fgt <- rep( argv$doit.fgt[1], M)
  if ( length( argv$prio.fgt) != M) 
    argv$prio.fgt <- rep( argv$prio.fgt[1], M)

  if ( length( argv$tpos.fgt) != (M*N)) 
    argv$tpos.fgt <- rep( argv$tpos.fgt[1], N*M)
  if ( length( argv$tneg.fgt) != (M*N))
    argv$tneg.fgt <- rep( argv$tneg.fgt[1], N*M)

  if ( length( argv$num_min_circle.fgt) != N)
    argv$num_min_circle.fgt <- rep( argv$num_min_circle.fgt[1], N)
  if ( length( argv$num_max_circle.fgt) != N)
    argv$num_max_circle.fgt <- rep( argv$num_max_circle.fgt[1], N)
  if ( length( argv$num_max_aggn.fgt) != N) 
    argv$num_max_aggn.fgt <- rep( argv$num_max_aggn.fgt[1], N)
  if ( length( argv$aggn_radius.fgt) != N) 
    argv$aggn_radius.fgt <- rep( argv$aggn_radius.fgt[1], N)

  #............................................................................
  # set doit/prio vectors
  doit <- vector( length=ndata, mode="numeric"); doit[]<-NA
  prio <- vector( length=ndata, mode="numeric"); prio[]<-NA
  for (f in 1:M) {
    ix <- which( data$prid == argv$prid[f])
    if ( length(ix) == 0) next
    doit[ix] <- argv$doit.fgt[f]
    prio[ix] <- argv$prio.fgt[f]
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

  # data transformation
  if (argv$transf.fgt) {
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

      if ( ( Bi <- length( ixfglab <- which( argv$fglab.fgt == ifg))) == 0 ) next

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
        if ( is.null(rfg)) boom( "ERROR in FGT")
        dfg <- getValues( rfg)
        if ( !is.null( rfgdem)) { dfgdem <- getValues( rfgdem) } else
                                { dfgdem <- rep( 0, length( dfg)) }

        # data transformation
        if (argv$transf.fgt) dfg <- boxcox( x=dfg, lambda=argv$boxcox.lambda)

        # get coordinates into CRS 
        ixx <- which( !is.na(dfg) & !is.na(dfgdem))
        if ( length( ixx) > 0) {
          fgxy         <- as.data.frame( xyFromCell( rfg, ixx))
          names( fgxy) <- c( "x", "y")
          coordinates( fgxy) <- c( "x", "y")
          proj4string( fgxy) <- CRS( fg_env$fg[[ifg]]$main.proj4)
          fgxy_transf <- as.data.frame( spTransform( fgxy, CRS=argv$proj4_where_dqc_is_done))
          fg_x        <- fgxy_transf[,1]
          fg_y        <- fgxy_transf[,2]
          fgll_transf <- as.data.frame( spTransform( fgxy, CRS="+proj=longlat +datum=WGS84"))
          fg_lon      <- fgll_transf[,1]
          fg_lat      <- fgll_transf[,2]
          fg_z        <- dfgdem[ixx]
          fg_val      <- dfg[ixx]
          rm( fgxy, fgxy_transf, fgll_transf)
        }

        # functions
        spatagg <- function(i) { val <- fg_val[nnix[i,]]; return( c( mean( val), sd( val))) }
        demspatagg <- function(i) { val <- fg_z[nnix[i,]]; mean( val) }

        # test
        for (i in 1:argv$i.fgt) {

          nsus[]<-0

          for (jj in 1:Bi) {

            j <- ixfglab[jj]

            tpos <- vector( length=ndata, mode="numeric"); tpos[]<-NA
            tneg <- vector( length=ndata, mode="numeric"); tneg[]<-NA
            for (f in 1:M) {
              ix <- which( data$prid == argv$prid[f])
              if ( length(ix) == 0) next
              tpos[ix] <- argv$tpos.fgt[(j-1)*M+f]
              tneg[ix] <- argv$tneg.fgt[(j-1)*M+f]
            }
            prio_unique <- sort( unique( prio, na.rm=T), decreasing=T)
            rm( ix)
            
            first_k <- T
            for (k in 1:length(prio_unique)) {
              # use only (probably) good observations with doit!=0
              flag_aux<-( (is.na(dqcflag_f) | dqcflag_f==argv$code.keep) &
                          !is.na(x) & !is.na(y) & !is.na(data$value) &
                          doit!=0 & !is.na(prio))
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
                obsToCheck_chk <- rep( 0, obsToCheck_n)
                # check only those observations with priorities geq than this
                obsToCheck_chk[prio[ix]>=prio_unique[k]] <- 1

                if (first_k) {
                  nn2 <- nn2( cbind( fg_x, fg_y), 
                              query = cbind( x[ix], y[ix]), 
                              k = argv$num_max_aggn.fgt[j], 
                              searchtype = "radius", radius = argv$aggn_radius.fgt[j])
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
                  fgstat_at_opoint <- array( data=NA, dim=c( ndata, 2))
                  fgstat_at_opoint[ix,1] <- aux[,1]
                  fgstat_at_opoint[ix,2] <- aux[,2]
                  rm( aux)

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
      
                background_values <- fgstat_at_opoint[ix,1]
                background_uncertainties <- fgstat_at_opoint[ix,2]

                if ( !is.null( rfgdem))
                  background_values <- background_values + 
                   argv$gamma.standard * ( z[ix] - demstat_at_opoint[ix])
 
                cat("\n") 
                res <- fgt( points = Points( obsToCheck_lat, 
                                             obsToCheck_lon, 
                                             obsToCheck_z),
                            obsToCheck_val,
                            obsToCheck_chk,
                            background_values,
                            background_uncertainties,
                            background_elab_type,
                            argv$num_min_circle.fgt[j],
                            argv$num_max_circle.fgt[j],
                            argv$circle_radius.fgt[j],
                            argv$circle_radius.fgt[j],
                            100,
                            num_min_prof,
                            min_elev_diff,
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
                            is.na(dqcflag_f[ix]) &
                            doit[ix]==1 )

                # set dqcflag
                if (length(sus)>0) dqcflag_f[ix[sus]] <- argv$code.fgt

              } else {
                cat( "no valid observations left, no FGT\n")
              }
              nsus[jj] <- ifelse( exists("sus"), length(sus), 0)
              t1a  <- Sys.time()
              str  <- " (#TOT "
              str1 <- ""
              for (f in 1:M) {
                if (f>1) str<-paste0(str,"; ")
                str <- paste0( str, "prid", argv$prid[f], "=", 
                        length( which( dqcflag_f==argv$code.fgt & 
                                       data$prid==argv$prid[f])))
                str1 <- paste0( str1, "::prid ", argv$prid[f]," ",
                                      "prio=", argv$prio.fgt[f],",",
                                      "tpos=", argv$tpos.fgt[(j-1)*M+f],",",
                                      "tneg=", argv$tneg.fgt[(j-1)*M+f])
              }
              str<-paste0(str,")")
              cat( paste0( "++++>> FGT-TMP ifg=",ifg,"of",B,
                           "/ens=", ens, "of", nens,
                           "/iteration=", i,
                           "/test=", j,
                           "/prio>=", prio_unique[k],
                           "/dqc param:",
                           "circle_rad=", argv$circle_radius.fgt[j],",",
                           str1,
                           "/time ",round(t1a-t0a,1),attr(t1a-t0a,"unit"),"\n"))
              cat( paste0( nsus[jj]," new suspect observations", str, "\n"))
              rm(str)
            } # end for k
          } # end for j
          if ( sum(nsus) <= argv$break.fgt) { cat("++++>> BREAK <<++++\n"); break}
        }   # end for i

        # accumulate dqcflag
        dqcflag_acc[which(dqcflag_f==argv$code.fgt)] <- dqcflag_acc[which(dqcflag_f==argv$code.fgt)] + 1

      } # end for ens

      # set dqcflag
      sus  <- which( dqcflag_acc >= (nens/2))
      nsus <- length( sus)
      if ( nsus > 0) dqcflag[sus] <- argv$code.fgt
      cat( paste0( "\n++++>> FGT-DEF ifg=", f, ", total number of suspect observations=", nsus, "\n"))

    } # end if
  } # end for first-guesses
  cat("+---------------------------------+\n")
  #
  return(dqcflag)
}

