#+ Read first-guess fields 
read_fgs <- function( argv, extent) {
#==============================================================================

  if (argv$verbose) cat("Read Background files\n")
  t0a<-Sys.time()

  if ( ( nfg <- length( fg_env$fg)) == 0) return( FALSE)

  for (f in 1:nfg) {

    # initializations

    fg_env$fg[[f]]$r_aux1 <- NULL
    fg_env$fg[[f]]$r_main <- NULL

    if ( is.null( fg_env$fg[[f]]$main.proj4)) 
      fg_env$fg[[f]]$main.proj4 <- ""
    if ( is.null( fg_env$fg[[f]]$main.proj4_var)) 
      fg_env$fg[[f]]$main.proj4_var <- ""
    if ( is.null( fg_env$fg[[f]]$main.proj4_att)) 
      fg_env$fg[[f]]$main.proj4_att <- ""
    if ( is.null( fg_env$fg[[f]]$main.t)) 
      fg_env$fg[[f]]$main.t <- NA
    if ( is.null( fg_env$fg[[f]]$main.acc)) 
      fg_env$fg[[f]]$main.acc <- F

    if ( is.null( fg_env$fg[[f]]$aux1.offset)) 
      fg_env$fg[[f]]$aux1.offset <- 0 
    if ( is.null( fg_env$fg[[f]]$aux1.cfact)) 
      fg_env$fg[[f]]$aux1.cfact <- 1 
    if ( is.null( fg_env$fg[[f]]$main.offset)) 
      fg_env$fg[[f]]$main.offset <- 0 
    if ( is.null( fg_env$fg[[f]]$main.cfact)) 
      fg_env$fg[[f]]$main.cfact <- 1 

    # read auxiliary file (e.g. dem for temperature)
    if ( !is.null( fg_env$fg[[f]]$aux1.file)) {
      if ( file.exists( fg_env$fg[[f]]$aux1.file)) {
        res <- get_data_from_ncfile( nc.file    = fg_env$fg[[f]]$aux1.file,
                                     nc.varname = fg_env$fg[[f]]$aux1.varname,
                                     nc.t       = fg_env$fg[[f]]$aux1.t,
                                     nc.e       = fg_env$fg[[f]]$aux1.e,
                                     topdown    = fg_env$fg[[f]]$aux1.topdown,
                                     var.dim    = list( ndim  = fg_env$fg[[f]]$aux1.ndim,
                                                        tpos  = fg_env$fg[[f]]$aux1.tpos,
                                                        epos  = fg_env$fg[[f]]$aux1.epos,
                                                        names = fg_env$fg[[f]]$aux1.dimnames),
                                     proj4      = fg_env$fg[[f]]$main.proj4,
                                     proj4_from_nc=list( var = fg_env$fg[[f]]$main.proj4_var, 
                                                         att = fg_env$fg[[f]]$main.proj4_att),
                                     xy_as_vars = F,
                                     return_raster = T,
                                     return_obsloc = F,
                                     extent      = extent,
                                     debug.file  = NA,
                                     nc.dqc_mode = "none")
        if (!is.null( res))
          fg_env$fg[[f]]$r_aux1 <- fg_env$fg[[f]]$aux1.offset + res$raster * fg_env$fg[[f]]$aux1.cfact 
        rm( res)
      }
    }

    # read the main file
    if ( !is.null( fg_env$fg[[f]]$main.file)) {
      if ( file.exists( fg_env$fg[[f]]$main.file)) {
        first <- T
        if ( is.null( fg_env$fg[[f]]$main.epos)) {
          ei <- 0
        } else {
          if ( is.null( fg_env$fg[[f]]$main.e)) {
            ei <- nc4.getDim( fg_env$fg[[f]]$main.file, 
                              varid = fg_env$fg[[f]]$main.dimnames[fg_env$fg[[f]]$main.epos])
          } else {
            ei <- fg_env$fg[[f]]$main.e
          }
        }
        for (ens in 1:length(ei)) {
          if( ei[ens] == 0) { nc_e <- NA} else { nc_e <- ei[ens]}
          res <- get_data_from_ncfile( nc.file    = fg_env$fg[[f]]$main.file,
                                       nc.varname = fg_env$fg[[f]]$main.varname,
                                       nc.t       = fg_env$fg[[f]]$main.t,
                                       nc.e       = nc_e,
                                       topdown    = fg_env$fg[[f]]$main.topdown,
                                       var.dim    = list( ndim  = fg_env$fg[[f]]$main.ndim,
                                                          tpos  = fg_env$fg[[f]]$main.tpos,
                                                          epos  = fg_env$fg[[f]]$main.epos,
                                                          names = fg_env$fg[[f]]$main.dimnames),
                                       var.de_acc = fg_env$fg[[f]]$main.acc,
                                       var.de_acc_by = "-1 hour",
                                       proj4         = fg_env$fg[[f]]$main.proj4,
                                       proj4_from_nc = list( var = fg_env$fg[[f]]$main.proj4_var, 
                                                             att = fg_env$fg[[f]]$main.proj4_att),
                                       xy_as_vars    = F,
                                       return_raster = T,
                                       return_obsloc = F,
                                       extent        = extent,
                                       debug.file    = NA,
                                       nc.dqc_mode   = "none")
          if ( !is.null( res)) {
            if ( first) {
              fg_env$fg[[f]]$r_main <- fg_env$fg[[f]]$main.offset + res$raster * fg_env$fg[[f]]$main.cfact 
            } else {
              fg_env$fg[[f]]$r_main <- stack( fg_env$fg[[f]]$r_main,
                                              fg_env$fg[[f]]$main.offset + res$raster * fg_env$fg[[f]]$main.cfact)
            }
            first <- F
          }
          rm(res)
        } # end loop over ensemble 
      } # end if main file exists 
    } # end read the main file
  } # end loop over fg files

  # print info
  if (argv$verbose) {
    t1a<-Sys.time()

    for (f in 1:nfg) {
      if (!is.null( fg_env$fg[[f]]$r_main)) {
        for (ens in 1:nlayers(fg_env$fg[[f]]$r_main)) {
          if ( class( fg_env$fg[[f]]$r_main) == "RasterLayer") {
            r <- fg_env$fg[[f]]$r_main
          } else {
            r <- raster( fg_env$fg[[f]]$r_main, ens)
          }
          if (ens == 1) {
            cat( paste( "FG", formatC( f, width=2, flag="0"),
                        "ensemble member", formatC( ens, width=2, flag="0"), ".",
                        "dim =", ncell( r), ".", 
                        "range =", round( min( getValues( r), na.rm=T), 2), 
                                   round( max( getValues( r), na.rm=T), 2), "\n"))
          } else {
            cat( paste( "      ensemble member", formatC( ens, width=2, flag="0"), ".",
                        "dim =", ncell( r), ".", 
                        "range =", round( min( getValues( r), na.rm=T), 2), 
                                   round( max( getValues( r), na.rm=T), 2), "\n"))
          }
        }
      }
      if (!is.null( fg_env$fg[[f]]$r_aux1)) {
        r <- fg_env$fg[[f]]$r_aux1
        cat( paste( "                   aux 1 .",
                    "dim =", ncell( r), ".", 
                    "range =", round( min( getValues( r), na.rm=T), 2), 
                               round( max( getValues( r), na.rm=T), 2), "\n"))
      }
    }
    cat( paste( "total time", round(t1a-t0a,1), attr(t1a-t0a,"unit"), "\n"))
    cat( "+---------------------------------+\n")
  }

  return( TRUE)
}
