#+ read input data (metadata and observations)
read_data_to_check <- function( argv)
#==============================================================================
# input arguments:
# argv. list with the input arguments (command line and/or config file)
# nfin. number of input files
# 
# output values:
# list of variables
# data. data frame: lat, lon, elev, value, prid
# dqcflag. numeric vector. data quality control flags
# z. numeric vector. elevation
# sctpog. numeric vector. spatial consistency test, probability of gross-error
# corep. numeric vector. coefficient of representativeness
# varidx. auxiliary variable used to write output
# varidx.opt. auxiliary variable used to write output
# dataopt. auxiliary variable used to write output
# extent. numeric vector. lonmin, lonmax, latmin, latmax
#==============================================================================
{

  nfin <- length( argv$input.files)

  first <- T

  # loop over input files (i.e. observation providers)
  for (f in 1:nfin) {
    
    if ( !file.exists(argv$input.files[f])) next 

    datain <- read.table( file   = argv$input.files[f],
                          header = T,
                          sep    = argv$separator[f],
                          stringsAsFactors = F,
                          strip.white      = T)

    # if elev is not present then create a fake one
    varidxtmp <- match( argv$varname.elev[f], names(datain))
    if (is.na(varidxtmp)) {
      argv$varname.elev[f] <- "elev"
      if (argv$elev_not_used) {
        datain$elev <- rep(  0, length = length( datain[[1]]))
      } else {
        datain$elev <- rep( NA, length = length( datain[[1]]))
      }
    }
    rm(varidxtmp)

    # varidx is used also in the output session
    varidxtmp<-match( c( argv$varname.lat[f],
                         argv$varname.lon[f],
                         argv$varname.elev[f],
                         argv$varname.value[f]),
                      names(datain) )
    if (any(is.na(varidxtmp))) {
      print("ERROR in the specification of the variable names")
      print(paste(" latitude=",argv$varname.lat[f]))
      print(paste("longitude=",argv$varname.lon[f]))
      print(paste("elevation=",argv$varname.elev[f]))
      print(paste("    value=",argv$varname.value[f]))
      print("header of input file:")
      print(argv$input.files[f])
      print(names(datain))
  #    boom()
      next
    }

    # build up the temporary data structure
    datatmp <- data.frame( datain[,varidxtmp])

    names(datatmp) <- c( "lat", "lon", "elev", "value")

    datatmp$lat <- suppressWarnings( as.numeric( datatmp$lat))
    datatmp$lon <- suppressWarnings( as.numeric( datatmp$lon))

    if (argv$elev_not_used) {
      datatmp$elev <- rep( 0, length(datatmp$lon))
    } else {
      datatmp$elev <- suppressWarnings( as.numeric( datatmp$elev))
    }
    auxz <- suppressWarnings( as.numeric( datatmp$elev))

    datatmp$value <- suppressWarnings(
      argv$input.offset[f] + argv$input.cfact[f] * as.numeric(datatmp$value))

    ndatatmp <- length( datatmp$lat)

    if (ndatatmp==0) next

    # set provider id
    datatmp$prid <- as.numeric( rep(argv$prid[f], ndatatmp))

    aux <- rep( NA, length=ndatatmp)
    # blacklist
    if ( any(!is.na(argv$blacklist.idx)) & 
         any(argv$blacklist.fidx==argv$prid[f])) {
      aux[argv$blacklist.idx[which(argv$blacklist.fidx==argv$prid[f])]] <- argv$code.black
    }
    if (any(!is.na(argv$blacklist.lat)) & 
        any(argv$blacklist.fll==argv$prid[f])) {
      eps <- 0.000001
      black_lat <- argv$blacklist.lat[which(argv$blacklist.fll==argv$prid[f])]
      black_lon <- argv$blacklist.lon[which(argv$blacklist.fll==argv$prid[f])]
      for (ii in 1:length(black_lat)) 
        aux[which( abs(datatmp$lon-black_lon[ii]) <= eps & abs(datatmp$lat-black_lat[ii])<eps)] <- argv$code.black
    }
    # keep-list
    if (any(!is.na(argv$keeplist.idx)) & 
        any(argv$keeplist.fidx==argv$prid[f])) {
      aux[argv$keeplist.idx[which(argv$keeplist.fidx==argv$prid[f])]]<-argv$code.keep  
    }
    if (any(!is.na(argv$keeplist.lat)) & 
        any(argv$keeplist.fll==argv$prid[f])) {
      eps <- 0.000001
      keep_lat <- argv$keeplist.lat[which(argv$keeplist.fll==argv$prid[f])]
      keep_lon <- argv$keeplist.lon[which(argv$keeplist.fll==argv$prid[f])]
      for (ii in 1:length(keep_lat)) 
        aux[which( abs(datatmp$lon-keep_lon[ii]) <= eps & abs(datatmp$lat-keep_lat[ii])<eps)] <- argv$code.keep
    }

    # ensure no duplicates
    if (argv$no_duplicates) {
      if (first) {
        dup_aux <- 0
        datacheck_lat  <- datatmp$lat
        datacheck_lon  <- datatmp$lon
        datacheck_elev <- datatmp$elev
      } else {
        dup_aux <- length(data$lat) 
        datacheck_lat  <- c( data$lat,  datatmp$lat)  
        datacheck_lon  <- c( data$lon,  datatmp$lon)
        datacheck_elev <- c( data$elev, datatmp$elev)
      }
      is_dup <- duplicate_check( points = Points( datacheck_lat,
                                                  datacheck_lon,
                                                  datacheck_elev),
                                 argv$no_duplicates_radius,
                                 argv$no_duplicates_vertical_range)
      ix_nodup <- which( is_dup == 0 & (1:length(is_dup)) > dup_aux) - dup_aux
      rm( is_dup, datacheck_lat, datacheck_lon, datacheck_elev, dup_aux)
    } else {
      ix_nodup <- 1:ndatatmp
    }

    # create and update the definitive data structure 
    ndatatmp<-length(ix_nodup)
    
    if ( ndatatmp > 0) {

      # datatmp$ lat lon elev value prid
      datatmp <- data.frame( lat   = suppressWarnings( as.numeric( datatmp$lat[ix_nodup])),
                             lon   = suppressWarnings( as.numeric( datatmp$lon[ix_nodup])),
                             elev  = suppressWarnings( as.numeric( datatmp$elev[ix_nodup])),
                             value = suppressWarnings( as.numeric( datatmp$value[ix_nodup])),
                             prid  = suppressWarnings( as.numeric( datatmp$prid[ix_nodup])))
      if ( first) {
        varidx  <- varidxtmp
        data    <- datatmp
        z       <- auxz[ix_nodup]
        dqcflag <- aux[ix_nodup]
        sctpog  <- rep(NA,length=ndatatmp)
        corep   <- rep(NA,length=ndatatmp)
        if ( any( !is.na( argv$varname.opt))) {
          # varidx.opt is used in the output session
          varidx.opt <- match( argv$varname.opt, names(datain))
          dataopt    <- as.data.frame( array(data=NA,
                                             dim=c(ndatatmp,length(argv$varname.opt))))
          names(dataopt) <- argv$varname.opt
          if ( any( !is.na( varidx.opt)))
            dataopt <- datain[ix_nodup,varidx.opt[which(!is.na(varidx.opt))],drop=F]
        }
        first <- F

      } else {

        data    <- rbind( data, datatmp)
        dqcflag <- c(  dqcflag, aux[ix_nodup])
        z       <- c(        z, auxz[ix_nodup])
        sctpog  <- c( sctpog, rep( NA, length=ndatatmp))
        corep   <- c( corep,  rep( NA, length=ndatatmp))

        # dataopt are auxiliary data that are not used in titanlib
        if ( any( !is.na( argv$varname.opt)) ) {
          varidx.opt.check <- match( argv$varname.opt, names(datain))
          if ( any(!is.na(varidx.opt.check) & is.na(varidx.opt)) ) {
            ixopt <- which( !is.na( varidx.opt.check) & is.na( varidx.opt))
            for (iopt in ixopt) {
              if ( varidx.opt.check[iopt] %in% varidx.opt) {
                varidx.opt[iopt] <- max(varidx.opt,na.rm=T)+1
              } else { 
                varidx.opt[iopt] <- varidx.opt.check[iopt]
              }
            }
            rm(ixopt,iopt)
          }
          dataopttmp <- as.data.frame( array(data=NA,
                                       dim=c(ndatatmp,length(argv$varname.opt))))
          names(dataopttmp)<-argv$varname.opt
          if (any(!is.na(varidx.opt.check)))
            dataopttmp<-datain[ix_nodup,varidx.opt.check[which(!is.na(varidx.opt.check))],
                               drop=F]
          dataopt<-rbind(dataopt,
                         dataopttmp)
          rm(dataopttmp)
        }
      }
    } # create and update the definitive data structure

    if (exists("ix_nodup")) rm(ix_nodup)
    if (exists("varidxtmp")) rm(varidxtmp)

    if ( any(!is.na( argv$varname.opt)) & any(!is.na(argv$blacklist.file_with_sourceIds))) {
      ix_souid <- which( names(datain) == argv$blacklist.file_sourceIds_varname)
      if ( length( ix_souid > 0)) {
        for (ii in 1:length(argv$blacklist.file_with_sourceIds)) {
          if (!file.exists(argv$blacklist.file_with_sourceIds[ii])) next
          ids_to_blacklist <- read.table( file=argv$blacklist.file_with_sourceIds[ii], header=F, sep=";", stringsAsFactors=F, strip.white=T)[,1]
          ixx <- which( dataopt[,ix_souid] %in% ids_to_blacklist)
          if ( length(ixx) > 0) dqcflag[ixx] <- argv$code.black 
        }
      }
    }

  } # END loop over input files (i.e. observation providers)

  rm( datatmp, datain, auxz, aux)

  ndata<-length(data$lat)
  if (ndata==0) {
    print("input file is empty")
    quit(status=0)
  }
  #
  # set domain extent for SCT and metadata tests
  argv$lonmin<-strings_to_numbers(strings=argv$lonmin)
  argv$lonmax<-strings_to_numbers(strings=argv$lonmax)
  argv$latmin<-strings_to_numbers(strings=argv$latmin)
  argv$latmax<-strings_to_numbers(strings=argv$latmax)
  if (argv$dqc_inbox_only) {
    extent_lonmin<-argv$lonmin
    extent_lonmax<-argv$lonmax
    extent_latmin<-argv$latmin
    extent_latmax<-argv$latmax
  } else {
    extent_lonmin<-min(data$lon,na.rm=T)
    extent_lonmax<-max(data$lon,na.rm=T)
    extent_latmin<-min(data$lat,na.rm=T)
    extent_latmax<-max(data$lat,na.rm=T)
  }

  if (argv$verbose | argv$debug) {
    print(paste("number of observations=",ndata))
    if ( any(!is.na(dqcflag) & dqcflag==argv$code.black))
      print(paste("number of blacklisted observations=",
            length(which(dqcflag==argv$code.black))) )
    if (any(!is.na(argv$keeplist.idx)) | any(!is.na(argv$keeplist.lat)))
      print(paste("number of keeplisted  observations=",
            length(which(dqcflag==argv$code.keep))) )
    if (nfin>1) {
      for (f in 1:nfin) { 
        print(paste("  number of observations provider",argv$prid[f],"=",
              length(which(data$prid==argv$prid[f]))))
        if (any(!is.na(argv$blacklist.idx)) | any(!is.na(argv$blacklist.lat)))
          print(paste("  number of blacklisted observations provider",
                argv$prid[f],"=",
                length(which(data$prid==argv$prid[f] & dqcflag==argv$code.black))) )
        if (any(!is.na(argv$keeplist.idx)) | any(!is.na(argv$keeplist.lat)))
          print(paste("  number of keeplisted  observations provider",
                argv$prid[f],"=",
                length(which(data$prid==argv$prid[f] & dqcflag==argv$code.keep))) )
      }
    }
    print(paste("extension of the domain considered (xmin,xmax,ymin,ymax)=",
                 round(extent_lonmin,6),",",round(extent_lonmax,6),",",
                 round(extent_latmin,6),",",round(extent_latmax,6)))
    print("+---------------------------------+")
  }
  #
  if (!exists("varidx.opt")) varidx.opt <- NULL
  if (!exists("dataopt"))    dataopt    <- NULL
  return( list( data    = data, 
                dqcflag = as.integer( dqcflag),
                z       = as.numeric( z),
                sctpog  = as.numeric( sctpog),
                corep   = as.numeric( corep),
                varidx  = varidx,
                varidx.opt  = varidx.opt,
                dataopt = dataopt,
                extent  = c( extent_lonmin, extent_lonmax, 
                             extent_latmin, extent_latmax)))
}
