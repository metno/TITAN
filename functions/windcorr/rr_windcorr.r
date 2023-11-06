#+ Correction for the wind-undercatch of precipitation
rr_windcorr <- function( argv, data, z, dqcflag, t2m=NULL) {
#==============================================================================
  cat( "Correction for the wind-undercatch of precipitation\n")

  if ( !file.exists(argv$t2m.file) | !file.exists(argv$t2m.demfile) | 
       !file.exists(argv$wind.file)) {
    cat( "File not found - Correction not possible\n")
    return( data)
  }

  if (is.null(t2m)) {
    #t2m
    t2m.offset   <- strings_to_numbers( strings = argv$t2m.offset,
                                        default = 0,
                                        neg     = argv$t2m.negoffset)
    t2m.cfact    <- strings_to_numbers( strings = argv$t2m.cfact,
                                        default = 0,
                                        neg     = argv$t2m.negcfact)
    t2m.demoffset <- strings_to_numbers( strings = argv$t2m.demoffset,
                                         default = 0,
                                         neg     = argv$t2m.demnegoffset)
    t2m.demcfact  <- strings_to_numbers( strings = argv$t2m.demcfact,
                                         default = 0,
                                         neg     = argv$t2m.demnegcfact)
    if (!is.na(argv$t2m.file)) {
      if ( argv$proj4t2m == "" & 
           argv$t2m.proj4_var == "" & 
           argv$t2m.proj4_att == "" ) {
        t2m.xy_as_vars<-T
        proj4t2m<-NULL
        proj4t2m_from_nc<-NULL
      } else {
        t2m.xy_as_vars<-F
        proj4t2m<-argv$proj4t2m
        proj4t2m_from_nc<-list(var=argv$t2m.proj4_var, att=argv$t2m.proj4_att)
      }
    }
    # read temperature from gridded field
    debug.file <- ifelse( argv$debug, file.path(argv$debug.dir,"input_data_rrwcor_t2m.RData"), NA)
    res<-get_data_from_ncfile(nc.file=argv$t2m.file,
                              nc.varname=argv$t2m.varname,
                              nc.t=argv$t2m.t,
                              nc.e=argv$t2m.e,
                              topdown=argv$t2m.topdown,
                              var.dim=list(ndim=argv$t2m.ndim,
                                           tpos=argv$t2m.tpos,
                                           epos=argv$t2m.epos,
                                           names=argv$t2m.dimnames),
                              var.de_acc=FALSE,
                              var.de_acc_by="",
                              proj4=proj4t2m,
                              proj4_from_nc=proj4t2m_from_nc,
                              xy_as_vars=t2m.xy_as_vars,
                              x_as_var.varname=argv$t2m.x_as_var.varname,
                              y_as_var.varname=argv$t2m.y_as_var.varname,
                              xy_as_var.dim=list(ndim=argv$t2m.xy_as_var.ndim,
                                                 tpos=argv$t2m.xy_as_var.tpos,
                                                 epos=NULL,
                                                 names=argv$t2m.xy_as_var.dimnames),
                              xy_as_var.dh_max=NA,
                              return_raster=F,
                              debug.file=debug.file,
                              nc.dqc_mode="none") 
    t2m<-t2m.offset+ res*t2m.cfact; rm(res)
    # read elevation from gridded field
    debug.file<-ifelse(argv$debug, file.path(argv$debug.dir,"input_data_rrwcor_t2mdem.RData"), NA)
    res<-get_data_from_ncfile(nc.file=argv$t2m.demfile,
                              nc.varname=argv$t2m.demvarname,
                              nc.t=argv$t2m.demt,
                              nc.e=argv$t2m.deme,
                              topdown=argv$t2m.demtopdown,
                              var.dim=list(ndim=argv$t2m.demndim,
                                           tpos=argv$t2m.demtpos,
                                           epos=argv$t2m.demepos,
                                           names=argv$t2m.demdimnames),
                              var.de_acc=FALSE,
                              var.de_acc_by="",
                              proj4=proj4t2m,
                              proj4_from_nc=proj4t2m_from_nc,
                              xy_as_vars=t2m.xy_as_vars,
                              x_as_var.varname=argv$t2m.x_as_var.varname,
                              y_as_var.varname=argv$t2m.y_as_var.varname,
                              xy_as_var.dim=list(ndim=argv$t2m.xy_as_var.ndim,
                                                 tpos=argv$t2m.xy_as_var.tpos,
                                                 epos=NULL,
                                                 names=argv$t2m.xy_as_var.dimnames),
                              xy_as_var.dh_max=NA,
                              return_raster=F,
                              debug.file=debug.file,
                              nc.dqc_mode="none") 
    zt2mdem<-t2m.demoffset+ res*t2m.demcfact; rm(res)
    t2m<-t2m+argv$gamma.standard*(z-zt2mdem)
  }
  # read windspeed from gridded field
  #  case of windspeed in the file
  if (argv$proj4wind=="" & argv$wind.proj4_var=="" & argv$wind.proj4_att=="" ) {
    wind.xy_as_vars<-T
    proj4wind<-NULL
    proj4wind_from_nc<-NULL
  } else {
    wind.xy_as_vars<-F
    proj4wind<-argv$proj4wind
    proj4wind_from_nc<-list(var=argv$wind.proj4_var, att=argv$wind.proj4_att)
  }
  if (!is.na(argv$windspeed.varname)) {
    debug.file<-ifelse(argv$debug, file.path(argv$debug.dir,"input_data_ws.RData"), NA)
    ws10m<-get_data_from_ncfile(nc.file=argv$wind.file,
                                nc.varname=argv$windspeed.varname,
                                nc.t=argv$wind.t,
                                nc.e=argv$wind.e,
                                topdown=argv$wind.topdown,
                                var.dim=list(ndim=argv$wind.ndim,
                                             tpos=argv$wind.tpos,
                                             epos=argv$wind.epos,
                                             names=argv$wind.dimnames),
                                var.de_acc=FALSE,
                                var.de_acc_by="",
                                proj4=proj4wind,
                                proj4_from_nc=proj4wind_from_nc,
                                xy_as_vars=wind.xy_as_vars,
                                x_as_var.varname=argv$wind.x_as_var.varname,
                                y_as_var.varname=argv$wind.y_as_var.varname,
                                xy_as_var.dim=list(ndim=argv$wind.xy_as_var.ndim,
                                                   tpos=argv$wind.xy_as_var.tpos,
                                                   epos=NULL,
                                                   names=argv$wind.xy_as_var.dimnames),
                                xy_as_var.dh_max=NA,
                                return_raster=F,
                                debug.file=debug.file,
                                nc.dqc_mode="none") 
  #  case of u,v in the file
  } else if (!is.na(argv$u.varname) & !is.na(argv$v.varname)) {
    debug.file<-ifelse(argv$debug, file.path(argv$debug.dir,"input_data_u10m.RData"), NA)
    u10m<-get_data_from_ncfile(nc.file=argv$wind.file,
                               nc.varname=argv$u.varname,
                               nc.t=argv$wind.t,
                               nc.e=argv$wind.e,
                               topdown=argv$wind.topdown,
                               var.dim=list(ndim=argv$wind.ndim,
                                            tpos=argv$wind.tpos,
                                            epos=argv$wind.epos,
                                            names=argv$wind.dimnames),
                               var.de_acc=FALSE,
                               var.de_acc_by="",
                               proj4=proj4wind,
                               proj4_from_nc=proj4wind_from_nc,
                               xy_as_vars=wind.xy_as_vars,
                               x_as_var.varname=argv$wind.x_as_var.varname,
                               y_as_var.varname=argv$wind.y_as_var.varname,
                               xy_as_var.dim=list(ndim=argv$wind.xy_as_var.ndim,
                                                  tpos=argv$wind.xy_as_var.tpos,
                                                  epos=NULL,
                                                  names=argv$wind.xy_as_var.dimnames),
                               xy_as_var.dh_max=NA,
                               return_raster=F,
                               debug.file=debug.file,
                               nc.dqc_mode="none") 
    debug.file<-ifelse(argv$debug, file.path(argv$debug.dir,"input_data_v10m.RData"), NA)
    v10m<-get_data_from_ncfile(nc.file=argv$wind.file,
                               nc.varname=argv$v.varname,
                               nc.t=argv$wind.t,
                               nc.e=argv$wind.e,
                               topdown=argv$wind.topdown,
                               var.dim=list(ndim=argv$wind.ndim,
                                            tpos=argv$wind.tpos,
                                            epos=argv$wind.epos,
                                            names=argv$wind.dimnames),
                               var.de_acc=FALSE,
                               var.de_acc_by="",
                               proj4=proj4wind,
                               proj4_from_nc=proj4wind_from_nc,
                               xy_as_vars=wind.xy_as_vars,
                               x_as_var.varname=argv$wind.x_as_var.varname,
                               y_as_var.varname=argv$wind.y_as_var.varname,
                               xy_as_var.dim=list(ndim=argv$wind.xy_as_var.ndim,
                                                  tpos=argv$wind.xy_as_var.tpos,
                                                  epos=NULL,
                                                  names=argv$wind.xy_as_var.dimnames),
                               xy_as_var.dh_max=NA,
                               return_raster=F,
                               debug.file=debug.file,
                               nc.dqc_mode="none") 

    ws10m<-sqrt(u10m*u10m+v10m*v10m)
    rm(u10m,v10m)
  #  case of wrong/no wind varname in the file
  } else {
    boom("ERROR precipitation correction for the wind undercatch: wind varname has not been specified")
  }
  if (argv$debug) save.image(file.path(argv$debug.dir,"rrcor_before.RData"))
  # precipitation data adjustment
  # 2023-11-06 CL begin
  if ( any( !is.na(argv$wind.prid))) {
    ix2corr <- which( data$prid %in% argv$wind.prid & !is.na(data$value) & is.na(dqcflag)) 
  } else {
    ix2corr <- which( !is.na(data$value) & is.na(dqcflag))
  }
  rawvalue <- data$value[ix2corr]
  res <- wolff_correction( par   = argv$rr.wcor.par,
                           t2m   = t2m[ix2corr],
                           ws10m = ws10m[ix2corr],
                           rr    = rawvalue)
  data$rawvalue <- data$value
  data$value[ix2corr]  <- res$rr.cor
  data$vsigma <- data$value; data$vsigma[] <- NA
  data$vsigma[ix2corr] <- res$sigma
  # 2023-11-06 CL end
  rm(res)
  if (argv$verbose) {
    print(paste0("# observations adjusted for wind-undercatch = ",
          length(ix2corr)))
    print(paste0("# not NAs-observations set to NAs after this correction = ",
          length(which( is.na(data$value[ix2corr]) & !is.na(data$rawvalue[ix2corr]) ))))
    print(paste0("# observations (ok-so-far) = ",
          length(which(!is.na(data$value) & is.na(dqcflag)))))
    print(paste0("# observations (>=0 & ok-so-far) = ",
          length(which(data$value>=0 & !is.na(data$value) & is.na(dqcflag)))))
    ix<-which(data$value[ix2corr]>=0 & !is.na(data$value[ix2corr]) & is.na(dqcflag[ix2corr]))
    if (length(ix)>0) {
      xx <- data$rawvalue[ix2corr][ix]; yy <- data$value[ix2corr][ix]
      lm<-lm(yy~xx+0)
      print(paste0("linear regression, obs_adjusted = ",
      round(as.numeric(lm$coefficients),3)," * obs_raw"))
    }
    print("+---------------------------------+")
  }
  if (argv$debug) save.image(file.path(argv$debug.dir,"rrcor.RData")) 
  #
  return( data)
}
