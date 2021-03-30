#+ parse the command line arguments
argparser <- function() {
# description:
# output:
#------------------------------------------------------------------------------
  # create parser object
  p <- arg_parser("titan")
  # specify our desired options
  #.............................................................................. 
  # INPUT/OUTPUT miscellaneous
  source( file.path( titan_fun_path, "argparsers/argparser_IO.r"), local=T)
  #.............................................................................. 
  # miscellaneous (debug, verbose, cores, paths)
  source( file.path( titan_fun_path, "argparsers/argparser_misc.r"), local=T)
  #.............................................................................. 
  # Quality control codes
  source( file.path( titan_fun_path, "argparsers/argparser_qccodes.r"), local=T)
  #.............................................................................. 
  # variable, data transformation and lapse rate  
  source( file.path( titan_fun_path, "argparsers/argparser_data_misc.r"), local=T)
  #.............................................................................. 
  # GEOGRAPHICAL domain definition  
  source( file.path( titan_fun_path, "argparsers/argparser_geopar.r"), local=T)
  #.............................................................................. 
  # duplicates
  source( file.path( titan_fun_path, "argparsers/argparser_duplicates.r"), local=T)
  #.............................................................................. 
  # metadata check
  source( file.path( titan_fun_path, "argparsers/argparser_metadata_check.r"), local=T)
  #.............................................................................. 
  # Plausibility check
  source( file.path( titan_fun_path, "argparsers/argparser_plausibility.r"), local=T)
  #.............................................................................. 
  # Climatological check
  source( file.path( titan_fun_path, "argparsers/argparser_climatcheck.r"), local=T)
  #.............................................................................. 
  # Buddy-check
  source( file.path( titan_fun_path, "argparsers/argparser_buddy.r"), local=T)
  #.............................................................................. 
  # isolated stations
  source( file.path( titan_fun_path, "argparsers/argparser_isolation.r"), local=T)
  #.............................................................................. 
  # first guess test
  source( file.path( titan_fun_path, "argparsers/argparser_fgt.r"), local=T)
  #.............................................................................. 
  # spatial consistency test
  source( file.path( titan_fun_path, "argparsers/argparser_sct.r"), local=T)
  #.............................................................................. 
  # spatial consistency test wit first guess
  source( file.path( titan_fun_path, "argparsers/argparser_sct_fg.r"), local=T)
  #.............................................................................. 
  # spatial consistency test dual
  source( file.path( titan_fun_path, "argparsers/argparser_sct_dual.r"), local=T)
  #.............................................................................. 
  # spatial consistency test dual wit first guess
  source( file.path( titan_fun_path, "argparsers/argparser_sct_fg_dual.r"), local=T)
  #.............................................................................. 
  # read dem file
  source( file.path( titan_fun_path, "argparsers/argparser_dem.r"), local=T)
  #.............................................................................. 
  # precipitation and temperature cross-check
  source( file.path( titan_fun_path, "argparsers/argparser_ccrrt.r"), local=T)
  #.............................................................................. 
  # blacklist
  source( file.path( titan_fun_path, "argparsers/argparser_blacklist.r"), local=T)
  #.............................................................................. 
  # keeplist
  source( file.path( titan_fun_path, "argparsers/argparser_keep.r"), local=T)
  #.............................................................................. 
  # doit flags
  source( file.path( titan_fun_path, "argparsers/argparser_doit.r"), local=T)
  #.............................................................................. 
  # precipitation correction for the wind-induced undercatch
  source( file.path( titan_fun_path, "argparsers/argparser_rr_wind_corr.r"), local=T)
  #.............................................................................. 

  #.............................................................................. 
  # PARSE arguments
  argv <- parse_args(p)

  #
  #-----------------------------------------------------------------------------
  # read configuration files

  if ( any( !is.na( argv$config.files))) {
    for (f in 1:length(argv$config.files)) {
      if ( file.exists( argv$config.files[f])) {
        source( argv$config.files[f], local=T)
        argv_tmp       <- append( argv, conf)
        names_argv_tmp <- names( argv_tmp)
        argv_def       <- list()
        names_argv_def <- integer(0)
        k <- 0
        for (i in 1:length(argv_tmp)) {
          if (names_argv_tmp[i] %in% names_argv_def) next
          k <- k + 1
          j <- which( names_argv_tmp == names_argv_tmp[i])
          argv_def[[k]]  <- argv_tmp[[j[length(j)]]]
          names_argv_def <- c( names_argv_def, names_argv_tmp[i])
        }
        names( argv_def) <- names_argv_def
        rm( argv, argv_tmp, names_argv_tmp, names_argv_def)
        argv <- argv_def
        rm( argv_def)
      } else {
        print( "WARNING: config file not found")
        print( argv$config.files[f])
      }
    }
  }

  #
  #-----------------------------------------------------------------------------
  # read fg config files
  if ( any( !is.na( argv$fg.files))) {
    for (f in 1:length(argv$fg.files)) {
      if ( file.exists( argv$fg.files[f])) {
        source( argv$fg.files[f], local=T)
        fg_env$fg[[f]] <- conf
        rm( conf)
      } else {
        print( "WARNING: config file not found")
        print( argv$fg.files[f])
      }
    }
  }

  # one can specify file names outside the fg config files
  if ( any( !is.na( argv$fg.filenames))) {
    if ( length( argv$fg.files) != length( argv$fg.filenames))
      boom( "since --fg.files and --fg.filenames are both defined, they must have the same lenght")
    for (f in 1:length(argv$fg.files)) {
      if ( file.exists( argv$fg.filenames[f])) {
        fg_env$fg[[f]]$main.file <- argv$fg.filenames[f]
      } else {
        boom( paste( "ERROR: file not foud", argv$fg.filenames[f]))
      }
    }
  }

  #
  #-----------------------------------------------------------------------------
  # CHECKS on input arguments
  # more than one input file
  if ( any( !is.na( argv$input.files))) {
    for (j in 1:length(argv$input.files)) {
      if ( !file.exists( argv$input.files[j])) 
        boom( paste( "ERROR: input file not found", argv$input.files[j]))
    }
  } else {
    boom( paste( "ERROR: specify at least one input file"))
  }
  nfin <- length(argv$input.files)

  # check consistency between number of files and provider ids
  if ( any( is.na( argv$prid))) {
    argv$prid <- 1:nfin
  } else {
    if ( length(argv$prid) != nfin) 
      boom( "ERROR: number of provider identifier is different from the number of input files")
  }

  #................................................................................
  # set input offsets and correction factors
  argv$input.offset <- strings_to_numbers( strings     = argv$input.offset, 
                                           default     = 0,
                                           strings_dim = nfin,
                                           neg         = argv$input.negoffset)
  argv$input.cfact <- strings_to_numbers( strings     = argv$input.cfact,
                                          default     = 1,
                                          strings_dim = nfin,
                                          neg         = argv$input.negcfact)
  #................................................................................
  # check variable
  if (!(argv$variable %in% c("T","RR"))) 
    boom("variable must be one of T, RR")

  # set proj4 variables (proj4from and proj4to are obsolete)
  if (argv$proj4from!=argv$proj4_input_obsfiles) {
    if (argv$proj4_input_obsfiles==proj4_input_obsfiles_default & 
        argv$proj4from!=proj4_input_obsfiles_default)
      argv$proj4_input_obsfiles<-argv$proj4from
  }
  if (argv$proj4to!=argv$proj4_where_dqc_is_done) {
    if (argv$proj4_where_dqc_is_done==proj4_where_dqc_is_done_default & 
        argv$proj4to!=proj4_where_dqc_is_done_default)
      argv$proj4_where_dqc_is_done<-argv$proj4to
  }
  #................................................................................
  # set variables to customize output
  if (argv$varname.lat.out!=argv$varname.y.out) {
    if (argv$varname.y.out==varname.y.out_default & 
        argv$varname.lat.out!=varname.y.out_default)
      argv$varname.y.out<-argv$varname.lat.out
  }
  if (argv$varname.lon.out!=argv$varname.x.out) {
    if (argv$varname.x.out==varname.x.out_default & 
        argv$varname.lon.out!=varname.x.out_default)
      argv$varname.x.out<-argv$varname.lon.out
  }
  if (argv$latlon.dig.out!=argv$xy.dig.out) {
    if (argv$xy.dig.out==xy.dig.out_default & 
        argv$latlon.dig.out!=xy.dig.out_default)
      argv$xy.dig.out<-argv$latlon.dig.out
  }
  #................................................................................
  # set the input arguments according to user specification

  # check external files
  if (argv$dem | argv$dem.fill) {
    if (!file.exists(argv$dem.file)) 
      boom(paste("ERROR: dem file not found",argv$dem.file))
  }

  # input column names
  if (any(is.na(argv$separator))) 
    argv$separator<-rep(";",nfin)
  if (length(argv$separator)==0) 
    argv$separator<-rep(";",nfin)
  if (length(argv$separator)!=nfin) 
    argv$separator<-rep(argv$separator[1],nfin)
  for (i in 1:length(argv$separator)) if (argv$separator[i]=="comma")  argv$separator[i]<-","

  if (any(is.na(argv$varname.lat))) 
    argv$varname.lat<-rep("lat",nfin)
  if (length(argv$varname.lat)==0) 
    argv$varname.lat<-rep("lat",nfin)
  if (length(argv$varname.lat)!=nfin) 
    argv$varname.lat<-rep(argv$varname.lat[1],nfin)

  if (any(is.na(argv$varname.lon))) 
    argv$varname.lon<-rep("lon",nfin)
  if (length(argv$varname.lon)==0) 
    argv$varname.lon<-rep("lon",nfin)
  if (length(argv$varname.lon)!=nfin) 
    argv$varname.lon<-rep(argv$varname.lon[1],nfin)

  if (any(is.na(argv$varname.elev))) 
    argv$varname.elev<-rep("elev",nfin)
  if (length(argv$varname.elev)==0) 
    argv$varname.elev<-rep("elev",nfin)
  if (length(argv$varname.elev)!=nfin) 
    argv$varname.elev<-rep(argv$varname.elev[1],nfin)

  if (any(is.na(argv$varname.value))) 
    argv$varname.value<-rep("value",nfin)
  if (length(argv$varname.value)==0) 
    argv$varname.value<-rep("value",nfin)
  if (length(argv$varname.value)!=nfin) 
    argv$varname.value<-rep(argv$varname.value[1],nfin)

  # check for spatial conversion flag
  if (!argv$spatconv) 
    boom(paste("ERROR: \"--spatconv\" (-c) option must be used on the command line. Default is: input is in lat-lon coordinates; some DQC tests takes place in kilometric coordinates specified by the user; output is in lat-lon coordinates"))

  # load netcdf library, if needed
  if ( argv$dem | argv$dem.fill |
      !is.na(argv$t2m.file))
    suppressPackageStartupMessages(library("ncdf4")) 
  #
  if (!is.na(argv$month.clim) & (argv$month.clim<1 | argv$month.clim>12)) {
    boom(paste("ERROR: month number is wrong. month number=",argv$month.clim))
  } else if (!is.na(argv$month.clim) & 
             (length(which(!is.na(argv$vmin.clim)))!=12 | 
              length(which(!is.na(argv$vmax.clim)))!=12) ) {
    boom("ERROR: climatological check, vmin.clim and/or vmax.clim vectors must have 12 arguments")
  }

  # blacklist
  if (any(!is.na(argv$blacklist.lat)) | 
      any(!is.na(argv$blacklist.lon)) |
      any(!is.na(argv$blacklist.fll)) ) {
    if ( (length(argv$blacklist.lat)!=length(argv$blacklist.lon))  |
         (length(argv$blacklist.lat)!=length(argv$blacklist.fll))  |
         (any(is.na(argv$blacklist.fll))) | 
         (any(is.na(argv$blacklist.lat))) | 
         (any(is.na(argv$blacklist.lon))) ) {
      boom(paste("ERROR in the blacklist definition, must have same number of lat,lon,IDprovider points. lat number=",argv$blacklist.lat,". lon number=",argv$blacklist.lon,". ID provider number=",argv$blacklist.fll))
    }
  }
  if (any(!is.na(argv$blacklist.idx)) | 
      any(!is.na(argv$blacklist.fidx)) ) {
    if ( (length(argv$blacklist.idx)!=length(argv$blacklist.fidx))  |
         (any(is.na(argv$blacklist.idx))) | 
         (any(is.na(argv$blacklist.fidx))) ) {
      print("ERROR in the blacklist definition, must have same number of index and IDprovider points")
      print(paste("index number=",argv$blacklist.idx))
      print(paste("ID provider number=",argv$blacklist.fidx))
      boom()
    }
  }
  # keeplist
  if (any(!is.na(argv$keeplist.lat)) | 
      any(!is.na(argv$keeplist.lon)) |
      any(!is.na(argv$keeplist.fll)) ) {
    if ( (length(argv$keeplist.lat)!=length(argv$keeplist.lon))  |
         (length(argv$keeplist.lat)!=length(argv$keeplist.fll))  |
         (any(is.na(argv$keeplist.fll))) | 
         (any(is.na(argv$keeplist.lat))) | 
         (any(is.na(argv$keeplist.lon))) ) {
      print("ERROR in the keeplist definition, must have same number of lat,lon,IDprovider points")
      print(paste("lat number=",argv$keeplist.lat))
      print(paste("lon number=",argv$keeplist.lon))
      print(paste("ID provider number=",argv$keeplist.fll))
      boom()
    }
  }
  if (any(!is.na(argv$keeplist.idx)) | 
      any(!is.na(argv$keeplist.fidx)) ) {
    if ( (length(argv$keeplist.idx)!=length(argv$keeplist.fidx))  |
         (any(is.na(argv$keeplist.idx))) | 
         (any(is.na(argv$keeplist.fidx))) ) {
      print("ERROR in the keeplist definition, must have same number of index and IDprovider points")
      print(paste("index number=",argv$keeplist.idx))
      print(paste("ID provider number=",argv$keeplist.fidx))
      boom()
    }
  }

  #
  # precip and temperature crosscheck
  argv$ccrrt.tmin <- strings_to_numbers( strings = argv$ccrrt.tmin,
                                         strings_dim = nfin)

  #
  # doit flags
  if ( any( is.na( argv$doit.buddy))) argv$doit.buddy <- rep( 1, length=nfin)
  if ( any( is.na( argv$doit.fgt)))   argv$doit.fgt   <- rep( 1, length=nfin)
  if ( any( is.na( argv$doit.sct)))   argv$doit.sct   <- rep( 1, length=nfin)
  if ( any( is.na( argv$doit.sct_fg)))   argv$doit.sct_fg   <- rep( 1, length=nfin)
  if ( any( is.na( argv$doit.sct_dual)))    argv$doit.sct_dual <- rep(1,length=nfin)
  if ( any( is.na( argv$doit.sct_fg_dual))) argv$doit.sct_fg_dual<-rep(1,length=nfin)
  if ( any( is.na( argv$doit.clim))) argv$doit.clim<-rep(1,length=nfin)
  if ( any( is.na( argv$doit.dem)))  argv$doit.dem<-rep(1,length=nfin)
  if ( any( is.na( argv$doit.iso))) argv$doit.iso<-rep(1,length=nfin)
  if ( any( !( argv$doit.buddy %in% c(0,1,2)))) boom("doit.buddy must contain only 0,1,2")
  if ( any( !( argv$doit.fgt %in% c(0,1,2)))) boom("doit.fgt must contain only 0,1,2")
  if ( any( !( argv$doit.sct %in% c(0,1,2)))) boom("doit.sct must contain only 0,1,2")
  if ( any( !( argv$doit.sct_fg %in% c(0,1,2)))) boom("doit.sct_fg must contain only 0,1,2")
  if ( any( !( argv$doit.sct_dual %in% c(0,1,2)))) boom("doit.sct_dual must contain only 0,1,2")
  if ( any( !( argv$doit.sct_fg_dual %in% c(0,1,2)))) boom("doit.sct_fg_dual must contain only 0,1,2")
  if ( any( !( argv$doit.clim %in% c(0,1,2)))) boom("doit.clim must contain only 0,1,2")
  if ( any( !( argv$doit.dem %in% c(0,1,2)))) boom("doit.dem must contain only 0,1,2")
  if ( any( !( argv$doit.iso %in% c(0,1,2)))) boom("doit.iso must contain only 0,1,2")

  #
  # set the thresholds for the plausibility check
  argv$vmin <- strings_to_numbers( strings=argv$vmin, neg=argv$vminsign)
  argv$vmax <- strings_to_numbers( strings=argv$vmax, neg=argv$vmaxsign)
  argv$vmin.clim <- strings_to_numbers( strings=argv$vmin.clim,
                                        strings_dim=12, neg=argv$vminsign.clim)
  argv$vmax.clim <- strings_to_numbers( strings=argv$vmax.clim,
                                        strings_dim=12, neg=argv$vmaxsign.clim)
  # priorities
  if ( is.null( argv$prio.buddy))        argv$prio.buddy <- rep( -1, length=nfin)
  if ( any( is.na( argv$prio.buddy)))    argv$prio.buddy <- rep( -1, length=nfin)
  if ( length( argv$prio.buddy) != nfin) argv$prio.buddy <- rep( -1, length=nfin)

  if ( is.null( argv$prio.sct_dual))        argv$prio.sct_dual <- rep( -1, length=nfin)
  if ( any( is.na( argv$prio.sct_dual)))    argv$prio.sct_dual <- rep( -1, length=nfin)
  if ( length( argv$prio.sct_dual) != nfin) argv$prio.sct_dual <- rep( -1, length=nfin)

  if ( is.null( argv$prio.sct_fg_dual))        argv$prio.sct_fg_dual <- rep( -1, length=nfin)
  if ( any( is.na( argv$prio.sct_fg_dual)))    argv$prio.sct_fg_dual <- rep( -1, length=nfin)
  if ( length( argv$prio.sct_fg_dual) != nfin) argv$prio.sct_fg_dual <- rep( -1, length=nfin)

  if ( is.null( argv$prio.sct))        argv$prio.sct <- rep( -1, length=nfin)
  if ( any( is.na( argv$prio.sct)))    argv$prio.sct <- rep( -1, length=nfin)
  if ( length( argv$prio.sct) != nfin) argv$prio.sct <- rep( -1, length=nfin)

  if ( is.null( argv$prio.sct_fg))        argv$prio.sct_fg <- rep( -1, length=nfin)
  if ( any( is.na( argv$prio.sct_fg)))    argv$prio.sct_fg <- rep( -1, length=nfin)
  if ( length( argv$prio.sct_fg) != nfin) argv$prio.sct_fg <- rep( -1, length=nfin)

  if ( is.null( argv$prio.fgt))        argv$prio.fgt <- rep( -1, length=nfin)
  if ( any( is.na( argv$prio.fgt)))    argv$prio.fgt <- rep( -1, length=nfin)
  if ( length( argv$prio.fgt) != nfin) argv$prio.fgt <- rep( -1, length=nfin)

  # wind-induced undercatch of precipitation, check consistency of inputs
  if (argv$rr.wcor & argv$variable!="RR") 
    boom("ERROR: wind-induced correction for undercatch is implemented for precipitation only")

  # proj4s
  if (!is.na(argv$wind.file)) {
    if (argv$proj4wind=="" & argv$wind.proj4_var=="" & argv$wind.proj4_att=="" ) {
      wind.xy_as_vars<-T
      proj4wind<-NULL
      proj4wind_from_nc<-NULL
    } else {
      wind.xy_as_vars<-F
      proj4wind<-argv$proj4wind
      proj4wind_from_nc<-list(var=argv$wind.proj4_var, att=argv$wind.proj4_att)
    }
  }

  # set the timestamp
  if (!is.na(argv$timestamp)) {
    if (is.na(argv$fg.t)) argv$fg.t<-argv$timestamp
    if (is.na(argv$fge.t)) argv$fge.t<-argv$timestamp
    if (is.na(argv$wind.t)) argv$wind.t<-argv$timestamp
    if (is.na(argv$t2m.t)) argv$t2m.t<-argv$timestamp
  }
  #----------------------------------------------------------------------------

  if ( argv$verbose) cat( ">> TITAN <<\n")

  #----------------------------------------------------------------------------

  return(argv)

}
