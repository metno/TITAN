  # ADDITIONAL input files / providers

  p <- add_argument(p, "--input.files",
                    help="input files (provider 1 provider2 provider3 ...)",
                    type="character",
                    default=NULL,
                    nargs=Inf)

  p <- add_argument(p, "--output.file",
                    help="output file",
                    type="character",
                    default="output.txt")

  # configuration file
  p <- add_argument(p, "--config.files",
                    help="configuration file(s) (use conf <- list( var1=..., var2=..., ... )",
                    type="character",
                    default=NULL,
                    nargs=Inf)

  #.............................................................................. 
  #
  # first-guess configuration files
  p <- add_argument(p, "--fg.files",
                    help="information used to read first-guess nc-file(s) (use conf <- list( var1=..., var2=..., ... )",
                    type="character",
                    default=NULL,
                    nargs=Inf)

  p <- add_argument(p, "--fg.filenames",
                    help="file names of the first-guess nc-files. It is an optional argument, if specified: i) must have the same length of --fg.files ii) override the main.file arguments in the fg.files",
                    type="character",
                    default=NULL,
                    nargs=Inf)


  #.............................................................................. 
  neg.str<-"Negative values can be specified in either one of these two ways: (1) by using the corresponding \"...neg...\" command line argument; (2) negative values start with \"_\" (e.g. _2=-2)"
  # ADDITIONAL input files / providers
  p <- add_argument(p, "--prid",
                    help="provider identifiers (provider1 provider2 provider3 ...)",
                    type="character",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--input.offset",
                    help=paste("offset applied to the input files (one for each provider, default=0).",neg.str),
                    type="character",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--input.negoffset",
                    help="sign for the offsets (1=negative; 0=positive, def=0)",
                    type="numeric",
                    default=NA,
                    nargs=Inf)
  #
  p <- add_argument(p, "--input.cfact",
                    help=paste("correction factor applied to the input files (one for each provider, default=1)",neg.str),
                    type="character",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--input.negcfact",
                    help="sign for the correction factors (1=negative; 0=positive, def=0)",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--separator",
                    help="character vector, input file(s) separator character(s) (default '';'')",
                    type="character",
                    nargs=Inf,
                    default=NA)

  p <- add_argument(p, "--separator.out",
                    help="separator character in the output file",
                    type="character",
                    default=";")

  p<- add_argument(p, "--latlon.dig.out",
                   help="number of decimal digits for latitude and longitude in the output file  (obsolete, use xy.dig.out)",
                   type="numeric",
                   default=5)

  p<- add_argument(p, "--xy.dig.out",
                   help="number of decimal digits for northing and easting coordinates in the output file",
                   type="numeric",
                   default=xy.dig.out_default)

  p<- add_argument(p, "--elev.dig.out",
                   help="number of decimal digits for elevation in the output file",
                   type="numeric",
                   default=1)

  p<- add_argument(p, "--value.dig.out",
                   help="number of decimal digits for the returned value in the output file",
                   type="numeric",
                   default=1)

  p <- add_argument(p, "--varname.lat",
                    help="character vector, latitude variable name(s) in the input file (default ''lat'')",
                    type="character",
                    nargs=Inf)

  p <- add_argument(p, "--varname.lat.out",
                    help="latitude variable name in the output file (obsolete, use varname.y.out)",
                    type="character",
                    default="lat")

  p <- add_argument(p, "--varname.y.out",
                    help="northing coordinate name in the output file",
                    type="character",
                    default=varname.y.out_default)

  p <- add_argument(p, "--varname.lon",
                    help="character vector, longitude variable name(s) in the input file (default ''lon'')",
                    type="character",
                    nargs=Inf)

  p <- add_argument(p, "--varname.lon.out",
                    help="longitude variable name in the output file (obsolete, use varname.x.out)",
                    type="character",
                    default="lon")

  p <- add_argument(p, "--varname.x.out",
                    help="easting coordinate name in the output file",
                    type="character",
                    default="lon")

  p <- add_argument(p, "--varname.elev",
                    help="character vector, elevation variable names(s) in the input file (default ''elev'')",
                    type="character",
                    nargs=Inf)

  p <- add_argument(p, "--varname.elev.out",
                    help="elevation variable name in the output file",
                    type="character",
                    default="elev")

  p <- add_argument(p, "--elev_not_used",
                    help="elevation is not used (will be set to zero)",
                    flag=T)

  p <- add_argument(p, "--varname.value",
                    help="character vector, variable name(s) in the input file (default ''value'')",
                    type="character",
                    nargs=Inf)

  p <- add_argument(p, "--varname.value.out",
                    help="name for the variable values (out)",
                    type="character",
                    default="value")

  # output file

  p <- add_argument(p, "--varname.opt",
       help="additional optional variables to be written on the output (out)",
                    type="character",
                    default=NA,
                    nargs=Inf)

  p<- add_argument(p, "--varname.prid",
                   help="name for the provider identifier (out)",
                   type="character",
                   default="prid")

  p<- add_argument(p, "--varname.dqc",
                   help="name for the data quality control flag (out)",
                   type="character",
                   default="dqc")

  p<- add_argument(p, "--varname.sct",
              help="name for the spatial consistency test returned value (out)",
                   type="character",
                   default="sct")

  p<- add_argument(p, "--varname.rep",
              help="name for the coefficient of representativeness (out)",
                   type="character",
                   default="rep")
#
  #.............................................................................. 
  # Timestamp valid for all the netcdf files
  p <- add_argument(p, "--timestamp",
                    help="timestamp, valid for all the netCDF file (YYYYMMDDHH00)",
                    type="character",
                    default=NA)

