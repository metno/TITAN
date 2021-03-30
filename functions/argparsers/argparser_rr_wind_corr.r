  #.............................................................................. 
  # precipitation correction for the wind-induced undercatch
  p <- add_argument(p, "--rr.wcor",
                help=" precipitation correction for the wind-induced undercatch",
                    flag=T)

  p <- add_argument(p,"--rr.wcor.filesetup",
                    help="predefined setup to read gridded fields: dem, t2m, wspeed. available options: meps",
                    type="character",
                    default=NULL)

  p <- add_argument(p,"--rr.wcor.filemeps",
                    help="meps netCDF file",
                    type="character",
                    default=NULL)

  # parameter 

  p <- add_argument(p, "--rr.wcor.par",
                    help=paste("Parameter used for correcting wind-induced loss",
                                " of solid precipitation (Wolff et al., 2015)"),
                    type="numeric",
                    default=c(4.24,1.81,0.18,0.99,0.66,1.07,0.18,0.11,2.35,0.12),
                    nargs=10)

  # temp 

  p <- add_argument(p,"--t2m.file",
                    help="air temperature netCDF file",
                    type="character",
                    default=NULL)

  p <- add_argument(p, "--t2m.offset",
                    help="air temperature offset",
                    type="character",
                    default="0")

  p <- add_argument(p, "--t2m.cfact",
                    help="air temperature correction factor",
                    type="character",
                    default="1")

  p <- add_argument(p, "--t2m.negoffset",
                    help="offset sign (1=neg, 0=pos)",
                    type="character",
                    default="0")

  p <- add_argument(p, "--t2m.negcfact",
                    help="correction factor sign (1=neg, 0=pos)",
                    type="character",
                    default="0")

  p <- add_argument(p, "--proj4t2m",
                    help="proj4 string for the air temperature file",
                    type="character",
                    default=proj4_where_dqc_is_done_default)

  p <- add_argument(p, "--t2m.varname",
                    help="air temperature variable name in the netCDF file",
                    type="character",
                    default=NULL)

  p <- add_argument(p, "--t2m.topdown",
                    help="logical, netCDF topdown parameter. If TRUE then turn the file upside down",
                    flag=T)

  p <- add_argument(p, "--t2m.ndim",
                    help="number of dimensions in the netCDF file",
                    type="numeric",
                    default=3)

  p <- add_argument(p, "--t2m.tpos",
                    help="position of the dimension ''time'' in the netCDF file",
                    type="numeric",
                    default=3)

  p <- add_argument(p, "--t2m.epos",
                    help="position of the dimension ''ensemble'' in the netCDF file",
                    type="numeric",
                    default=3)

  p <- add_argument(p, "--t2m.dimnames",
                    help="dimension names in the netCDF file",
                    type="character",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--t2m.proj4_var",
                    help="variable that include the specification of the proj4 string",
                    type="character",
                    default="projection_lambert")

  p <- add_argument(p, "--t2m.proj4_att",
                    help="attribute with the specification of the proj4 string",
                    type="character",
                    default="proj4")

  p <- add_argument(p, "--t2m.t",
                    help="timestamp to read in the air temperature netCDF file (YYYYMMDDHH00)",
                    type="character",
                    default=NA)

  p <- add_argument(p, "--t2m.e",
                    help="ensemble member to read in the air temperature netCDF file",
                    type="numeric",
                    default=NA)

  p <- add_argument(p, "--t2m.x_as_var.varname",
                    help="easting coordinate, variable name (used when proj4 is not specified)",
                    type="character",
                    default=NA)

  p <- add_argument(p, "--t2m.y_as_var.varname",
                    help="northing coordinate, variable name (used when proj4 is not specified)",
                    type="character",
                    default=NA)

  p <- add_argument(p, "--t2m.xy_as_var.ndim",
                    help="easting/northing coordinates, number of dimensions",
                    type="numeric",
                    default=NA)

  p <- add_argument(p, "--t2m.xy_as_var.tpos",
                    help="easting/northing coordinates, position of time dimension",
                    type="numeric",
                    default=NA)

  p <- add_argument(p, "--t2m.xy_as_var.dimnames",
                    help="easting/northing coordinates, dimension names in the netCDF file",
                    type="character",
                    default=NA,
                    nargs=Inf)

  # dem for temperature adjustments

  p <- add_argument(p,"--t2m.demfile",
                    help="dem file associated to the first-guess or background file",
                    type="character",
                    default=NULL)

  p <- add_argument(p, "--t2m.demoffset",
                    help="offset",
                    type="character",
                    default="0")

  p <- add_argument(p, "--t2m.demcfact",
                    help="correction factor",
                    type="character",
                    default="1")

  p <- add_argument(p, "--t2m.demnegoffset",
                    help="offset sign (1=neg, 0=pos)",
                    type="numeric",
                    default=0)

  p <- add_argument(p, "--t2m.demnegcfact",
                    help="correction factor sign (1=neg, 0=pos)",
                    type="numeric",
                    default=0)

  p <- add_argument(p, "--t2m.demt",
                    help="timestamp to read in the netCDF file (YYYYMMDDHH00)",
                    type="character",
                    default=NA)

  p <- add_argument(p, "--t2m.demvarname",
                    help="variable name in the netCDF file (dem associated to the first-guess)",
                    type="character",
                    default="none")

  p <- add_argument(p, "--t2m.demtopdown",
                    help="logical, netCDF topdown parameter. If TRUE then turn the field upside down",
                    flag=T)

  p <- add_argument(p, "--t2m.demndim",
                    help="number of dimensions in the netCDF file",
                    type="numeric",
                    default=3)

  p <- add_argument(p, "--t2m.demepos",
                    help="position of the dimension ''ensemble'' in the netCDF file",
                    type="numeric",
                    default=3)

  p <- add_argument(p, "--t2m.demtpos",
                    help="position of the dimension ''time'' in the netCDF file",
                    type="numeric",
                    default=3)

  p <- add_argument(p, "--t2m.deme",
                    help="ensemble member to read in the netCDF file",
                    type="numeric",
                    default=NA)

  p <- add_argument(p, "--t2m.demdimnames",
                    help="dimension names in the netCDF file",
                    type="character",
                    default=NA,
                    nargs=Inf)

  # wind

  p <- add_argument(p,"--wind.file",
                    help="air temperature netCDF file",
                    type="character",
                    default=NULL)

  p <- add_argument(p, "--proj4wind",
                    help="proj4 string for the air temperature file",
                    type="character",
                    default=proj4_where_dqc_is_done_default)

  p <- add_argument(p, "--windspeed.varname",
                    help="air temperature variable name in the netCDF file",
                    type="character",
                    default=NULL)

  p <- add_argument(p, "--u.varname",
                    help="air temperature variable name in the netCDF file",
                    type="character",
                    default=NULL)

  p <- add_argument(p, "--v.varname",
                    help="air temperature variable name in the netCDF file",
                    type="character",
                    default=NULL)

  p <- add_argument(p, "--wind.topdown",
                    help="logical, netCDF topdown parameter. If TRUE then turn the file upside down",
                    flag=T)

  p <- add_argument(p, "--wind.ndim",
                    help="number of dimensions in the netCDF file",
                    type="numeric",
                    default=3)

  p <- add_argument(p, "--wind.tpos",
                    help="position of the dimension ''time'' in the netCDF file",
                    type="numeric",
                    default=3)

  p <- add_argument(p, "--wind.epos",
                    help="position of the dimension ''ensemble'' in the netCDF file",
                    type="numeric",
                    default=3)

  p <- add_argument(p, "--wind.dimnames",
                    help="dimension names in the netCDF file",
                    type="character",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--wind.proj4_var",
                    help="variable that include the specification of the proj4 string",
                    type="character",
                    default="projection_lambert")

  p <- add_argument(p, "--wind.proj4_att",
                    help="attribute with the specification of the proj4 string",
                    type="character",
                    default="proj4")

  p <- add_argument(p, "--wind.t",
                    help="timestamp to read in the air temperature netCDF file (YYYYMMDDHH00)",
                    type="character",
                    default=NA)

  p <- add_argument(p, "--wind.e",
                    help="ensemble member to read in the air temperature netCDF file",
                    type="numeric",
                    default=NA)

  p <- add_argument(p, "--wind.x_as_var.varname",
                    help="easting coordinate, variable name (used when proj4 is not specified)",
                    type="character",
                    default=NA)

  p <- add_argument(p, "--wind.y_as_var.varname",
                    help="northing coordinate, variable name (used when proj4 is not specified)",
                    type="character",
                    default=NA)

  p <- add_argument(p, "--wind.xy_as_var.ndim",
                    help="easting/northing coordinates, number of dimensions",
                    type="numeric",
                    default=NA)

  p <- add_argument(p, "--wind.xy_as_var.tpos",
                    help="easting/northing coordinates, position of time dimension",
                    type="numeric",
                    default=NA)

  p <- add_argument(p, "--wind.xy_as_var.dimnames",
                    help="easting/northing coordinates, dimension names in the netCDF file",
                    type="character",
                    default=NA,
                    nargs=Inf)
