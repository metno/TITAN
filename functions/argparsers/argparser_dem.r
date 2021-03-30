  #.............................................................................. 
  # check elevation against dem
  p <- add_argument(p, "--dem",
       help="check elevation against digital elevation model (dem)",
                    flag=T,
                    short="-dm")
  p <- add_argument(p, "--dz.dem",
       help="maximum allowed deviation between observation and dem elevations [m]",
                    type="numeric",
                    default=500,
                    short="-zD")
  p <- add_argument(p, "--dem.fill",
                    help="fill missing elevation with data from dem",
                    flag=T,
                    short="-df")
  p <- add_argument(p, "--dem.file",
       help="land area fraction file (netCDF in kilometric coordinates)",
                    type="character",
                    default=NULL,
                    short="-dmf")
  p <- add_argument(p, "--proj4dem",
                    help="proj4 string for the dem",
                    type="character",
                    default=proj4_where_dqc_is_done_default,
                    short="-pd")
  p <- add_argument(p, "--dem.varname",
                    help="variable name in the netCDF file",
                    type="character",
                    default="altitude",
                    short="-dmv")
  p <- add_argument(p, "--dem.topdown",
                    help="logical, netCDF topdown parameter. If TRUE then turn the dem upside down",
                    flag=T,
                    short="-dmtd")
  p <- add_argument(p, "--dem.ndim",
                    help="number of dimensions in the netCDF file",
                    type="numeric",
                    default=3,
                    short="-dmnd")
  p <- add_argument(p, "--dem.tpos",
                    help="position of the dimension ''time'' in the netCDF file",
                    type="numeric",
                    default=3,
                    short="-dmti")
  p <- add_argument(p, "--dem.dimnames",
                    help="dimension names in the netCDF file",
                    type="character",
                    default=c("x","y","time"),
                    short="-dmna",
                    nargs=Inf)
  p <- add_argument(p, "--dem.proj4_var",
                    help="variable that include the specification of the proj4 string",
                    type="character",
                    default="projection_lambert",
                    short="-dmp4v")
  p <- add_argument(p, "--dem.proj4_att",
                    help="attribute with the specification of the proj4 string",
                    type="character",
                    default="proj4",
                    short="-dmp4a")
  p <- add_argument(p, "--dem.x_as_var.varname",
                    help="easting coordinate, variable name (used when proj4 is not specified)",
                    type="character",
                    default=NA)
  p <- add_argument(p, "--dem.y_as_var.varname",
                    help="northing coordinate, variable name (used when proj4 is not specified)",
                    type="character",
                    default=NA)
  p <- add_argument(p, "--dem.xy_as_var.ndim",
                    help="easting/northing coordinates, number of dimensions",
                    type="numeric",
                    default=NA)
  p <- add_argument(p, "--dem.xy_as_var.tpos",
                    help="easting/northing coordinates, position of time dimension",
                    type="numeric",
                    default=NA)
  p <- add_argument(p, "--dem.xy_as_var.dimnames",
                    help="easting/northing coordinates, dimension names in the netCDF file",
                    type="character",
                    default=NA,
                    nargs=Inf)
