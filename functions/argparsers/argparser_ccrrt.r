  #.............................................................................. 
  # precipitation and temperature cross-check
  p <- add_argument(p, "--ccrrt",
             help="precipitation (in-situ) and temperature (field) cross-check",
                    flag=T,
                    short="-ccrrtf")
  p <- add_argument(p, "--ccrrt.tmin",
                    help="temperature thresholds (vector, negative values start with \"_\" (e.g. _2=-2)",
                    type="character",
                    default=NA,
                    nargs=Inf,
                    short="-ccrrtt")
  p <- add_argument(p,"--ccrrt.filesetup",
                    help="predefined setup to read gridded fields: dem, t2m. available options: meps",
                    type="character",
                    default=NULL,
                    short="-ccrrtfs")
  p <- add_argument(p,"--ccrrt.filemeps",
                    help="meps netCDF file",
                    type="character",
                    default=NULL,
                    short="-ccrrtfm")
