  #.............................................................................. 
  # GEOGRAPHICAL domain definition  
  # NOTE: lat-lon setup to have Oslo in a single box
  p <- add_argument(p, "--lonmin",
                    help="longitude of south-eastern domain corner",
                    type="character",
                    default="3",
                    short="-lon")
  p <- add_argument(p, "--lonmax",
                    help="longitude of south-western domain corner",
                    type="character",
                    default="33",
                    short="-lox")
  p <- add_argument(p, "--latmin",
                    help="latitude of south-eastern domain corner",
                    type="character",
                    default="53.25",
                    short="-lan")
  p <- add_argument(p, "--latmax",
                    help="latitude of north-western domain corner",
                    type="character",
                    default="71.8",
                    short="-lax")
  p <- add_argument(p, "--dqc_inbox_only",
                    help="perform dqc only in the defined box (lonmin,lonmax,latmin,latmax)",
                    flag=T)
  
  # transformation between coordinate reference systems
  p <- add_argument(p, "--spatconv",
                    help="flag for conversion of spatial coordinates before running the data quality checks",
                    flag=T,
                    short="-c")
  p <- add_argument(p, "--proj4from",
                    help="proj4 string for the original coordinate reference system (obsolete, use \"proj4_input_obsfiles\")",
                    type="character",
                    default=proj4_input_obsfiles_default,
                    short="-pf")
  p <- add_argument(p, "--proj4_input_obsfiles",
                    help="proj4 string for the original coordinate reference system",
                    type="character",
                    default=proj4_input_obsfiles_default)
  p <- add_argument(p, "--proj4to",
                    help="proj4 string for the coordinate reference system where the DQC is performed (obsolete, use \"proj4_where_dqc_is_done\"",
                    type="character",
                    default=proj4_where_dqc_is_done_default,
                    short="-pt")
  p <- add_argument(p, "--proj4_where_dqc_is_done",
                    help="proj4 string for the coordinate reference system where the DQC is performed",
                    type="character",
                    default=proj4_where_dqc_is_done_default)
  p <- add_argument(p, "--proj4_output_files",
                    help="proj4 string for the output coordinate reference system",
                    type="character",
                    default=proj4_input_obsfiles_default)
