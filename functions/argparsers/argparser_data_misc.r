  #.............................................................................. 
  # VARIABLE definition 
  p<- add_argument(p, "--variable",
                   help=paste("meteorological variable (T temperature,",
                              " RR precipitation, RH relative humidity,",
                              " SD surface_snow_thickness)"),
                   type="character",
                   default="T",
                   short="-var")

  # parameter for the Box-Cox transformation (rquired for var=RR)
  p <- add_argument(p, "--boxcox.lambda",
                    help="parameter used in the Box-Cox transformation (var RR)",
                    type="numeric",default=0.5,short="-l")
  #.............................................................................. 
  # standard value for the moist adiabatic lapse rate
  p <- add_argument(p, "--gamma.standard",
    help="standard value for the moist adiabatic temperature lapse rate dT/dz",
                    type="numeric",
                    default=-0.0065)

