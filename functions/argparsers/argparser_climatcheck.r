  # default based on Norwegian hourly temperature from 2010-2017
  p <- add_argument(p, "--vmin.clim",
                    help=paste("minimum allowed value [units of the variable specified]",neg.str),
                    type="character",
                    nargs=12,
                    default=c("_45","_45","_40","_35","_20","_15","_10","_15","_15","_20","_35","_45"))
  p <- add_argument(p, "--vmax.clim",
                    help=paste("maximum allowed value [units of the variable specified]",neg.str),
                    type="character",
                    nargs=12,
                    default=c("20","20","25","25","35","35","40","40","35","30","25","20"))
  p <- add_argument(p, "--month.clim",
                    help="month (number 1-12)",
                    type="numeric",
                    short="-mC",
                    default=NA)
  p <- add_argument(p, "--vminsign.clim",
                    help="minimum allowed value, sign [1=neg, 0=pos]",
                    type="numeric",
                    nargs=12,
                    default=c(0,0,0,0,0,0,0,0,0,0,0,0))
  p <- add_argument(p, "--vmaxsign.clim",
                    help="maximum allowed value, sign [1=neg, 0=pos]",
                    type="numeric",
                    nargs=12,
                    default=c(0,0,0,0,0,0,0,0,0,0,0,0))

