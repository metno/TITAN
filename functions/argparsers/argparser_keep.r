  #.............................................................................. 
  # keep (keep them)
  # specified by triple/pairs of numbers: either (lat,lon,IDprovider) OR (index,IDprovider)
  p <- add_argument(p, "--keeplist.lat",
                    help="observation keeplist (latitude)",
                    type="numeric",
                    default=NA,
                    nargs=Inf,
                    short="-kla")
  p <- add_argument(p, "--keeplist.lon",
                    help="observation keeplist (longitude)",
                    type="numeric",
                    default=NA,
                    nargs=Inf,
                    short="-klo")
  p <- add_argument(p, "--keeplist.fll",
                    help="observation keeplist (ID provider)",
                    type="numeric",
                    default=NA,
                    nargs=Inf,
                    short="-kfll")
  p <- add_argument(p, "--keeplist.idx",
                    help="observation keeplist (position in input file)",
                    type="numeric",
                    default=NA,
                    nargs=Inf, 
                    short="-kix")
  p <- add_argument(p, "--keeplist.fidx",
                    help="observation keeplist (ID provider)",
                    type="numeric",
                    default=NA,
                    nargs=Inf,
                    short="-kfix")
