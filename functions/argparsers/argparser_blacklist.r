  #.............................................................................. 
  # blacklist
  # specified by triple/pairs of numbers: either (lat,lon,IDprovider) OR (index,IDprovider)
  p <- add_argument(p, "--blacklist.lat",
                    help="observation blacklist (latitude)",
                    type="numeric",default=NA,nargs=Inf,short="-bla")
  p <- add_argument(p, "--blacklist.lon",
                    help="observation blacklist (longitude)",
                    type="numeric",default=NA,nargs=Inf,short="-blo")
  p <- add_argument(p, "--blacklist.fll",
                    help="observation blacklist (ID provider)",
                    type="numeric",default=NA,nargs=Inf,short="-bfll")
  p <- add_argument(p, "--blacklist.idx",
                    help="observation blacklist (position in input file)",
                    type="numeric",default=NA,nargs=Inf,short="-bix")
  p <- add_argument(p, "--blacklist.fidx",
                    help="observation blacklist (ID provider)",
                    type="numeric",default=NA,nargs=Inf,short="-bfix")
  p <- add_argument(p, "--blacklist.file_with_sourceIds",
                    help="blacklist based on sourceId (list of source Ids)",
                    type="character",default=NA,nargs=Inf)
  p <- add_argument(p, "--blacklist.file_sourceIds_varname",
                    help="label for the sourceId variable in the input file",
                    type="character",default="sourceId")

