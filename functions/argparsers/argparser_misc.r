  #.............................................................................. 
  # titan path (use for the sct with smart boxes)
  p <- add_argument(p, "--titan_path",
                    help="path to the directory where the TITAN code is",
                    type="character",
                    default=NULL,
                    short="-tip")
  #.............................................................................. 
  # titanlib path 
  p <- add_argument(p, "--titanlib_path",
                    help="path to the directory where the TITAN code is",
                    type="character",
                    default=NULL,
                    short="-tlp")
  #.............................................................................. 
  # DEBUG
  p <- add_argument(p, "--debug",
                    help="debug mode",
                    flag=T,
                    short="-dbg")
  p <- add_argument(p, "--debug.dir",
                    help="directory for debug output",
                    type="character",
                    default=".",
                    short="-dbgd")
  p <- add_argument(p, "--verbose",
                    help="debug mode",
                    flag=T,
                    short="-v")
  #.............................................................................. 
  # run on several cores 
  p <- add_argument(p, "--cores",
                    help="set the number of cores for parallel runs. Rpackage \"parallel\" required. 0 stands for \"use detectCores\". Default do not use it.",
                    type="numeric",
                    default=NA)

