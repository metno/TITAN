  #.............................................................................. 
  # doit flags

  comstr<-" Decide if the test should be applied to all, none or only to a selection of observations based on the provider. Possible values are 0, 1, 2. It is possible to specify either one global value or one value for each provider. Legend: 1 => the observations will be used in the elaboration and they will be tested; 0 => the observations will not be used and they will not be tested; 2 => the observations will be used but they will not be tested."

  p <- add_argument(p, "--doit.buddy",
                    help=paste("customize the buddy-test application.",comstr),
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--doit.sct",
                    help=paste("customize the application of SCT.",comstr),
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--doit.sct_fg",
                    help=paste("customize the sct_fg check application.",comstr),
                    type="integer",
                    default=NA,
                    nargs=Inf)


  p <- add_argument(p, "--doit.clim",
                    help=paste("customize the application of the climatological check.",comstr),
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--doit.dem",
                    help=paste("customize the application of the test of observation elevation against the digital elevation model.",comstr),
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--doit.iso",
                    help=paste("customize the application of the isolation test.",comstr),
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--doit.sct_dual",
                    help=paste("customize the sct_dual check application.",comstr),
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--doit.sct_fg_dual",
                    help=paste("customize the sct_fg_dual check application.",comstr),
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--doit.fgt",
                    help=paste("customize the FGT application.",comstr),
                    type="integer",
                    default=NA,
                    nargs=Inf)
