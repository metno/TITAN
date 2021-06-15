# not really spatial data quality control

  p <- add_argument(p, "--code.nometa",
                    help="quality code returned in case of missing metadata",
                    type="integer",
                    default=901)

  p <- add_argument(p, "--code.dem",
                    help="quality code returned in case of dem check fails",
                    type="integer",
                    default=902)

  p <- add_argument(p, "--code.ccrt",
                    help=paste("quality code returned in case of precipitation",
                               "and temperature crosscheck fails"),
                    type="integer",
                    default=301)

  p <- add_argument(p, "--code.p",
    help="quality code returned in case of the check on plausible values fails",
                    type="integer",
                    default=501)

  p <- add_argument(p, "--code.clim",
   help="quality code returned in case of the check on climatological values fails",
                    type="integer",
                    default=502)

# statistics of deviations

  p <- add_argument(p, "--code.buddy",
    help="quality code returned in case of the buddy check fails",
                    type="integer",
                    default=10)

  p <- add_argument(p, "--code.fgt",
   help="quality code returned in case of check against a first-guess field fails",
                    type="integer",
                    default=11)

# statistics of deviations, spatial analysis

  p <- add_argument(p, "--code.sct",
                    help="quality code returned in case of SCT fails",
                    type="integer",
                    default=1)

  p <- add_argument(p, "--code.sct_fg",
                    help="quality code returned in case of SCT with background fails",
                    type="integer",
                    default=2)

  p <- add_argument(p, "--code.sct_dual",
                    help="quality code returned in case of SCT dual fails",
                    type="integer",
                    default=3)

  p <- add_argument(p, "--code.sct_fg_dual",
                    help="quality code returned in case of SCT dual with background fails",
                    type="integer",
                    default=4)

# isolation

  p <- add_argument(p, "--code.iso",
          help="quality code returned in case of isolation check fails",
                    type="integer",
                    default=800)

# blacklist/ keeplist

  p <- add_argument(p, "--code.black",
      help="quality code assigned to observations listed in the blacklist",
                    type="integer",
                    default=999)

  p <- add_argument(p, "--code.keep",
      help="quality code assigned to observations listed in the keep-list",
                    type="integer",
                    default=990)

