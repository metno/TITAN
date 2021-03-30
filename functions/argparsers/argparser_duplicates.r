  #.............................................................................. 
  # duplicates
  p <- add_argument(p, "--no_duplicates",
                    help="remove duplicates from input data",
                    flag=T)
  p <- add_argument(p, "--dup.match_tol_x",
                    help="remove duplicates, matching tolerance for lat and lon (degrees)",
                    type="numeric",
                    default=0.00001)
  p <- add_argument(p, "--dup.match_tol_z",
                    help="remove duplicates, matching tolerance for elevation (m)",
                    type="numeric",
                    default=1)

