  #.............................................................................. 
  # duplicates
  p <- add_argument(p, "--no_duplicates",
                    help="remove duplicates from input data",
                    flag=T)
  p <- add_argument(p, "--no_duplicates_radius",
                    help="remove duplicates, matching tolerance in the horizontal (meters)",
                    type="numeric",
                    default=500)
  p <- add_argument(p, "--no_duplicates_vertical_range",
                    help="remove duplicates, matching tolerance for elevation (m)",
                    type="numeric",
                    default=1)

