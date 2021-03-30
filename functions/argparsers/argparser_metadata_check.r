  #.............................................................................. 
  # metadata check
  p <- add_argument(p, "--zmin",
                    help="minimum allowed elevation in the domain [m amsl]",
                    type="numeric",
                    default=0,
                    short="-z")
  p <- add_argument(p, "--zmax",
                    help="maximum allowed elevation in the domain [m amsl]",
                    type="numeric",
                    default=2500,
                    short="-Z")
