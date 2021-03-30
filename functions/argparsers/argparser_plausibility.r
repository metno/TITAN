  #.............................................................................. 
  # Plausibility check

  p <- add_argument(p, "--vmin",
                    help=paste("minimum allowed value [units of the variable specified]",neg.str),
                    type="character",
                    default="_50")

  p <- add_argument(p, "--vmax",
                    help=paste("maximum allowed value [units of the variable specified]",neg.str),
                    type="character",
                    default="40")

  p <- add_argument(p, "--vminsign",
                    help="minimum allowed value, sign [1=neg, 0=pos]",
                    type="numeric",
                    default=0)

  p <- add_argument(p, "--vmaxsign",
                    help="maximum allowed value, sign [1=neg, 0=pos]",
                    type="numeric",
                    default=0)

