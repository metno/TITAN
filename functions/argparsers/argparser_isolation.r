
# Isolation check

  p <- add_argument(p, "iso",
                    help="do the isolation check",
                    flag=T)

  p <- add_argument(p, "--dr.iso",
                    help="check for the number of observation in a dr-by-dr square-box around each observation [m]",
                    type="numeric",
                    default=25000,
                    short="-dI")

  p <- add_argument(p, "--n.iso",
                    help="threshold (number of neighbouring observations) for the identification of isolated observations.",
                    type="integer",
                    default=10,
                    short="-nI")

