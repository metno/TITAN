  #.............................................................................. 
  # spatial consistency test

  p <- add_argument(p, "--sct_dual",
                    help="do SCT dual",
                    flag=T)

  p <- add_argument(p, "--prio.sct_dual",
                    help="specify priorities. This is a numeric vector, with the dimension of the number of providers. One positive value for each provider. The smaller the value, the higher the priority.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--i.sct_dual",
                    help="number of SCT dual iterations. In case of more tests (N>1), then the whole sequence is repeated several times.",
                    type="integer",
                    default=1)

  p <- add_argument(p, "--inner_radius.sct_dual",
                    help="radius (m) of the inner circle. Can be a N-vector when a sequence of N tests is needed. Note: this is the parameter defining N.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--outer_radius.sct_dual",
                    help="radius (m) of the outer circle. Either one value valid for all N or N values.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_min_outer.sct_dual",
                    help="minimum number of neighbouring observations (within the outer circle) required. Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_max_outer.sct_dual",
                    help="maximum number of neighbouring observations (within the outer circle). Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--min_horizontal_scale.sct_dual",
                    help=paste("OI, minimum allowed value for the horizontal de-correlation lenght (of the background error correlation) [m]. Either one value valid for all N or N values."),
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--max_horizontal_scale.sct_dual",
                    help=paste("OI, maximum allowed value for the horizontal de-correlation lenght (of the background error correlation) [m]. Either one value valid for all N or N values."),
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--kth_closest_obs_horizontal_scale.sct_dual",
                    help=paste("OI, horizontal de-correlation lenght (of the background error correlation) is computed adaptively based on the distance to the closest observations, as specified by this value. Either one value valid for all N or N values."),
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--vertical_scale.sct_dual",
                    help="OI, vertical de-correlation lenght  (of the background error correlation) [m]. Either one value valid for all N or N values.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--thr_relinfo.sct_dual",
                    help="SCT dual threshold for the relative information content (typically a number from 0 to 1). Usuful to avoid flagging observations in areas where there is a transition between \"yes\" and \"no\". Either a single value or one value for each test AND provider (e.g. N tests, M providers; then N*M values; order is: val-1 (test-1,provider-1), val-2 (test-1,provider-2), val-3 (test-1,provider-3), ... ).",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--event_thresholds.sct_dual",
                    help="event threshold used to divide the dataset into \"yes, the event happend\" or \"no, the event did not happen\". Either a single value or one value for each test AND provider (e.g. N tests, M providers; then N*M values; order is: val-1 (test-1,provider-1), val-2 (test-1,provider-2), val-3 (test-1,provider-3), ... ).",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--conditions.sct_dual",
                    help="condition used to divide the dataset into \"yes, the event happend\" or \"no, the event did not happen\". One of: \"Eq\", \"Gt\", \"Geq\", \"Lt\", \"Leq\". Either one value valid for all N or N values.",
                    type="character",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--break.sct_dual",
                    help="break the loop if the number of flagged observations in the last iretation (by considering al the test) is euqual to or less than this value.",
                    type="numeric",
                    default=0)
