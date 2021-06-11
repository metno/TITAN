  #.............................................................................. 
  # spatial consistency test with background

  p <- add_argument(p, "--sct_fg",
                    help="do the sct_fg",
                    flag=T)

  p <- add_argument(p, "--i.sct_fg",
                    help="number of SCT_fg iterations",
                    type="integer",
                    default=1)

  p <- add_argument(p, "--fglab.sct_fg",
                    help="labels identifying which of the background fields have to be used in the test. The labels refers to the order of the first-guess files. Note: this is the parameter defining N.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--inner_radius.sct_fg",
                    help="radius (m) of the inner circle. Can be a N-vector when a sequence of N tests is needed. Note: this is the parameter defining N.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--outer_radius.sct_fg",
                    help="radius (m) of the outer circle. Either one value valid for all N or N values.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--prio.sct_fg",
                    help="specify priorities. This is a numeric vector, with the dimension of the number of providers. One positive value for each provider. The smaller the value, the higher the priority.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--tpos.sct_fg",
                    help="threshold when the observed value is greater than the spatial trend. Either a single value or one value for each test AND provider (e.g. N tests, M providers; then N*M values; order is: val-1 (test-1,provider-1), val-2 (test-1,provider-2), val-3 (test-1,provider-3), ... ).",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--tneg.sct_fg",
                    help="threshold when the observed value is smaller than the spatial trend. Either a single value or one value for each test AND provider (e.g. N tests, M providers; then N*M values; order is: val-1 (test-1,provider-1), val-2 (test-1,provider-2), val-3 (test-1,provider-3), ... ).",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--eps2.sct_fg",
                    help="relative precision (0-1). Either a single value or one value for each test AND provider (e.g. N tests, M providers; then N*M values; order is: val-1 (test-1,provider-1), val-2 (test-1,provider-2), val-3 (test-1,provider-3), ... ).",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_min_outer.sct_fg",
                    help="minimum number of neighbouring observations (within the outer circle) required. Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_max_outer.sct_fg",
                    help="maximum number of neighbouring observations (within the outer circle). Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--min_horizontal_scale.sct_fg",
                    help=paste("OI, minimum allowed value for the horizontal de-correlation lenght (of the background error correlation) [m]. Either one value valid for all N or N values."),
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--max_horizontal_scale.sct_fg",
                    help=paste("OI, maximum allowed value for the horizontal de-correlation lenght (of the background error correlation) [m]. Either one value valid for all N or N values."),
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--kth_closest_obs_horizontal_scale.sct_fg",
                    help=paste("OI, horizontal de-correlation lenght (of the background error correlation) is computed adaptively based on the distance to the closest observations, as specified by this value. Either one value valid for all N or N values."),
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--vertical_scale.sct_fg",
                    help="OI, vertical de-correlation lenght  (of the background error correlation) [m]. Either one value valid for all N or N values.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--transf.sct_fg",
                    help="apply Box-Cox transformation",
                    flag=T)

  p <- add_argument(p, "--break.sct_fg",
                    help="break the loop if the number of flagged observations in the last iretation (by considering al the test) is euqual to or less than this value.",
                    type="numeric",
                    default=0)

  p <- add_argument(p, "--aggn_radius.sct_fg",
                    help="radius defining the background aggregation area. Either one value valid for all N or N values.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_max_aggn.sct_fg",
                    help="maximum number of points considered inside the background aggregation area. Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--a_delta.sct_fg",
                    help="range of admissibile values is plus/minus this number",
                    type="numeric",
                    default=15)

  p <- add_argument(p, "--v_delta.sct_fg",
                    help="range of valid values is plus/minus this number",
                    type="numeric",
                    default=0.5)

  p <- add_argument(p, "--basic.sct_fg",
                    help="should we use the \"basic\" mode or the \"resistant\" mode?",
                    flag=T)

