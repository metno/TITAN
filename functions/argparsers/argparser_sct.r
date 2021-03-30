  #.............................................................................. 
  # spatial consistency test

  p <- add_argument(p, "--sct",
                    help="do the sct",
                    flag=T)

  p <- add_argument(p, "--i.sct",
                    help="number of SCT iterations",
                    type="integer",
                    default=1)

  p <- add_argument(p, "--background_elab_type.sct",
                    help="background elaboration type (\"VerticalProfile\", \"VerticalProfileTheilSen\", \"MeanOuterCircle\", \"MedianOuterCircle\" and \"External\")",
                    type="character",
                    default="VerticalProfileTheilSen")

  p <- add_argument(p, "--inner_radius.sct",
                    help="radius (m) of the inner circle. Can be a N-vector when a sequence of N tests is needed. Note: this is the parameter defining N.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--outer_radius.sct",
                    help="radius (m) of the outer circle. Either one value valid for all N or N values.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--prio.sct",
                    help="specify priorities. This is a numeric vector, with the dimension of the number of providers. One positive value for each provider. The smaller the value, the higher the priority.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--tpos.sct",
                    help="threshold when the observed value is greater than the spatial trend. Either a single value or one value for each test AND provider (e.g. N tests, M providers; then N*M values; order is: val-1 (test-1,provider-1), val-2 (test-1,provider-2), val-3 (test-1,provider-3), ... ).",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--tneg.sct",
                    help="threshold when the observed value is smaller than the spatial trend. Either a single value or one value for each test AND provider (e.g. N tests, M providers; then N*M values; order is: val-1 (test-1,provider-1), val-2 (test-1,provider-2), val-3 (test-1,provider-3), ... ).",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--eps2.sct",
                    help="relative precision (0-1). Either a single value or one value for each test AND provider (e.g. N tests, M providers; then N*M values; order is: val-1 (test-1,provider-1), val-2 (test-1,provider-2), val-3 (test-1,provider-3), ... ).",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_min_outer.sct",
                    help="minimum number of neighbouring observations (within the outer circle) required. Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_max_outer.sct",
                    help="maximum number of neighbouring observations (within the outer circle). Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_min_prof.sct",
                    help="background elaboration. Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--min_elev_diff.sct",
                    help="background elaboration. Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--min_horizontal_scale.sct",
                    help=paste("OI, minimum allowed value for the horizontal de-correlation lenght (of the background error correlation) [m]. Either one value valid for all N or N values."),
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--max_horizontal_scale.sct",
                    help=paste("OI, maximum allowed value for the horizontal de-correlation lenght (of the background error correlation) [m]. Either one value valid for all N or N values."),
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--kth_closest_obs_horizontal_scale.sct",
                    help=paste("OI, horizontal de-correlation lenght (of the background error correlation) is computed adaptively based on the distance to the closest observations, as specified by this value. Either one value valid for all N or N values."),
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--vertical_scale.sct",
                    help="OI, vertical de-correlation lenght  (of the background error correlation) [m]. Either one value valid for all N or N values.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--transf.sct",
                    help="apply Box-Cox transformation",
                    flag=T)

  p <- add_argument(p, "--break.sct",
                    help="break the loop if the number of flagged observations in the last iretation (by considering al the test) is euqual to or less than this value.",
                    type="numeric",
                    default=0)
