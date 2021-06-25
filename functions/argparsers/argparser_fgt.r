  #.............................................................................. 
  # FGT

  p <- add_argument(p, "--fgt",
                    help="do the first-guess test",
                    flag=T)

  p <- add_argument(p, "--prio.fgt",
                    help="specify priorities. This is a numeric vector, with the dimension of the number of providers. One positive value for each provider. The smaller the value, the higher the priority.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--fgt_fglab.fgt",
                    help="labels identifying which of the background fields have to be used in the test. The labels refers to the order of the first-guess files. Note: this is the parameter defining N.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--circle_radius.fgt",
                    help="radius (m) of the circle. Can be a N-vector when a sequence of N tests is needed.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--i.fgt",
                    help="number of FGT iterations",
                    type="integer",
                    default=1)

  p <- add_argument(p, "--tpos.fgt",
                    help="threshold when the observed value is greater than the spatial trend. Either a single value or one value for each test AND provider (e.g. N tests, M providers; then N*M values; order is: val-1 (test-1,provider-1), val-2 (test-1,provider-2), val-3 (test-1,provider-3), ... ).",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--tneg.fgt",
                    help="threshold when the observed value is smaller than the spatial trend. Either a single value or one value for each test AND provider (e.g. N tests, M providers; then N*M values; order is: val-1 (test-1,provider-1), val-2 (test-1,provider-2), val-3 (test-1,provider-3), ... ).",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_min_circle.fgt",
                    help="minimum number of neighbouring observations (within the outer circle) required. Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_max_circle.fgt",
                    help="maximum number of neighbouring observations (within the outer circle). Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--aggn_radius.fgt",
                    help="radius defining the background aggregation area. Either one value valid for all N or N values.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_max_aggn.fgt",
                    help="maximum number of points considered inside the background aggregation area. Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--transf.fgt",
                    help="apply Box-Cox transformation",
                    flag=T)

  p <- add_argument(p, "--break.fgt",
                    help="break the loop if the number of flagged observations in the last iretation (by considering al the test) is euqual to or less than this value.",
                    type="numeric",
                    default=0)

  p <- add_argument(p, "--a_delta.fgt",
                    help="range of admissibile values is plus/minus this number",
                    type="numeric",
                    default=15)

  p <- add_argument(p, "--v_delta.fgt",
                    help="range of valid values is plus/minus this number",
                    type="numeric",
                    default=0.5)

  p <- add_argument(p, "--a_fact.fgt",
                    help="factor used in the definition of the range of admissibile values",
                    type="numeric",
                    default=0.5)

  p <- add_argument(p, "--v_fact.fgt",
                    help="factor used in the definition of the range of admissibile values",
                    type="numeric",
                    default=0.1)

  p <- add_argument(p, "--basic.fgt",
                    help="should we use the \"basic\" mode or the \"resistant\" mode?",
                    flag=T)
