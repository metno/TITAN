  #.............................................................................. 
  # Buddy-check

  p <- add_argument(p, "--buddy",
                    help="do the buddy check",
                    flag=T)

  p <- add_argument(p, "--i.buddy",
                    help="number of iterations",
                    type="integer",
                    default=1)

  p <- add_argument(p, "--prio.buddy",
                    help="specify priorities. This is a numeric vector, with the dimension of the number of providers. One positive value for each provider. The smaller the value, the greater the priority. Priorities are used only in the first round of the check, when each observation is compared only against other observations having a priority equal to or grater than its own.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--break.buddy",
                    help="break the loop if the number of flagged observations in the last iretation (by considering al the test) is euqual to or less than this value.",
                    type="integer",
                    default=0)

  p <- add_argument(p, "--transf.buddy",
                    help="apply Box-Cox transformation",
                    flag=T)

  p <- add_argument(p, "--inner_radius.buddy",
                    help="radius (m) of the inner circle. Can be a N-vector when a sequence of N tests is needed. Note: this is the parameter defining N.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--outer_radius.buddy",
                    help="radius (m) of the outer circle. Either one value valid for all N or N values.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_min_outer.buddy",
                    help="minimum number of neighbouring observations (within the outer circle) required. Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_max_outer.buddy",
                    help="maximum number of neighbouring observations (within the outer circle). Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--num_min_prof.buddy",
                    help="minimum number of observations in the outer circle to compute a vertical profile. Either one value valid for all N or N values.",
                    type="integer",
                    default=NA,
                    nargs=Inf)

  p <- add_argument(p, "--min_elev_diff.buddy",
                    help="minimum elevation difference in the outer circle to compute a vertical profile. Either one value valid for all N or N values.",
                    type="numeric",
                    default=NA,
                    nargs=Inf)


