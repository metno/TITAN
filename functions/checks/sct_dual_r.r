#+ spatial consistency test for dichotomous (yes/no) variables
sct_dual_r <- function( argv, 
                        data,
                        z, 
                        dqcflag) {
#------------------------------------------------------------------------------
# split observations into yes/no and perform spatial consistency check 
# See the wiki 
#==============================================================================

  cat( paste0( "sct_dual (code=", argv$code.sct_dual, ")\n"))

  ndata <- length( data$lat)
  
  nfin  <- length( argv$input.files)

  # number of observation providers
  M <- nfin

  # number of tests
  if ( ( N <- length( argv$inner_radius.sct_dual)) == 0) return( FALSE)

  nsus <- vector( mode="numeric", length=N)

  debug <- F

  if ( length( argv$doit.sct_dual) != M) 
    argv$doit.sct_dual <- rep( argv$doit.sct_dual[1], M)
  if ( length( argv$prio.sct_dual) != M) 
    argv$prio.sct_dual <- rep( argv$prio.sct_dual[1], M)

  if ( length( argv$event_thresholds.sct_dual) != (M*N)) 
    argv$event_thresholds.sct_dual <- rep( argv$event_thresholds.sct_dual[1], N*M)
  if ( length( argv$conditions.sct_dual) != N) 
    argv$conditions.sct_dual <- rep( argv$conditions.sct_dual[1], N)

  if ( length( argv$thr_relinfo.sct_dual) != (M*N)) 
    argv$thr_relinfo.sct_dual <- rep( argv$thr_relinfo.sct_dual[1], N*M)

  if ( length( argv$num_min_outer.sct_dual) != N) 
    argv$num_min_outer.sct_dual <- rep( argv$num_min_outer.sct_dual[1], N)
  if ( length( argv$num_max_outer.sct_dual) != N) 
    argv$num_max_outer.sct_dual <- rep( argv$num_max_outer.sct_dual[1], N)
  if ( length( argv$outer_radius.sct_dual) != N) 
    argv$outer_radius.sct_dual <- rep( argv$outer_radius.sct_dual[1], N)
  if ( length( argv$min_horizontal_scale.sct_dual) != N) 
    argv$min_horizontal_scale.sct_dual <- rep( argv$min_horizontal_scale.sct_dual[1], N)
  if ( length( argv$max_horizontal_scale.sct_dual) != N) 
    argv$max_horizontal_scale.sct_dual <- rep( argv$max_horizontal_scale.sct_dual[1], N)
  if ( length( argv$kth_closest_obs_horizontal_scale.sct_dual) != N) 
    argv$kth_closest_obs_horizontal_scale.sct_dual <- rep( argv$kth_closest_obs_horizontal_scale.sct_dual[1], N)
  if ( length( argv$vertical_scale.sct_dual) != N) 
    argv$vertical_scale.sct_dual <- rep( argv$vertical_scale.sct_dual[1], N)

  #............................................................................
  # set doit/prio vectors
  doit <- vector( length=ndata, mode="numeric"); doit[]<-NA
  prio <- vector( length=ndata, mode="numeric"); prio[]<-NA
  for (f in 1:M) {
    ix <- which( data$prid == argv$prid[f])
    if ( length(ix) == 0) next
    doit[ix] <- argv$doit.sct_dual[f]
    prio[ix] <- argv$prio.sct_dual[f]
  }
  prio_unique <- sort( unique( prio, na.rm=T), decreasing=T)
  rm( ix)

  #............................................................................
  # test

  flagok_pre <- !is.na( data$lat) & !is.na( data$lon) & !is.na( data$value) &
              doit != 0 & !is.na( prio)

  for (i in 1:argv$i.sct_dual) {

    nsus[]<-0

    for (j in 1:N) {

      t     <- vector( length=ndata, mode="numeric"); t[]<-NA
      r_eve <- vector( length=ndata, mode="numeric"); r_eve[]<-NA

      for (f in 1:M) {
        ix <- which( data$prid == argv$prid[f])
        if ( length(ix) == 0) next
        t[ix]     <- argv$thr_relinfo.sct_dual[(j-1)*M+f]
        r_eve[ix] <- argv$event_thresholds.sct_dual[(j-1)*M+f]
      }
      rm( ix)

      for (k in 1:length(prio_unique)) {

        # use only (probably) good observations with doit!=0
        ix <- which( ( ( is.na( dqcflag) | dqcflag == argv$code.keep) &
                         flagok_pre))
        t0a <- Sys.time()

        obsToCheck_n<-length(ix)

        if (obsToCheck_n>0) {
          # define global 1D vector used in statSpat (1D for fast access)
          obsToCheck_lon  <- data$lon[ix]
          obsToCheck_lat  <- data$lat[ix]
          obsToCheck_z    <- z[ix]
          obsToCheck_val  <- data$value[ix]
          obsToCheck_t    <- t[ix]
          obsToCheck_r    <- r_eve[ix]
          obsToCheck_chk  <- rep( 0, obsToCheck_n)
          # check only those observations with priorities geq than this
          obsToCheck_chk[prio[ix]>=prio_unique[k]] <- 1
          cat("\n")
          res <- sct_dual(
                      points = Points( obsToCheck_lat, 
                                       obsToCheck_lon, 
                                       obsToCheck_z),
                      obsToCheck_val,
                      obsToCheck_chk,
                      obsToCheck_r,
                      argv$conditions.sct_dual[j],
                      argv$num_min_outer.sct_dual[j],
                      argv$num_max_outer.sct_dual[j],
                      argv$inner_radius.sct_dual[j],
                      argv$outer_radius.sct_dual[j],
                      100,
                      argv$min_horizontal_scale.sct_dual[j],
                      argv$max_horizontal_scale.sct_dual[j],
                      argv$kth_closest_obs_horizontal_scale.sct_dual[j],
                      argv$vertical_scale.sct_dual[j],
                      obsToCheck_t,
                      debug)
 
          flag <- res

          # suspect if: 
          sus<-which( flag[1:obsToCheck_n] == 1 &
                      is.na(dqcflag[ix]) &
                      doit[ix]==1 )

          # set dqcflag
          if (length(sus)>0) dqcflag[ix[sus]] <- argv$code.sct_dual

        } else {
          cat( "no valid observations left, no sct check\n")
        }
        nsus[j] <- ifelse( exists("sus"), length(sus), 0)
        t1a  <- Sys.time()
        str  <- " (#TOT "
        str1 <- ""
        for (f in 1:M) {
          if (f>1) str<-paste0(str,"; ")
          str <- paste0( str, "prid", argv$prid[f], "=", 
                  length( which( dqcflag==argv$code.sct_dual & 
                                 data$prid==argv$prid[f])))
          str1 <- paste0( str1, "::prid ", argv$prid[f]," ",
                                "prio=", argv$prio.sct_dual[f],",",
                                "t=", argv$thr_relinfo.sct_dual[(j-1)*M+f],",",
                                "r=", argv$event_thresholds.sct_dual[(j-1)*M+f])
        }
        str<-paste0(str,")")
        cat( paste0( "iteration=", i,
                     "/test=", j,
                     "/prio>=", prio_unique[k],
                     "/dqc param:",
                     "inner_rad=", argv$inner_radius.sct_dual[j],",",
                     "outer_rad=", argv$outer_radius.sct_dual[j],",",
                     "nmin_outer=", argv$num_min_outer.sct_dual[j],",",
                     "nmax_outer=", argv$num_max_outer.sct_dual[j],",",
                     "min_hscale=", argv$min_horizontal_scale.sct_dual[j],",",
                     "max_hscale=", argv$max_horizontal_scale.sct_dual[j],",",
                     "kth_dh=", argv$kth_closest_obs_horizontal_scale.sct_dual[j],",",
                     "dz=", argv$vertical_scale.sct_dual[j],",",
                     "cond=", argv$conditions.sct_dual[j],",",
                     str1,
                     "/time ",round(t1a-t0a,1),attr(t1a-t0a,"unit"),"\n"))
        cat( paste0(nsus[j]," new suspect observations",str,"\n"))
        rm(str)
      } # end for k
    } # end for j
    if ( sum(nsus) <= argv$break.sct_dual) break
  }  # end for i
  cat("+---------------------------------+\n")
  #
  return(dqcflag)
}

