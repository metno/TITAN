#+ number of stations within "drmin" units from location "xy" 
nstat<-function(xy,drmin) {
# input
# xy=vector. 2=x;3=y
# output
# 1=nobs
#------------------------------------------------------------------------------
  if (any(is.na(xy))) return(NA)
  i<-which(abs(xy[1]-xtot)<=drmin & 
           abs(xy[2]-ytot)<=drmin)
  return(length(i))
}

#+ summary statistics centered on a point 
statSpat_mapply<-function(i,                           # index over tot-vectors
                          dr,                          # m
                          gamma=-0.0065,               # degC/m
                          adjust_for_elev_diff=F,
                          pmax=200,                    # max number of stations
                          priority=F,
                          statistics="buddy_standard", # OR buddy_event
                          event_threshold=NULL,
                          event_def=NULL) {
#------------------------------------------------------------------------------
# The point has coordinates (obsToCheck_x[i],obsToCheck_y[i],obsToCheck_z[i]). 
# The neighbourhood considered includes the (max pmax) closest dataToUse 
# within a radius dr. The so-called buddies.
# The i-th point is excluded from the reference statistics.
#
# input
# 6 vectors. itot=index; xtot=x; ytot=y; ztot=z; ttot=variable; priotot=priority
#
# -- output "buddy_standard"
# 1 = nobs 
# 2 = maxVertDist[m]
# 3 = mean(buddies)
# 4 = standard_dev(buddies)
#
# -- output "buddy_event"
# 1 = nobs
# 2 = maxVertDist[m]
# 3 = event yes/no at the i-th point (1=yes,0=no)
# 4 = percentage of event=yes among the buddies
#------------------------------------------------------------------------------
  deltax<-abs(obsToCheck_x[i]-dataToUse_x)
  deltay<-abs(obsToCheck_y[i]-dataToUse_y)
  if (priority) {
    jx<-which(deltax<=dr & 
              deltay<=dr &
              dataToUse_i!=obsToCheck_i[i]            &
              (dataToUse_prio<obsToCheck_prio[i] | dataToUse_prio<0) )
  } else {
    jx<-which(deltax<=dr & 
              deltay<=dr &
              dataToUse_i!=obsToCheck_i[i])
  }
  njx<-length(jx)
  if (njx==0) return(c(0,NA,NA,NA))
  if (njx>pmax) {
    disth2<-deltax[jx]*deltax[jx]+deltay[jx]*deltay[jx]
    jx<-jx[order(disth2, decreasing=F)[1:pmax]]
    njx<-length(jx)
  }
  if (njx==0) return(c(0,NA,NA,NA))
  dz<-dataToUse_z[jx]-obsToCheck_z[i]
  dz_mx<-max(abs(dz))
  if (adjust_for_elev_diff) {
    buddies_val<-dataToUse_val[jx]+gamma*dz
  } else {
    buddies_val<-dataToUse_val[jx]
  }
  # buddy check, standard
  if (statistics=="buddy_standard") {
    if (njx==1) return(c(1,dz,buddies_val,NA))
    return(c(njx,dz_mx,mean(buddies_val),sd(buddies_val)))
  # buddy check, based on the definition of a binary event
  } else if (statistics=="buddy_event") {
    if (event_def=="gt") {
      i_eve<-as.integer(obsToCheck_val[i]>event_threshold)
      buddies_eve_yes.prob<-length(which(buddies_val>event_threshold))/njx
    } else if (event_def=="ge") {
      i_eve<-as.integer(obsToCheck_val[i]>=event_threshold)
      buddies_eve_yes.prob<-length(which(buddies_val>=event_threshold))/njx
    } else if (event_def=="lt") {
      i_eve<-as.integer(obsToCheck_val[i]<event_threshold)
      buddies_eve_yes.prob<-length(which(buddies_val<event_threshold))/njx
    } else if (event_def=="le") {
      i_eve<-as.integer(obsToCheck_val[i]<=event_threshold)
      buddies_eve_yes.prob<-length(which(buddies_val<=event_threshold))/njx
    } else {
      i_eve<-NA
      buddies_eve_yes.prob<-NA
    }
    return(c(njx,dz_mx,i_eve,buddies_eve_yes.prob))
  } else {
    return(c(NA,NA,NA,NA))
  }
}

#+ Box-Cox transformation
boxcox<-function(x,lambda) {
  if (lambda==0) {
    return(log(x))
  } else {
    return((x**lambda-1)/lambda)
  }
}

