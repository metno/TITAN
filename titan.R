#!/usr/bin/env Rscript
# + TITAN - spatial data quality control for in-situ observations
# mailto: cristianl@met.no
# # https://github.com/metno/TITAN
#
# command line:
#  >titan.R input_file output_file [options]
# to list available options:
#  >titan.R --help 
#-----------------------------------------------------------------------------
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#-----------------------------------------------------------------------------
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("sp"))
suppressPackageStartupMessages(library("raster"))
suppressPackageStartupMessages(library("rgdal"))
options(warn = 2, scipen = 999)
#------------------------------------------------------------------------------
# FUNCTIONS

# + manage fatal error
boom<-function(str=NULL) {
  print("Fatal Error:")
  if (!is.null(str)) print(str)
  quit(status=1)
}

# auxiliary function to keep/blacklist observations
setCode_lonlat<-function(lonlat,code) {
# lonlat. vector. 1=lon; 2=lat
  ix<-which(datatmp$lon==lonlat[1] & datatmp$lat==lonlat[2])
  if (length(ix)>0)  {
    aux[ix]<-code
    assign("aux",aux,envir=.GlobalEnv)
  }
  return(length(ix))
}

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

#+ summary statistics in a box of /pm "drmin" centered on location "ixyz" 
statSpat_mapply<-function(i,
                          drmin,
                          gamma=-0.0065,
                          priority=F,
                          statistics="buddy_standard", #buddy_event
                          event_threshold=NULL,
                          event_def=NULL) {
# input
# 6 vectors. itot=index; xtot=x; ytot=y; ztot=z; ttot=variable; priotot=priority
# output
# "buddy_standard"
# 1=nobs;2=maxVertDist[m];3=mean(temp);4=sd(temp)
# "buddy_event"
# 1=nobs;2=maxVertDist[m];
# 3=event occurrence (1=yes,0=no);
# 4=frequency of event=yes in the surroundings
# NOTE: temperature is adjusted for elevation differences (gamma=-0.0065 K/m)
#------------------------------------------------------------------------------
#    ixyztp_tot<-cbind(itot,xtot,ytot,ztot,ttot,priotot)
  if ( any(is.na(itot)) | any(is.na(xtot)) | any(is.na(ytot)) | 
       any(is.na(ztot)) | any(is.na(ttot)) | any(is.na(priotot)) ) {
    return(c(NA,NA,NA,NA))
  }
  deltax<-abs(xtot[i]-xtot)
  deltay<-abs(ytot[i]-ytot)
  if (priority) {
    ixx0<-which(deltax<=drmin & 
                deltay<=drmin &
                itot!=itot[i]            &
                (priotot<priotot[i] | priotot<0) )
  } else {
    ixx0<-which(deltax<=drmin & 
                deltay<=drmin &
                itot!=itot[i])
  }
  if (length(ixx0)==0) return(c(0,NA,NA,NA))
  ixx<-ixx0[which(sqrt(deltax[ixx0]*deltax[ixx0]+
                       deltay[ixx0]*deltay[ixx0])<=drmin)]
  rm(ixx0)
  if (length(ixx)==0) return(c(0,NA,NA,NA))
  dz<-ztot[ixx]-ztot[i]
  if (argv$variable=="T") {
    t2use<-ttot[ixx]+gamma*dz
  } else {
    t2use<-ttot[ixx]
  }
  if (length(ixx)==1) {
    if (statistics=="buddy_standard") {
      return(c(1,dz,t2use,NA))
    } else {
      return(c(1,NA,NA,NA))
    }
  }
  if (statistics=="buddy_standard") {
    tmean<-mean(t2use) 
    tsd<-sd(t2use) 
    dz_mx<-max(abs(dz))
    return(c(length(ixx),round(dz_mx,0),round(tmean,1),round(tsd,3)))
  } else if (statistics=="buddy_event") {
    if (event_def=="gt") {
      eve<-as.integer(ifelse(ttot[i]>event_threshold,1,0))
      eve_yes.prob<-length(which(t2use>event_threshold))/length(ixx)
    } else if (event_def=="ge") {
      eve<-as.integer(ifelse(ttot[i]>=event_threshold,1,0))
      eve_yes.prob<-length(which(t2use>=event_threshold))/length(ixx)
    } else if (event_def=="lt") {
      eve<-as.integer(ifelse(ttot[i]<event_threshold,1,0))
      eve_yes.prob<-length(which(t2use<event_threshold))/length(ixx)
    } else if (event_def=="le") {
      eve<-as.integer(ifelse(ttot[i]<=event_threshold,1,0))
      eve_yes.prob<-length(which(t2use<=event_threshold))/length(ixx)
    } else {
      eve<-NA
      eve_yes.prob<-NA
    }
    dz_mx<-max(abs(dz))
    return(c(length(ixx),round(dz_mx,0),eve,eve_yes.prob))
  } else {
    return(c(NA,NA,NA,NA))
  }
}

#+ Correction for wind-induced undercatch of precipitation (Wolff et al., 2015)
wolff_correction<-function(par,rr,t2m,ws10m) {
# Wolff MA, Isaksen K, Petersen-Øverleir A, Ødemark K, Reitan T, Brækkan R. 
# Derivation of a new continuous adjustment function for correcting 
# wind-induced loss of solid precipitation: results of a Norwegian field study.
# Hydrology and Earth System Sciences. 2015 Feb 1;19(2):951.
#------------------------------------------------------------------------------
  theta<-par[1]
  beta<-par[2]
  tau1<-par[3]
  tau2<-par[4]
  Ttau<-par[5]
  stau<-par[6]
  sig1<-par[7]
  sig2<-par[8]
  Tsig<-par[9]
  ssig<-par[10]
  n<-length(t2m)
  aux0<-exp((t2m-Ttau)/stau)
  aux1<-aux0/(1+aux0)
  f<- ( 1 - tau1 - (tau2-tau1) * aux1 ) * exp(-(ws10m/theta)**beta) + 
     tau1 + (tau2-tau1) * aux1
  aux2<-exp((t2m-Tsig)/ssig)
  s<-sig1+(sig2-sig1)*(aux2)/(1+aux2)
  rr.cor<-rr*1./f
  return(list(f=f,sigma=s,rr.cor=rr.cor))
}

#+ Box-Cox transformation
boxcox<-function(x,lambda) {
  if (lambda==0) {
    return(log(x))
  } else {
    return((x**lambda-1)/lambda)
  }
}

#+ vertical profile of temperature (Frei, 2014)
tvertprof<-function(z,t0,gamma,a,h0,h1i) {
# ref:
# Frei, C. (2014). Interpolation of temperature in a mountainous region 
#  using nonlinear profiles and non‐Euclidean distances.
#  International Journal of Climatology, 34(5), 1585-1605.
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
#  a= numeric. inversion (spatial) length
#  h0= numeric. z where inversion starts [m]
#  h1i= numeric. h0+h1i is z where inversion stops [m]
#       (Frei uses h1 directly, I use an increment to h0 so to avoid ending
#        up with h1<=h0 during the optimization)
# Output
#  t= array. temperature [K or degC]
#------------------------------------------------------------------------------
  t<-z
  t[]<-NA
  h1<-h0+abs(h1i)
  z.le.h0<-which(z<=h0)
  z.ge.h1<-which(z>=h1)
  z.in<-which(z>h0 & z<h1)
  if (length(z.le.h0)>0)
   t[z.le.h0]<-t0-gamma*z[z.le.h0]-a 
  if (length(z.ge.h1)>0)
   t[z.ge.h1]<-t0-gamma*z[z.ge.h1] 
  if (length(z.in)>0)
   t[z.in]<-t0-gamma*z[z.in]-a/2*(1+cos(pi*(z[z.in]-h0)/(h1-h0)))
  return(t)
}

#+ cost function used for optimization of tvertprof parameter
tvertprof2opt<-function(par) {
  te<-tvertprof(z=zopt,t0=par[1],gamma=par[2],a=par[3],h0=par[4],h1i=par[5])
  return(log((mean((te-topt)**2))**0.5))
}

#+ vertical profile of temperature (linear)
tvertprof_basic<-function(z,t0,gamma) {
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
# Output
#  t= array. temperature [K or degC]
#------------------------------------------------------------------------------
  return(t0+gamma*z)
}

#+ cost function used for optimization of tvertprof parameter
tvertprofbasic2opt<-function(par) {
  te<-tvertprof_basic(z=zopt,t0=par[1],gamma=argv$gamma.standard)
  return(log((mean((te-topt)**2))**0.5))
}

#+ cost function used for optimization of vertprof parameter
vertprofbasic2opt<-function(par,vert_coord,gamma,obs) {
  pred<-tvertprof_basic(z=vert_coord,t0=par[1],gamma=gamma)
  return(log((mean((pred-obs)**2))**0.5))
}

#+ cost function used for optimization of vertprof parameter
vertprof2opt<-function(par,vert_coord,obs) {
  pred<-tvertprof(z=vert_coord,
                  t0=par[1],
                  gamma=par[2],
                  a=par[3],
                  h0=par[4],
                  h1i=par[5])
  return(log((mean((pred-obs)**2))**0.5))
}


#+ function to plot points
plotp<-function(x,y,val,br,col,
                map=NULL,map.br=NULL,map.col=NULL,
                xl=NULL,yl=NULL) {
  if (is.null(xl)) {
    plot(x,y)
  } else {
    plot(x,y,xlim=xl,ylim=yl)
  }
  if (!is.null(map)) {
    if (!is.null(map.br)) {
      image(map,add=T,breaks=map.br,col=map.col)
    } else {
      image(map,add=T)
    }
  } 
  for (i in 1:length(col)) {
    ix<-which(val>=br[i] & val<br[i+1])
    if (length(ix)>0) 
      points(x[ix],y[ix],col=col[i],pch=19)
    if (i==1) {
      legstr<-paste("<=",br[1],sep="")
    } else if (i==length(col)) {
      legstr<-c(legstr,paste(">=",br[length(col)],sep=""))
    } else {
      legstr<-c(legstr,br[i])
    }
  }
  legend(x="bottomright",fill=rev(col),legend=rev(legstr))
}

#+ SCT - spatial consistency test
sct<-function(ixynp,
              nmin=50,
              dzmin=30,
              Dhmin=10000,
              Dz=200,
              Dz.bg=1500,
              eps2.bg=0.5,
              eps2=NA,
              T2=NA,
              T2pos=NA,
              T2neg=NA,
              sus.code=4,
              faster=F) {
# ref:
#  Lussana, C., Uboldi, F., & Salvati, M. R. (2010). A spatial consistency 
#   test for surface observations from mesoscale meteorological networks.
#   Quarterly Journal of the Royal Meteorological Society, 136(649), 1075-1088.
# input
#  ixynp= vector(4). 1=box identifier on the grid; 
#                    2/3=easting/northing coord (center of the box);
#                    4=number of stations within the box
#  NOTE: stations are associated to a box before running this function
#  nmin= numeric. minimum number of stations to fit a vertical profile
#  dzmin= numeric. minimum elevation range to fit a vertical profile [m]
#  Dhmin= numeric. minimum value for OI horizontal decorellation length [m]
#  Dz= numeric. OI vertical decorellation length [m]
#  Dz.bg= numeric. OI vertical decorellation length 
#         for the construction of the background (RH,RR,SD) [m]
#  eps2.bg= numeric. OI ratio between obs_err_variance/backg_err_variance
#  eps2= numeric. OI ratio between obs_err_variance/backg_err_variance
#  T2=numeric. SCT threshold. (obs-pred)^2/(varObs+varPred)^2 > T2, suspect!
#  T2pos=numeric. SCT threshold. (obs-pred)>=0 AND 
#                                (obs-pred)^2/(varObs+varPred)^2 > T2, suspect!
#  T2neg=numeric. SCT threshold. (obs-pred) <0 AND
#                                (obs-pred)^2/(varObs+varPred)^2 > T2, suspect!
#  NOTE: if T2pos and T2neg are specified then they are used, otherwise use T2
#  sus.code=numeric. identifier code for suspect observation
# output
#  number of rejected stations. (special cases: (i) NA if the function has not
#   been applied; (ii) -1 if just one station in the domain
#  
#  NOTE: "dqcflag" global variable is also updated
#------------------------------------------------------------------------------
  # check for something strange with the number of stations
  if (is.na(ixynp[4]) | is.null(ixynp[4]) | !is.finite(ixynp[4]) ) return(NA)
  # j, index for the stations in the box
  if (argv$usefge.sct) {
    j<-which(itot==ixynp[1] & !is.na(fge.mu[ix]))
  } else if (argv$usefg.sct) {
    j<-which(itot==ixynp[1] & !is.na(fg[ix]))
  } else {
    j<-which(itot==ixynp[1])
  }
  if (length(j)==0) {
    if (argv$verbose) print(paste("warning: no stations in box",ixynp[1]))
    return(-1)
  }
  # case of just one station
  if (ixynp[4]==1) {
    sctpog[ix[j]]<--1
    assign("sctpog",sctpog,envir=.GlobalEnv)
    if (argv$verbose) print(paste("warning: just one station in box",ixynp[1]))
    return(-1)
  }
  # "zopt" and "topt" are set as global variables, so that they can be used by 
  #   other functions in the optimization procedure
  assign("zopt",ztot[j],envir=.GlobalEnv)
  dz<-as.numeric(quantile(zopt,probs=0.95))-
      as.numeric(quantile(zopt,probs=0.05))
  assign("topt",ttot[j],envir=.GlobalEnv)
  # distance matrices
  disth<-(outer(xtot[j],xtot[j],FUN="-")**2.+
          outer(ytot[j],ytot[j],FUN="-")**2.)**0.5
  distz<-abs(outer(ztot[j],ztot[j],FUN="-"))
  # if the background is derived from a first-guess field
  if (argv$usefge.sct) {
    tb<-fge.mu[ix[j]]
  } else if (argv$usefg.sct) {
    tb<-fg[ix[j]]
  # no external first-guess AND variable is temperature
  } else if (argv$variable=="T") {
    # if region is too flat or not enough obs, then go for a basic profile
    if (dz<dzmin | ixynp[4]<nmin) {
      par<-c(mean(topt))
      opt<-optimize(f=tvertprofbasic2opt,interval=c(argv$vmin,argv$vmax))
      tb<-tvertprof_basic(zopt,
                          t0=opt$minimum,
                          gamma=argv$gamma.standard)
    # otherwise, fit a vertical profile of temperature
    } else {
      par<-c(mean(topt),argv$gamma.standard,5,
             as.numeric(quantile(zopt,probs=0.1)),
             as.numeric(quantile(zopt,probs=0.9)))
      opt<-optim(par,tvertprof2opt)
      tb<-tvertprof(zopt,
                    t0=opt$par[1],
                    gamma=opt$par[2],
                    a=opt$par[3],
                    h0=opt$par[4],
                    h1i=opt$par[5])
    }
  # no external first-guess AND variable is relative humidity, prec, snow-depth
  } else if (argv$variable %in% c("RH","RR","SD")) {
    tb<-vector(mode="numeric",length=ixynp[4])
    # if less than five observations in the box... not that much one can do 
    if (ixynp[4]<5) {
      tb[]<-mean(topt)
    } else {
      # set the reference length for the box-dependent large-scale 
      Dh.bg<-max((3*Dhmin),
                 mean(apply(cbind(1:nrow(disth),disth),
                     MARGIN=1,
                     FUN=function(x)
          {as.numeric(quantile(x[which(!((1:length(x))%in%c(1,(x[1]+1))))],
                      probs=0.75))})))
      # background error correlation matrix
      S.bg<-exp(-0.5*(disth/Dh.bg)**2.-0.5*(distz/Dz.bg)**2.)
      if (argv$laf.sct) {
        S.bg<-S.bg * 
              (1-(1-argv$lafmin.sct)*abs(outer(laftot[j],laftot[j],FUN="-")))
      }
      # background=OI(obs-mean(obs)) with Dh.bg,Dz.bg,eps2.bg
      tb[]<-mean(topt)+
            crossprod(S.bg,
                      crossprod(chol2inv(chol((S.bg+eps2.bg*diag(ncol(S.bg))))),
                      (topt-mean(topt))))
      rm(S.bg,Dh.bg)
    }
  }
  if (any(tb<sctvmin)) tb[which(tb<sctvmin)]<-sctvmin
  if (any(tb>sctvmax)) tb[which(tb>sctvmax)]<-sctvmax
  # OI for SCT (Lussana et al., 2010)
  # initialize variables
  # expand eps2 to eps2.vec
  eps2.vec<-vector(length=length(j))
  t2.vec<-vector(length=length(j))
  t2pos.vec<-vector(length=length(j))
  t2neg.vec<-vector(length=length(j))
  t2.vec[]<-NA
  t2pos.vec[]<-NA
  t2neg.vec[]<-NA
  for (f in 1:nfin) { 
    if (any(pridtot[j]==argv$prid[f])) {
      ipx<-which(pridtot[j]==argv$prid[f])
      eps2.vec[ipx]<-eps2[f]
      if (!any(is.na(T2))) 
        t2.vec[ipx]<-T2[f]
      if (!any(is.na(T2pos))) 
        t2pos.vec[ipx]<-T2pos[f]
      if (!any(is.na(T2neg))) 
        t2neg.vec[ipx]<-T2neg[f]
    }
  }

  # Use the smart-box approach implemented in C
  if(argv$smartbox.sct & argv$variable == "T") {
    ns = length(j)
    # TODO: Deal with sel vs sel2check in sct_wrapper
    out<-.C("sct_smart_boxes",ns=as.integer(ns),
                                       x=as.double(xtot[j]),
                                       y=as.double(ytot[j]),
                                       z=as.double(ztot[j]),
                                       t=as.double(ttot[j]),
                                       nmax=as.integer(1000),
                                       nmin=as.integer(100),
                                       nminprof=as.integer(nmin),
                                       dzmin=as.double(dzmin),
                                       dhmin=as.double(Dhmin),
                                       dz=as.double(Dz),
                                       t2pos=as.double(t2pos.vec),
                                       t2neg=as.double(t2neg.vec),
                                       eps2=as.double(eps2.vec),
                                       flags=as.integer(1:ns),
                                       sct=as.double(1:ns),
                                       rep=as.double(1:ns),
                                       boxids=as.integer(1:ns))

    keeping = which(out$flags == 0)
    remove = which(out$flags == 1)
    dqcflag[ix[j[remove]]] <- sus.code

    assign("dqcflag",dqcflag,envir=.GlobalEnv)
    corep[ix[j]] <- out$rep
    corep[ix[j[remove]]] <- NA
    assign("corep",corep,envir=.GlobalEnv)
    sctpog[ix[j]] <- out$sct
    assign("sctpog",sctpog,envir=.GlobalEnv)

    return(length(remove))
  }
  else {
    #  dqc flags
    dqctmp<-dqcflag[ix[j]]
    doittmp<-doit[ix[j]]
    #  probability of gross error (pog)
    pog<-dqctmp
    pog[]<-NA
    # set to optimal Dh for the SCT
    # either Dhmin or the average 10-percentile of the set of distances between a 
    # station and all the others
    Dh<-max(Dhmin,
            mean(apply(cbind(1:nrow(disth),disth),
                       MARGIN=1,
                       FUN=function(x){
         as.numeric(quantile(x[which(!((1:length(x))%in%c(1,(x[1]+1))))],
                            probs=0.1))})))
    # background error correlation matrix
    S<-exp(-0.5*(disth/Dh)**2.-0.5*(distz/Dz)**2.)
    if (argv$laf.sct) {
      S<-S * (1-(1-argv$lafmin.sct)*abs(outer(laftot[j],laftot[j],FUN="-")))
    }
    # S+eps2I
    diag(S)<-diag(S)+eps2.vec
    # innovation
    d<-topt-tb
    # select observations used in the test
    sel<-which(is.na(dqctmp) | dqctmp==argv$keep.code)
    # select observations to test 
    sel2check<-which(is.na(dqctmp) & doittmp==1)
    first<-T
    # loop over SCT iterations 
    # NOTE: SCT flags at most one observation, iterate until all observations pass
    while (length(sel)>1) { 
      # first iteration, inver the matrix
      if (first) {
        SRinv<-chol2inv(chol(S))
        # from S+R go back to S
        diag(S)<-diag(S)-eps2.vec
        first<-F
      } else if (length(indx)>1) {
        S<-S[-indx,-indx]
        eps2.vec<-eps2.vec[-indx]
        diag(S)<-diag(S)+eps2.vec
        SRinv<-chol2inv(chol(S))
        # from S+R go back to S
        diag(S)<-diag(S)-eps2.vec
      } else {
        # Update inverse matrix (Uboldi et al 2008, Appendix AND erratum!)
        aux<-SRinv
        SRinv<-aux[-indx,-indx]-
               (tcrossprod(aux[indx,-indx],aux[-indx,indx]))*Zinv[indx]
        S<-S[-indx,-indx]
        eps2.vec<-eps2.vec[-indx]
        rm(aux)
      }
      # next tree lines: compute cvres=(obs - CrossValidation_prediction)
      Zinv<-1/diag(SRinv)
      SRinv.d<-crossprod(SRinv,d[sel])
      ares<-crossprod(S,SRinv.d)-d[sel] #   a-Obs
      cvres<--Zinv*SRinv.d              # CVa-Obs, Lussana et al 2010 Eq(A13)
      sig2o<-mean(d[sel]*(-ares))       # Lussana et al 2010, Eq(32)
      if (sig2o<0.01) sig2o<-0.01       # better not use sig2o too small
      # pog=cvres/(sig2obs+sig2CVpred), Lussana et al 2010 Eq(20)
      pog[sel]<-(ares*cvres)/sig2o
      sctpog[ix[j[sel]]]<-pog[sel]
      assign("sctpog",sctpog,envir=.GlobalEnv)
      if (length(sel2check)==0) break
      # define the threshold for each observation
      if (any(!is.na(t2pos.vec))) {
        cvres<-(-cvres) # Obs-CVa 
        # cvres is a sel-vector, while T2vec is a sel2check-vector
        # sel includes all the sel2check values, plus possibly others
        T2vec<-vector(length=length(sel2check),mode="numeric")
        match<-match(sel2check,sel)
        ipos2check<-which(cvres[match]>=0 & is.na(dqctmp[sel2check]))
        if (length(ipos2check)>0) 
          T2vec[ipos2check]<- t2pos.vec[sel2check[ipos2check]]
        ipos2check<-which(cvres[match]<0 & is.na(dqctmp[sel2check]))
        if (length(ipos2check)>0)
          T2vec[ipos2check]<- t2neg.vec[sel2check[ipos2check]]
        rm(ipos2check)
      } else {
        T2vec<-t2.vec[sel2check]
      }
      # check if any obs fails the test
      if (any(pog[sel2check]>T2vec)) {
        # flag as suspect only the observation with the largest cvres 
        # allow for more than obs being flagged at the same time
        if (faster) {
          if (any(pog[sel2check]>(2*T2vec))) {
            indx<-which(pog[sel2check]>(2*T2vec))
          } else {
            ixprobsus<-which(pog[sel2check]>T2vec)
            indx<-ixprobsus[which.max(pog[ixprobsus])]
          }
        } else {
          ixprobsus<-which(pog[sel2check]>T2vec)
          indx<-ixprobsus[which.max(pog[ixprobsus])]
        }
        dqctmp[sel2check[indx]]<-sus.code
        # update global variable with flags
        dqcflag[ix[j[sel2check[indx]]]]<-sus.code
        assign("dqcflag",dqcflag,envir=.GlobalEnv)
        # update corep (useful from 2nd iteration onwards, if an obs fails the
        # sct, then no need for representativeness
        corep[ix[j[sel2check[indx]]]]<-NA
        assign("corep",corep,envir=.GlobalEnv)
        # update selection
        sel<-which(is.na(dqctmp) | dqctmp==argv$keep.code)
        sel2check<-which(is.na(dqctmp) & doittmp==1)
      } else {
        break
      }
    } # end cycle SCT model
    # coefficient of observation representativeness (more than 1 obs left)
    if (length(sel)>1) { 
      corep[ix[j[sel]]]<-(d[sel]*(-ares))/sig2o
      assign("corep",corep,envir=.GlobalEnv)
    }
    return(length(which(dqctmp==sus.code)))
  }
}

`nc4.getTime` <-
function(file.name,varid="time") {
    # open netCDF file
    defs<-nc_open(file.name)
    tcors <- ncvar_get(defs, varid=varid)
    tunits <- ncatt_get(defs, varid,"units")$value
    tt <- nc4t2str(nct=tcors,nct.unit=tunits,format="%Y%m%d%H%M")
    # close file
    nc_close(defs)
    # return the result
    tt    
}

`nc4.getDim` <-
function(file.name,varid="ensemble_member") {
    # open netCDF file
    defs<-nc_open(file.name)
    i<-0
    for (j in 1:defs$ndims) if (defs$dim[[j]]$name==varid) i<-j
    if (i==0) {
      for (j in 1:defs$nvars) if (defs$var[[j]]$name==varid) i<-j
      if (i==0) { nc_close(defs); return(NULL) }
      vals<-ncvar_get(defs, varid=varid)
    } else {
      vals<-defs$dim[[i]]$vals
    }
    nc_close(defs)
    # close file
    # return the result
    vals    
}


`str2Rdate` <-
function(ts,format="%Y-%m-%d %H:%M:%S") {
# Author: Christoph Frei
# ===========================================================================
# converts a string into an R (POSIXt,POSIXct) date object
# date objects can be used with arithmetic operations +/-
# ts is a character or a vector of characters with date information
# in the format specified in format
# Output is a date object

     # the lengthy bunch of testing is necessary because strptime needs
     # explicit specifications for month and day, otherwise it returns NA.
     # this extension allows inputs of the format "%Y-%m" in which case the
     # the first day of the month is taken as a reference.

     #Êcheck if year/month/day is specified
     ysp <- length(c(grep("%Y",format,fixed=TRUE),
                     grep("%y",format,fixed=TRUE)))
     msp <- length(c(grep("%m",format,fixed=TRUE),
                     grep("%b",format,fixed=TRUE),
                     grep("%B",format,fixed=TRUE)))
     jsp <- length(c(grep("%j",format,fixed=TRUE)))
     dsp <- length(c(grep("%d",format,fixed=TRUE)))
     if (ysp > 1) { stop("ERROR: Multiple specification of year in 
                         date format.") }
     if (ysp == 0) { stop("ERROR: No year specification in 
                         date format.") }
     if (msp > 1) { stop("ERROR: Multiple specification of month in 
                         date format.") }
     if (dsp > 1) { stop("ERROR: Multiple specification of day in 
                         date format.") }

     # append month or day if not specified
     tss <- ts
     formati <- format
     if (jsp == 0) {
     if (msp == 0) { 
        tss <- paste(tss,"01",sep="")
        formati <- paste(formati,"%m",sep="")
     }
     if (dsp == 0) { 
        tss <- paste(tss,"01",sep="") 
        formati <- paste(formati,"%d",sep="")
     }
     }

     # this is necessary because strptime() returns NA otherwise
     as.POSIXct(strptime(tss,format=formati),tz="GMT")
}

`Rdate2str` <-
function(date,format="%Y-%m-%d %H:%M:%S") {
#Author: Christoph Frei
# ===========================================================================
# converts a R (POSIXt,POSIXct) date object into a string
# date is one or a vector of R date-time objects
# format specifies the desired character format analogous to str2Rdate
# Output is one or a vector of characters
     format(date,format=format,tz="GMT")
}


`nc4t2str` <- 
function (nct, nct.unit, format = "%Y%m%d%H%M") {
#Author: Christoph Frei
    fmt.nc <- "%Y-%m-%d %H:%M:00"
    rd.expl <- str2Rdate("190001010000", format = "%Y%m%d%H%M")
    ss.expl <- Rdate2str(rd.expl, format = fmt.nc)
    ll <- nchar(ss.expl)
    reference <- nct.unit
    t.unit <- switch(substr(reference, 1, 4), 
                     minu = "T",
                     hour = "H", 
                     days = "D",
                     mont = "M",
                     year = "Y",
                     seco = "S")
    c.start <- switch(t.unit, T = 15, H = 13, D = 12, M = 14, Y = 13, S=15)
    ref.date.str <- substr(reference, c.start, stop = c.start + ll - 1)
    if (nchar(ref.date.str)<nchar(ss.expl)) {
      ll1 <- nchar(ref.date.str)
      aux<- substr(ss.expl, (ll1+1), stop = ll)
      ref.date.str<-paste(ref.date.str,":00",sep="")
    }
#           ref.date <- str2Rdate(ref.date.str, format = fmt.nc)
    ref.date <- str2Rdate(ref.date.str)
    if (!(t.unit %in% c("Y", "M", "T", "H", "D", "S"))) {
        stop("Time conversion not available")
    }
    if (t.unit %in% c("T", "H", "D", "S")) {
        t.scal <- switch(t.unit, T = 60, H = 3600, D = 86400, S=1)
        secs.elapsed <- nct * t.scal
        r.tims <- ref.date + secs.elapsed
    }
    if (t.unit == "Y") {
        yy.ref <- as.numeric(Rdate2str(ref.date, format = "%Y"))
        rr.ref <- Rdate2str(ref.date, format = "%m%d%H%M")
        yy.ch <- sprintf("%03d", nct + yy.ref)
        dd.ch <- sapply(yy.ch, FUN = function(x) paste(x, rr.ref, 
            sep = ""))
        r.tims <- str2Rdate(dd.ch, format = "%Y%m%d%H%M")
    }
    if (t.unit == "M") {
        yy.ref <- as.numeric(Rdate2str(ref.date, format = "%Y"))
        mm.ref <- as.numeric(Rdate2str(ref.date, format = "%m"))
        rr.ref <- Rdate2str(ref.date, format = "%d%H%M")
        mm <- nct + (yy.ref * 12 + mm.ref)
        hhy <- trunc((mm - 1)/12)
        hhm <- mm - hhy * 12
        dd.ch <- sapply(1:length(hhy), 
                 FUN = function(ii) paste(sprintf("%04d", 
                     hhy[ii]), sprintf("%02d", hhm[ii]), rr.ref, sep = ""))
        r.tims <- str2Rdate(dd.ch, format = "%Y%m%d%H%M")
    }
    Rdate2str(r.tims, format = format)
}


# read netCDF file format
nc4in<-function(nc.file,
                nc.varname,
                topdown=F,
                out.dim=NULL,
                proj4=NULL,
                nc.proj4=NULL,
                selection=NULL,
                verbose=F) {
#-------------------------------------------------------------------------------
# Read data from ncdf file.
# Parameters:
# + nc.file. character. netCDF file name
# + nc.varname. character. variable name in the netCDF file
# + out.dim. list. variable dimension parameters 
#   out.dim$
#    ndim. integer. number of dimensions
#    tpos. integer. position of the dimension identifying "time"
#    epos. integer. position of the dimension identifying "ensemble member"
#    zpos. integer. position of the dimension identifying "vertical level"
#    names. character vector (1...ndim). variable-specific dimension names
#           names[1]=x must be the Easting coordinate
#           names[2]=y must be the Northing coordinate
#           ...
# + selection. list. select gridded field to read
#   selection$
#    t. character. timestamp to read from the netCDF file (YYYYMMDDHH00)
#    z. numeric. hight level to read from the netCDF file
#    e. numeric. ensemble member to read from the netCDF file
# + proj4 string can be obtained either: 
#     proj4. character. proj4 string set by the user (override nc.proj4)
#   OR
#     nc.proj4. list. read proj4 string from file
#     nc.proj4$
#      var. variable which contains proj4 attribute
#      att. attribute of var which contains proj4 string
#
# example:
# file<-"/lustre/storeB/project/metproduction/products/meps/meps_extracted_2_5km_20161230T06Z.nc"
# nc.varname<-"precipitation_amount_acc"
# ndim<-4
# tpos<-4
# epos<-3
# names<-vector(length=4,mode="character")
# names<-c("x","y","ensemble_member","time") # same names as the nc dimensions
#   # morevoer, these must be all the dimensions of the nc.varname in the nc file
#   # (apart from dimensions that have just 1-field, i.e. they are not really 
#   #  useful dimensions)
#===============================================================================
  # open/read/close netcdf file
  if (!file.exists(nc.file)) return(NULL)
  nc <- nc_open(nc.file)
  nc.varnames<-attr(nc$var,"names")
  nc.dimnames<-attr(nc$dim,"names")
  if (is.null(out.dim)) {
    ndim<-2
    names<-vector(length=ndim,mode="character")
    names<-c(nc.dimnames[1],nc.dimnames[2])
    out.dim<-list(ndim=ndim,names=names)
    rm(ndim,names)
  }
  out.dimids4nc<-match(out.dim$names,nc.dimnames)
  if (any(is.na(out.dimids4nc))) return(NULL)
  if (!(nc.varname%in%nc.varnames)) return(NULL)
  # proj4 string is needed either (1) set by the user or (2) read from nc
  if (is.null(proj4) & is.null(nc.proj4)) {
    proj4<-NA
  } else if (is.null(proj4)) {
    aux<-ncatt_get(nc,nc.proj4$var,nc.proj4$att)
    if (!aux$hasatt) return(NULL)
    proj4<-aux$value
  }
  idx<-which(nc.varnames==nc.varname)
  var<-nc$var[[idx]]
  var.size  <- var$varsize
  var.ndims <- var$ndims
  var.dimids<- var$dimids + 1
  out.dimids4var<-match(out.dimids4nc,var.dimids)
  if (any(is.na(out.dimids4var))) return(NULL)
  # define raster
  text.xdim<-paste("nc$dim$",out.dim$names[1],"$vals",sep="")
  text.ydim<-paste("nc$dim$",out.dim$names[2],"$vals",sep="")
  dx<-abs(eval(parse(text=paste(text.xdim,"[1]",sep="")))-
          eval(parse(text=paste(text.xdim,"[2]",sep=""))))
  ex.xmin<-min(eval(parse(text=text.xdim)))-dx/2
  ex.xmax<-max(eval(parse(text=text.xdim)))+dx/2
  dy<-abs(eval(parse(text=paste(text.ydim,"[1]",sep="")))-
          eval(parse(text=paste(text.ydim,"[2]",sep=""))))
  ex.ymin<-min(eval(parse(text=text.ydim)))-dy/2
  ex.ymax<-max(eval(parse(text=text.ydim)))+dy/2
  nx<-eval(parse(text=paste("nc$dim$",out.dim$names[1],"$len",sep="")))
  ny<-eval(parse(text=paste("nc$dim$",out.dim$names[2],"$len",sep="")))
  if (nx>1 & ny>1) {
    is.raster<-T
    r <-raster(ncol=nx, nrow=ny,
               xmn=ex.xmin, xmx=ex.xmax,
               ymn=ex.ymin, ymx=ex.ymax,
               crs=proj4)
    r[]<-NA
  } else {
    is.raster<-F
    s<-NULL
  }
  # time coordinate
  if (!is.null(out.dim$tpos)) {
    if (!is.na(out.dim$tpos) & !is.nan(out.dim$tpos)) {
      tcors  <- ncvar_get(nc, varid=out.dim$names[out.dim$tpos])
      tunits <- ncatt_get(nc, out.dim$names[out.dim$tpos],"units")$value
      tall <- nc4t2str(nct=tcors,nct.unit=tunits,format="%Y%m%d%H%M")
      ntall<-length(tall)
      tsel<-1:ntall
      if (!is.null(selection)) {
        if (!is.null(selection$t)) {
          if (any(selection$t %in% tall)) {
            tsel<-which(tall %in% selection$t)
          }
        }
      }
      ntsel<-length(tsel)
    }
  } else {
    ntsel<-1
  }
  # ensemble member
  esel<-NULL
  nesel<-0
  if (!is.null(out.dim$epos)) {
    ecors<-ncvar_get(nc, varid=out.dim$names[out.dim$epos])
    eall<-ecors
    neall<-length(ecors)
    esel<-1:neall
    if (!is.null(selection)) {
      if (!is.null(selection$e)) {
        if (any(selection$e %in% eall)) {
          esel<-which(eall %in% selection$e)
          nesel<-length(esel)
        }
      }
    }
  }
  # elevation
  zsel<-NULL
  nzsel<-0
  if (!is.null(out.dim$zpos)) {
    zcors<-ncvar_get(nc, varid=out.dim$names[out.dim$zpos])
    zall<-zcors
    nzall<-length(zcors)
    zsel<-1:nzall
    if (!is.null(selection)) {
      if (!is.null(selection$z)) {
        if (any(selection$z %in% zall)) {
          zsel<-which(zall %in% selection$z)
          nzsel<-length(zsel)
        }
      }
    }
  }
  # labels
  jtot<-ntsel
  if (nzsel>0) jtot<-jtot*nzsel
  if (nesel>0) jtot<-jtot*nesel
  j<-0
  tlab<-vector(mode="character",length=jtot)
  if (nesel>0) {
    elab<-vector(mode="character",length=jtot)
  } else {
    elab<-NULL
  }
  zlab<-vector(mode="character",length=jtot)
  if (nesel>0) {
    zlab<-vector(mode="character",length=jtot)
  } else {
    zlab<-NULL
  }
  labels<-list(t=tlab,e=elab,z=zlab)
  data.vec<-array(data=NA,dim=c(nx*ny,jtot))
  if (out.dim$ndim==2) {
    j<-1
    start <- rep(1,var.ndims) # begin with start=(1,1,1,...,1)
    count <- var.size # begin w/count=(nx,ny,nz,...,nt), reads entire var
    data <- ncvar_get( nc, var, start=start, count=count )
    if (length(dim(data))!=2) return(NULL)
    if (topdown) for (i in 1:nx) data[i,1:ny]<-data[i,ny:1]
    if (is.raster) {
      r[]<-t(data)
      data.vec[,j]<-getValues(r)
      if (exists("s")) s<-stack(s,r)
      if (j==1) {
        s<-r
      } else {
        s<-stack(s,r)
      }
    }
    labels$t[j]<-NA
  } else {
    for (t in tsel) {
      # == cycle for data with 5 dimension (x,y,z,e,t)
      if ((nzsel>0) & (nesel>0)) {
        for (e in esel) {
          for (z in zsel) {
            j<-j+1
            start <- rep(1,var.ndims) # begin with start=(1,1,1,...,1)
            start[out.dimids4var[out.dim$tpos]] <- t # change to start=(1,1,1,...,t) 
                                                     # to read timestep t
            start[out.dimids4var[out.dim$epos]] <- e # change to start=(1,1,1,...,t) 
                                                     # to read timestep t, ens e
            start[out.dimids4var[out.dim$zpos]] <- z # change to start=(1,1,1,...,t) 
                                                     # to read timestep t, ens e, lev z
            count <- var.size # begin w/count=(nx,ny,nz,...,nt), reads entire var
            count[out.dimids4var[out.dim$tpos]] <- 1 # change to count=(nx,ny,nz,...,1)
                                                     # to read 1 tstep
            count[out.dimids4var[out.dim$epos]] <- 1 # change to count=(nx,ny,nz,...,1)
                                                     # to read 1 tstep, 1ens
            count[out.dimids4var[out.dim$zpos]] <- 1 # change to count=(nx,ny,nz,...,1)
                                                     # to read 1 tstep, 1ens, 1z
            data <- ncvar_get( nc, var, start=start, count=count )
            if (length(dim(data))!=2) return(NULL)
            if (topdown) for (i in 1:nx) data[i,1:ny]<-data[i,ny:1]
            if (is.raster) {
              r[]<-t(data)
              data.vec[,j]<-getValues(r)
              if (j==1) {
                s<-r
              } else {
                s<-stack(s,r)
              }
            }
            labels$t[j]<-tall[t]
            labels$e[j]<-eall[e]
            labels$z[j]<-zall[z]
            if (verbose) print(paste("Data for variable",nc.varname,
                                     "at timestep",labels$t[j],
                                     "for ensemble_member",labels$e[j],
                                     "at elevation",labels$z[j]))
          }   # end for z
        } # end for e
      # == cycle for data with 4 dimension (x,y,z,t)
      } else if (nzsel>0) {
        for (z in zsel) {
          j<-j+1
          start <- rep(1,var.ndims) # begin with start=(1,1,1,...,1)
          start[out.dimids4var[out.dim$tpos]] <- t # change to start=(1,1,1,...,t) 
                                                   # to read timestep t
          start[out.dimids4var[out.dim$zpos]] <- z # change to start=(1,1,1,...,t) 
                                                   # to read timestep t, ens e, lev z
          count <- var.size # begin w/count=(nx,ny,nz,...,nt), reads entire var
          count[out.dimids4var[out.dim$tpos]] <- 1 # change to count=(nx,ny,nz,...,1)
                                                   # to read 1 tstep
          count[out.dimids4var[out.dim$zpos]] <- 1 # change to count=(nx,ny,nz,...,1) 
                                                   # to read 1 tstep, 1ens, 1z
          data <- ncvar_get( nc, var, start=start, count=count )
          if (length(dim(data))!=2) return(NULL)
          if (topdown) for (i in 1:nx) data[i,1:ny]<-data[i,ny:1]
          if (is.raster) {
            r[]<-t(data)
            data.vec[,j]<-getValues(r)
            if (j==1) {
              s<-r
            } else {
              s<-stack(s,r)
            }
          }
          labels$t[j]<-tall[t]
          labels$z[j]<-zall[z]
          if (verbose) print(paste("Data for variable",nc.varname,
                                   "at timestep",labels$t[j],
                                   "at elevation",labels$z[j]))
        } # end for z
      # == cycle for data with 4 dimension (x,y,e,t)
      } else if (nesel>0) {
        for (e in esel) {
          j<-j+1
          start <- rep(1,var.ndims) # begin with start=(1,1,1,...,1)
          start[out.dimids4var[out.dim$tpos]] <- t # change to start=(1,1,1,...,t) 
                                                   # to read timestep t
          start[out.dimids4var[out.dim$epos]] <- e # change to start=(1,1,1,...,t) 
                                                 # to read timestep t, ens e, lev z
          count <- var.size # begin w/count=(nx,ny,nz,...,nt), reads entire var
          count[out.dimids4var[out.dim$tpos]] <- 1 # change to count=(nx,ny,nz,...,1)
                                                   # to read 1 tstep
          count[out.dimids4var[out.dim$epos]] <- 1 # change to count=(nx,ny,nz,...,1) 
                                                   # to read 1 tstep, 1ens, 1z
          data <- ncvar_get( nc, var, start=start, count=count )
          if (is.raster) {
            if (length(dim(data))!=2) return(NULL)
            if (topdown) for (i in 1:nx) data[i,1:ny]<-data[i,ny:1]
            r[]<-t(data)
            data.vec[,j]<-getValues(r)
            if (j==1) {
              s<-r
            } else {
              s<-stack(s,r)
            }
          } else {
            if (topdown) data[1:length(data)]<-data[length(data):1]
            data.vec[,j]<-data
          }
          labels$t[j]<-tall[t]
          labels$e[j]<-eall[e]
          if (verbose) print(paste("Data for variable",nc.varname,
                                   "at timestep",labels$t[j],
                                   "for ensemble_member",labels$e[j]))
        } # end for e
      # == cycle for data with 3 dimension (x,y,t)
      } else if (out.dim$ndim>2) {
        j<-j+1
        start <- rep(1,var.ndims) # begin with start=(1,1,1,...,1)
        start[out.dimids4var[out.dim$tpos]] <- t # change to start=(1,1,1,...,t) 
                                                 # to read timestep t
        count <- var.size # begin w/count=(nx,ny,nz,...,nt), reads entire var
        count[out.dimids4var[out.dim$tpos]] <- 1 # change to count=(nx,ny,nz,...,1) 
                                                 # to read 1 tstep
        data <- ncvar_get( nc, var, start=start, count=count )
        if (length(dim(data))!=2) return(NULL)
        if (topdown) for (i in 1:nx) data[i,1:ny]<-data[i,ny:1]
        if (is.raster) {
          r[]<-t(data)
          data.vec[,j]<-getValues(r)
          if (exists("s")) s<-stack(s,r)
          if (j==1) {
            s<-r
          } else {
            s<-stack(s,r)
          }
        }
        labels$t[j]<-tall[t]
        if (verbose) print(paste("Data for variable",nc.varname,
                                 "at timestep",labels$t[j]))
      } else {
        print("stranger things are going on right now!")
      } 
    } # end for "time"
  }
  nc_close(nc)
  return(list(data=data.vec,stack=s,labels=labels))
}

#+ plot (used for debugging)
plotSCTgrid<-function() {
  r<-raster(e,ncol=argv$grid.sct[2],nrow=argv$grid.sct[1])
  xy<-xyFromCell(r,1:ncell(r))
  xr<-xy[,1]
  yr<-xy[,2]
  ir<-1:ncell(r)
  r[]<-1:ncell(r)
  par(mar=c(1,1,1,1))
  plot(e,add=T)
  #points(xr,yr,pch=19,col="black")
  xrv<-unique(xr)
  yrv<-unique(yr)
  yrmn<-min(yr,na.rm=T)-res(r)[2]/2
  yrmx<-max(yr,na.rm=T)+res(r)[2]/2
  xrmn<-min(xr,na.rm=T)-res(r)[1]/2
  xrmx<-max(xr,na.rm=T)+res(r)[1]/2
  for (iii in 1:(length(yrv))) {
    lines(c(xrmn,xrmx),c(yrv[iii]+res(r)[2]/2,yrv[iii]+res(r)[2]/2),lwd=2.5)
  }
  for (iii in 1:(length(xrv))) {
    lines(c(xrv[iii]+res(r)[1]/2,xrv[iii]+res(r)[1]/2),c(yrmn,yrmx),lwd=2.5)
  }
  box()
}

#+ plot (used for debugging)
plot_debug<-function(ff,
                     r,
                     r1=NULL,
                     lbr=20,
                     x,
                     y,
                     proj4,
                     proj4plot=NULL) {
  rmn<-range(getValues(r),na.rm=T)[1]
  rmx<-range(getValues(r),na.rm=T)[2]
  rbr<-seq(rmn,rmx,length=lbr)
  col<-c(rev(rainbow((lbr-1))))
  png(file=ff,width=800,height=800)
  image(r,breaks=rbr,col=col)
  if (!is.null(r1)) contour(r1,levels=c(0,1),add=T)
  xy<-as.data.frame(cbind(x,y))
  coordinates(xy)<-c("x","y")
  proj4string(xy)<-CRS(proj4)
  if (!is.null(proj4plot)) xy<-spTransform(xy,CRS(proj4plot))
  points(xy,cex=0.8,pch=19)
  dev.off()
}

#------------------------------------------------------------------------------
# optimal interpolation 
oi_var_gridpoint_by_gridpoint<-function(i,
                                        dh=10000, #m
                                        box_o_nearest_halfwidth=100000, #m
                                        dz=NA,
                                        lafmin=NA,
                                        dh_adaptive=F,
                                        corr="soar",
                                        pmax,
                                        fg=NA,
                                        fg_gamma=NA,
                                        fg_min=NA,
                                        fg_max=NA,
                                        succ_corr=F,
                                        y_elab=F,
                                        loocv=F,
                                        o_errvar_min=0.001,
                                        o_errvar_max=4,
                                        xa_errvar_min=0.001,
                                        xa_errvar_max=4) {
# global variables: xgrid_spint, ygrid_spint, zgrid_spint, lafgrid_spint,
#                   xobs_spint, yobs_spint, zobs_spint, lafobs_spint,
#                   yo_spint, yb_spint, xb_spint, eps2_spint
#                   nobs
# Description:
# OI analysis at the i-th point (at xgrid_spint[i],ygrid_spint[i],...)
# given the set of observations yo_spint (at xobs_spint,yobs_spint,...)
# with or without a backgorund
# 
# Input arguments
# i: gridpoint index (refers to vectors Xgrid_spint,...,xb_spint
# dh: horizontal de-correlation length (m)
# box_o_nearest_halfwidth: half-width of the square box used to select the 
#                          nearest observations
# dz: vertical de-correlation length (m, NA if z is not considered)
# lafmin: land-area fraction minimum value for the de-correlation factor 
#         (0-1, NA if laf is not considered)
# dh_adaptive: logical. F if dh is the actual horizontal de-correlation length
#                       T if the actual horizontal de-correlation length is 
#                         obtained as the 10-th percentile of the distances
#                         between the i-th gridpoint and the nearest observations
#                         (up to the nearest pmax observations). In this case dh
#                         is used as a lower limit.
# corr: model for the correlation functions. 
#       "soar" second order auto-regressive (only horizontal distance is used).
#       "gaussian" gaussian.
# pmax: mximum number of nearest observations to consider
# fg: method used to compute the background (first-guess). 
#  NA (default), background available in yb_spint, xb_spint
#  "linear", background as a linear function of the vertical coordinate  
#  "Frei", background as a non-linear function of the vertical coordinate  
#  "mean", background set to the mean value of the nearest observations
# fg_gamma: vertical lapse rate used for fg="linear" and optimized for fg="Frei"
# fg_min: minimum allowed value for the background
# fg_max: maximum allowed value for the background
# succ_corr: logical.
#  F (default). the background is used as it is.
#  T. three steps of successive corrections are applied to the background.
# y_elab: logical.
#  F (default). grid.. vectors and obs.. vectors refer to the same locations
#  T. grid.. and obs.. vectors refer to different locations
# loocv. logical. Leave one out cross-validation
#  F (default). standard run, use all the observations
#  T. run in loocv mode, without considering the observations at the i-th gridpoint
# o_errvar_min. minimum allowed value for the observation error variance.
# o_errvar_max. maximum allowed value for the observation error variance.
# xa_errvar_min. minimum allowed value for the analysis error variance.
# xa_errvar_max. maximum allowed value for the analysis error variance.
#
# return(c(xa,xa_errvar,o_errvar,xidi,idiv,av,dh))
# NOTE: av is the leave-one-out CV. However, its errvar is not returned.
#       todo: figure out how to compute the leave-ione-out errvar.
#------------------------------------------------------------------------------
# select the p_i observations nearest to the i-th gridpoint
  deltax<-abs(xgrid_spint[i]-xobs_spint)
  deltay<-abs(ygrid_spint[i]-yobs_spint)
  if (!is.na(dz)) {deltaz<-abs(zgrid_spint[i]-zobs_spint);dz2<-dz*dz}
  if (!is.na(lafmin)) deltalaf<-abs(lafgrid_spint[i]-lafobs_spint)
  consider_obs<-rep(T,nobs)
  if (y_elab) {
    res<-ifelse(exists("yb_spint"),yb_spint[i],NA)
    av<-NA
    idiv<-0
    xidi<-1/(1+eps2_spint[i])
    if (loocv) consider_obs[i]<-F
    if (length(which(deltax<(box_o_nearest_halfwidth)))==1) 
      return(c(res,NA,NA,xidi,idiv,av,dh))
    if (length(which(deltay<(box_o_nearest_halfwidth)))==1) 
      return(c(res,NA,NA,xidi,idiv,av,dh))
  } else {
    res<-ifelse(exists("xb_spint"),xb_spint[i],NA)
    av<-NA
    idiv<-NA
    xidi<-0
    if (loocv) consider_obs[which(deltax<1 & deltay<1)]<-F
    if (!any(deltax<(box_o_nearest_halfwidth))) 
      return(c(res,NA,NA,xidi,idiv,av,dh))
    if (!any(deltay<(box_o_nearest_halfwidth))) 
      return(c(res,NA,NA,xidi,idiv,av,dh))
  }
  ixa<-which( deltax<(box_o_nearest_halfwidth) & 
              deltay<(box_o_nearest_halfwidth) & 
              consider_obs )
  if (length(ixa)==0) return(c(res,NA,NA,xidi,idiv,av,dh))
  disth2<-deltax[ixa]*deltax[ixa]+deltay[ixa]*deltay[ixa]
  if (length(ixa)>pmax) {
    ixb<-order(disth2, decreasing=F)[1:pmax]
    ixa<-ixa[ixb]
    disth2<-disth2[ixb]
    rm(ixb)
  }
  p_i<-length(ixa)
  # correlation matrices
  if (dh_adaptive) {
    dh<-max(dh,as.numeric(quantile(sqrt(disth2),probs=0.1)))
  }
  dh2<-dh*dh
  if (corr=="gaussian") {
    rloc<-exp( -0.5* disth2 / dh2 )
  } else if (corr=="soar")  {
    distnorm_loc<-sqrt(disth2) / dh
    rloc<-(1+distnorm_loc)*exp(-distnorm_loc)
    if (!succ_corr) rm(distnorm_loc)
  }
  if (corr=="gaussian") {
    S<-exp(-0.5*(outer(yobs_spint[ixa],yobs_spint[ixa],FUN="-")**2. + 
                 outer(xobs_spint[ixa],xobs_spint[ixa],FUN="-")**2)/dh2)
  } else if (corr=="soar")  {
    distnorm<-sqrt(outer(yobs_spint[ixa],yobs_spint[ixa],FUN="-")**2. + 
                   outer(xobs_spint[ixa],xobs_spint[ixa],FUN="-")**2) / dh 
    S<-(1+distnorm)*exp(-distnorm)
    if (!succ_corr) rm(distnorm)
  }
  # first-guess from nearest observations
  if (!is.na(fg)) {
    yo_mean<-mean(yo_spint[ixa])
    # Frei profile may provide unrealistic values if 
    #  it is required to extrapole a value for elevations 
    #  far from the ones used to optimize the parameters
    #  OR if the elevations used are within a narrow layer
    if (fg=="Frei") {
      zmin<-sort(zobs_spint[ixa])[min(c(length(ixa),2))]
      zmax<-sort(zobs_spint[ixa])[max(c(1,(p_i-1)))]
      if ( (zmin-zgrid_spint[i])>50 |
           (zgrid_spint[i]-zmax)>50 |
           (zmax-zmin)<25 ) fg<-"linear"
    } 
    if (fg=="linear") {
      par<-c(yo_mean)
      opt<-optimize(f=vertprofbasic2opt,
                    interval=c(argv$vmin,argv$vmax),
                    vert_coord=zobs_spint[ixa],
                    gamma=fg_gamma,
                    obs=yo_spint[ixa])
#      aaa<-tvertprof_basic(1:2000,
#                                  t0=opt$minimum,
#                                  gamma=fg_gamma)
      yb_spint_i<-tvertprof_basic(zobs_spint[ixa],
                                  t0=opt$minimum,
                                  gamma=fg_gamma)
      xb_i<-tvertprof_basic(zgrid_spint[i],
                            t0=opt$minimum,
                            gamma=fg_gamma)
    } else if (fg=="Frei") {
      par<-c(yo_mean,
             fg_gamma,
             5,
             zmin,
             zmax)
      opt<-optim(par,vertprof2opt,vert_coord=zobs_spint[ixa],obs=yo_spint[ixa])
      yb_spint_i<-tvertprof(z=zobs_spint[ixa],
                            t0=opt$par[1],
                            gamma=opt$par[2],
                            a=opt$par[3],
                            h0=opt$par[4],
                            h1i=opt$par[5])
      xb_i<-tvertprof(z=zgrid_spint[i],
                      t0=opt$par[1],
                      gamma=opt$par[2],
                      a=opt$par[3],
                      h0=opt$par[4],
                      h1i=opt$par[5])
#      aaa<-tvertprof(z=1:2000,
#                      t0=opt$par[1],
#                      gamma=opt$par[2],
#                      a=opt$par[3],
#                      h0=opt$par[4],
#                      h1i=opt$par[5])
    } else if (fg=="mean") {
      yb_spint_i<-rep(yo_mean,length=p_i)
      xb_i<-yo_mean
    }
  } else {
    yb_spint_i<-yb_spint[ixa]
    xb_i<-xb_spint[i]
  } # end compute first-guess
  if (!is.null(fg_min)) {
    if (!is.na(fg_min) & !is.nan(fg_min)) {
      yb_spint_i[which(yb_spint_i<fg_min)]<-fg_min
      xb_i<-max(xb_i,fg_min)
    }
  }
  if (!is.null(fg_max)) {
    if (!is.na(fg_max) & !is.nan(fg_max)) {
      yb_spint_i[which(yb_spint_i>fg_max)]<-fg_max
      xb_i<-min(xb_i,fg_max)
    }
  }
  # successive corrections (Barnes scheme) step
  if (succ_corr) {
    for (sc in 3:1) {
      if (corr=="gaussian") {
        S1<-S**(1/sc**2)
        rloc1<-rloc**(1/sc**2)
      } else if (corr=="soar")  {
        distnorm1<-distnorm/sc
        S1<-(1+distnorm1)*exp(-distnorm1)
        distnorm1_loc<-distnorm_loc/sc
        rloc1<-(1+distnorm1_loc)*exp(-distnorm1_loc)
      }
      yb_spint_i<-yb_spint_i+crossprod(S1,(yo_spint[ixa]-yb_spint_i)) / 
                  (rowSums(S1)+eps2_spint[ixa])
      xb_i<-xb_i+sum(rloc1*(yo_spint[ixa]-yb_spint_i))/sum(rloc1+eps2_spint[ixa])
    }
    rm(S1,rloc1)
    if (corr=="soar") rm(distnorm1,distnorm,distnorm_loc,distnorm1_loc)
  }
  # adjust gaussian correlations by taking into account more geo-parameters
  if (corr=="gaussian") {
    if (!is.na(dz)) {
      S<-S*exp(-0.5*abs(outer(zobs_spint[ixa],zobs_spint[ixa],FUN="-"))/dz2) 
      rloc<-rloc*exp(-0.5*deltaz[ixa]/dz2)
    }
    if (!is.na(lafmin)) {
      S<-S*(1-(1-lafmin)*
         abs(outer(lafobs_spint[ixa],lafobs_spint[ixa],FUN="-")))
      rloc<-rloc*(1-(1-lafmin)*deltalaf[ixa])
    }
  }
  # innovation
  di<-yo_spint[ixa]-yb_spint_i
  # OI analysis
  SRinv<-chol2inv(chol( (S+diag(x=eps2_spint[ixa],length(ixa))) ))
  xidi<-sum(rloc*as.vector(rowSums(SRinv)))
  SRinv_di<-crossprod(SRinv,di)       
  o_errvar<-min(c(o_errvar_max,
                  max(c(o_errvar_min,
                        mean(di*(di-crossprod(S,SRinv_di)))))))
  rm(S)
  xa_errvar<-min(c(xa_errvar_max,
                   max(c(xa_errvar_min,
                         (o_errvar/ mean(eps2_spint[ixa])) * 
                         (1-sum(as.vector(crossprod(rloc,SRinv))*rloc))))))
  xa<-xb_i+sum(rloc*as.vector(SRinv_di))
  if (y_elab & !loocv) {
    ii<-which(ixa==i)
    Wii<-sum(rloc*SRinv[ii,])
    idiv<-(xidi-Wii)/(1-Wii)
    av<-(xa-Wii*yo_spint[i])/(1-Wii)
  }
  # debug
  #if (yo_to_check[i]>20) {
  #png(file=paste0("png/tvert_",i,".png"),width=800,height=800)
  #plot(aaa,1:2000,ylim=range(c(zobs_spint[ixa],zgrid_spint[i])),xlim=range(c(yo_spint[ixa],yb_spint_i,xb_i,yo_to_check[i])),
  #main=paste(round(xa,1),round(xa_errvar,3),round(o_errvar,3),round(xidi,3),round(dh,1),round((yo_to_check[i]-xa)**2/(xa_errvar+o_errvar),2)))
  #points(yo_spint[ixa],zobs_spint[ixa],pch=21,bg="blue")
  #points(yb_spint_i,zobs_spint[ixa],pch=21,bg="cyan")
  #points(xb_i,zgrid_spint[i],pch=21,bg="red")
  #points(yo_to_check[i],zgrid_spint[i],pch=21,bg="red")
  #dev.off()
  #}
  return(c(xa,xa_errvar,o_errvar,xidi,idiv,av,dh))
}

#+ spatial interpolation for cool test 
spint_cool<-function(i,
                     thres,  # same unit as yo_cool
                     dh_max, # max distance to be considered a nn
                     condition="lt",
                     mode="nearest") {
#------------------------------------------------------------------------------
# spatial interpolation for the cool test. Defualt is nearest neighbour.
#------------------------------------------------------------------------------
  if (mode=="nearest") { 
    deltax<-abs(xgrid_cool[i]-xobs_cool)
    deltay<-abs(ygrid_cool[i]-yobs_cool)
    if (!any(deltax<dh_max & deltay<dh_max)) return(NA)
    if (any(deltax<(dh_max/20) & deltay<(dh_max/20))) {
      ixnear<-which(deltax<(dh_max/20) & deltay<(dh_max/20))
    } else if (any(deltax<(dh_max/10) & deltay<(dh_max/10))) {
      ixnear<-which(deltax<(dh_max/10) & deltay<(dh_max/10))
    } else {
      ixnear<-which(deltax<(dh_max) & deltay<(dh_max))
    }
    ixnearest<-which.min(deltax[ixnear]*deltax[ixnear]+
                         deltay[ixnear]*deltay[ixnear])
#    if (sqrt((xgrid_cool[i]-xobs_cool[ixnear])**2+
#             (ygrid_cool[i]-yobs_cool[ixnear])**2)>dh_max) 
#      return(NA)
    if (condition=="lt") {
      ifelse(yo_cool[ixnear[ixnearest]]< thres,return(0),return(1))
    } else if (condition=="le") {
      ifelse(yo_cool[ixnear[ixnearest]]<=thres,return(0),return(1))
    } else if (condition=="gt") {
      ifelse(yo_cool[ixnear[ixnearest]]> thres,return(0),return(1))
    } else if (condition=="ge") {
      ifelse(yo_cool[ixnear[ixnearest]]>=thres,return(0),return(1))
    }
  } else {
    return(NA)
  }
}

#==============================================================================
#  MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN
#==============================================================================
t0<-Sys.time()
#
# default values that need to be specified
proj4_input_obsfiles_default<-"+proj=longlat +datum=WGS84"
proj4_where_dqc_is_done_default<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
xy.dig.out_default<-5
varname.y.out_default<-"lat"
varname.x.out_default<-"lon"
#
# create parser object
p <- arg_parser("titan")
# specify our desired options
#.............................................................................. 
# REQUIRED input 
p <- add_argument(p, "input",
                  help="input file",
                  type="character")
p <- add_argument(p, "output",
                  help="output file",
                  type="character",
                  default="output.txt")
# configuration file
p <- add_argument(p, "--config.file",
                  help="configuration file",
                  type="character",
                  default=NULL,
                  short="cf")
#.............................................................................. 
# VARIABLE definition 
p<- add_argument(p, "--variable",
                 help=paste("meteorological variable (T temperature,",
                            " RR precipitation, RH relative humidity,",
                            " SD surface_snow_thickness)"),
                 type="character",
                 default="T",
                 short="-var")
# parameter for the Box-Cox transformation (rquired for var=RR)
p <- add_argument(p, "--boxcox.lambda",
                  help="parameter used in the Box-Cox transformation (var RR)",
                  type="numeric",default=0.5,short="-l")
#.............................................................................. 
# titan path (use for the sct with smart boxes)
p <- add_argument(p, "--titan_path",
                  help="path to the directory where the TITAN code is",
                  type="character",
                  default=NULL,
                  short="-tip")
#.............................................................................. 
neg.str<-"Negative values can be specified in either one of these two ways: (1) by using the corresponding \"...neg...\" command line argument; (2) negative values start with \"_\" (e.g. _2=-2)"
# ADDITIONAL input files / providers
p <- add_argument(p, "--input.files",
                  help="additional input files (provider2 provider3 ...)",
                  type="character",
                  default=NULL,
                  nargs=Inf,
                  short="-i")
p <- add_argument(p, "--prid",
                  help="provider identifiers (provider1 provider2 provider3 ...)",
                  type="character",
                  default=NA,
                  nargs=Inf,
                  short="-pr")
p <- add_argument(p, "--input.offset",
                  help=paste("offset applied to the input files (one for each provider, default=0).",neg.str),
                  type="character",
                  default=NA,
                  nargs=Inf,
                  short="-io")
p <- add_argument(p, "--input.negoffset",
                  help="sign for the offsets (1=negative; 0=positive, def=0)",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-ion")
#
p <- add_argument(p, "--input.cfact",
                  help=paste("correction factor applied to the input files (one for each provider, default=1)",neg.str),
                  type="character",
                  default=NA,
                  nargs=Inf, 
                  short="-icf")
p <- add_argument(p, "--input.negcfact",
                  help="sign for the correction factors (1=negative; 0=positive, def=0)",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-icfn")
#.............................................................................. 
# DEBUG
p <- add_argument(p, "--debug",
                  help="debug mode",
                  flag=T,
                  short="-dbg")
p <- add_argument(p, "--debug.dir",
                  help="directory for debug output",
                  type="character",
                  default=".",
                  short="-dbgd")
p <- add_argument(p, "--verbose",
                  help="debug mode",
                  flag=T,
                  short="-v")
#.............................................................................. 
# Quality control codes
p <- add_argument(p, "--nometa.code",
                  help="quality code returned in case of missing metadata",
                  type="numeric",
                  default=1,
                  short="-nometac")
p <- add_argument(p, "--p.code",
  help="quality code returned in case of the check on plausible values fails",
                  type="numeric",
                  default=2,
                  short="-pcodec")
p <- add_argument(p, "--clim.code",
 help="quality code returned in case of the check on climatological values fails",
                  type="numeric",
                  default=3,
                  short="-climc")
p <- add_argument(p, "--buddy.code",
  help="quality code returned in case of the buddy check fails",
                  type="numeric",
                  default=4,
                  short="-buddyc")
p <- add_argument(p, "--sct.code",
                  help="quality code returned in case of SCT fails",
                  type="numeric",
                  default=5,
                  short="-sctc")
p <- add_argument(p, "--dem.code",
                  help="quality code returned in case of SCT fails",
                  type="numeric",
                  default=6,
                  short="-demc")
p <- add_argument(p, "--isol.code",
        help="quality code returned in case of isolation check fails",
                  type="numeric",
                  default=7,
                  short="-isolc")
p <- add_argument(p, "--fg.code",
 help="quality code returned in case of check against a first-guess field (deterministic) fails",
                  type="numeric",
                  default=8,
                  short="-fgc")
p <- add_argument(p, "--fge.code",
 help="quality code returned in case of check against a first-guess field (ensemble) fails",
                  type="numeric",
                  default=10,
                  short="-fgc")
p <- add_argument(p, "--ccrrt.code",
                  help=paste("quality code returned in case of precipitation",
                             "and temperature crosscheck fails"),
                  type="numeric",
                  default=11,
                  short="-ccrrtc")
p <- add_argument(p, "--cool.code",
             help=paste("quality code returned in case of cool check fails"),
                  type="numeric",
                  default=15,
                  short="-coolc")
p <- add_argument(p, "--buddy_eve.code",
  help="quality code returned in case of the buddy check event-based fails",
                  type="numeric",
                  default=13,
                  short="-buddyec")
p <- add_argument(p, "--black.code",
    help="quality code assigned to observations listed in the blacklist",
                  type="numeric",
                  default=100,
                  short="-blackc")
p <- add_argument(p, "--keep.code",
    help="quality code assigned to observations listed in the keep-list",
                  type="numeric",
                  default=100,
                  short="-keepc")
#.............................................................................. 
# standard value for the moist adiabatic lapse rate
p <- add_argument(p, "--gamma.standard",
  help="standard value for the moist adiabatic temperature lapse rate dT/dz",
                  type="numeric",
                  default=-0.0065)
#.............................................................................. 
# GEOGRAPHICAL domain definition  
# NOTE: lat-lon setup to have Oslo in a single box
p <- add_argument(p, "--lonmin",
                  help="longitude of south-eastern domain corner",
                  type="numeric",
                  default=5,
                  short="-lon")
p <- add_argument(p, "--lonmax",
                  help="longitude of south-western domain corner",
                  type="numeric",
                  default=28,
                  short="-lox")
p <- add_argument(p, "--latmin",
                  help="latitude of south-eastern domain corner",
                  type="numeric",
                  default=53.25,
                  short="-lan")
p <- add_argument(p, "--latmax",
                  help="latitude of north-western domain corner",
                  type="numeric",
                  default=71.8,
                  short="-lax")
p <- add_argument(p, "--dqc_inbox_only",
                  help="perform dqc only in the defined box (lonmin,lonmax,latmin,latmax)",
                  type="logical",
                  default=F)

# transformation between coordinate reference systems
p <- add_argument(p, "--spatconv",
                  help="flag for conversion of spatial coordinates before running the data quality checks",
                  flag=T,
                  short="-c")
p <- add_argument(p, "--proj4from",
                  help="proj4 string for the original coordinate reference system (obsolete, use \"proj4_input_obsfiles\")",
                  type="character",
                  default=proj4_input_obsfiles_default,
                  short="-pf")
p <- add_argument(p, "--proj4_input_obsfiles",
                  help="proj4 string for the original coordinate reference system",
                  type="character",
                  default=proj4_input_obsfiles_default)
p <- add_argument(p, "--proj4to",
                  help="proj4 string for the coordinate reference system where the DQC is performed (obsolete, use \"proj4_where_dqc_is_done\"",
                  type="character",
                  default=proj4_where_dqc_is_done_default,
                  short="-pt")
p <- add_argument(p, "--proj4_where_dqc_is_done",
                  help="proj4 string for the coordinate reference system where the DQC is performed",
                  type="character",
                  default=proj4_where_dqc_is_done_default)
p <- add_argument(p, "--proj4_output_files",
                  help="proj4 string for the output coordinate reference system",
                  type="character",
                  default=proj4_input_obsfiles_default)
#.............................................................................. 
# INPUT/OUTPUT names
p <- add_argument(p, "--separator",
                  help="character vector, input file(s) separator character(s) (default '';'')",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--separator.out",
                  help="separator character in the output file",
                  type="character",
                  default=";")
p<- add_argument(p, "--latlon.dig.out",
                 help="number of decimal digits for latitude and longitude in the output file  (obsolete, use xy.dig.out)",
                 type="numeric",
                 default=5,
                 short="-lldo")
p<- add_argument(p, "--xy.dig.out",
                 help="number of decimal digits for northing and easting coordinates in the output file",
                 type="numeric",
                 default=xy.dig.out_default)
p<- add_argument(p, "--elev.dig.out",
                 help="number of decimal digits for elevation in the output file",
                 type="numeric",
                 default=1,
                 short="-edo")
p<- add_argument(p, "--value.dig.out",
                 help="number of decimal digits for the returned value in the output file",
                 type="numeric",
                 default=1,
                 short="-vdo")
p <- add_argument(p, "--varname.lat",
                  help="character vector, latitude variable name(s) in the input file (default ''lat'')",
                  type="character",
                  nargs=Inf,
                  short="-vlat")
p <- add_argument(p, "--varname.lat.out",
                  help="latitude variable name in the output file (obsolete, use varname.y.out)",
                  type="character",
                  default="lat")
p <- add_argument(p, "--varname.y.out",
                  help="northing coordinate name in the output file",
                  type="character",
                  default=varname.y.out_default)
p <- add_argument(p, "--varname.lon",
                  help="character vector, longitude variable name(s) in the input file (default ''lon'')",
                  type="character",
                  nargs=Inf,
                  short="-vlon")
p <- add_argument(p, "--varname.lon.out",
                  help="longitude variable name in the output file (obsolete, use varname.x.out)",
                  type="character",
                  default="lon")
p <- add_argument(p, "--varname.x.out",
                  help="easting coordinate name in the output file",
                  type="character",
                  default="lon")
p <- add_argument(p, "--varname.elev",
                  help="character vector, elevation variable names(s) in the input file (default ''elev'')",
                  type="character",
                  nargs=Inf,
                  short="-vele")
p <- add_argument(p, "--varname.elev.out",
                  help="elevation variable name in the output file",
                  type="character",
                  default="elev")
p <- add_argument(p, "--elev_not_used",
                  help="elevation is not used (will be set to zero)",
                  type="logical",
                  default=F)
p <- add_argument(p, "--varname.value",
                  help="character vector, variable name(s) in the input file (default ''value'')",
                  type="character",
                  nargs=Inf,
                  short="-vval")
p <- add_argument(p, "--varname.value.out",
                  help="name for the variable values (in/out)",
                  type="character",
                  default="value")
# output file
p <- add_argument(p, "--varname.opt",
     help="additional optional variables to be written on the output (out)",
                  type="character",
                  default=NA,
                  nargs=Inf,
                  short="-vopt")
p<- add_argument(p, "--varname.prid",
                 help="name for the provider identifier (out)",
                 type="character",
                 default="prid",
                 short="-vprid")
p<- add_argument(p, "--varname.dqc",
                 help="name for the data quality control flag (out)",
                 type="character",
                 default="dqc",
                 short="-vdqc")
p<- add_argument(p, "--varname.sct",
            help="name for the spatial consistency test returned value (out)",
                 type="character",
                 default="sct",
                 short="-vsct")
p<- add_argument(p, "--varname.rep",
            help="name for the coefficient of representativeness (out)",
                 type="character",
                 default="rep",
                 short="-vrep")
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
#.............................................................................. 
# Plausibility check
p <- add_argument(p, "--tmin",
                  help="minimum allowed temperature [K or degC] (deprecated, use vmin instead)",
                  type="numeric",
                  default=NA,
                  short="-tP")
p <- add_argument(p, "--tmax",
                  help="maximum allowed temperature [K or degC] (deprecated, use vmax instead)",
                  type="numeric",
                  default=NA,
                  short="-TP")
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
#.............................................................................. 
# Climatological check
# default based on Norwegian hourly temperature from 2010-2017
p <- add_argument(p, "--tmin.clim",
                  help="minimum allowed temperature [K or degC] (deprecated)",
                  type="numeric",
                  nargs=12,
                  short="-tC",
                  default=rep(NA,12))
p <- add_argument(p, "--tmax.clim",
                  help="maximum allowed temperature [K or degC] (deprecated)",
                  type="numeric",
                  nargs=12,
                  short="-TC",
                  default=rep(NA,12))
p <- add_argument(p, "--vmin.clim",
                  help=paste("minimum allowed value [units of the variable specified]",neg.str),
                  type="character",
                  nargs=12,
                  default=c("_45","_45","_40","_35","_20","_15","_10","_15","_15","_20","_35","_45"))
p <- add_argument(p, "--vmax.clim",
                  help=paste("maximum allowed value [units of the variable specified]",neg.str),
                  type="character",
                  nargs=12,
                  default=c("20","20","25","25","35","35","40","40","35","30","25","20"))
p <- add_argument(p, "--month.clim",
                  help="month (number 1-12)",
                  type="numeric",
                  short="-mC",
                  default=NA)
p <- add_argument(p, "--vminsign.clim",
                  help="minimum allowed value, sign [1=neg, 0=pos]",
                  type="numeric",
                  nargs=12,
                  default=c(0,0,0,0,0,0,0,0,0,0,0,0))
p <- add_argument(p, "--vmaxsign.clim",
                  help="maximum allowed value, sign [1=neg, 0=pos]",
                  type="numeric",
                  nargs=12,
                  default=c(0,0,0,0,0,0,0,0,0,0,0,0))
#.............................................................................. 
# Buddy-check (event-based)
p <- add_argument(p, "--buddy_eve",
                  help="do the buddy check event-based",
                  flag=T,
                  short="-Be")
p <- add_argument(p, "--thr_eve.buddy_eve",
                  help=paste("numeric vector with the thresholds used to define events (same units of the specified variable). Each threshold defines two events: (i) less than the threshold and (ii) greater or equal to it."),
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--dr.buddy_eve",
                  help="numeric vector of the same dimension of thr_eve.buddy_eve. perform the buddy_eve-check in a dr-by-dr square-box around each observation [m] (default is 3000m)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--i.buddy_eve",
                  help="number of buddy_eve-check iterations",
                  type="integer",
                  default=1,
                  short="-iBe")
p <- add_argument(p, "--thr.buddy_eve",
                  help="buddy_eve-check threshold (same dimension of thr_eve.buddy_eve). flag observation as suspect if: is event(no) and less than a fraction of thr.buddy_eve of the neighbouring observations are event(no) OR is event(yes) and less than a fraction of thr.buddy_eve of the neighbouring observations are event(yes). (default is 0.05, i.e. 5% of the neighbouring observations)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--n.buddy_eve",
                  help="numeric vector of the same dimension of thr_eve.buddy_eve. minimum number of neighbouring observations to perform the buddy-check (dafualt is set to 5)",
                  type="integer",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--dz.buddy_eve",
                  help="numeric vector of the same dimension of thr_eve.buddy_eve. maximum allowed range of elevation in a square-box to perform the buddy_eve-check (i.e. no check if elevation > dz.buddy_eve) (default is 30m)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--prio.buddy_eve",
                  help="priorities used in the buddy_eve check. One positive value for each provider. The lower the number the more priority that provider gets. In the first round of the buddy_eve check, high priority observations are not compared against lower priority observations. From the second round onwards, all the providers are set to the same priority value.",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-prBe")
#.............................................................................. 
# Buddy-check
p <- add_argument(p, "--dr.buddy",
                  help="consider only the observations within a radius of \"dr\" meters from each observation [m]",
                  type="numeric",
                  default=3000,
                  short="-dB")
p <- add_argument(p, "--i.buddy",
                  help="number of buddy-check iterations",
                  type="integer",
                  default=1,
                  short="-iB")
p <- add_argument(p, "--thr.buddy",
                  help="buddy-check threshold. flag observation if: abs(obs-pred)/st_dev > thr.buddy",
                  type="numeric",
                  default=3,
                  short="-thB")
p <- add_argument(p, "--sdmin.buddy",
                  help="minimum allowed value for the standard deviation",
                  type="numeric",
                  default=0.5)
p <- add_argument(p, "--n.buddy",
                  help="minimum number of neighbouring observations to perform the buddy-check",
                  type="integer",
                  default=5,
                  short="-nB")
p <- add_argument(p, "--dz.buddy",
                  help="maximum allowed range of elevation in a square-box to perform the buddy-check (i.e. no check if elevation > dz.buddy)",
                  type="numeric",
                  default=30,
                  short="-zB")
p <- add_argument(p, "--prio.buddy",
                  help="priorities used in the buddy check. One positive value for each provider. Lower numbers get higher priorities. In the first round of the buddy check, high priority observations are not compared against lower priority observations. From the second round onwards, all the providers are set to the same priority value.",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-prB")
#.............................................................................. 
# isolated stations
p <- add_argument(p, "--dr.isol",
                  help="check for the number of observation in a dr-by-dr square-box around each observation [m]",
                  type="numeric",
                  default=25000,
                  short="-dI")
p <- add_argument(p, "--n.isol",
                  help="threshold (number of neighbouring observations) for the identification of isolated observations.",
                  type="integer",
                  default=10,
                  short="-nI")
#.............................................................................. 
# spatial consistency test
p <- add_argument(p, "--grid.sct",
                  help="nrow ncol (i.e. number_of_rows number_of_columns). SCT is performed independently over several boxes. The regular nrow-by-ncol grid is used to define those rectangular boxes where the SCT is performed.",
                  type="integer",
                  nargs=2,
                  default=c(20,20),
                  short="-gS")
p <- add_argument(p, "--i.sct",
                  help="number of SCT iterations",
                  type="integer",
                  default=1,
                  short="-iS")
p <- add_argument(p, "--n.sct",
                  help="minimum number of stations in a box to run SCT",
                  type="integer",default=50,short="-nS")
p <- add_argument(p, "--dz.sct",
                  help="minimum range of elevation in a box to run SCT [m]",
                  type="numeric",
                  default=30,
                  short="-zS")
p <- add_argument(p, "--DhorMin.sct",
                  help=paste("OI, minimum allowed value for the horizontal de-correlation",
                  "lenght (of the background error correlation) [m]"),
                  type="numeric",
                  default=10000,
                  short="-hS")
p <- add_argument(p, "--Dver.sct",
                  help="OI, vertical de-correlation lenght  (of the background error correlation) [m]",
                  type="numeric",
                  default=200,
                  short="-vS")
p <- add_argument(p, "--eps2.sct",
                  help="OI, ratio between observation error variance and background error variance. eps2.sct is a vector of positive values (not NAs). If eps2.sct has the same length of the number of input files, then a provider dependent eps2 will be used in the SCT. Otherwise, the value of eps2.sct[1] will be used for all providers and any other eps2.sct[2:n] value will be ignored",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-eS")
p <- add_argument(p, "--thr.sct",
                  help="SCT threshold. flag observation if: (obs-Cross_Validation_pred)^2/(varObs+varCVpred) > thr.sct. thr.sct is a vector of positive values (not NAs). If thr.sct has the same length of the number of input files, then a provider dependent threshold will be used in the SCT. Otherwise, the value of thr.sct[1] will be used for all providers and any other thr.sct[2:n] value will be ignored ",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-tS")
p <- add_argument(p, "--thrpos.sct",
                  help="SCT threshold. flag observation if: (obs-Cross_Validation_pred)^2/(varObs+varCVpred) > thrpos.sct AND (obs-Cross_Validation_pred)>=0. thrpos.sct is a vector of positive values (not NAs). If thrpos.sct has the same length of the number of input files, then a provider dependent threshold will be used in the SCT. Otherwise, the value of thrpos.sct[1] will be used for all providers and any other thrpos.sct[2:n] value will be ignored ",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-tpS")
p <- add_argument(p, "--thrneg.sct",
                  help="SCT threshold. flag observation if: (obs-Cross_Validation_pred)^2/(varObs+varCVpred) > thrneg.sct AND (obs-Cross_Validation_pred)<0.  thrneg.sct is a vector of positive values (not NAs). If thrneg.sct has the same length of the number of input files, then a provider dependent threshold will be used in the SCT. Otherwise, the value of thrneg.sct[1] will be used for all providers and any other thrneg.sct[2:n] value will be ignored ",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-tnS")
p <- add_argument(p, "--laf.sct",
                  help="use land area fraction in the OI correlation function (0-100%)",
                  flag=T,
                  short="-lS")
p <- add_argument(p, "--lafmin.sct",
                  help="land area fraction influence in the OI correlation function",
                  type="numeric",
                  default=0.5,
                  short="-lmS")
p <- add_argument(p, "--fast.sct",
                  help="faster spatial consistency test. Allow for flagging more than one observation simulataneously",
                  flag=T,
                  short="-fS")
p <- add_argument(p, "--smartbox.sct",
                  help="use smart boxes in the spatial consistency test for temperature",
                  flag=T)
p <- add_argument(p, "--stn_by_stn.sct",
                  help="spatial consistency test, station by station mode",
                  flag=T)
p <- add_argument(p, "--corr.sct",
                  help="correlation function to use (''gaussian'',''soar''",
                  type="character",
                  default="gaussian")
p <- add_argument(p, "--box_o_nearest_halfwidth.sct",
                  help="half-width of the square box used to select the nearest observations",
                  type="numeric",
                  default=100000)
p <- add_argument(p, "--pmax.sct",
                  help="maximum number of observations to use in the neighbourhood of each observations",
                  type="integer",
                  default=50)
p <- add_argument(p, "--succ_corr.sct",
                  help="successive correction step (yes/no)",
                  flag=T)
p <- add_argument(p, "--o_errvar_min.sct",
                  help="minimum allowed value for the observation error variance",
                  type="numeric",
                  default=0.001)
p <- add_argument(p, "--o_errvar_max.sct",
                  help="maximum allowed value for the observation error variance",
                  type="numeric",
                  default=4)
p <- add_argument(p, "--xa_errvar_min.sct",
                  help="minimum allowed value for the analysis error variance",
                  type="numeric",
                  default=0.001)
p <- add_argument(p, "--xa_errvar_max.sct",
                  help="maximum allowed value for the analysis error variance",
                  type="numeric",
                  default=4)
p <- add_argument(p, "--fg_lab.sct",
                  help="method used to create the first-guess (\"linear\",\"Frei\",NA used together with usefg(e).sct )",
                  type="character",
                  default="Frei")
p <- add_argument(p, "--fg_gamma.sct",
                  help="lapse rate value",
                  type="numeric",
                  default=-0.0065)
#.............................................................................. 
# observation representativeness
p <- add_argument(p, "--mean.corep",
                  help="average coefficient for the observation representativeness. mean.corep is a vector of positive values (not NAs). If mean.corep has the same length of the number of input files, then a provider dependent value will be used. Otherwise, the value of mean.corep[1] will be used for all providers and any other mean.corep[2:n] value will be ignored",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-avC")
p <- add_argument(p, "--min.corep",
     help="minimum value for the coefficient for the observation representativeness. If min.corep has the same length of the number of input files, then a provider dependent value will be used. Otherwise, the value of min.corep[1] will be used for all providers and any other min.corep[2:n] value will be ignored",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-mnC")
p <- add_argument(p, "--max.corep",
     help="maximum value for the coefficient for the observation representativeness. If max.corep has the same length of the number of input files, then a provider dependent value will be used. Otherwise, the value of max.corep[1] will be used for all providers and any other max.corep[2:n] value will be ignored",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-mxC")
p <- add_argument(p, "--const.corep",
     help="value assigned to the coefficient for the observation representativeness. If const.corep has the same length of the number of input files, then a provider dependent value will be used. Otherwise, the value of const.corep[1] will be used for all providers and any other const.corep[2:n] value will be ignored. If specified, const.corep has priority over the other corep parameters.",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-coC")
#.............................................................................. 
# laf
p <- add_argument(p, "--laf.file",
                  help="land area fraction file (netCDF in kilometric coordinates)",
                  type="character",
                  default=NULL,
                  short="-lfS")
p <- add_argument(p, "--proj4laf",
                  help="proj4 string for the laf",
                  type="character",
                  default=proj4_where_dqc_is_done_default,
                  short="-pl")
p <- add_argument(p, "--laf.varname",
                  help="variable name in the netCDF file",
                  type="character",
                  default="land_area_fraction",
                  short="-lfv")
p <- add_argument(p, "--laf.topdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the laf upside down",
                  flag=T,
                  short="-lftd")
p <- add_argument(p, "--laf.ndim",
                  help="number of dimensions in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-lfnd")
p <- add_argument(p, "--laf.tpos",
                  help="position of the dimension ''time'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-lfti")
p <- add_argument(p, "--laf.dimnames",
                  help="dimension names in the netCDF file",
                  type="character",
                  default=c("x","y","time"),
                  short="-lfna",
                  nargs=Inf)
p <- add_argument(p, "--laf.proj4_var",
                  help="variable that include the specification of the proj4 string",
                  type="character",
                  default="projection_lambert",
                  short="-lfp4v")
p <- add_argument(p, "--laf.proj4_att",
                  help="attribute with the specification of the proj4 string",
                  type="character",
                  default="proj4",
                  short="-lfp4a")
#.............................................................................. 
# check elevation against dem
p <- add_argument(p, "--dem",
     help="check elevation against digital elevation model (dem)",
                  flag=T,
                  short="-dm")
p <- add_argument(p, "--dz.dem",
     help="maximum allowed deviation between observation and dem elevations [m]",
                  type="numeric",
                  default=500,
                  short="-zD")
p <- add_argument(p, "--dem.fill",
                  help="fill missing elevation with data from dem",
                  flag=T,
                  short="-df")
p <- add_argument(p, "--dem.file",
     help="land area fraction file (netCDF in kilometric coordinates)",
                  type="character",
                  default=NULL,
                  short="-dmf")
p <- add_argument(p, "--proj4dem",
                  help="proj4 string for the dem",
                  type="character",
                  default=proj4_where_dqc_is_done_default,
                  short="-pd")
p <- add_argument(p, "--dem.varname",
                  help="variable name in the netCDF file",
                  type="character",
                  default="altitude",
                  short="-dmv")
p <- add_argument(p, "--dem.topdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the dem upside down",
                  flag=T,
                  short="-dmtd")
p <- add_argument(p, "--dem.ndim",
                  help="number of dimensions in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-dmnd")
p <- add_argument(p, "--dem.tpos",
                  help="position of the dimension ''time'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-dmti")
p <- add_argument(p, "--dem.dimnames",
                  help="dimension names in the netCDF file",
                  type="character",
                  default=c("x","y","time"),
                  short="-dmna",
                  nargs=Inf)
p <- add_argument(p, "--dem.proj4_var",
                  help="variable that include the specification of the proj4 string",
                  type="character",
                  default="projection_lambert",
                  short="-dmp4v")
p <- add_argument(p, "--dem.proj4_att",
                  help="attribute with the specification of the proj4 string",
                  type="character",
                  default="proj4",
                  short="-dmp4a")
#.............................................................................. 
# precipitation and temperature cross-check
p <- add_argument(p, "--ccrrt",
           help="precipitation (in-situ) and temperature (field) cross-check",
                  flag=T,
                  short="-ccrrtf")
p <- add_argument(p, "--ccrrt.tmin",
                  help="temperature thresholds (vector, negative values start with \"_\" (e.g. _2=-2)",
#                  type="numeric",
                  type="character",
                  default=NA,
                  nargs=Inf,
                  short="-ccrrtt")
p <- add_argument(p,"--ccrrt.filesetup",
                  help="predefined setup to read gridded fields: dem, t2m. available options: meps",
                  type="character",
                  default=NULL,
                  short="-ccrrtfs")
p <- add_argument(p,"--ccrrt.filemeps",
                  help="meps netCDF file",
                  type="character",
                  default=NULL,
                  short="-ccrrtfm")
#.............................................................................. 
# cool check
p <- add_argument(p, "--cool",
                  help="cool (Check fOr hOLes in the field) test",
                  flag=T)
p <- add_argument(p, "--i.cool",
                  help="number of cool-test iterations",
                  type="integer",
                  default=1)
p <- add_argument(p, "--thres.cool",
                  help="numeric vector with the thresholds used to define events (same units of the specified variable). A threshold transforms an observation into a binary event: observation is less than the threshold OR observations is greater or equal to it.",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--condition.cool",
                  help="character vector specifying the conditions to apply at each of the thresholds (\"lt\"=less than; \"le\"=less or equal than; \"gt\"=greater than; \"ge\"=greater or equal than).",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--n.cool",
                  help="minimum acceptable number of observations producing a ''hole'' in the field (specified as a function of the provider and the threshold). If a clump of connected cells originates from a small number of observations, then it cannot be properly resolved by the observational network. As a result, the observations causing those small-scale events are assumed to be affected by large representativeness errors and flagged as ''suspect''. Format: array (ntot[thres1],nprid1[thres1],nprid2[thres1],...,ntot[thres2],nprid1[thres2],nprid2[thres2],...)",
                  type="integer",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--grid_res.cool",
                  help="resolution of the grid used to compute the field (same units as ).",
                  type="integer",
                  default=1000)
p <- add_argument(p, "--dh_max.cool",
                  help="gridpoints having the nearest observation more than dh_max units apart are set to NAs.",
                  type="integer",
                  default=100000)
#.............................................................................. 
# blacklist
# specified by triple/pairs of numbers: either (lat,lon,IDprovider) OR (index,IDprovider)
p <- add_argument(p, "--blacklist.lat",
                  help="observation blacklist (latitude)",
                  type="numeric",default=NA,nargs=Inf,short="-bla")
p <- add_argument(p, "--blacklist.lon",
                  help="observation blacklist (longitude)",
                  type="numeric",default=NA,nargs=Inf,short="-blo")
p <- add_argument(p, "--blacklist.fll",
                  help="observation blacklist (ID provider)",
                  type="numeric",default=NA,nargs=Inf,short="-bfll")
p <- add_argument(p, "--blacklist.idx",
                  help="observation blacklist (position in input file)",
                  type="numeric",default=NA,nargs=Inf,short="-bix")
p <- add_argument(p, "--blacklist.fidx",
                  help="observation blacklist (ID provider)",
                  type="numeric",default=NA,nargs=Inf,short="-bfix")
#.............................................................................. 
# keep (keep them)
# specified by triple/pairs of numbers: either (lat,lon,IDprovider) OR (index,IDprovider)
p <- add_argument(p, "--keeplist.lat",
                  help="observation keeplist (latitude)",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-kla")
p <- add_argument(p, "--keeplist.lon",
                  help="observation keeplist (longitude)",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-klo")
p <- add_argument(p, "--keeplist.fll",
                  help="observation keeplist (ID provider)",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-kfll")
p <- add_argument(p, "--keeplist.idx",
                  help="observation keeplist (position in input file)",
                  type="numeric",
                  default=NA,
                  nargs=Inf, 
                  short="-kix")
p <- add_argument(p, "--keeplist.fidx",
                  help="observation keeplist (ID provider)",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-kfix")
#.............................................................................. 
# radar-derived precipitation as output
p <- add_argument(p, "--radarout",
                  help="include the radar-derived precipitation as output.",
                  flag=T,
                  short="-rado")
p <- add_argument(p, "--radarout.prid",
                  help="provider identifier for the radar data",
                  type="numeric",
                  default=100,
                  short="-radop")
p <- add_argument(p, "--radarout.aggfact",
                  help="aggregation factor for the radar-derived precipitation",
                  type="numeric",
                  default=3,
                  short="-radop")
#.............................................................................. 
# doit flags
comstr<-" Decide if the test should be applied to all, none or only to a selection of observations based on the provider. Possible values are 0, 1, 2. It is possible to specify either one global value or one value for each provider. Legend: 1 => the observations will be used in the elaboration and they will be tested; 0 => the observations will not be used and they will not be tested; 2 => the observations will be used but they will not be tested."
p <- add_argument(p, "--doit.buddy",
                  help=paste("customize the buddy-test application.",comstr),
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-dob")
p <- add_argument(p, "--doit.buddy_eve",
                  help=paste("customize the buddy_eve-test application.",comstr),
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-dobe")
p <- add_argument(p, "--doit.sct",
                  help=paste("customize the application of SCT.",comstr),
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-dos")
p <- add_argument(p, "--doit.clim",
                  help=paste("customize the application of the climatological check.",comstr),
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-doc")
p <- add_argument(p, "--doit.dem",
                  help=paste("customize the application of the test of observation elevation against the digital elevation model.",comstr),
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-dod")
p <- add_argument(p, "--doit.isol",
                  help=paste("customize the application of the isolation test.",comstr),
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-doi")
p <- add_argument(p, "--doit.fg",
 help=paste("customize the application of the check against a deterministic first-guess field.",comstr),
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-dofg")
p <- add_argument(p, "--doit.fge",
 help=paste("customize the application of the check against an ensemble of first-guess fields.",comstr),
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-dofge")
p <- add_argument(p, "--doit.cool",
                  help=paste("customize the cool check application.",comstr),
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-dopud")
#.............................................................................. 
# precipitation correction for the wind-induced undercatch
p <- add_argument(p, "--rr.wcor",
              help=" precipitation correction for the wind-induced undercatch",
                  flag=T)
p <- add_argument(p,"--rr.wcor.filesetup",
                  help="predefined setup to read gridded fields: dem, t2m, wspeed. available options: meps",
                  type="character",
                  default=NULL,
                  short="-rrwf")
p <- add_argument(p,"--rr.wcor.filemeps",
                  help="meps netCDF file",
                  type="character",
                  default=NULL,
                  short="-rrwfm")
# parameter 
p <- add_argument(p, "--rr.wcor.par",
                  help=paste("Parameter used for correcting wind-induced loss",
                              " of solid precipitation (Wolff et al., 2015)"),
                  type="numeric",
                  default=c(4.24,1.81,0.18,0.99,0.66,1.07,0.18,0.11,2.35,0.12),
                  nargs=10,
                  short="-rrwpar")
# temp 
p <- add_argument(p,"--t2m.file",
                  help="air temperature netCDF file",
                  type="character",
                  default=NULL,
                  short="-rrwt")
p <- add_argument(p, "--t2m.offset",
                  help="air temperature offset",
                  type="numeric",
                  default=0,
                  short="-rrwto")
p <- add_argument(p, "--t2m.cfact",
                  help="air temperature correction factor",
                  type="numeric",
                  default=1,
                  short="-rrwtc")
p <- add_argument(p, "--t2m.negoffset",
                  help="offset sign (1=neg, 0=pos)",
                  type="numeric",
                  default=0,
                  short="-rrwtos")
p <- add_argument(p, "--t2m.negcfact",
                  help="correction factor sign (1=neg, 0=pos)",
                  type="numeric",
                  default=0,
                  short="-rrwton")
p <- add_argument(p, "--proj4t2m",
                  help="proj4 string for the air temperature file",
                  type="character",
                  default=proj4_where_dqc_is_done_default,
                  short="-rrwtp")
p <- add_argument(p, "--t2m.varname",
                  help="air temperature variable name in the netCDF file",
                  type="character",
                  default=NULL,
                  short="-rrwtv")
p <- add_argument(p, "--t2m.topdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the file upside down",
                  flag=T,
                  short="-rrwtt")
p <- add_argument(p, "--t2m.ndim",
                  help="number of dimensions in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-rrwtn")
p <- add_argument(p, "--t2m.tpos",
                  help="position of the dimension ''time'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-rrwtti")
p <- add_argument(p, "--t2m.epos",
                  help="position of the dimension ''ensemble'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-rrwtei")
p <- add_argument(p, "--t2m.dimnames",
                  help="dimension names in the netCDF file",
                  type="character",
                  default=NA,
                  short="-rrwtna",
                  nargs=Inf)
p <- add_argument(p, "--t2m.proj4_var",
                  help="variable that include the specification of the proj4 string",
                  type="character",
                  default="projection_lambert",
                  short="-rrwtp4v")
p <- add_argument(p, "--t2m.proj4_att",
                  help="attribute with the specification of the proj4 string",
                  type="character",
                  default="proj4",
                  short="-rrwtp4a")
p <- add_argument(p, "--t2m.t",
                  help="timestamp to read in the air temperature netCDF file (YYYYMMDDHH00)",
                  type="character",
                  default=NA,
                  short="-rrwttt")
p <- add_argument(p, "--t2m.e",
                  help="ensemble member to read in the air temperature netCDF file",
                  type="numeric",
                  default=NA,
                  short="-rrwtee")
# dem for temperature adjustments
p <- add_argument(p,"--t2m.demfile",
                  help="dem file associated to the first-guess or background file",
                  type="character",
                  default=NULL,
                  short="-rrwdf")
p <- add_argument(p, "--t2m.demoffset",
                  help="offset",
                  type="numeric",
                  default=0,
                  short="-rrwdoff")
p <- add_argument(p, "--t2m.demcfact",
                  help="correction factor",
                  type="numeric",
                  default=1,
                  short="-rrwdcf")
p <- add_argument(p, "--t2m.demnegoffset",
                  help="offset sign (1=neg, 0=pos)",
                  type="numeric",
                  default=0,
                  short="-rrwdnoff")
p <- add_argument(p, "--t2m.demnegcfact",
                  help="correction factor sign (1=neg, 0=pos)",
                  type="numeric",
                  default=0,
                  short="-rrwdncf")
p <- add_argument(p, "--t2m.demt",
                  help="timestamp to read in the netCDF file (YYYYMMDDHH00)",
                  type="character",
                  default=NA,
                  short="-rrwdt")
p <- add_argument(p, "--t2m.demvarname",
                  help="variable name in the netCDF file (dem associated to the first-guess)",
                  type="character",
                  default="none",
                  short="-rrwdv")
p <- add_argument(p, "--t2m.demtopdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the field upside down",
                  flag=T,
                  short="-rrwdtd")
p <- add_argument(p, "--t2m.demndim",
                  help="number of dimensions in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-rrwdnd")
p <- add_argument(p, "--t2m.demepos",
                  help="position of the dimension ''ensemble'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-rrwdei")
p <- add_argument(p, "--t2m.demtpos",
                  help="position of the dimension ''time'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-rrwdti")
p <- add_argument(p, "--t2m.deme",
                  help="ensemble member to read in the netCDF file",
                  type="numeric",
                  default=NA,
                  short="-rrwdee")
p <- add_argument(p, "--t2m.demdimnames",
                  help="dimension names in the netCDF file",
                  type="character",
                  default=NA,
                  short="-rrwdna",
                  nargs=Inf)
# wind
p <- add_argument(p,"--wind.file",
                  help="air temperature netCDF file",
                  type="character",
                  default=NULL,
                  short="-rrww")
p <- add_argument(p, "--proj4wind",
                  help="proj4 string for the air temperature file",
                  type="character",
                  default=proj4_where_dqc_is_done_default,
                  short="-rrwwp")
p <- add_argument(p, "--windspeed.varname",
                  help="air temperature variable name in the netCDF file",
                  type="character",
                  default=NULL,
                  short="-rrwwv")
p <- add_argument(p, "--u.varname",
                  help="air temperature variable name in the netCDF file",
                  type="character",
                  default=NULL,
                  short="-rrwwv")
p <- add_argument(p, "--v.varname",
                  help="air temperature variable name in the netCDF file",
                  type="character",
                  default=NULL,
                  short="-rrwwv")
p <- add_argument(p, "--wind.topdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the file upside down",
                  flag=T,
                  short="-rrwwt")
p <- add_argument(p, "--wind.ndim",
                  help="number of dimensions in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-rrwwn")
p <- add_argument(p, "--wind.tpos",
                  help="position of the dimension ''time'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-rrwwti")
p <- add_argument(p, "--wind.epos",
                  help="position of the dimension ''ensemble'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-rrwwei")
p <- add_argument(p, "--wind.dimnames",
                  help="dimension names in the netCDF file",
                  type="character",
                  default=NA,
                  short="-rrwwna",
                  nargs=Inf)
p <- add_argument(p, "--wind.proj4_var",
                  help="variable that include the specification of the proj4 string",
                  type="character",
                  default="projection_lambert",
                  short="-rrwwp4v")
p <- add_argument(p, "--wind.proj4_att",
                  help="attribute with the specification of the proj4 string",
                  type="character",
                  default="proj4",
                  short="-rrwwp4a")
p <- add_argument(p, "--wind.t",
                  help="timestamp to read in the air temperature netCDF file (YYYYMMDDHH00)",
                  type="character",
                  default=NA,
                  short="-rrwwtt")
p <- add_argument(p, "--wind.e",
                  help="ensemble member to read in the air temperature netCDF file",
                  type="numeric",
                  default=NA,
                  short="-rrwwee")
#.............................................................................. 
# first-guess or background file / deterministic
p <- add_argument(p, "--fg",
        help="check against a deteministic first-guess (fg) field on a regular grid",
                  flag=T)
p <- add_argument(p, "--thr.fg",
     help="maximum allowed deviation between observation and fg (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--thrpos.fg",
     help="maximum allowed deviation between observation and fg (if obs>fg) (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--thrneg.fg",
     help="maximum allowed deviation between observation and fg (if obs<fg) (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--perc.fg_minval",
     help="do the prec.fg test only for values greater than the specified threshold (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--thrperc.fg",
     help="maximum allowed deviation between observation and fg (as %, e.g. 0.1=10%) (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--thrposperc.fg",
     help="maximum allowed deviation between observation and fg (if obs>fg, as %, e.g. 0.1=10%) (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--thrnegperc.fg",
     help="maximum allowed deviation between observation and fg (if obs<fg, as %, e.g. 0.1=10%) (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--fg.dodqc",
                  help="check the first-guess field for weird values",
                  type="logical",
                  default=T)
p <- add_argument(p,"--fg.type",
                  help="file type for the first-guess (e.g., meps)",
                  type="character",
                  default=NULL,
                  short="-fgty")
p <- add_argument(p,"--fg.file",
                  help="first-guess or background file",
                  type="character",
                  default=NULL,
                  short="-fgf")
p <- add_argument(p, "--fg.offset",
                  help="offset",
                  type="numeric",
                  default=0,
                  short="-fgoff")
p <- add_argument(p, "--fg.cfact",
                  help="correction factor",
                  type="numeric",
                  default=1,
                  short="-fgcf")
p <- add_argument(p, "--fg.negoffset",
                  help="offset sign (1=neg, 0=pos)",
                  type="numeric",
                  default=0,
                  short="-fgnoff")
p <- add_argument(p, "--fg.negcfact",
                  help="correction factor sign (1=neg, 0=pos)",
                  type="numeric",
                  default=0,
                  short="-fgncf")
p <- add_argument(p, "--proj4fg",
                  help="proj4 string for the first-guess file",
                  type="character",
                  default=proj4_where_dqc_is_done_default,
                  short="-pfg")
p <- add_argument(p, "--usefg.sct",
         help="use the first-guess field provided as the SCT-background",
                  flag=T)
p <- add_argument(p, "--usefg.cool",
         help="use the first-guess field provided in the cool test",
                  flag=T)
p <- add_argument(p, "--fg.varname",
                  help="variable name in the netCDF file",
                  type="character",
                  default="land_area_fraction",
                  short="-fgv")
p <- add_argument(p, "--fg.topdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the fg upside down",
                  flag=T,
                  short="-fgtd")
p <- add_argument(p, "--fg.ndim",
                  help="number of dimensions in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-fgnd")
p <- add_argument(p, "--fg.tpos",
                  help="position of the dimension ''time'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-fgti")
p <- add_argument(p, "--fg.epos",
                  help="position of the dimension ''ensemble'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-fgti")
p <- add_argument(p, "--fg.dimnames",
                  help="dimension names in the netCDF file",
                  type="character",
                  default=NA,
                  short="-fgna",
                  nargs=Inf)
p <- add_argument(p, "--fg.proj4_var",
                  help="variable that include the specification of the proj4 string",
                  type="character",
                  default="projection_lambert",
                  short="-fgp4v")
p <- add_argument(p, "--fg.proj4_att",
                  help="attribute with the specification of the proj4 string",
                  type="character",
                  default="proj4",
                  short="-fgp4a")
p <- add_argument(p, "--fg.t",
                  help="timestamp to read from the netCDF file (YYYYMMDDHH00)",
                  type="character",
                  default=NA,
                  short="-fgtt")
p <- add_argument(p, "--fg.e",
                  help="ensemble member to read in the netCDF file",
                  type="numeric",
                  default=NA,
                  short="-fgee")
p <- add_argument(p, "--fg.acc",
                  help="first-guess field is accumulated",
                  flag=T,
                  short="-fgacc")
p <- add_argument(p,"--fg.demfile",
                  help="dem file associated to the first-guess or background file",
                  type="character",
                  default=NULL,
                  short="-fgdf")
p <- add_argument(p, "--fg.demoffset",
                  help="offset",
                  type="numeric",
                  default=0,
                  short="-fgdoff")
p <- add_argument(p, "--fg.demcfact",
                  help="correction factor",
                  type="numeric",
                  default=1,
                  short="-fgdcf")
p <- add_argument(p, "--fg.demnegoffset",
                  help="offset sign (1=neg, 0=pos)",
                  type="numeric",
                  default=0,
                  short="-fgdnoff")
p <- add_argument(p, "--fg.demnegcfact",
                  help="correction factor sign (1=neg, 0=pos)",
                  type="numeric",
                  default=0,
                  short="-fgdncf")
p <- add_argument(p, "--fg.demt",
                  help="timestamp to read in the netCDF file (YYYYMMDDHH00)",
                  type="character",
                  default=NA,
                  short="-fgdt")
p <- add_argument(p, "--fg.demvarname",
                  help="variable name in the netCDF file (dem associated to the first-guess)",
                  type="character",
                  default="none",
                  short="-fgdv")
p <- add_argument(p, "--fg.demtopdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the field upside down",
                  flag=T,
                  short="-fgdtd")
p <- add_argument(p, "--fg.demndim",
                  help="number of dimensions in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-fgdnd")
p <- add_argument(p, "--fg.demepos",
                  help="position of the dimension ''ensemble'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-fgdei")
p <- add_argument(p, "--fg.demtpos",
                  help="position of the dimension ''time'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-fgdti")
p <- add_argument(p, "--fg.deme",
                  help="ensemble member to read in the netCDF file",
                  type="numeric",
                  default=NA,
                  short="-fgdee")
p <- add_argument(p, "--fg.demdimnames",
                  help="dimension names in the netCDF file",
                  type="character",
                  default=NA,
                  short="-fgdna",
                  nargs=Inf)
#.............................................................................. 
# first-guess or background file / ensemble
p <- add_argument(p, "--fge",
        help="check against an ensemble of first-guess (fge) fields on a regular grid",
                  flag=T)
p <- add_argument(p, "--thr.fge",
     help="maximum allowed deviation between observation and fg (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--thrpos.fge",
     help="maximum allowed deviation between observation and fg (if obs>fg) (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--thrneg.fge",
     help="maximum allowed deviation between observation and fg (if obs<fg) (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--perc.fge_minval",
     help="do the prec.fg test only for values greater than the specified threshold (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--thrperc.fge",
     help="maximum allowed deviation between observation and fg (as %, e.g. 0.1=10%) (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--thrposperc.fge",
     help="maximum allowed deviation between observation and fg (if obs>fg, as %, e.g. 0.1=10%) (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--thrnegperc.fge",
     help="maximum allowed deviation between observation and fg (if obs<fg, as %, e.g. 0.1=10%) (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--throut.fge",
     help="observation is an outlier if abs(obs-mean_ens)/sd_ens>threshold (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--throutpos.fge",
     help="same as throut.fge, used only if obs>fg (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--throutneg.fge",
     help="same as throut.fge, used only if obs<fg (provider dependent)",
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--sdmin.fge",
                  help="ensemble standard deviation, minimum value",
                  type="numeric",
                  default=NULL,
                  short="-fgesdmn")
p <- add_argument(p,"--fge.type",
                  help="file type for the ensemble first-guess (e.g., meps)",
                  type="character",
                  default=NULL,
                  short="-fgety")
p <- add_argument(p, "--fge.offset",
                  help="offset",
                  type="numeric",
                  default=0,
                  short="-fgedoff")
p <- add_argument(p, "--fge.cfact",
                  help="correction factor",
                  type="numeric",
                  default=1,
                  short="-fgedcf")
p <- add_argument(p, "--fge.negoffset",
                  help="offset sign (1=neg, 0=pos)",
                  type="numeric",
                  default=0,
                  short="-fgednoff")
p <- add_argument(p, "--fge.negcfact",
                  help="correction factor sign (1=neg, 0=pos)",
                  type="numeric",
                  default=0,
                  short="-fgedncf")
p <- add_argument(p,"--fge.file",
                  help="file with the ensemble members",
                  type="character",
                  default=NULL,
                  short="-fgef")
p <- add_argument(p, "--proj4fge",
                  help="proj4 string for the first-guess file",
                  type="character",
                  default=proj4_where_dqc_is_done_default,
                  short="-pfge")
p <- add_argument(p, "--usefge.sct",
         help="use the ensemble mean as the SCT-background",
                  flag=T)
p <- add_argument(p, "--fge.varname",
                  help="variable name in the netCDF file",
                  type="character",
                  default="land_area_fraction",
                  short="-fgev")
p <- add_argument(p, "--fge.topdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the fge upside down",
                  flag=T,
                  short="-fgetd")
p <- add_argument(p, "--fge.acc",
                  help="ensemble first-guess field is accumulated",
                  flag=T,
                  short="-fgacc")
p <- add_argument(p, "--fge.ndim",
                  help="number of dimensions in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-fgend")
p <- add_argument(p, "--fge.tpos",
                  help="position of the dimension ''time'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-fgeti")
p <- add_argument(p, "--fge.epos",
                  help="position of the dimension ''ensemble'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-fgeti")
p <- add_argument(p, "--fge.dimnames",
                  help="dimension names in the netCDF file",
                  type="character",
                  default=NA,
                  short="-fgena",
                  nargs=Inf)
p <- add_argument(p, "--fge.proj4_var",
                  help="variable that include the specification of the proj4 string",
                  type="character",
                  default="projection_lambert",
                  short="-fgep4v")
p <- add_argument(p, "--fge.proj4_att",
                  help="attribute with the specification of the proj4 string",
                  type="character",
                  default="proj4",
                  short="-fgep4a")
p <- add_argument(p, "--fge.t",
                  help="timestamp to read in the netCDF file (YYYYMMDDHH00)",
                  type="character",
                  default=NA,
                  short="-fgett")
p <- add_argument(p, "--fge.e",
                  help="ensemble member to read in the netCDF file",
                  type="numeric",
                  default=NA,
                  short="-fgeee")
p <- add_argument(p,"--fge.demfile",
                  help="dem file associated to the first-guess or background file",
                  type="character",
                  default=NULL,
                  short="-fgedf")
p <- add_argument(p, "--fge.demoffset",
                  help="offset",
                  type="numeric",
                  default=0,
                  short="-fgedoff")
p <- add_argument(p, "--fge.demcfact",
                  help="correction factor",
                  type="numeric",
                  default=1,
                  short="-fgedcf")
p <- add_argument(p, "--fge.demnegoffset",
                  help="offset sign (1=neg, 0=pos)",
                  type="numeric",
                  default=0,
                  short="-fgednoff")
p <- add_argument(p, "--fge.demnegcfact",
                  help="correction factor sign (1=neg, 0=pos)",
                  type="numeric",
                  default=0,
                  short="-fgedncf")
p <- add_argument(p, "--fge.demt",
                  help="timestamp to read in the netCDF file (YYYYMMDDHH00)",
                  type="character",
                  default=NA,
                  short="-fgedt")
p <- add_argument(p, "--fge.demvarname",
                  help="variable name in the netCDF file (dem associated to the first-guess)",
                  type="character",
                  default="none",
                  short="-fgedv")
p <- add_argument(p, "--fge.demtopdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the field upside down",
                  flag=T,
                  short="-fgedtd")
p <- add_argument(p, "--fge.demndim",
                  help="number of dimensions in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-fgednd")
p <- add_argument(p, "--fge.demepos",
                  help="position of the dimension ''ensemble'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-fgedei")
p <- add_argument(p, "--fge.demtpos",
                  help="position of the dimension ''time'' in the netCDF file",
                  type="numeric",
                  default=3,
                  short="-fgedti")
p <- add_argument(p, "--fge.deme",
                  help="ensemble member to read in the netCDF file",
                  type="numeric",
                  default=NA,
                  short="-fgedee")
p <- add_argument(p, "--fge.demdimnames",
                  help="dimension names in the netCDF file",
                  type="character",
                  default=NA,
                  short="-fgedna",
                  nargs=Inf)
#.............................................................................. 
# Timestamp valid for all the netcdf files
p <- add_argument(p, "--timestamp",
                  help="timestamp, valid for all the netCDF file (YYYYMMDDHH00)",
                  type="character",
                  default=NA,
                  short="-t")
#.............................................................................. 
# run on several cores 
p <- add_argument(p, "--cores",
                  help="set the number of cores for parallel runs. Rpackage \"parallel\" required. 0 stands for \"use detectCores\". Default do not use it.",
                  type="numeric",
                  default=NA)
#.............................................................................. 
# PARSE arguments
argv <- parse_args(p)
#
#-----------------------------------------------------------------------------
# read configuration file
if (!is.na(argv$config.file)) {
  if (file.exists(argv$config.file)) {
    source(argv$config.file)
    argv_tmp<-append(argv,conf)
    names_argv_tmp<-names(argv_tmp)
    argv_def<-list()
    names_argv_def<-integer(0)
    k<-0
    for (i in 1:length(argv_tmp)) {
      if (names_argv_tmp[i] %in% names_argv_def) next
      k<-k+1
      j<-which(names_argv_tmp==names_argv_tmp[i])
      argv_def[[k]]<-argv_tmp[[j[length(j)]]]
      names_argv_def<-c(names_argv_def,names_argv_tmp[i])
    }
    names(argv_def)<-names_argv_def
    rm(argv_tmp,names_argv_tmp,names_argv_def)
    rm(argv)
    argv<-argv_def
    rm(argv_def)
  } else {
    print("WARNING: config file not found")
    print(argv$config.file)
  }
}
#
#-----------------------------------------------------------------------------
if (!is.na(argv$cores)) {
  suppressPackageStartupMessages(library("parallel"))
  if (argv$cores==0) argv$cores <- detectCores()
  print(paste("--> multi-core run, cores=",argv$cores))
}
#
#-----------------------------------------------------------------------------
# CHECKS on input arguments
if (!file.exists(argv$input)) {
  print("ERROR: input file not found")
  print(argv$input)
  quit(status=1)
}
# more than one input file
if (any(!is.na(argv$input.files))) {
  for (j in 1:length(argv$input.files)) {
    if (!file.exists(argv$input.files[j])) {
      print("ERROR: input file not found")
      print(argv$input.files[j])
      quit(status=1)
    }
  }
  argv$input.files<-c(argv$input,argv$input.files)
} else {
  argv$input.files<-argv$input
}
nfin<-length(argv$input.files)
# check consistency between number of files and provider ids
if (any(is.na(argv$prid))) {
  argv$prid<-1:nfin
} else {
  if (length(argv$prid)!=nfin) {
    print("ERROR: number of provider identifier is different from the number of input files")
    quit(status=1)
  }
}
# set input offsets and correction factors
if (any(is.na(argv$input.offset))) {
  argv$input.offset<-rep(0,length=nfin)
} else {
  if (length(argv$input.offset)!=nfin) 
    argv$input.offset<-rep(argv$input.offset[1],length=nfin)
  aux<-vector(length=nfin,mode="numeric")
  for (i in 1:nfin) aux[i]<-as.numeric(gsub("_","-",argv$input.offset[i]))
  argv$input.offset<-aux
  rm(aux)
  if (!any(is.na(argv$input.negoffset))) {
    if (length(argv$input.negoffset)!=nfin) 
      argv$input.negoffset<-rep(argv$input.negoffset[1],length=nfin)
    argv$input.offset<-argv$input.offset*(-1)**argv$input.negoffset
  }
}
if (any(is.na(argv$input.cfact))) {
  argv$input.cfact<-rep(1,length=nfin)
} else {
  if (length(argv$input.cfact)!=nfin) 
    argv$input.cfact<-rep(argv$input.cfact[1],length=nfin)
  aux<-vector(length=nfin,mode="numeric")
  for (i in 1:nfin) aux[i]<-as.numeric(gsub("_","-",argv$input.cfact[i]))
  argv$input.cfact<-aux
  rm(aux)
  if (!any(is.na(argv$input.negcfact))) {
    if (length(argv$input.negcfact)!=nfin) 
      argv$input.negcfact<-rep(argv$input.negcfact[1],length=nfin)
    argv$input.cfact<-argv$input.cfact*(-1)**argv$input.negcfact
  }
}
# set offsets and correction factors
#if (any(is.na(argv$input.offset))) argv$input.offset<-rep(0,nfin)
#if (any(is.na(argv$input.negoffset))) argv$input.negoffset<-rep(0,nfin)
#if (any(is.na(argv$input.cfact))) argv$input.cfact<-rep(1,nfin)
#if (any(is.na(argv$input.negcfact))) argv$input.negcfact<-rep(0,nfin)
#argv$input.offset<-argv$input.offset*(-1)**argv$input.negoffset
#argv$input.cfact<-argv$input.cfact*(-1)**argv$input.negcfact
# check variable
if (!(argv$variable %in% c("T","RH","RR","SD"))) {
  print("variable must be one of T, RH, RR, SD")
  quit(status=1)
}
# set proj4 variables (proj4from and proj4to are obsolete)
if (argv$proj4from!=argv$proj4_input_obsfiles) {
  if (argv$proj4_input_obsfiles==proj4_input_obsfiles_default & 
      argv$proj4from!=proj4_input_obsfiles_default)
    argv$proj4_input_obsfiles<-argv$proj4from
}
if (argv$proj4to!=argv$proj4_where_dqc_is_done) {
  if (argv$proj4_where_dqc_is_done==proj4_where_dqc_is_done_default & 
      argv$proj4to!=proj4_where_dqc_is_done_default)
    argv$proj4_where_dqc_is_done<-argv$proj4to
}
# set variables to customize output
if (argv$varname.lat.out!=argv$varname.y.out) {
  if (argv$varname.y.out==varname.y.out_default & 
      argv$varname.lat.out!=varname.y.out_default)
    argv$varname.y.out<-argv$varname.lat.out
}
if (argv$varname.lon.out!=argv$varname.x.out) {
  if (argv$varname.x.out==varname.x.out_default & 
      argv$varname.lon.out!=varname.x.out_default)
    argv$varname.x.out<-argv$varname.lon.out
}
if (argv$latlon.dig.out!=argv$xy.dig.out) {
  if (argv$xy.dig.out==xy.dig.out_default & 
      argv$latlon.dig.out!=xy.dig.out_default)
    argv$xy.dig.out<-argv$latlon.dig.out
}
# set the input arguments according to user specification
if (!is.na(argv$fg.type)) {
  if (argv$fg.type=="meps") {
    if (argv$variable=="T") {
      argv$fg.epos<-5
      argv$fg.e<-0
      argv$fg.varname<-"air_temperature_2m"
      argv$fg.ndim<-5 
      argv$fg.tpos<-3
      if (any(is.na(argv$fg.dimnames))) {
        argv$fg.dimnames<-c("x","y","time","height1","ensemble_member")
      } else {
        argv$fg.ndim<-length(argv$fg.dimnames)
      }
      argv$proj4fg<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
      argv$fg.offset<-273.15
      argv$fg.negoffset<-1
      argv$fg.demvarname<-"surface_geopotential" 
      argv$fg.demndim<-5
      argv$fg.demtpos<-3
      argv$fg.demepos<-5
      argv$fg.deme<-0
      argv$fg.demdimnames<-c("x","y","time","height0","ensemble_member")
      # divide geopotential by g=9.80665. This calculates geopotential height (above mean sea level)
      argv$fg.demcfact<-0.0980665 
      argv$fg.topdown<-TRUE
      argv$fg.demtopdown<-TRUE
    } else if (argv$variable=="RR") {
      argv$fg.epos<-5
      argv$fg.e<-0
      argv$fg.varname<-"precipitation_amount_acc"
      argv$fg.ndim<-5 
      argv$fg.tpos<-3
      if (any(is.na(argv$fg.dimnames))) {
        argv$fg.dimnames<-c("x","y","time","height0","ensemble_member")
      } else {
        argv$fg.ndim<-length(argv$fg.dimnames)
      }
      argv$proj4fg<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
      argv$fg.acc<-TRUE
      argv$fg.topdown<-TRUE
    } else if (argv$variable=="RH") {
      argv$fg.epos<-5
      argv$fg.e<-0
      argv$fg.varname<-"relative_humidity_2m"
      argv$fg.ndim<-5 
      argv$fg.tpos<-3
      if (any(is.na(argv$fg.dimnames))) {
        argv$fg.dimnames<-c("x","y","time","height1","ensemble_member")
      } else {
        argv$fg.ndim<-length(argv$fg.dimnames)
      }
      argv$proj4fg<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
      argv$fg.cfact<-100.
      argv$fg.topdown<-TRUE
    } else {
      print("ERROR in --fg.type, combination of type/variable not available")
      quit(status=1)
    } 
  } else if (argv$fg.type=="radar") {
    if (argv$variable=="RR") {
      argv$fg.epos<-NA
      argv$fg.e<-NULL
      argv$fg.varname<-"lwe_precipitation_rate"
      argv$fg.ndim<-3 
      argv$fg.tpos<-3
      if (any(is.na(argv$fg.dimnames))) {
        argv$fg.dimnames<-c("Xc","Yc","time")
      } else {
        argv$fg.ndim<-length(argv$fg.dimnames)
      }
      argv$proj4fg<-"+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
      argv$fg.topdown<-FALSE
    } else {
      print("ERROR in --fg.type, combination of type/variable not available")
      quit(status=1)
    }
  } else if (argv$fg.type=="surfex_T") {
    if (argv$variable=="RH") {
      argv$fg.epos<-NA
      argv$fg.e<-NULL
      argv$fg.varname<-"relative_humidity_2m"
      argv$fg.ndim<-3 
      argv$fg.tpos<-3
      argv$fg.dimnames<-c("x","y","time")
      argv$proj4fg<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
      argv$fg.topdown<-TRUE
    } else if (argv$variable=="T") {
      argv$fg.epos<-NA
      argv$fg.e<-NULL
      argv$fg.varname<-"air_temperature_2m"
      argv$fg.ndim<-3 
      argv$fg.tpos<-3
      argv$fg.dimnames<-c("x","y","time")
      argv$proj4fg<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
      argv$fg.topdown<-TRUE
    } else if (argv$variable=="SD") {
      argv$fg.epos<-NA
      argv$fg.e<-NULL
      argv$fg.varname<-"surface_snow_thickness"
      argv$fg.ndim<-3 
      argv$fg.tpos<-3
      argv$fg.dimnames<-c("x","y","time")
      argv$proj4fg<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
      argv$fg.topdown<-TRUE
    } else {
      print("ERROR in --fg.type, combination of type/variable not available")
      quit(status=1)
    }
  } else {
    print("ERROR in --fg.type, type not recognized")
    quit(status=1)
  }
}
if (!is.na(argv$fge.type)) {
  if (argv$fge.type=="meps") {
    if (argv$variable=="T") {
      argv$fge.epos<-5
      argv$fge.varname<-"air_temperature_2m"
      argv$fge.ndim<-5 
      argv$fge.tpos<-3
      if (is.na(argv$fge.dimnames)) {
        argv$fge.dimnames<-c("x","y","time","height1","ensemble_member")
      } else {
        argv$fge.ndim<-length(argv$fge.dimnames)
      }
      argv$proj4fge<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
      argv$fge.offset<-273.15
      argv$fge.negoffset<-1
      argv$fge.demvarname<-"surface_geopotential" 
      argv$fge.demndim<-5
      argv$fge.demtpos<-3
      argv$fge.demepos<-5
      argv$fge.deme<-0
      if (is.na(argv$fge.demdimnames)) {
        argv$fge.demdimnames<-c("x","y","time","height0","ensemble_member")
      } else {
        argv$fge.demndim<-length(argv$fge.demdimnames)
      }
      # divide geopotential by g=9.80665. This calculates geopotential height (above mean sea level)
      argv$fge.demcfact<-0.0980665 
      argv$fge.topdown<-TRUE
      argv$fge.demtopdown<-TRUE
    } else if (argv$variable=="RR") {
      argv$fge.epos<-5
      argv$fge.varname<-"precipitation_amount_acc"
      argv$fge.ndim<-5 
      argv$fge.tpos<-3
      if (is.na(argv$fge.dimnames)) {
        argv$fge.dimnames<-c("x","y","time","height0","ensemble_member")
      } else {
        argv$fge.ndim<-length(argv$fge.dimnames)
      }
      argv$proj4fge<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
      argv$fge.acc<-TRUE
      argv$fge.topdown<-TRUE
    } else if (argv$variable=="RH") {
      argv$fge.epos<-5
      argv$fge.varname<-"relative_humidity_2m"
      argv$fge.ndim<-5 
      argv$fge.tpos<-3
      if (is.na(argv$fge.dimnames)) {
        argv$fge.dimnames<-c("x","y","time","height1","ensemble_member")
      } else {
        argv$fge.ndim<-length(argv$fge.dimnames)
      }
      argv$proj4fge<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
      argv$fge.cfact<-100.
      argv$fge.topdown<-TRUE
    } else {
      print("ERROR in --fge.type, combination of type/variable not available")
      quit(status=1)
    } 
  } else {
    print("ERROR in --fge.type, type not recognized")
    quit(status=1)
  }
}
# shortcut for meps file in the precip correction for wind undercatch
if (!is.na(argv$rr.wcor.filesetup)) {
  if (argv$rr.wcor.filesetup=="meps") {
    argv$t2m.file<-argv$rr.wcor.filemeps
    argv$t2m.epos<-5
    argv$t2m.e<-0
    argv$t2m.varname<-"air_temperature_2m"
    argv$t2m.ndim<-5 
    argv$t2m.tpos<-3
    argv$t2m.dimnames<-c("x","y","time","height1","ensemble_member")
    argv$proj4t2m<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
    argv$t2m.offset<-273.15
    argv$t2m.negoffset<-1
    argv$t2m.demfile<-argv$rr.wcor.filemeps
    argv$t2m.demvarname<-"surface_geopotential" 
    argv$t2m.demndim<-5
    argv$t2m.demtpos<-3
    argv$t2m.demepos<-5
    argv$t2m.deme<-0
    argv$t2m.demdimnames<-c("x","y","time","height0","ensemble_member")
    # divide geopotential by g=9.80665. This calculates geopotential height (above mean sea level)
    argv$t2m.demcfact<-0.0980665 
    argv$t2m.topdown<-TRUE
    argv$t2m.demtopdown<-TRUE
    argv$wind.file<-argv$rr.wcor.filemeps
    argv$wind.epos<-5
    argv$wind.e<-0
    argv$u.varname<-"x_wind_10m" 
    argv$v.varname<-"y_wind_10m" 
    argv$wind.ndim<-5 
    argv$wind.tpos<-3
    argv$wind.dimnames<-c("x","y","time","height3","ensemble_member")
    argv$proj4wind<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
    argv$wind.topdown<-TRUE
  } else {
    print("ERROR --rr.wcor.filesetup argument not recognized")
    quit(status=1)
  }
}
# shortcut for meps file in the precip-temp cross-check
if (!is.na(argv$ccrrt.filesetup)) {
  if (argv$ccrrt.filesetup=="meps") {
    argv$t2m.file<-argv$ccrrt.filemeps
    argv$t2m.epos<-5
    argv$t2m.e<-0
    argv$t2m.varname<-"air_temperature_2m"
    argv$t2m.ndim<-5 
    argv$t2m.tpos<-3
    argv$t2m.dimnames<-c("x","y","time","height1","ensemble_member")
    argv$proj4t2m<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
    argv$t2m.offset<-273.15
    argv$t2m.negoffset<-1
    argv$t2m.demfile<-argv$ccrrt.filemeps
    argv$t2m.demvarname<-"surface_geopotential" 
    argv$t2m.demndim<-5
    argv$t2m.demtpos<-3
    argv$t2m.demepos<-5
    argv$t2m.deme<-0
    argv$t2m.demdimnames<-c("x","y","time","height0","ensemble_member")
    # divide geopotential by g=9.80665. This calculates geopotential height (above mean sea level)
    argv$t2m.demcfact<-0.0980665 
    argv$t2m.topdown<-TRUE
    argv$t2m.demtopdown<-TRUE
  } else {
    print("ERROR --ccrrt.filesetup argument not recognized")
    quit(status=1)
  }
}
# check external files
if (argv$dem | argv$dem.fill) {
  if (!file.exists(argv$dem.file)) {
    print("ERROR: dem file not found")
    print(argv$dem.file)
    quit(status=1)
  }
}
if (argv$laf.sct) {
  if (!file.exists(argv$laf.file)) {
    print("ERROR: laf file not found")
    print(argv$laf.file)
    quit(status=1)
  }
}
if (!is.na(argv$fg.file)) {
  if (!file.exists(argv$fg.file)) {
    print("ERROR: first-guess file not found")
    print(argv$fg.file)
    quit(status=1)
  }
}
# input column names
if (any(is.na(argv$separator))) 
  argv$separator<-rep(";",nfin)
if (length(argv$separator)==0) 
  argv$separator<-rep(";",nfin)
if (length(argv$separator)!=nfin) 
  argv$separator<-rep(argv$separator[1],nfin)
for (i in 1:length(argv$separator)) if (argv$separator[i]=="comma")  argv$separator[i]<-","
if (any(is.na(argv$varname.lat))) 
  argv$varname.lat<-rep("lat",nfin)
if (length(argv$varname.lat)==0) 
  argv$varname.lat<-rep("lat",nfin)
if (length(argv$varname.lat)!=nfin) 
  argv$varname.lat<-rep(argv$varname.lat[1],nfin)
if (any(is.na(argv$varname.lon))) 
  argv$varname.lon<-rep("lon",nfin)
if (length(argv$varname.lon)==0) 
  argv$varname.lon<-rep("lon",nfin)
if (length(argv$varname.lon)!=nfin) 
  argv$varname.lon<-rep(argv$varname.lon[1],nfin)
if (any(is.na(argv$varname.elev))) 
  argv$varname.elev<-rep("elev",nfin)
if (length(argv$varname.elev)==0) 
  argv$varname.elev<-rep("elev",nfin)
if (length(argv$varname.elev)!=nfin) 
  argv$varname.elev<-rep(argv$varname.elev[1],nfin)
if (any(is.na(argv$varname.value))) 
  argv$varname.value<-rep("value",nfin)
if (length(argv$varname.value)==0) 
  argv$varname.value<-rep("value",nfin)
if (length(argv$varname.value)!=nfin) 
  argv$varname.value<-rep(argv$varname.value[1],nfin)
# TODO: adapt the procedure for input data others than lat-lon
if (!argv$spatconv) {
  print("ERROR: \"--spatconv\" (-c) option must be used on the command line")
  print("input metadata are expected to be lat-lon coordinates")
  print(" some DQC tests takes place in kilometric coordinates specified by the user")
  print("output is in lat-lon coordinates")
  quit(status=1)
}
if (argv$laf.sct | argv$dem | argv$dem.fill |
    !is.na(argv$fg.file) | !is.na(argv$fge.file) |
    !is.na(argv$t2m.file))
  suppressPackageStartupMessages(library("ncdf4")) 
# first-guess file
if (!is.na(argv$fg.file)) {
  if (!file.exists(argv$fg.file)) {
    print("ERROR file not found")
    print(argv$fg.file)
    quit(status=1)
  }
  if ( argv$variable=="T" &
       !file.exists(argv$fg.demfile)) {
    print("ERROR: for temperature, a digital elevation model must be specified together with a first-guess file (det)")
    quit(status=1)
  }
  suppressPackageStartupMessages(library("ncdf4")) 
} else if (argv$fg) {
  print("ERROR no first-guess file provided for the first-guess test")
  quit(status=1)
} else if (argv$usefg.sct) {
  print("ERROR: SCT requested with the background derived from a first-guess field that is not being provided as input")
  quit(status=1)
}
#
if (!is.na(argv$month.clim) & (argv$month.clim<1 | argv$month.clim>12)) {
  print("ERROR: month number is wrong:")
  print(paste("month number=",argv$month.clim))
  quit(status=1)
} else if (!is.na(argv$month.clim) & 
           (length(which(!is.na(argv$vmin.clim)))!=12 | 
            length(which(!is.na(argv$vmax.clim)))!=12) ) {
  print("ERROR: climatological check, vmin.clim and/or vmax.clim vectors must have 12 arguments")
  quit(status=1)
}
# blacklist
if (any(!is.na(argv$blacklist.lat)) | 
    any(!is.na(argv$blacklist.lon)) |
    any(!is.na(argv$blacklist.fll)) ) {
  if ( (length(argv$blacklist.lat)!=length(argv$blacklist.lon))  |
       (length(argv$blacklist.lat)!=length(argv$blacklist.fll))  |
       (any(is.na(argv$blacklist.fll))) | 
       (any(is.na(argv$blacklist.lat))) | 
       (any(is.na(argv$blacklist.lon))) ) {
    print("ERROR in the blacklist definition, must have same number of lat,lon,IDprovider points")
    print(paste("lat number=",argv$blacklist.lat))
    print(paste("lon number=",argv$blacklist.lon))
    print(paste("ID provider number=",argv$blacklist.fll))
    quit(status=1)
  }
}
if (any(!is.na(argv$blacklist.idx)) | 
    any(!is.na(argv$blacklist.fidx)) ) {
  if ( (length(argv$blacklist.idx)!=length(argv$blacklist.fidx))  |
       (any(is.na(argv$blacklist.idx))) | 
       (any(is.na(argv$blacklist.fidx))) ) {
    print("ERROR in the blacklist definition, must have same number of index and IDprovider points")
    print(paste("index number=",argv$blacklist.idx))
    print(paste("ID provider number=",argv$blacklist.fidx))
    quit(status=1)
  }
}
# keeplist
if (any(!is.na(argv$keeplist.lat)) | 
    any(!is.na(argv$keeplist.lon)) |
    any(!is.na(argv$keeplist.fll)) ) {
  if ( (length(argv$keeplist.lat)!=length(argv$keeplist.lon))  |
       (length(argv$keeplist.lat)!=length(argv$keeplist.fll))  |
       (any(is.na(argv$keeplist.fll))) | 
       (any(is.na(argv$keeplist.lat))) | 
       (any(is.na(argv$keeplist.lon))) ) {
    print("ERROR in the keeplist definition, must have same number of lat,lon,IDprovider points")
    print(paste("lat number=",argv$keeplist.lat))
    print(paste("lon number=",argv$keeplist.lon))
    print(paste("ID provider number=",argv$keeplist.fll))
    quit(status=1)
  }
}
if (any(!is.na(argv$keeplist.idx)) | 
    any(!is.na(argv$keeplist.fidx)) ) {
  if ( (length(argv$keeplist.idx)!=length(argv$keeplist.fidx))  |
       (any(is.na(argv$keeplist.idx))) | 
       (any(is.na(argv$keeplist.fidx))) ) {
    print("ERROR in the keeplist definition, must have same number of index and IDprovider points")
    print(paste("index number=",argv$keeplist.idx))
    print(paste("ID provider number=",argv$keeplist.fidx))
    quit(status=1)
  }
}
#
# observation representativeness
if ( (any(is.na(argv$mean.corep)) | 
      any(is.na(argv$min.corep))  | 
      any(is.na(argv$max.corep)) ) & 
      any(is.na(argv$const.corep)) ) {
  print("++WARNING")
  print("parameters related to the coefficient of observation representativeness are not properly specified")
  print("--mean.corep --min.corep and --max.corep or --const.corep should be specified")
  print("Because they are not specified, it is assumed that the coefficient of observation representativeness is not considered an interesting output. As a  consequence, the corep parameters are set to default values (min.corep=0.9 mean.corep=1 max.corep=1.1")
  argv$min.corep<-0.9
  argv$mean.corep<-1
  argv$max.corep<-1.1
}
if (length(argv$min.corep)!=nfin) 
  argv$min.corep<-rep(argv$min.corep[1],length=nfin)
if (length(argv$mean.corep)!=nfin) 
  argv$mean.corep<-rep(argv$mean.corep[1],length=nfin)
if (length(argv$max.corep)!=nfin) 
  argv$max.corep<-rep(argv$max.corep[1],length=nfin)
if (length(argv$const.corep)!=nfin) 
  argv$const.corep<-rep(argv$const.corep[1],length=nfin)
#
# precip and temperature crosscheck
if (length(argv$ccrrt.tmin)!=nfin) 
  argv$ccrrt.tmin<-rep(argv$ccrrt.tmin[1],length=nfin)
aux<-vector(length=nfin,mode="numeric")
for (i in 1:nfin) aux[i]<-as.numeric(gsub("_","-",argv$ccrrt.tmin[i]))
argv$ccrrt.tmin<-aux
rm(aux)
#
# fg
if (argv$fg) {
  if (length(argv$thrpos.fg)!=nfin) 
    argv$thrpos.fg<-rep(argv$thrpos.fg[1],length=nfin)
  if (length(argv$thrneg.fg)!=nfin)
    argv$thrneg.fg<-rep(argv$thrneg.fg[1],length=nfin)
  if (length(argv$thr.fg)!=nfin)
    argv$thr.fg<-rep(argv$thr.fg[1],length=nfin)
  if (length(argv$thrposperc.fg)!=nfin) 
    argv$thrposperc.fg<-rep(argv$thrposperc.fg[1],length=nfin)
  if (length(argv$thrnegperc.fg)!=nfin)
    argv$thrnegperc.fg<-rep(argv$thrnegperc.fg[1],length=nfin)
  if (length(argv$thrperc.fg)!=nfin)
    argv$thrperc.fg<-rep(argv$thrperc.fg[1],length=nfin)
  if (length(argv$perc.fg_minval)!=nfin) 
    argv$perc.fg_minval<-rep(argv$perc.fg_minval[1],length=nfin)
  if ( !any(!is.na(argv$thrpos.fg)) &
       !any(!is.na(argv$thrneg.fg)) &
       !any(!is.na(argv$thr.fg)) &
       !any(!is.na(argv$thrposperc.fg)) &
       !any(!is.na(argv$thrnegperc.fg)) &
       !any(!is.na(argv$thrperc.fg)) &
       !any(!is.na(argv$perc.fg_minval)) ) {
    print("Error in specification of fg-thresholds")
    quit(status=1)
  }
}
#
# fge
if (argv$fge) {
  if (length(argv$thrpos.fge)!=nfin) 
    argv$thrpos.fge<-rep(argv$thrpos.fge[1],length=nfin)
  if (length(argv$thrneg.fge)!=nfin)
    argv$thrneg.fge<-rep(argv$thrneg.fge[1],length=nfin)
  if (length(argv$thr.fge)!=nfin) 
    argv$thr.fge<-rep(argv$thr.fge[1],length=nfin)
  if (length(argv$thrposperc.fge)!=nfin) 
    argv$thrposperc.fge<-rep(argv$thrposperc.fge[1],length=nfin)
  if (length(argv$thrnegperc.fge)!=nfin)
    argv$thrnegperc.fge<-rep(argv$thrnegperc.fge[1],length=nfin)
  if (length(argv$thrperc.fge)!=nfin) 
    argv$thrperc.fge<-rep(argv$thrperc.fge[1],length=nfin)
  if (length(argv$perc.fge_minval)!=nfin) 
    argv$perc.fge_minval<-rep(argv$perc.fge_minval[1],length=nfin)
  if (length(argv$thrposout.fge)!=nfin) 
    argv$thrposout.fge<-rep(argv$thrposout.fge[1],length=nfin)
  if (length(argv$thrnegout.fge)!=nfin)
    argv$thrnegout.fge<-rep(argv$thrnegout.fge[1],length=nfin)
  if (length(argv$throut.fge)!=nfin) 
    argv$throut.fge<-rep(argv$throut.fge[1],length=nfin)
  if ( !any(!is.na(argv$thrpos.fge)) &
       !any(!is.na(argv$thrneg.fge)) &
       !any(!is.na(argv$thr.fge)) &
       !any(!is.na(argv$thrposperc.fge)) &
       !any(!is.na(argv$thrnegperc.fge)) &
       !any(!is.na(argv$thrperc.fge)) &
       !any(!is.na(argv$perc.fge_minval)) &
       !any(!is.na(argv$thrposout.fge)) &
       !any(!is.na(argv$thrnegout.fge)) &
       !any(!is.na(argv$throut.fge)) ) {
    print("Error in specification of fge-thresholds")
    quit(status=1)
  }
}
#
# SCT
# if defined, thrpos.sct and thrneg.sct have the priority on thr.sct
if ( (any(is.na(argv$thrpos.sct)) & any(!is.na(argv$thrpos.sct))) |
     (any(is.na(argv$thrneg.sct)) & any(!is.na(argv$thrneg.sct))) ) {
  print("SCT thresholds for positive and negative deviations are not properly specified")
  print(paste("threshold(s) when (Obs-CVpred) <0 (thrneg.sct)",argv$thrneg.sct))
  print(paste("threshold(s) when (Obs-CVpred)>=0 (thrpos.sct)",argv$thrpos.sct))
  quit(status=1)
}
if (length(argv$thrpos.sct)!=nfin) {
  argv$thrpos.sct<-rep(argv$thrpos.sct[1],length=nfin)
}
if (length(argv$thrneg.sct)!=nfin) {
  argv$thrneg.sct<-rep(argv$thrneg.sct[1],length=nfin)
}
if ( any(is.na(argv$thr.sct)) & 
     is.na(argv$thrneg.sct[1]) ) {
  print("++ WARNING")
  print("thr.sct should be specified and it must not contain NAs")
  print(" because either it has not been specified or it has been set to NA, ")
  print(" then TITAN will use the default value of 16")
  argv$thr.sct<-vector();argv$thr.sct[1]<-16
}
if (length(argv$thr.sct)!=nfin) 
  argv$thr.sct<-rep(argv$thr.sct[1],length=nfin)
#
# eps2
if (any(is.na(argv$eps2.sct))) {
  print("++ WARNING")
  print("eps2.sct should be specified and it must not contain NAs")
  print(" because either it has not been specified or it has been set to NA, ")
  print(" then TITAN will use the default value of 0.5")
  argv$eps2.sct<-vector();argv$eps2.sct[1]<-0.5
}
if (length(argv$eps2.sct)!=nfin) 
  argv$eps2.sct<-rep(argv$eps2.sct[1],length=nfin)
#
# cool test
if (argv$cool) {
  if (any(is.na(argv$thres.cool))) {
    print("++ WARNING")
    print("COOL test has NAs among the specified thresholds")
    print(" TITAN will use only the default threshold of 0.1")
    argv$thres.cool<-vector();argv$thres.cool[1]<-0.1
  }
  if ( any( is.na(argv$condition.cool) | 
       !(argv$condition.cool %in% c("lt","le","gt","ge")) )) {
    print("++ WARNING")
    print("COOL test has NAs and/or not allowed strings among the specified conditions")
    print(" TITAN will use only the default condition \"lt\"")
    argv$condition.cool<-vector();argv$condition.cool[1]<-"lt"
  }
  if (length(argv$condition.cool)!=length(argv$thres.cool)) {
    print("++ WARNING")
    print("COOL test, different lengths for vectors thres.cool and condition.cool")
    print(" TITAN will use the default condition \"lt\" for all thresholds")
    argv$condition.cool<-vector(mode="character",length=length(argv$thres.cool))
    argv$condition.cool[]<-"lt"
  }
  n.cool<-array(data=NA,dim=c(length(argv$thres.cool),(nfin+1))) 
  for (i in 1:length(argv$thres.cool)) {
    for (j in 1:(nfin+1)) { 
      n.cool[i,j]<-argv$n.cool[((i-1)*(nfin+1)+j)]
    }
  }
  if (any(is.na(n.cool))) {
    print("++ ERROR")
    print("COOL test. something wrong in the specification of the n.cool argument")
    print(n.cool)
    quit(status=1)
  }
}
#
# doit flags
if (any(is.na(argv$doit.buddy))) argv$doit.buddy<-rep(1,length=nfin)
if (any(is.na(argv$doit.buddy_eve))) argv$doit.buddy_eve<-rep(1,length=nfin)
if (any(is.na(argv$doit.sct))) argv$doit.sct<-rep(1,length=nfin)
if (any(is.na(argv$doit.clim))) argv$doit.clim<-rep(1,length=nfin)
if (any(is.na(argv$doit.dem))) argv$doit.dem<-rep(1,length=nfin)
if (any(is.na(argv$doit.isol))) argv$doit.isol<-rep(1,length=nfin)
if (any(is.na(argv$doit.fg))) argv$doit.fg<-rep(1,length=nfin)
if (any(is.na(argv$doit.fge))) argv$doit.fge<-rep(1,length=nfin)
if (any(is.na(argv$doit.cool))) argv$doit.cool<-rep(1,length=nfin)
if (any(!(argv$doit.buddy %in% c(0,1,2)))) {
  print("doit.buddy must contain only 0,1,2")
  quit(status=1)
}
if (any(!(argv$doit.buddy_eve %in% c(0,1,2)))) {
  print("doit.buddy_eve must contain only 0,1,2")
  quit(status=1)
}
if (any(!(argv$doit.sct %in% c(0,1,2)))) {
  print("doit.sct must contain only 0,1,2")
  quit(status=1)
}
if (any(!(argv$doit.clim %in% c(0,1,2)))) {
  print("doit.clim must contain only 0,1,2")
  quit(status=1)
}
if (any(!(argv$doit.dem %in% c(0,1,2)))) {
  print("doit.dem must contain only 0,1,2")
  quit(status=1)
}
if (any(!(argv$doit.isol %in% c(0,1,2)))) {
  print("doit.isol must contain only 0,1,2")
  quit(status=1)
}
if (any(!(argv$doit.fg %in% c(0,1,2)))) {
  print("doit.fg must contain only 0,1,2")
  quit(status=1)
}
if (any(!(argv$doit.cool %in% c(0,1,2)))) {
  print("doit.cool must contain only 0,1,2")
  quit(status=1)
}
#
# set the thresholds for the plausibility check
if (!is.na(argv$tmin) & is.na(argv$vmin)) argv$vmin<-argv$tmin
if (!is.na(argv$tmax) & is.na(argv$vmax)) argv$vmax<-argv$tmax
if (any(is.na(argv$vmin.clim)) & !any(is.na(argv$tmin.clim))) 
  argv$vmin.clim<-argv$tmin.clim
if (any(is.na(argv$vmax.clim)) & !any(is.na(argv$tmax.clim))) 
  argv$vmax.clim<-argv$tmax.clim
argv$vmin<-as.numeric(gsub("_","-",argv$vmin))
argv$vmin<-argv$vmin*(-1)**argv$vminsign
argv$vmax<-as.numeric(gsub("_","-",argv$vmax))
argv$vmax<-argv$vmax*(-1)**argv$vmaxsign
argv$vmin.clim<-as.numeric(gsub("_","-",argv$vmin.clim))
argv$vmin.clim<-argv$vmin.clim*(-1)**argv$vminsign.clim
argv$vmax.clim<-as.numeric(gsub("_","-",argv$vmax.clim))
argv$vmax.clim<-argv$vmax.clim*(-1)**argv$vmaxsign.clim
# buddy priorities
if (is.null(argv$prio.buddy)) argv$prio.buddy<-rep(-1,length=nfin)
if (any(is.na(argv$prio.buddy))) argv$prio.buddy<-rep(-1,length=nfin)
if (length(argv$prio.buddy)!=nfin) argv$prio.buddy<-rep(-1,length=nfin)
if (is.null(argv$prio.buddy_eve)) argv$prio.buddy_eve<-rep(-1,length=nfin)
if (any(is.na(argv$prio.buddy_eve))) argv$prio.buddy_eve<-rep(-1,length=nfin)
if (length(argv$prio.buddy_eve)!=nfin) argv$prio.buddy_eve<-rep(-1,length=nfin)
# SCT with smart boxes
if (argv$smartbox.sct & argv$variable == "T") {
  if (!file.exists(file.path(argv$titan_path,"sct","sct_smart_boxes.so"))) {
    print("ERROR: file not found")
    print(file.path(argv$titan_path,"sct","sct_smart_boxes.so"))
    quit(status=1)
  }
  dyn.load(file.path(argv$titan_path,"sct","sct_smart_boxes.so"))
}
# buddy_eve checks
if (any(is.na(argv$thr_eve.buddy_eve)))
  argv$thr_eve.buddy_eve<-c(0.1,1,10)
if (any(is.na(argv$dr.buddy_eve)))
  argv$dr.buddy_eve<-c(3000,3000,3000)
if (length(argv$dr.buddy_eve)!=length(argv$thr_eve.buddy_eve))
  argv$dr.buddy_eve<-c(3000,3000,3000)
if (any(is.na(argv$thr.buddy_eve)))
  argv$thr.buddy_eve<-c(0.05,0.05,1)
if (length(argv$thr.buddy_eve)!=length(argv$thr_eve.buddy_eve))
  argv$thr.buddy_eve<-c(0.05,0.05,1)
if (any(is.na(argv$n.buddy_eve)))
  argv$n.buddy_eve<-c(5,5,5)
if (length(argv$n.buddy_eve)!=length(argv$n_eve.buddy_eve))
  argv$n.buddy_eve<-c(5,5,5)
if (any(is.na(argv$dz.buddy_eve)))
  argv$dz.buddy_eve<-c(1500,1500,1500)
if (length(argv$dz.buddy_eve)!=length(argv$n_eve.buddy_eve))
  argv$dz.buddy_eve<-c(1500,1500,1500)
# wind-induced undercatch of precipitation, check consistency of inputs
if (argv$rr.wcor & argv$variable!="RR") {
  print(paste("ERROR: wind-induced correction for undercatch is implemented",
              "for precipitation only"))
  quit(status=1)
}
# proj4s
if (argv$dem | argv$dem.fill) {
  if (argv$proj4dem=="" & argv$dem.proj4_var=="" & argv$dem.proj4_att=="" ) {
    dem.xy_as_vars<-T
    proj4dem<-NULL
    proj4dem_from_nc<-NULL
  } else {
    dem.xy_as_vars<-F
    proj4dem<-argv$proj4dem
    proj4dem_from_nc<-list(var=argv$dem.proj4_var, att=argv$dem.proj4_att)
  }
}
# set the timestamp
if (!is.na(argv$timestamp)) {
  if (is.na(argv$fg.t)) argv$fg.t<-argv$timestamp
  if (is.na(argv$fge.t)) argv$fge.t<-argv$timestamp
  if (is.na(argv$wind.t)) argv$wind.t<-argv$timestamp
  if (is.na(argv$t2m.t)) argv$t2m.t<-argv$timestamp
}
#
#-----------------------------------------------------------------------------
if (argv$verbose | argv$debug) print(">> TITAN <<")
if (argv$debug) 
  capture.output(print(argv), file=file.path(argv$debug.dir,"argv.txt"))
#
#-----------------------------------------------------------------------------
# read data
first<-T
for (f in 1:nfin) {
  datain<-read.table(file=argv$input.files[f],
                     header=T,
                     sep=argv$separator[f],
                     stringsAsFactors=F,
                     strip.white=T)
  # if elev is not present that create a fake one
  varidxtmp<-match(argv$varname.elev[f],names(datain))
  if (is.na(varidxtmp)) {
    argv$varname.elev[f]<-"elev"
    if (argv$elev_not_used) {
      datain$elev<-rep(0,length=length(datain[[1]]))
    } else {
      datain$elev<-rep(NA,length=length(datain[[1]]))
    }
  }
  rm(varidxtmp)
  # varidx is used also in the output session
  varidxtmp<-match(c(argv$varname.lat[f],
                     argv$varname.lon[f],
                     argv$varname.elev[f],
                     argv$varname.value[f]),
                   names(datain))
  if (any(is.na(varidxtmp))) {
    print("ERROR in the specification of the variable names")
    print(paste("latitutde=",argv$varname.lat[f]))
    print(paste("longitude=",argv$varname.lon[f]))
    print(paste("elevation=",argv$varname.elev[f]))
    print(paste("    value=",argv$varname.value[f]))
    print("header of input file:")
    print(argv$input.files[f])
    print(names(datain))
    quit(status=1)
  }
  datatmp<-data.frame(datain[,varidxtmp])
  names(datatmp)<-c("lat","lon","elev","value")
  datatmp$lat<-suppressWarnings(as.numeric(datatmp$lat))
  datatmp$lon<-suppressWarnings(as.numeric(datatmp$lon))
  datatmp$elev<-suppressWarnings(as.numeric(datatmp$elev))
  auxz<-datatmp$elev
  datatmp$value<-suppressWarnings(
    argv$input.offset[f]+
    argv$input.cfact[f]*as.numeric(datatmp$value))
  ndatatmp<-length(datatmp$lat)
  if (ndatatmp==0) next
  # set provider id
  datatmp$prid<-rep(argv$prid[f],ndatatmp)
  aux<-rep(NA,length=ndatatmp)
  if (any(!is.na(argv$blacklist.idx)) & 
      any(argv$blacklist.fidx==argv$prid[f])) {
    aux[argv$blacklist.idx[which(argv$blacklist.fidx==argv$prid[f])]]<-argv$black.code  
  }
  if (any(!is.na(argv$blacklist.lat)) & 
      any(argv$blacklist.fll==argv$prid[f])) {
    out<-apply(cbind(argv$blacklist.lon[which(argv$blacklist.fll==argv$prid[f])],
                     argv$blacklist.lat[which(argv$blacklist.fll==argv$prid[f])])
               ,FUN=setCode_lonlat,MARGIN=1,code=argv$black.code)
    rm(out)
  }
  if (any(!is.na(argv$keeplist.idx)) & 
      any(argv$keeplist.fidx==argv$prid[f])) {
    aux[argv$keeplist.idx[which(argv$keeplist.fidx==argv$prid[f])]]<-argv$keep.code  
  }
  if (any(!is.na(argv$keeplist.lat)) & 
      any(argv$keeplist.fll==argv$prid[f])) {
    out<-apply(cbind(argv$keeplist.lon[which(argv$keeplist.fll==argv$prid[f])],
                     argv$keeplist.lat[which(argv$keeplist.fll==argv$prid[f])])
               ,FUN=setCode_lonlat,MARGIN=1,code=argv$keep.code)
    rm(out)
  }
  if (first) {
    varidx<-varidxtmp
    data<-datatmp
    first<-F
    z<-auxz
    dqcflag<-aux
    sctpog<-rep(NA,length=ndatatmp)
    corep<-rep(NA,length=ndatatmp)
    if (any(!is.na(argv$varname.opt))) {
      # varidx.opt is used in the output session
      varidx.opt<-match(argv$varname.opt,
                        names(datain))
      dataopt<-as.data.frame(array(data=NA,
                             dim=c(ndatatmp,length(argv$varname.opt))))
      names(dataopt)<-argv$varname.opt
      if (any(!is.na(varidx.opt)))
        dataopt<-datain[,varidx.opt[which(!is.na(varidx.opt))],drop=F]
    }
  } else {
    data<-rbind(data,datatmp)
    dqcflag<-c(dqcflag,aux)
    z<-c(z,auxz)
    sctpog<-c(sctpog,rep(NA,length=ndatatmp))
    corep<-c(corep,rep(NA,length=ndatatmp))
    if (any(!is.na(argv$varname.opt))) {
      varidx.opt.check<-match(argv$varname.opt,
                        names(datain))
      if (any(!is.na(varidx.opt.check) & is.na(varidx.opt))) {
        ixopt<-which(!is.na(varidx.opt.check) & is.na(varidx.opt))
        for (iopt in ixopt) {
          if (varidx.opt.check[iopt] %in% varidx.opt) {
            varidx.opt[iopt]<-max(varidx.opt,na.rm=T)+1
          } else { 
            varidx.opt[iopt]<-varidx.opt.check[iopt]
          }
        }
        rm(ixopt,iopt)
      }
      dataopttmp<-as.data.frame(array(data=NA,
                                dim=c(ndatatmp,length(argv$varname.opt))))
      names(dataopttmp)<-argv$varname.opt
      if (any(!is.na(varidx.opt.check)))
        dataopttmp<-datain[,varidx.opt.check[which(!is.na(varidx.opt.check))],
                           drop=F]
      dataopt<-rbind(dataopt,
                     dataopttmp)
      rm(dataopttmp)
    }
  }
  rm(varidxtmp)
}
rm(datatmp,datain,auxz,aux)
ndata<-length(data$lat)
if (ndata==0) {
  print("input file is empty")
  quit(status=0)
}
if (argv$verbose | argv$debug) {
  print(paste("number of observations=",ndata))
  if (any(!is.na(argv$blacklist.idx)) | any(!is.na(argv$blacklist.lat)))
    print(paste("number of blacklisted observations=",
          length(which(dqcflag==argv$black.code))) )
  if (any(!is.na(argv$keeplist.idx)) | any(!is.na(argv$keeplist.lat)))
    print(paste("number of keeplisted  observations=",
          length(which(dqcflag==argv$keep.code))) )
  if (nfin>1) {
    for (f in 1:nfin) { 
      print(paste("  number of observations provider",argv$prid[f],"=",
            length(which(data$prid==argv$prid[f]))))
      if (any(!is.na(argv$blacklist.idx)) | any(!is.na(argv$blacklist.lat)))
        print(paste("  number of blacklisted observations provider",
              argv$prid[f],"=",
              length(which(data$prid==argv$prid[f] & dqcflag==argv$black.code))) )
      if (any(!is.na(argv$keeplist.idx)) | any(!is.na(argv$keeplist.lat)))
        print(paste("  number of keeplisted  observations provider",
              argv$prid[f],"=",
              length(which(data$prid==argv$prid[f] & dqcflag==argv$keep.code))) )
    }
  }
  print("+---------------------------------+")
}
#
#-----------------------------------------------------------------------------
# test for no metadata (1st round) 
# NOTE: keep-listed stations could be flagged here
# Two rounds required because we need to extract info from the first guesses 
# that are in-line with the dem-filled elevation
dqcflag.bak<-dqcflag # .bak, if dem.fill to distinguish NAs / keep.code
ix<-which(is.na(dqcflag) | dqcflag==argv$keep.code)
if (length(ix)>0) {
  meta<-is.na(data$lat[ix]) | 
        is.na(data$lon[ix]) |
        is.na(z[ix]) | 
        z[ix]<argv$zmin | 
        z[ix]>argv$zmax |
        is.na(data$value[ix]) 
  if (argv$dqc_inbox_only) 
    meta<- meta | ( data$lat[ix] < argv$latmin | 
                    data$lat[ix] > argv$latmax |
                    data$lon[ix] < argv$lonmin | 
                    data$lon[ix] > argv$lonmax )
  if (any(meta)) dqcflag[ix[which(meta)]]<-argv$nometa.code
} else {
  print("no valid observations left, no metadata check")
}
if (argv$verbose | argv$debug) {
  flagaux<-dqcflag==argv$nometa.code & !is.na(dqcflag)
  print("test for no metdata, statistics over the whole dataset")
  print(paste("# observations lacking metadata and/or NAs=",
        length(which(flagaux))))
  print(paste("  # NAs                 =",
        length(which(flagaux & is.na(data$value))))) # coincides with all the NAs
  print(paste("  # lon-lat missing (*) =",
        length(which(flagaux & (is.na(data$lat) | is.na(data$lon) | 
                     data$lat < argv$latmin | data$lat > argv$latmax | 
                     data$lon < argv$lonmin | data$lon > argv$lonmax ) ))))
  print(paste("  # z missing           =",
        length(which(flagaux & is.na(z)))))
  print(paste("  # z out of range      =",
        length(which(flagaux & !is.na(z) & 
                     (z<argv$zmin | z>argv$zmax) ))))
  print("(*) or outside the specified box")
  rm(flagaux)
  print("+---------------------------------+")
}
rm(ix)
if (exists("meta")) rm(meta)
#
#-----------------------------------------------------------------------------
# coordinate transformation
if (argv$spatconv) {
  if (argv$debug) print("conversion of spatial coordinates")
  # initialization
  x<-data$lon
  y<-data$lon
  x[]<-NA
  y[]<-NA
  # do it
  coord<-SpatialPoints(cbind(data$lon,data$lat),
                       proj4string=CRS(argv$proj4_input_obsfiles))
  coord.new<-spTransform(coord,CRS(argv$proj4_where_dqc_is_done))
  xy.new<-coordinates(coord.new)
  x<-round(xy.new[,1],0)
  y<-round(xy.new[,2],0)
  xp<-expand.grid(c(argv$lonmin,argv$lonmax),c(argv$latmin,argv$latmax))
  coord<-SpatialPoints(xp,
                       proj4string=CRS(argv$proj4_input_obsfiles))
  coord.new<-spTransform(coord,CRS(argv$proj4_where_dqc_is_done))
  # define the extent for the SCT grid
  e<-extent(coord.new)
  xl<-e[1:2]
  yl<-e[3:4]
  rm(coord,coord.new,xy.new,xp)
} else {
  x<-data$lon
  y<-data$lat
  xl<-c(argv$lonmin,argv$lonmax)
  yl<-c(argv$latmin,argv$latmax)
  e<-extent(c(xl,yl))
}
if (argv$debug) {
  save.image(file.path(argv$debug.dir,"input_data.RData")) 
  print("+---------------------------------+")
}
#
#-----------------------------------------------------------------------------
# Read geographical information (optional) 
if (argv$dem | argv$dem.fill) {
  if (argv$verbose | argv$debug) print("read digital elevation model")

#+
get_data_from_ncfile<-function(nc.file,
                               nc.varname,
                               topdown,
                               var.dim,
                               proj4,
                               proj4_from_nc,
                               selection,
                               return_just_raster=T) {
#------------------------------------------------------------------------------
  ti<-nc4.getTime(argv$dem.file)
  raux<-try(nc4in(nc.file=nc.file,
                  nc.varname=nc.varname,
                  topdown=topdown,
                  out.dim=var.dim,
                  proj4=proj4,
                  nc.proj4=proj4_from_nc,
                  selection=selection))
  if (is.null(raux)) boom(paste("ERROR while reading file:",argv$nc.file))
  if (return_just_raster) {
    return(raux$stack)
  } else {
    return(raux)
  }
  if (dem.xy_as_vars) {
    raux<-try(nc4in(nc.file=argv$dem.file,
                    nc.varname=argv$dem.x_as_var.varname,
                    topdown=argv$dem.topdown,
                    out.dim=list(ndim=argv$dem.xy_as_var.ndim,
                                 tpos=argv$dem.xy_as_var.tpos,
                                 epos=NULL,
                                 names=argv$dem.xy_as_var.dimnames),
                    proj4=proj4dem,
                    nc.proj4=proj4dem_from_nc,
                    selection=list(t=ti[1],e=NULL)))
    if (is.null(raux)) boom(paste("ERROR while reading file:",argv$dem.file))
    rx<-raux$stack
    rm(raux)
    
  } else if (argv$proj4dem!=argv$proj4_where_dqc_is_done) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4_input_obsfiles))
    coord.new<-spTransform(coord,CRS(argv$proj4dem))
    xy.tmp<-coordinates(coord.new)
    zdem<-extract(rdem,xy.tmp)
    rm(coord,coord.new,xy.tmp)
  } else {
    zdem<-extract(rdem,cbind(x,y))
  }

}

  ti<-nc4.getTime(argv$dem.file)
  raux<-try(nc4in(nc.file=argv$dem.file,
                  nc.varname=argv$dem.varname,
                  topdown=argv$dem.topdown,
                  out.dim=list(ndim=argv$dem.ndim,
                               tpos=argv$dem.tpos,
                               epos=NULL,
                               names=argv$dem.dimnames),
                  proj4=proj4dem,
                  nc.proj4=proj4dem_from_nc,
                  selection=list(t=ti[1],e=NULL)))
  if (is.null(raux)) boom(paste("ERROR while reading file:",argv$dem.file))
  rdem<-raux$stack
  rm(raux)
  if (dem.xy_as_vars) {
    raux<-try(nc4in(nc.file=argv$dem.file,
                    nc.varname=argv$dem.x_as_var.varname,
                    topdown=argv$dem.topdown,
                    out.dim=list(ndim=argv$dem.xy_as_var.ndim,
                                 tpos=argv$dem.xy_as_var.tpos,
                                 epos=NULL,
                                 names=argv$dem.xy_as_var.dimnames),
                    proj4=proj4dem,
                    nc.proj4=proj4dem_from_nc,
                    selection=list(t=ti[1],e=NULL)))
    if (is.null(raux)) boom(paste("ERROR while reading file:",argv$dem.file))
    rx<-raux$stack
    rm(raux)
    
  } else if (argv$proj4dem!=argv$proj4_where_dqc_is_done) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4_input_obsfiles))
    coord.new<-spTransform(coord,CRS(argv$proj4dem))
    xy.tmp<-coordinates(coord.new)
    zdem<-extract(rdem,xy.tmp)
    rm(coord,coord.new,xy.tmp)
  } else {
    zdem<-extract(rdem,cbind(x,y))
  }
  # fill missing elevation with dem
  if (argv$dem.fill) {
    iz<-which( !is.na(data$value) &
               !is.na(data$lat)   & 
               !is.na(data$lon)   &
               !is.na(zdem) &
               !(zdem<argv$zmin | zdem>argv$zmax) &
               (is.na(z) | is.nan(z) | (z<argv$zmin | z>argv$zmax)) )
    z[iz]<-zdem[iz]
    dqcflag[iz]<-dqcflag.bak[iz]
    rm(dqcflag.bak)
    if (argv$verbose | argv$debug) {
      print(paste("# stations with elevation derived from DEM=",length(iz)))
      print("+---------------------------------+")
    }
    rm(iz)
  }  
  if (argv$debug) {
    png(file=file.path(argv$debug.dir,"dem.png"),width=800,height=800)
    image(rdem,
          breaks=c(-1000,seq(0,1500,length=10),seq(1501,2500,length=5),8000),
          col=c("azure",rev(gray.colors(15))))
    plotSCTgrid()
    dev.off()
  }
  rm(rdem)
  if (argv$debug) save.image(file.path(argv$debug.dir,"dem.RData")) 
} # END read DEM
if (argv$laf.sct) {
  if (argv$verbose | argv$debug)
    print("read land area fraction")
  ti<-nc4.getTime(argv$laf.file)
  raux<-try(nc4in(nc.file=argv$laf.file,
                  nc.varname=argv$laf.varname,
                  topdown=argv$laf.topdown,
                  out.dim=list(ndim=argv$laf.ndim,
                               tpos=argv$laf.tpos,
                               epos=NULL,
                               names=argv$laf.dimnames),
                  proj4=argv$proj4laf,
                  nc.proj4=list(var=argv$laf.proj4_var,
                                att=argv$laf.proj4_att),
                  selection=list(t=ti[1],e=NULL)))
  if (is.null(raux)) {
    print("ERROR while reading file:")
    print(argv$laf.file)
    quit(status=1)
  }
  rlaf<-raux$stack
  rm(raux,ti)
  if (argv$proj4laf!=argv$proj4_where_dqc_is_done) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4_input_obsfiles))
    coord.new<-spTransform(coord,CRS(argv$proj4laf))
    xy.tmp<-coordinates(coord.new)
    laf<-extract(rlaf,xy.tmp)/100.
    rm(coord,coord.new,xy.tmp)
  } else {
    laf<-extract(rlaf,cbind(x,y))/100.
  }
  if (any(is.na(laf))) laf[which(is.na(laf))]<-1
  if (argv$debug) {
    png(file=file.path(argv$debug.dir,"laf.png"),width=800,height=800)
    image(rlaf,breaks=c(seq(0,1,length=20)),col=c(rev(rainbow(19))))
    plotSCTgrid()
    dev.off()
  }
  if (!argv$debug) rm(rlaf)
} else {
  # use a fake laf
  laf<-rep(1,ndata)
} # END read LAF
if (argv$laf.sct & argv$debug) {
  save.image(file.path(argv$debug.dir,"input_data_laf.RData")) 
  print("+---------------------------------+")
}
#
#-----------------------------------------------------------------------------
# precipitation (in-situ) and temperature (field) cross-check (optional)
if (argv$ccrrt) {
  if (argv$debug | argv$verbose) 
    print(paste0("precipitation (in-situ) and temperature (field) cross-check (",argv$ccrrt.code,")"))
  # read temperature from gridded field
  ti<-nc4.getTime(argv$t2m.file)
  if (is.na(argv$t2m.t)) argv$t2m.t<-ti[1]
  if (!(argv$t2m.t %in% ti)) {
    print("ERROR timestamp requested is not in the file:")
    print(argv$t2m.t)
    print(ti)
    quit(status=1)
  }
  t2m.epos<-argv$t2m.epos
  if (is.na(argv$t2m.epos)) t2m.epos<-NULL
  t2m.e<-argv$t2m.e
  if (is.na(argv$t2m.e)) t2m.e<-NULL
  raux<-try(nc4in(nc.file=argv$t2m.file,
                  nc.varname=argv$t2m.varname,
                  topdown=argv$t2m.topdown,
                  out.dim=list(ndim=argv$t2m.ndim,
                               tpos=argv$t2m.tpos,
                               epos=t2m.epos,
                               names=argv$t2m.dimnames),
                  proj4=argv$proj4t2m,
                  nc.proj4=list(var=argv$t2m.proj4_var,
                                att=argv$t2m.proj4_att),
                  selection=list(t=argv$t2m.t,e=t2m.e)))
  rt2m<-raux$stack
  rm(raux)
  if (argv$debug)
    plot_debug(ff=file.path(argv$debug.dir,"rrwcor_t2m.png"),
               r=rt2m,x=x,y=y,proj4=argv$proj4_where_dqc_is_done,proj4plot=argv$proj4t2m)
  coord<-SpatialPoints(cbind(data$lon,data$lat),
                       proj4string=CRS(argv$proj4_input_obsfiles))
  coord.new<-spTransform(coord,CRS(argv$proj4t2m))
  xy.tmp<-coordinates(coord.new)
  t2m<-extract(rt2m,xy.tmp,method="bilinear")
  t2m<-argv$t2m.offset*(-1)**(argv$t2m.negoffset)+
       t2m*argv$t2m.cfact*(-1)**(argv$t2m.negcfact)
  rm(coord,coord.new,xy.tmp)
  # read elevation from gridded field
  if (!is.null(argv$t2m.demtpos)) {
    t2m.demtpos<-argv$t2m.demtpos
    ti<-nc4.getTime(argv$t2m.demfile)
    if (is.na(argv$t2m.demt)) argv$t2m.demt<-ti[1]
    if (!(argv$t2m.demt %in% ti)) {
      print("ERROR timestamp requested is not in the file:")
      print(argv$t2m.demt)
      print(ti)
      quit(status=1)
    }
  } else {
    t2m.demtpos<-argv$t2m.demtpos
  }
  if (!is.null(argv$t2m.demepos)) {
    t2m.demepos<-argv$t2m.demepos
    if (is.na(argv$t2m.demepos)) t2m.demepos<-NULL
  } else {
    t2m.demepos<-NULL
  }
  t2m.deme<-argv$t2m.deme
  if (is.na(argv$t2m.deme)) t2m.deme<-NULL
  raux<-try(nc4in(nc.file=argv$t2m.demfile,
                  nc.varname=argv$t2m.demvarname,
                  topdown=argv$t2m.demtopdown,
                  out.dim=list(ndim=argv$t2m.demndim,
                               tpos=t2m.demtpos,
                               epos=t2m.demepos,
                               names=argv$t2m.demdimnames),
                  proj4=argv$proj4t2m,
                  nc.proj4=list(var=argv$t2m.proj4_var,
                                att=argv$t2m.proj4_att),
                  selection=list(t=argv$t2m.demt,e=t2m.deme)))
  if (is.null(raux)) {
    print("ERROR while reading file:")
    print(argv$t2m.demfile)
    quit(status=1)
  }
  rt2mdem<-raux$stack
  rm(raux,ti)
  if (argv$proj4t2m!=argv$proj4_where_dqc_is_done) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4_input_obsfiles))
    coord.new<-spTransform(coord,CRS(argv$proj4t2m))
    xy.tmp<-coordinates(coord.new)
    zt2mdem<-extract(rt2mdem,xy.tmp,method="bilinear")
    rm(coord,coord.new,xy.tmp)
  } else {
    zt2mdem<-extract(rt2mdem,cbind(x,y),method="bilinear")
  }
  zt2mdem<-argv$t2m.demoffset*(-1)**(argv$t2m.demnegoffset)+
          zt2mdem*argv$t2m.demcfact*(-1)**(argv$t2m.demnegcfact)
  if (argv$debug)
    plot_debug(ff=file.path(argv$debug.dir,"rrwcor_dem.png"),
               r=rt2mdem,x=x,y=y,proj4=argv$proj4_where_dqc_is_done,proj4plot=argv$proj4t2m)
  t2m<-t2m+argv$gamma.standard*(z-zt2mdem)
  # cross-check
  for (f in 1:nfin) 
    dqcflag[which(data$prid==argv$prid[f] & 
                  t2m<argv$ccrrt.tmin[f] &
                  is.na(dqcflag))]<-argv$ccrrt.code
  #
  if (argv$verbose | argv$debug) {
    print("precipitaton and temperature  crosscheck")
    print(paste("temp thresholds =",toString(argv$ccrrt.tmin)))
    print(paste("# suspect observations=",
          length(which(dqcflag==argv$ccrrt.code & !is.na(dqcflag)))))
    print("+---------------------------------+")
  }
  if (argv$debug)
    save.image(file.path(argv$debug.dir,"dqcres_ccrrt.RData")) 
  rm(t2m,zt2mdem)
}
#
#-----------------------------------------------------------------------------
# Correction for the wind-undercatch of precipitation (optional)
if (argv$rr.wcor) {
  if (argv$debug | argv$verbose)
    print("Correction for the wind-undercatch of precipitation")
  # read temperature from gridded field
  ti<-nc4.getTime(argv$t2m.file)
  if (is.na(argv$t2m.t)) argv$t2m.t<-ti[1]
  if (!(argv$t2m.t %in% ti)) {
    print("ERROR timestamp requested is not in the file:")
    print(argv$t2m.t)
    print(ti)
    quit(status=1)
  }
  t2m.epos<-argv$t2m.epos
  if (is.na(argv$t2m.epos)) t2m.epos<-NULL
  t2m.e<-argv$t2m.e
  if (is.na(argv$t2m.e)) t2m.e<-NULL
  raux<-try(nc4in(nc.file=argv$t2m.file,
                  nc.varname=argv$t2m.varname,
                  topdown=argv$t2m.topdown,
                  out.dim=list(ndim=argv$t2m.ndim,
                               tpos=argv$t2m.tpos,
                               epos=t2m.epos,
                               names=argv$t2m.dimnames),
                  proj4=argv$proj4t2m,
                  nc.proj4=list(var=argv$t2m.proj4_var,
                                att=argv$t2m.proj4_att),
                  selection=list(t=argv$t2m.t,e=t2m.e)))
  rt2m<-raux$stack
  rm(raux)
  coord<-SpatialPoints(cbind(data$lon,data$lat),
                       proj4string=CRS(argv$proj4_input_obsfiles))
  coord.new<-spTransform(coord,CRS(argv$proj4t2m))
  xy.tmp<-coordinates(coord.new)
  t2m<-extract(rt2m,xy.tmp,method="bilinear")
  t2m<-argv$t2m.offset*(-1)**(argv$t2m.negoffset)+
       t2m*argv$t2m.cfact*(-1)**(argv$t2m.negcfact)
  if (argv$debug)
    plot_debug(ff=file.path(argv$debug.dir,"rrwcor_t2m.png"),
               r=rt2m,x=x,y=y,proj4=argv$proj4_where_dqc_is_done,proj4plot=argv$proj4t2m)
  rm(coord,coord.new,xy.tmp)
  # read elevation from gridded field
  if (!is.null(argv$t2m.demtpos)) {
    t2m.demtpos<-argv$t2m.demtpos
    ti<-nc4.getTime(argv$t2m.demfile)
    if (is.na(argv$t2m.demt)) argv$t2m.demt<-ti[1]
    if (!(argv$t2m.demt %in% ti)) {
      print("ERROR timestamp requested is not in the file:")
      print(argv$t2m.demt)
      print(ti)
      quit(status=1)
    }
  } else {
    t2m.demtpos<-argv$t2m.demtpos
  }
  if (!is.null(argv$t2m.demepos)) {
    t2m.demepos<-argv$t2m.demepos
    if (is.na(argv$t2m.demepos)) t2m.demepos<-NULL
  } else {
    t2m.demepos<-NULL
  }
  t2m.deme<-argv$t2m.deme
  if (is.na(argv$t2m.deme)) t2m.deme<-NULL
  raux<-try(nc4in(nc.file=argv$t2m.demfile,
                  nc.varname=argv$t2m.demvarname,
                  topdown=argv$t2m.demtopdown,
                  out.dim=list(ndim=argv$t2m.demndim,
                               tpos=t2m.demtpos,
                               epos=t2m.demepos,
                               names=argv$t2m.demdimnames),
                  proj4=argv$proj4t2m,
                  nc.proj4=list(var=argv$t2m.proj4_var,
                                att=argv$t2m.proj4_att),
                  selection=list(t=argv$t2m.demt,e=t2m.deme)))
  if (is.null(raux)) {
    print("ERROR while reading file:")
    print(argv$t2m.demfile)
    quit(status=1)
  }
  rt2mdem<-raux$stack
  rm(raux,ti)
  if (argv$proj4t2m!=argv$proj4_where_dqc_is_done) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4_input_obsfiles))
    coord.new<-spTransform(coord,CRS(argv$proj4t2m))
    xy.tmp<-coordinates(coord.new)
    zt2mdem<-extract(rt2mdem,xy.tmp,method="bilinear")
    rm(coord,coord.new,xy.tmp)
  } else {
    zt2mdem<-extract(rt2mdem,cbind(x,y),method="bilinear")
  }
  zt2mdem<-argv$t2m.demoffset*(-1)**(argv$t2m.demnegoffset)+
          zt2mdem*argv$t2m.demcfact*(-1)**(argv$t2m.demnegcfact)
  if (argv$debug)
    plot_debug(ff=file.path(argv$debug.dir,"rrwcor_dem.png"),
               r=rt2mdem,x=x,y=y,proj4=argv$proj4_where_dqc_is_done,proj4plot=argv$proj4t2m)
  t2m<-t2m+argv$gamma.standard*(z-zt2mdem)
  # read windspeed from gridded field
  #  case of windspeed in the file
  if (!is.na(argv$windspeed.varname)) {
    ti<-nc4.getTime(argv$wind.file)
    if (is.na(argv$wind.t)) argv$wind.t<-ti[1]
    if (!(argv$wind.t %in% ti)) {
      print("ERROR timestamp requested is not in the file:")
      print(argv$wind.t)
      print(ti)
      quit(status=1)
    }
    wind.epos<-argv$wind.epos
    if (is.na(argv$wind.epos)) wind.epos<-NULL
    wind.e<-argv$wind.e
    if (is.na(argv$wind.e)) wind.e<-NULL
    raux<-try(nc4in(nc.file=argv$wind.file,
                    nc.varname=argv$windspeed.varname,
                    topdown=argv$wind.topdown,
                    out.dim=list(ndim=argv$wind.ndim,
                                 tpos=argv$wind.tpos,
                                 epos=wind.epos,
                                 names=argv$wind.dimnames),
                    proj4=argv$proj4wind,
                    nc.proj4=list(var=argv$wind.proj4_var,
                                  att=argv$wind.proj4_att),
                    selection=list(t=argv$wind.t,e=wind.e)))
    if (is.null(raux)) {
      print("ERROR while reading file:")
      print(argv$wind.file)
      quit(status=1)
    }
    rwind<-raux$stack
    rm(raux,ti)
  #  case of u,v in the file
  } else if (!is.na(argv$u.varname) & !is.na(argv$v.varname)) {
    ti<-nc4.getTime(argv$wind.file)
    if (is.na(argv$wind.t)) argv$wind.t<-ti[1]
    if (!(argv$wind.t %in% ti)) {
      print("ERROR timestamp requested is not in the file:")
      print(argv$wind.t)
      print(ti)
      quit(status=1)
    }
    wind.epos<-argv$wind.epos
    if (is.na(argv$wind.epos)) wind.epos<-NULL
    wind.e<-argv$wind.e
    if (is.na(argv$wind.e)) wind.e<-NULL
    raux<-try(nc4in(nc.file=argv$wind.file,
                    nc.varname=argv$u.varname,
                    topdown=argv$wind.topdown,
                    out.dim=list(ndim=argv$wind.ndim,
                                 tpos=argv$wind.tpos,
                                 epos=wind.epos,
                                 names=argv$wind.dimnames),
                    proj4=argv$proj4wind,
                    nc.proj4=list(var=argv$wind.proj4_var,
                                  att=argv$wind.proj4_att),
                    selection=list(t=argv$wind.t,e=wind.e)))
    if (is.null(raux)) {
      print("ERROR while reading file:")
      print(argv$wind.file)
      quit(status=1)
    }
    ru<-raux$stack
    rm(raux,ti)
    if (argv$debug)
      plot_debug(ff=file.path(argv$debug.dir,"rrwcor_u.png"),
                 r=ru,x=x,y=y,proj4=argv$proj4_where_dqc_is_done,proj4plot=argv$proj4wind)
    raux<-try(nc4in(nc.file=argv$wind.file,
                    nc.varname=argv$v.varname,
                    topdown=argv$wind.topdown,
                    out.dim=list(ndim=argv$wind.ndim,
                                 tpos=argv$wind.tpos,
                                 epos=wind.epos,
                                 names=argv$wind.dimnames),
                    proj4=argv$proj4wind,
                    nc.proj4=list(var=argv$wind.proj4_var,
                                  att=argv$wind.proj4_att),
                    selection=list(t=argv$wind.t,e=wind.e)))
    if (is.null(raux)) {
      print("ERROR while reading file:")
      print(argv$wind.file)
      quit(status=1)
    }
    rv<-raux$stack
    rm(raux)
    if (argv$debug)
      plot_debug(ff=file.path(argv$debug.dir,"rrwcor_v.png"),
                 r=rv,x=x,y=y,proj4=argv$proj4_where_dqc_is_done,proj4plot=argv$proj4wind)
    rwind<-rv
    rwind[]<-sqrt(getValues(ru)**2+getValues(rv)**2)
    rm(ru,rv)
  #  case of wrong/no wind varname in the file
  } else {
    print(paste("ERROR precipitation correction for the wind undercatch:",
                " wind varname has not been specified"))
    quit(status=1)
  }
  if (argv$proj4wind!=argv$proj4_where_dqc_is_done) {
    ws10m<-extract(rwind, 
           coordinates( spTransform( SpatialPoints(cbind(data$lon,data$lat), 
                                           proj4string=CRS(argv$proj4_input_obsfiles)),
                                    CRS(argv$proj4wind)) ), method="bilinear")
  } else {
    ws10m<-extract(rwind,cbind(x,y),method="bilinear")
  }
  if (argv$debug)
    plot_debug(ff=file.path(argv$debug.dir,"rrwcor_ws.png"),
               r=rwind,x=x,y=y,proj4=argv$proj4_where_dqc_is_done,proj4plot=argv$proj4wind)
  rm(rwind) 
  # precipitation data adjustment
  data$rawvalue<-data$value
  res<-wolff_correction(par=argv$rr.wcor.par,
                        t2m=t2m,
                        ws10m=ws10m,
                        rr=data$rawvalue)
  data$value<-res$rr.cor
  data$vsigma<-res$sigma
  rm(res)
  if (argv$debug | argv$verbose) {
    print(paste0("# observations (ok-so-far) = ",
          length(which(!is.na(data$value) & is.na(dqcflag)))))
    print(paste0("# observations (>=0 & ok-so-far) = ",
          length(which(data$value>=0 & !is.na(data$value) & is.na(dqcflag)))))
    print(paste0("# not NAs-observations set to NAs after this correction = ",
          length(which( is.na(data$value) & !is.na(data$rawvalue) ))))
    ix<-which(data$value>=0 & !is.na(data$value) & is.na(dqcflag))
    if (length(ix)>0) {
      lm<-lm(data$value[ix]~data$rawvalue[ix]+0)
      print(paste0("linear regression, obs_adjusted = ",
      round(as.numeric(lm$coefficients),3)," * obs_raw"))
    }
    print("+---------------------------------+")
  }
  if (argv$debug) {
    save.image(file.path(argv$debug.dir,"rrcor.RData")) 
  }
}
#
#-----------------------------------------------------------------------------
# Read deterministic first guess (optional)
if (!is.na(argv$fg.file)) {
  if (argv$verbose | argv$debug)
    print("Read deterministic first-guess")
  ti<-nc4.getTime(argv$fg.file)
  if (is.na(argv$fg.t)) argv$fg.t<-ti[1]
  if (!(argv$fg.t %in% ti)) {
    print("ERROR timestamp requested is not in the file:")
    print(argv$fg.t)
    print(ti)
    quit(status=1)
  }
  fg.epos<-argv$fg.epos
  if (is.na(argv$fg.epos)) fg.epos<-NULL
  fg.e<-argv$fg.e
  if (is.na(argv$fg.e)) fg.e<-NULL
  if (argv$fg.acc) {
    tminus1h<-format(as.POSIXlt(
      seq(as.POSIXlt(strptime(argv$fg.t,"%Y%m%d%H%M",tz="UTC")),
          length=2,by="-1 hour"),"UTC")[2],"%Y%m%d%H%M",tz="UTC")
    raux<-try(nc4in(nc.file=argv$fg.file,
                    nc.varname=argv$fg.varname,
                    topdown=argv$fg.topdown,
                    out.dim=list(ndim=argv$fg.ndim,
                                 tpos=argv$fg.tpos,
                                 epos=fg.epos,
                                 names=argv$fg.dimnames),
                    proj4=argv$proj4fg,
                    nc.proj4=list(var=argv$fg.proj4_var,
                                  att=argv$fg.proj4_att),
                    selection=list(t=c(tminus1h,argv$fg.t),e=fg.e)))
    if (is.null(raux)) {
      print("ERROR while reading file:")
      print(argv$fg.file)
      quit(status=1)
    }
    rfg<-raster(raux$stack,"layer.2")-raster(raux$stack,"layer.1")
  } else {
    raux<-try(nc4in(nc.file=argv$fg.file,
                    nc.varname=argv$fg.varname,
                    topdown=argv$fg.topdown,
                    out.dim=list(ndim=argv$fg.ndim,
                                 tpos=argv$fg.tpos,
                                 epos=fg.epos,
                                 names=argv$fg.dimnames),
                    proj4=argv$proj4fg,
                    nc.proj4=list(var=argv$fg.proj4_var,
                                  att=argv$fg.proj4_att),
                    selection=list(t=argv$fg.t,e=fg.e)))
    if (is.null(raux)) {
      print("ERROR while reading file (var):")
      print(argv$fg.file)
      quit(status=1)
    }
    rfg<-raux$stack
  }
  rm(raux,ti,fg.e,fg.epos)
  # radar fg, data quality control 
  if (!is.na(argv$fg.type)) {
    if (argv$fg.type=="radar" & argv$fg.dodqc) {
      if (argv$verbose | argv$debug)
        print("Read radar and do the radar-DQC")
      t0a<-Sys.time()
      suppressPackageStartupMessages(library("igraph"))
      dfg<-getValues(rfg)
      # a. remove not plausible values
      if (argv$verbose | argv$debug) print(" remove not plausible values")
      radardqc.min<-0
      radardqc.max<-300
      ix<-which( !is.na(dfg) & (dfg<radardqc.min | dfg>radardqc.max) )
      if (length(ix)>0) {
        dfg[ix]<-NA
        rfg[]<-dfg
      }
      # b. remove patches of connected cells that are too small
      #  check for small and isolated clumps (patches) of connected cells with 
      #  precipitation greater than a predefined threshold
      #   threshold 0 mm/h. remove all the clumps made of less than 100 cells
      #   threshold 1 mm/h. remove all the clumps made of less than 50 cells
      if (argv$verbose | argv$debug) print(" remove small clumps")
      radardqc.clump.thr<-c(0,1)
      radardqc.clump.n<-c(100,50)
      for (i in 1:length(radardqc.clump.thr)) {
        raux<-rfg
        if (any(dfg<=radardqc.clump.thr[i])) 
          raux[which(dfg<=radardqc.clump.thr[i])]<-NA
        rclump<-clump(raux)
        fr<-freq(rclump)
        ix<-which(!is.na(fr[,2]) & fr[,2]<=radardqc.clump.n[i])
        dfg[getValues(rclump) %in% fr[ix,1]]<-NA
        rfg[]<-dfg
        rm(raux,fr,ix,rclump)
      }
      # c. remove outliers. Check for outliers in square boxes of 51km by 51km
#      if (argv$verbose | argv$debug) print(" remove outliers")
      raux<-rfg
      daux<-boxcox(x=dfg,lambda=0.5)
      raux[]<-daux
#      radardqc.outl.fact<-ceiling(51000/xres(rfg))
#      if ((radardqc.outl.fact%%2)==0) radardqc.outl.fact<-radardqc.outl.fact+1
#      if (argv$verbose | argv$debug) print(" remove outliers -mean")
      # aggregate over boxes of fact x fact cells, take the mean
#      avg<-getValues(
#            crop(
#             focal(
#              extend(raux,c(radardqc.outl.fact,radardqc.outl.fact)),
#              w=matrix(1,radardqc.outl.fact,radardqc.outl.fact),fun=mean,na.rm=T),
#            raux) )
#      if (argv$verbose | argv$debug) print(" remove outliers -sd")
#      stdev<-getValues(
#              crop(
#               focal(
#                extend(raux,c(radardqc.outl.fact,radardqc.outl.fact)),
#                w=matrix(1,radardqc.outl.fact,radardqc.outl.fact),fun=sd,na.rm=T),
#              raux) )
#      ix<-which(stdev>0)
      # compute mean and sd
      raux_agg<-aggregate(raux,fact=5,fun=mean,na.rm=T)
      daux_agg<-getValues(raux_agg)
      ix_aux<-which(!is.na(daux_agg))
      xyaux<-xyFromCell(raux_agg,ix_aux)
      xrad_aux<-xyaux[,1]
      yrad_aux<-xyaux[,2]
      vrad_aux<-daux_agg[ix_aux]
      get_rad_stat<-function(i,dh_ref=25000) { 
        deltax<-abs(xrad_aux[i]-xrad_aux)
        deltay<-abs(yrad_aux[i]-yrad_aux)
        ix<-which( deltax<dh_ref & deltay<dh_ref )
        dist<-deltax; dist[]<-NA
        dist[ix]<-sqrt(deltax[ix]*deltax[ix]+deltay[ix]*deltay[ix])
        ix<-which( !is.na(dist) & dist<dh_ref )
        return(c(mean(vrad_aux[ix]),sd(vrad_aux[ix])))
      }
      if (!is.na(argv$cores)) {
        arr<-mcmapply(get_rad_stat,
                      1:length(ix_aux),
                      mc.cores=argv$cores,
                      SIMPLIFY=T)
      # no-multicores
      } else {
        arr<-mapply(get_rad_stat,
                    1:length(ix_aux),
                    SIMPLIFY=T)
      }
      raux_agg[]<-NA; raux_agg[ix_aux]<-arr[1,]
      raux<-disaggregate(raux_agg,fact=5)
      if (ncell(raux)>ncell(rfg)) {
        raux<-crop(raux,rfg)
      } else if (ncell(raux)<ncell(rfg)) {
        raux<-extend(raux,rfg)
      }
      avg<-getValues(raux)
      raux_agg[]<-NA; raux_agg[ix_aux]<-arr[2,]
      raux<-disaggregate(raux_agg,fact=5,method="bilinear",na.rm=T)
      if (ncell(raux)>ncell(rfg)) {
        raux<-crop(raux,rfg)
      } else if (ncell(raux)<ncell(rfg)) {
        raux<-extend(raux,rfg)
      }
      stdev<-getValues(raux)
      ix<-which(stdev>0 & !is.na(daux) & !is.na(avg) & !is.na(stdev))
      rm(arr,raux_agg,ix_aux,xrad_aux,yrad_aux,vrad_aux,daux_agg,xyaux)
      # outliers are defined as in Lanzante,1997: abs(value-mean)/st.dev > 5
      suspect<-which((abs(daux[ix]-avg[ix])/stdev[ix])>5) 
      if (length(suspect)>0) dfg[ix[suspect]]<-NA
      rfg[]<-dfg
      rm(raux,daux,avg,stdev,ix,suspect,dfg)
      if (argv$radarout) rrad<-rfg
      t1a<-Sys.time()
      if (argv$verbose | argv$debug) {
        print(paste(" remove outliers - time",round(t1a-t0a,1),
                                              attr(t1a-t0a,"unit")))
      }
    }
  }
  if (argv$proj4fg!=argv$proj4_where_dqc_is_done) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4_input_obsfiles))
    coord.new<-spTransform(coord,CRS(argv$proj4fg))
    xy.tmp<-coordinates(coord.new)
    fg<-extract(rfg,xy.tmp,method="bilinear")
    rm(coord,coord.new,xy.tmp)
  } else {
    fg<-extract(rfg,cbind(x,y),method="bilinear")
  }
  fg<-argv$fg.offset*(-1)**(argv$fg.negoffset)+
      fg*argv$fg.cfact*(-1)**(argv$fg.negcfact)
  if (argv$cool & argv$usefg.cool) {
    aux<-argv$fg.offset*(-1)**(argv$fg.negoffset)+
         getValues(rfg)*argv$fg.cfact*(-1)**(argv$fg.negcfact)
    ix_aux<-which(!is.na(aux)& !is.nan(aux) & is.finite(aux) & is.numeric(aux))
    if (!exists("xobs_cool_aux")) xobs_cool_aux<-integer(0)
    if (!exists("yobs_cool_aux")) yobs_cool_aux<-integer(0)
    if (!exists("pridobs_cool_aux")) pridobs_cool_aux<-integer(0)
    if (!exists("yo_cool_aux")) yo_cool_aux<-integer(0)
    if (length(ix_aux)>0) {
      yo_cool_aux<-c(yo_cool_aux,aux[ix_aux])
      xygrid_cool<-xyFromCell(rfg,ix_aux)
      xobs_cool_aux<-c(xobs_cool_aux,xygrid_cool[,1])
      yobs_cool_aux<-c(yobs_cool_aux,xygrid_cool[,2])
      pridobs_cool_aux<-c(pridobs_cool_aux,rep(NA,length(ix_aux)))
      rm(xygrid_cool)
    }
    rm(aux,ix_aux)
  }
  if (!argv$debug) rm(rfg)
  # temperature: adjust for elevation differences
  if (argv$variable=="T") {
    ti<-nc4.getTime(argv$fg.demfile)
    if (is.na(argv$fg.demt)) argv$fg.demt<-ti[1]
    if (!(argv$fg.demt %in% ti)) {
      print("ERROR timestamp requested is not in the file:")
      print(argv$fg.demt)
      print(ti)
      quit(status=1)
    }
    fg.demepos<-argv$fg.demepos
    if (is.na(argv$fg.demepos)) fg.demepos<-NULL
    fg.deme<-argv$fg.deme
    if (is.na(argv$fg.deme)) fg.deme<-NULL
    raux<-try(nc4in(nc.file=argv$fg.demfile,
                    nc.varname=argv$fg.demvarname,
                    topdown=argv$fg.demtopdown,
                    out.dim=list(ndim=argv$fg.demndim,
                                 tpos=argv$fg.demtpos,
                                 epos=fg.demepos,
                                 names=argv$fg.demdimnames),
                    proj4=argv$proj4fg,
                    nc.proj4=list(var=argv$fg.proj4_var,
                                  att=argv$fg.proj4_att),
                    selection=list(t=argv$fg.demt,e=fg.deme)))
    if (is.null(raux)) {
      print("ERROR while reading file:")
      print(argv$fg.demfile)
      quit(status=1)
    }
    rfgdem<-raux$stack
    rm(raux,ti)
    if (argv$proj4fg!=argv$proj4_where_dqc_is_done) {
      coord<-SpatialPoints(cbind(data$lon,data$lat),
                           proj4string=CRS(argv$proj4_input_obsfiles))
      coord.new<-spTransform(coord,CRS(argv$proj4fg))
      xy.tmp<-coordinates(coord.new)
      zfgdem<-extract(rfgdem,xy.tmp,method="bilinear")
      rm(coord,coord.new,xy.tmp)
    } else {
      zfgdem<-extract(rfgdem,cbind(x,y),method="bilinear")
    }
    zfgdem<-argv$fg.demoffset*(-1)**(argv$fg.demnegoffset)+
            zfgdem*argv$fg.demcfact*(-1)**(argv$fg.demnegcfact)
    if (!argv$debug) rm(rfgdem)
    fg<-fg+argv$gamma.standard*(z-zfgdem)
    if (!argv$debug) rm(zfgdem)
  }
  if (argv$debug) {
    png(file=file.path(argv$debug.dir,"fg.png"),width=800,height=800)
    image(rfg,
          breaks=seq(range(getValues(rfg),na.rm=T)[1],
                     range(getValues(rfg),na.rm=T)[2],length=20),
          col=c(rev(rainbow(19))))
    if (exists("rlaf")) contour(rlaf,levels=c(0,1),add=T)
    xy.tmp<-as.data.frame(cbind(x,y))
    coordinates(xy.tmp)<-c("x","y")
    proj4string(xy.tmp)<-CRS(argv$proj4_where_dqc_is_done)
    xy.tmp.fg<-spTransform(xy.tmp,CRS(argv$proj4fg))
    points(xy.tmp.fg,cex=0.8,pch=19)
    dev.off()
    rm(rfg,xy.tmp,xy.tmp.fg)
    if (argv$variable=="T") {
      png(file=file.path(argv$debug.dir,"fgdem.png"),width=800,height=800)
      image(rfgdem,
            breaks=seq(range(getValues(rfgdem),na.rm=T)[1],
                       range(getValues(rfgdem),na.rm=T)[2],length=20),
            col=c(rev(rainbow(19))))
      if (exists("rlaf")) contour(rlaf,levels=c(0,1),add=T)
      points(x,y,cex=0.8,pch=19)
      dev.off()
    }
    save.image(file.path(argv$debug.dir,"input_data_fg.RData")) 
    if (exists("rfgdem")) rm(rfgdem)
  }
  if (argv$verbose | argv$debug)
    print("+---------------------------------+")
}
#
#-----------------------------------------------------------------------------
# Read first guess ensemble (optional)
if (!is.na(argv$fge.file)) {
  if (argv$verbose | argv$debug)
    print("Read ensemble first-guess")
  t0a<-Sys.time()
  ti<-nc4.getTime(argv$fge.file)
  if (is.na(argv$fge.t)) argv$fge.t<-ti[1]
  if (!(argv$fge.t %in% ti)) {
    print("ERROR timestamp requested is not in the file:")
    print(argv$fge.t)
    print(ti)
    quit(status=1)
  }
  # temperature: adjust for elevation differences
  if (argv$variable=="T") {
    ti<-nc4.getTime(argv$fge.demfile)
    if (is.na(argv$fge.demt)) argv$fge.demt<-ti[1]
    if (!(argv$fge.demt %in% ti)) {
      print("ERROR timestamp requested is not in the file:")
      print(argv$fge.demt)
      print(ti)
      quit(status=1)
    }
    fge.demepos<-argv$fge.demepos
    if (is.na(argv$fge.demepos)) fge.demepos<-NULL
    fge.deme<-argv$fge.deme
    if (is.na(argv$fge.deme)) fge.deme<-NULL
    raux<-nc4in(nc.file=argv$fge.demfile,
                nc.varname=argv$fge.demvarname,
                topdown=argv$fge.demtopdown,
                out.dim=list(ndim=argv$fge.demndim,
                             tpos=argv$fge.demtpos,
                             epos=fge.demepos,
                             names=argv$fge.demdimnames),
                proj4=argv$proj4fge,
                nc.proj4=list(var=argv$fge.proj4_var,
                              att=argv$fge.proj4_att),
                selection=list(t=argv$fge.demt,e=fge.deme))
    rfgedem<-raux$stack
    if (argv$debug) {
      png(file=file.path(argv$debug.dir,"fgedem.png"),width=800,height=800)
      image(rfgedem,
            breaks=seq(range(getValues(rfgedem),na.rm=T)[1],
                       range(getValues(rfgedem),na.rm=T)[2],length=20),
            col=c(rev(rainbow(19))))
      if (exists("rlaf")) contour(rlaf,levels=c(0,1),add=T)
      points(x,y,cex=0.8,pch=19)
      dev.off()
    }
    rm(raux,ti)
    if (argv$proj4fge!=argv$proj4_where_dqc_is_done) {
      coord<-SpatialPoints(cbind(data$lon,data$lat),
                           proj4string=CRS(argv$proj4_input_obsfiles))
      coord.new<-spTransform(coord,CRS(argv$proj4fge))
      xy.tmp<-coordinates(coord.new)
      zfgedem<-extract(rfgedem,xy.tmp,method="bilinear")
      rm(coord,coord.new,xy.tmp)
    } else {
      zfgedem<-extract(rfgedem,cbind(x,y),method="bilinear")
    }
    zfgedem<-argv$fg.demoffset*(-1)**(argv$fg.demnegoffset)+
             zfgedem*argv$fg.demcfact*(-1)**(argv$fg.demnegcfact)
    rm(rfgedem)
  }
  #
  if (is.na(argv$fge.epos)) {
    print("ERROR fge.epos must have a valid value")
    quit(status=1)
  }
  ei<-nc4.getDim(argv$fge.file,
                 varid=argv$fge.dimnames[argv$fge.epos])
  edata<-array(data=NA,dim=c(ndata,length(ei)))
  first<-T
  for (ens in 1:length(ei)) {
    if (argv$fge.acc) {
      tminus1h<-format(as.POSIXlt(
        seq(as.POSIXlt(strptime(argv$fge.t,"%Y%m%d%H%M",tz="UTC")),
            length=2,by="-1 hour"),"UTC")[2],"%Y%m%d%H%M",tz="UTC")
      raux<-try(nc4in(nc.file=argv$fge.file,
                      nc.varname=argv$fge.varname,
                      topdown=argv$fge.topdown,
                      out.dim=list(ndim=argv$fge.ndim,
                                   tpos=argv$fge.tpos,
                                   epos=argv$fge.epos,
                                   names=argv$fge.dimnames),
                      proj4=argv$proj4fge,
                      nc.proj4=list(var=argv$fge.proj4_var,
                                    att=argv$fge.proj4_att),
                      selection=list(t=c(tminus1h,argv$fge.t),
                                     e=ei[ens])))
      if (is.null(raux)) {
        print("ERROR while reading file:")
        print(argv$fge.file)
        quit(status=1)
      }
      rfge<-raster(raux$stack,"layer.2")-
            raster(raux$stack,"layer.1")
    } else {
      raux<-nc4in(nc.file=argv$fge.file,
                  nc.varname=argv$fge.varname,
                  topdown=argv$fge.topdown,
                  out.dim=list(ndim=argv$fge.ndim,
                               tpos=argv$fge.tpos,
                               epos=argv$fge.epos,
                               names=argv$fge.dimnames),
                  proj4=argv$proj4fge,
                  nc.proj4=list(var=argv$fge.proj4_var,
                                att=argv$fge.proj4_att),
#                  selection=list(t=argv$fge.t,e=ei[ens],z=4))
                  selection=list(t=argv$fge.t,e=ei[ens]))
      rfge<-raux$stack
    }
    # debug
    if (argv$debug) {
      png(file=file.path(argv$debug.dir,
                         paste("fge_",formatC(ens,width=2,flag="0"),".png",sep="")),
          width=800,height=800)
      image(rfge,
            breaks=seq(range(getValues(rfge),na.rm=T)[1],
                       range(getValues(rfge),na.rm=T)[2],length=20),
            col=c(rev(rainbow(19))))
      if (exists("rlaf")) contour(rlaf,levels=1,add=T)
      if (exists("rlaf")) rm(rlaf)
      points(x,y,cex=0.8,pch=19)
      dev.off()
    }
    rm(raux)
    # extract data at observation locations (coord conversion)
    if (argv$proj4fge!=argv$proj4_where_dqc_is_done) {
      coord<-SpatialPoints(cbind(data$lon,data$lat),
                           proj4string=CRS(argv$proj4_input_obsfiles))
      coord.new<-spTransform(coord,CRS(argv$proj4fge))
      xy.tmp<-coordinates(coord.new)
      edata[,ens]<-extract(rfge,xy.tmp,method="bilinear")
      rm(coord,coord.new,xy.tmp)
    } else {
      edata[,ens]<-extract(rfge,cbind(x,y),method="bilinear")
    }
    edata[,ens]<-argv$fge.offset*(-1)**(argv$fge.negoffset)+
               edata[,ens]*argv$fge.cfact*(-1)**(argv$fge.negcfact)
    if (argv$variable=="T") { 
      edata[,ens]<-edata[,ens]+argv$gamma.standard*(z-zfgedem)
    }
    rm(rfge)
  } # end of cycle over ensemble members
  # compute mean and sd 
  fge.mu<-rowMeans(edata,na.rm=T)
  fge_sd<-function(i) { sd(edata[i,],na.rm=T) }
  if (!is.na(argv$cores)) {
    fge.sd<-mcmapply(fge_sd,
                     1:dim(edata)[1],
                     mc.cores=argv$cores,
                     SIMPLIFY=T)
  # no-multicores
  } else {
    fge.sd<-mapply(fge_sd,
                   1:dim(edata)[1],
                   SIMPLIFY=T)
  }
  fge.sd<-pmax(fge.sd,argv$sdmin.fge,na.rm=T)
  if (argv$debug)
    save.image(file.path(argv$debug.dir,"input_data_fge.RData")) 
  rm(edata)
  # debug
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    print(paste("fge sd(5,25,50,75,95-th)=",
          toString(round(as.vector(quantile(fge.sd,
                             probs=c(0.05,0.25,0.5,0.75,0.95),
                             na.rm=T)),2))))
    print(paste("time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    print("+---------------------------------+")
  }
}
#
#-----------------------------------------------------------------------------
# test for no metadata (final) 
# use only (probably) good observations
# NOTE: keep-listed stations could be flagged here
ix<-which(is.na(dqcflag) | dqcflag==argv$keep.code)
if (length(ix)>0) {
  meta<-is.na(data$lat[ix]) | 
        is.na(data$lon[ix]) |
        data$lat[ix] < argv$latmin | 
        data$lat[ix] > argv$latmax | 
        data$lon[ix] < argv$lonmin | 
        data$lon[ix] > argv$lonmax |
        is.na(z[ix]) | 
        z[ix]<argv$zmin | 
        z[ix]>argv$zmax |
        is.na(data$value[ix]) 
  if (any(meta)) dqcflag[ix[which(meta)]]<-argv$nometa.code
} else {
  print("no valid observations left, no metadata check")
}
if (argv$verbose | argv$debug) {
  flagaux<-dqcflag==argv$nometa.code & !is.na(dqcflag)
  print("test for no metdata (2nd round), statistics over the whole dataset")
  print(paste("# observations lacking metadata and/or NAs=",
        length(which(flagaux))))
  print(paste("  # NAs                 =",
        length(which(flagaux & is.na(data$value))))) # coincides with all the NAs
  print(paste("  # lon-lat missing (*) =",
        length(which(flagaux & (is.na(data$lat) | is.na(data$lon) | 
                     data$lat < argv$latmin | data$lat > argv$latmax | 
                     data$lon < argv$lonmin | data$lon > argv$lonmax ) ))))
  print(paste("  # z missing           =",
        length(which(flagaux & is.na(z)))))
  print(paste("  # z out of range      =",
        length(which(flagaux & !is.na(z) & 
                     (z<argv$zmin | z>argv$zmax) ))))
  print("(*) or outside the specified box")
  rm(flagaux)
  print("+---------------------------------+")
}
rm(ix)
if (exists("meta")) rm(meta)
if (argv$debug) 
  save.image(file.path(argv$debug.dir,"dqcres_meta.RData")) 
#
#-----------------------------------------------------------------------------
# check elevation against dem 
# NOTE: keep-listed stations canNOT be flagged here
if (argv$dem) {
  if (argv$verbose | argv$debug) 
    print(paste0("check station elevations against digital elevation model (",argv$dem.code,")"))
  # set doit vector
  doit<-vector(length=ndata,mode="numeric")
  doit[]<-NA
  for (f in 1:nfin) doit[data$prid==argv$prid[f]]<-argv$doit.dem[f]
  # use only (probably) good observations
  ix<-which(is.na(dqcflag))
  if (length(ix)>0) {
    ixna<-which(!is.na(z) & !is.na(zdem) & is.na(dqcflag))
    sus<-which( abs(z[ixna]-zdem[ixna])>argv$dz.dem &
                doit[ixna]==1 )
    # set dqcflag
    if (length(sus)>0) dqcflag[ixna[sus]]<-argv$dem.code
  }  else {
    print("no valid observations left, no dem check")
  }
  if (argv$verbose | argv$debug) {
    print(paste("#stations with elevations too different from digital elevation model =",
                length(which(dqcflag==argv$dem.code))))
    print("+---------------------------------+")
  }
  rm(doit)
}
if (argv$debug) 
  save.image(file.path(argv$debug.dir,"dqcres_demcheck.RData")) 
#
#-----------------------------------------------------------------------------
# plausibility test
# NOTE: keep-listed stations could be flagged here
ix<-which( (is.na(dqcflag) | dqcflag==argv$keep.code) &
           (data$value<argv$vmin | data$value>argv$vmax))
if (length(ix)>0) dqcflag[ix]<-argv$p.code
if (argv$verbose | argv$debug) {
  print(paste0("plausibility test (",argv$p.code,")"))
  print(paste("min/max thresholds =",argv$vmin,argv$vmax))
  print(paste("# <min=",length(which(dqcflag==argv$p.code & 
                                     data$value<argv$vmin))))
  print(paste("# >max=",length(which(dqcflag==argv$p.code & 
                                     data$value>argv$vmax))))
  print(paste("# suspect observations=",length(which(dqcflag==argv$p.code & 
                                                     !is.na(dqcflag)))))
  print("+---------------------------------+")
}
if (argv$debug) 
  save.image(file.path(argv$debug.dir,"dqcres_plausibility.RData")) 
#
#-----------------------------------------------------------------------------
# climatological check 
# NOTE: keep-listed stations canNOT be flagged here
# use only (probably) good observations
if (!is.na(argv$month.clim)) {
  # set doit vector
  doit<-vector(length=ndata,mode="numeric")
  doit[]<-NA
  for (f in 1:nfin)
    doit[data$prid==argv$prid[f]]<-argv$doit.clim[f]
  # apply the test on all the observations except blacklist/keeplist 
  ix<-which(is.na(dqcflag))
  if (length(ix)>0) {
    # flag only observations that are suspect and have doit==1
    sus<-which( (data$value[ix]<argv$vmin.clim[argv$month.clim] | 
                 data$value[ix]>argv$vmax.clim[argv$month.clim]) &
                 doit[ix]==1)
    # set dqcflag
    if (length(sus)>0) dqcflag[ix[sus]]<-argv$clim.code
  } else {
    print("no valid observations left, no climatological check")
  }
  if (argv$verbose | argv$debug) {
    print(paste("climatological test (month=",argv$month.clim,")",sep=""))
    print(paste("# suspect observations=",length(which(dqcflag==argv$clim.code))))
    print("+---------------------------------+")
  }
  rm(doit)
  if (argv$debug) 
    save.image(file.path(argv$debug.dir,"dqcres_monthclim.RData")) 
}
#
#-----------------------------------------------------------------------------
# buddy check (event-based)
#  Define an event compare each observation against the average of neighbouring observations 
# NOTE: keep-listed stations are used but they canNOT be flagged here
if (argv$buddy_eve) {
  nsus<-vector(mode="numeric",length=length(argv$thr_eve.buddy_eve))
  # set doit/prio vectors
  doit<-vector(length=ndata,mode="numeric"); doit[]<-NA
  prio<-vector(length=ndata,mode="numeric"); prio[]<-NA
  for (f in 1:nfin) {
    aux<-which(data$prid==argv$prid[f])
    if (length(aux)==0) next
    doit[aux]<-argv$doit.buddy_eve[f]
    prio[aux]<-argv$prio.buddy_eve[f]
  }
  rm(aux)
  # test
  print(paste0("buddy_eve-check (",argv$buddy_eve.code,")"))
  print(paste0("priorities ",toString(argv$prio.buddy_eve)))
  for (i in 1:argv$i.buddy_eve) {
    priority<-ifelse((i==1 & any(prio!=(-1))),T,F)
    nsus[]<-0
    for (j in 1:length(argv$thr_eve.buddy_eve)) {
      # use only (probably) good observations with doit!=0
      ix<-which( (is.na(dqcflag) | dqcflag==argv$keep.code) &
                 doit!=0 )
      t0a<-Sys.time()
      if (length(ix)>0) {
        # define global 1D vector used in statSpat (1D for fast access)
        itot<-1:length(ix)
        xtot<-x[ix]
        ytot<-y[ix]
        ztot<-as.numeric(z[ix])
        priotot<-as.numeric(prio[ix])
        ttot<-data$value[ix]
        if (!is.na(argv$cores)) {
          stSp_buddy_eve<-mcmapply(statSpat_mapply,
                                   1:length(itot),
                                   mc.cores=argv$cores,
                                   SIMPLIFY=T,
                                   drmin=argv$dr.buddy_eve[j],
                                   priority=priority,
                                   statistics="buddy_event", #buddy_event
                                   event_threshold=argv$thr_eve.buddy_eve[j],
                                   event_def="lt")
        # no-multicores
        } else {
          stSp_buddy_eve<-mapply(statSpat_mapply,
                                 1:length(itot),
                                 SIMPLIFY=T,
                                 drmin=argv$dr.buddy,
                                 priority=priority,
                                 statistics="buddy_event", #buddy_event
                                 event_threshold=argv$thr_eve.buddy_eve[j],
                                 event_def="lt")
        }
        n.buddy_eve<-ifelse(priority,0,argv$n.buddy_eve)
        # suspect if:
        if (argv$thr.buddy_eve[j]<1) { 
          sus<-which( 
            (stSp_buddy_eve[1,]>n.buddy_eve[j] & 
             stSp_buddy_eve[2,]<argv$dz.buddy_eve[j] &
             is.na(dqcflag[ix])) &
             doit[ix]==1 & 
            ( ( stSp_buddy_eve[3,]==0 & 
               ((1-stSp_buddy_eve[4,])<=argv$thr.buddy_eve[j])) |
              ( stSp_buddy_eve[3,]==1 & 
               (    stSp_buddy_eve[4,]<=argv$thr.buddy_eve[j])) ))
        } else if (argv$thr.buddy_eve[j]>=1) {
          nyes<-round(stSp_buddy_eve[1,]*stSp_buddy_eve[4,],0)
          nno<-stSp_buddy_eve[1,]-nyes
          sus<-which( 
            (stSp_buddy_eve[1,]>n.buddy_eve[j] & 
            stSp_buddy_eve[2,]<argv$dz.buddy_eve[j] &
            is.na(dqcflag[ix])) &
            doit[ix]==1 & 
            ( (stSp_buddy_eve[3,]==0 & nno<argv$thr.buddy_eve[j]) |
            (stSp_buddy_eve[3,]==1 & nyes<argv$thr.buddy_eve[j]) ))
          rm(nyes,nno)
        } else {
          sus<-integer(0)
        }
        # set dqcflag
        if (length(sus)>0) dqcflag[ix[sus]]<-argv$buddy_eve.code
      } else {
        print("no valid observations left, no buddy_eve check")
      }
      nsus[j]<-ifelse(exists("sus"),length(sus),0)
      if (argv$verbose | argv$debug) {
        t1a<-Sys.time()
        str<-" (#TOT "
        for (f in 1:nfin) {
          if (f>1) str<-paste0(str,"; ") 
          aux<-which(dqcflag==argv$buddy_eve.code & data$prid==argv$prid[f])
          str<-paste0(str,"prid",argv$prid[f],"=",length(aux))
          rm(aux)
        }
        str<-paste0(str,")")
        prestr<-ifelse(j==1,"","  ")
        print(paste0(prestr,"buddy_eve-check, iteration=",i,
                     "/event def: value less than ",argv$thr_eve.buddy_eve[j],
                     "/dqc param: thres=",argv$thr.buddy_eve[j],
                     "rad=",argv$dr.buddy_eve[j],
                     "/time ",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
        print(paste0(prestr,nsus[j]," new suspect observations",str))
        rm(prestr,str)
      }
    } # end for j
    if (sum(nsus)==0) break
  }  # end for i
  rm(doit)
  if (argv$debug) 
    save.image(file.path(argv$debug.dir,"dqcres_buddy_eve.RData")) 
  if (argv$verbose | argv$debug) 
    print("+---------------------------------+")
  if (exists("stSp_buddy_eve")) rm(stSp_buddy_eve)
  if (exists("sus")) rm(sus)
  if (exists("xtot")) rm(xtot)
  if (exists("ytot")) rm(ytot)
  if (exists("ztot")) rm(ztot)
  if (exists("itot")) rm(itot)
  if (exists("ixyztp_tot")) rm(ixyztp_tot)
}
#
#-----------------------------------------------------------------------------
# buddy check (standard)
#  compare each observation against the average of neighbouring observations 
# NOTE: keep-listed stations are used but they canNOT be flagged here
nprev<-0
# set doit/prio vectors
doit<-vector(length=ndata,mode="numeric"); doit[]<-NA
prio<-vector(length=ndata,mode="numeric"); prio[]<-NA
for (f in 1:nfin) {
  aux<-which(data$prid==argv$prid[f])
  if (length(aux)==0) next
  doit[aux]<-argv$doit.buddy[f]
  prio[aux]<-argv$prio.buddy[f]
}
rm(aux)
# test
print(paste0("buddy-check (",argv$buddy.code,")"))
print(paste0("priorities ",toString(argv$prio.buddy)))
for (i in 1:argv$i.buddy) {
  priority<-ifelse((i==1 & any(prio!=(-1))),T,F)
  # use only (probably) good observations with doit!=0
  ix<-which( (is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0 )
  t0a<-Sys.time()
  if (length(ix)>0) {
    # define global 1D vector used in statSpat (1D for fast access)
    itot<-1:length(ix)
    xtot<-x[ix]
    ytot<-y[ix]
    ztot<-as.numeric(z[ix])
    priotot<-as.numeric(prio[ix])
    if (argv$variable=="RR") {
      ttot<-boxcox(x=data$value[ix],lambda=argv$boxcox.lambda)
    } else {
      ttot<-as.numeric(data$value[ix])
    }
    if (!is.na(argv$cores)) {
      stSp_buddy<-mcmapply(statSpat_mapply,
                           1:length(itot),
                           mc.cores=argv$cores,
                           SIMPLIFY=T,
                           drmin=argv$dr.buddy,
                           priority=priority)
    # no-multicores
    } else {
      stSp_buddy<-mapply(statSpat_mapply,
                         1:length(itot),
                         SIMPLIFY=T,
                         drmin=argv$dr.buddy,
                         priority=priority)
    }
    # probability of gross error
    stSp_buddy[4,]<-pmax(argv$sdmin.buddy,stSp_buddy[4,])
    stSp_buddy[4,which(stSp_buddy[1,]==1)]<-argv$sdmin.buddy
    pog<-abs(ttot-stSp_buddy[3,])/stSp_buddy[4,]
    n.buddy<-ifelse(priority,0,argv$n.buddy) 
    # suspect if: 
    sus<-which( (pog>argv$thr.buddy & 
                 stSp_buddy[1,]>n.buddy & 
                 stSp_buddy[2,]<argv$dz.buddy &
                 is.na(dqcflag[ix])) &
                 doit[ix]==1 )
    # set dqcflag
    if (length(sus)>0) dqcflag[ix[sus]]<-argv$buddy.code
  } else {
    print("no valid observations left, no buddy check")
  }
  ncur<-length(which(dqcflag==argv$buddy.code))
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    str<-"(#TOT "
    for (f in 1:nfin) {
      if (f>1) str<-paste0(str,"; ") 
      aux<-which(dqcflag==argv$buddy.code & data$prid==argv$prid[f])
      str<-paste0(str,"prid",argv$prid[f],"=",length(aux))
      rm(aux)
    }
    str<-paste0(str,")")
    print(paste("buddy-check, iteration=",i,
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    print(paste("#new suspect observations=",ncur-nprev,str))
    rm(str)
  }
#  if (argv$debug) {
#    auxok<-which(is.na(dqcflag))
#    aux<-which(dqcflag==argv$buddy.code & !is.na(dqcflag))
#    aux9<-which(data$prid==9 & is.na(dqcflag))
#    for (ii in aux) {
#      png(file=paste0(argv$debug.dir,"/buddy_",i,"_",ii,".png"),width=800,height=800)
#      plot(1,1,xlim=c(x[ii]-argv$dr.buddy,x[ii]+argv$dr.buddy),
#           ylim=c(y[ii]-argv$dr.buddy,y[ii]+argv$dr.buddy))
#      text(x[auxok],y[auxok],data$value[auxok])
#      text(x[aux9],y[aux9],data$value[aux9],col="gray")
#      text(x[aux],y[aux],data$value[aux],col="red")
#      dev.off()
#    }
#    rm(aux)
#  }
  if ((ncur-nprev)==0 & i>2) break
  nprev<-ncur
}
rm(doit,prio)
if (argv$debug) 
  save.image(file.path(argv$debug.dir,"dqcres_buddy.RData")) 
if (argv$verbose | argv$debug) 
  print("+---------------------------------+")
if (exists("stSp_buddy")) rm(stSp_buddy)
if (exists("sus")) rm(sus)
if (exists("pog")) rm(pog)
if (exists("xtot")) rm(xtot)
if (exists("ytot")) rm(ytot)
if (exists("ztot")) rm(ztot)
if (exists("itot")) rm(itot)
if (exists("ixyztp_tot")) rm(ixyztp_tot)
#
#-----------------------------------------------------------------------------
# check against a first-guess (deterministic)
if (argv$fg) {
  if (argv$verbose | argv$debug) {
    print(paste0("first-guess check det (",argv$fg.code,")"))
  }
  # set doit vector
  doit<-vector(length=ndata,mode="numeric"); doit[]<-NA
  thrvec<-vector(length=ndata,mode="numeric"); thrvec[]<-NA
  thrposvec<-vector(length=ndata,mode="numeric"); thrposvec[]<-NA
  thrnegvec<-vector(length=ndata,mode="numeric"); thrnegvec[]<-NA
  perc_minvalvec<-vector(length=ndata,mode="numeric"); perc_minvalvec[]<-NA
  thrpercvec<-vector(length=ndata,mode="numeric"); thrpercvec[]<-NA
  thrpospercvec<-vector(length=ndata,mode="numeric"); thrpospercvec[]<-NA
  thrnegpercvec<-vector(length=ndata,mode="numeric"); thrnegpercvec[]<-NA
  for (f in 1:nfin) {
    if (!any(data$prid==argv$prid[f])) next
    aux<-which(data$prid==argv$prid[f])
    doit[aux]<-argv$doit.fg[f]
    thrvec[aux]<-argv$thr.fg[f]
    thrposvec[aux]<-argv$thrpos.fg[f]
    thrnegvec[aux]<-argv$thrneg.fg[f]
    perc_minvalvec[aux]<-argv$perc.fg_minval[f]
    thrpercvec[aux]<-argv$thrperc.fg[f]
    thrpospercvec[aux]<-argv$thrposperc.fg[f]
    thrnegpercvec[aux]<-argv$thrnegperc.fg[f]
    rm(aux)
  }
  # use only (probably) good observations
  ix<-which(is.na(dqcflag) & doit!=0)
  if (length(ix)>0) {
    dev<-data$value-fg
    devperc<-dev/fg
    flag_sus<-rep(F,ndata)
    flag_to_check<-is.na(dqcflag) & doit==1 &
                   !is.na(data$value) & 
                   !is.nan(data$value) & 
                   is.finite(data$value) &
                   !is.na(fg) & !is.nan(fg) & is.finite(fg)
    if (any(!is.na(thrvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrvec) & flag_to_check & abs(dev)>thrvec)
    if (any(!is.na(thrposvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrposvec) & flag_to_check & dev>thrposvec)
    if (any(!is.na(thrnegvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrnegvec) & flag_to_check & dev<0 & abs(dev)>thrnegvec)
    flag_to_check<-flag_to_check & fg>=perc_minvalvec
    if (any(!is.na(thrpercvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrpercvec) & flag_to_check & abs(devperc)>thrpercvec)
    if (any(!is.na(thrpospercvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrpospercvec) & flag_to_check & 
        dev>0 & abs(devperc)>thrpospercvec)
    if (any(!is.na(thrnegpercvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrnegpercvec) & flag_to_check & 
        dev<0 & abs(devperc)>thrnegpercvec)
    ix_sus<-which(flag_sus)
    rm(flag_sus,flag_to_check,dev,devperc)
    rm(doit,thrvec,thrposvec,thrnegvec,perc_minvalvec,thrpercvec)
    rm(thrpospercvec,thrnegpercvec)
    # set dqcflag
    if (length(ix_sus)>0) dqcflag[ix_sus]<-argv$fg.code
    rm(ix_sus)
  }  else {
    print("no valid observations left, no first-guess check")
  }
  if (argv$verbose | argv$debug) {
    print(paste("# observations that fail the first-guess check (det)=",
                length(which(dqcflag==argv$fg.code))))
    print("+---------------------------------+")
  }
  if (argv$debug) 
    save.image(file.path(argv$debug.dir,"dqcres_fg.RData")) 
  if (exists("doit")) rm(doit)
  if (exists("thrvec")) rm(thrvec)
  if (exists("thrposvec")) rm(thrposvec)
  if (exists("thrnegvec")) rm(thrnegvec)  
  if (exists("perc_minvalvec")) rm(perc_minvalvec)
  if (exists("thrpercvec")) rm(thrpercvec)
  if (exists("thrpospercvec")) rm(thrpospercvec)
  if (exists("thrnegpercvec")) rm(thrnegpercvec)
}
#
#-----------------------------------------------------------------------------
# check against a first-guess (ensemble)
if (argv$fge) {
  if (argv$verbose | argv$debug) {
    print(paste0("first-guess check ens (",argv$fge.code,")"))
  }
  # set doit vector
  doit<-vector(length=ndata,mode="numeric"); doit[]<-NA
  thrvec<-vector(length=ndata,mode="numeric"); thrvec[]<-NA
  thrposvec<-vector(length=ndata,mode="numeric"); thrposvec[]<-NA
  thrnegvec<-vector(length=ndata,mode="numeric"); thrnegvec[]<-NA
  thrpercvec<-vector(length=ndata,mode="numeric"); thrpercvec[]<-NA
  thrpospercvec<-vector(length=ndata,mode="numeric"); thrpospercvec[]<-NA
  thrnegpercvec<-vector(length=ndata,mode="numeric"); thrnegpercvec[]<-NA
  throutvec<-vector(length=ndata,mode="numeric"); throutvec[]<-NA
  thrposoutvec<-vector(length=ndata,mode="numeric"); thrposoutvec[]<-NA
  thrnegoutvec<-vector(length=ndata,mode="numeric"); thrnegoutvec[]<-NA
  for (f in 1:nfin) {
    if (!any(data$prid==argv$prid[f])) next
    aux<-which(data$prid==argv$prid[f])
    doit[aux]<-argv$doit.fge[f]
    thrvec[aux]<-argv$thr.fge[f]
    thrposvec[aux]<-argv$thrpos.fge[f]
    thrnegvec[aux]<-argv$thrneg.fge[f]
    perc_minvalvec[aux]<-argv$perc.fge_minval[f]
    thrpercvec[aux]<-argv$thrperc.fge[f]
    thrpospercvec[aux]<-argv$thrposperc.fge[f]
    thrnegpercvec[aux]<-argv$thrnegperc.fge[f]
    throutvec[aux]<-argv$throut.fge[f]
    thrposoutvec[aux]<-argv$thrposout.fge[f]
    thrnegoutvec[aux]<-argv$thrnegout.fge[f]
    rm(aux)
  }
  # use only (probably) good observations
  ix<-which(is.na(dqcflag) & doit!=0)
  if (length(ix)>0) {
    dev<-data$value-fge.mu
    devperc<-dev/fge.mu
    devout<-dev/fge.sd
    flag_sus<-rep(F,ndata)
    flag_to_check<-is.na(dqcflag) & doit==1 &
                   !is.na(data$value) & 
                   !is.nan(data$value) & 
                   is.finite(data$value) &
                   !is.na(fge.mu) & !is.nan(fge.mu) & is.finite(fge.mu)
    if (any(!is.na(thrvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrvec) & flag_to_check & abs(dev)>thrvec)
    if (any(!is.na(thrposvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrposvec) & flag_to_check & dev>thrposvec)
    if (any(!is.na(thrnegvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrnegvec) & flag_to_check & dev<0 & abs(dev)>thrnegvec)
    flag_to_check<-flag_to_check & fge.mu>=perc_minvalvec
    if (any(!is.na(thrpercvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrpercvec) & flag_to_check & abs(devperc)>thrpercvec)
    if (any(!is.na(thrpospercvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrpospercvec) & flag_to_check & 
        dev>0 & abs(devperc)>thrpospercvec)
    if (any(!is.na(thrnegpercvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrnegpercvec) & flag_to_check & 
        dev<0 & abs(devperc)>thrnegpercvec)
    flag_to_check<-flag_to_check & 
                   !is.na(fge.sd) & !is.nan(fge.sd) & is.finite(fge.sd) 
    if (any(!is.na(throutvec))) 
      flag_sus<-flag_sus | 
       (!is.na(throutvec) & flag_to_check & abs(devout)>throutvec)
    if (any(!is.na(thrposoutvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrposoutvec) & flag_to_check & 
        dev>0 & abs(devout)>thrposoutvec)
    if (any(!is.na(thrnegoutvec))) 
      flag_sus<-flag_sus | 
       (!is.na(thrnegoutvec) & flag_to_check & 
        dev<0 & abs(devout)>thrnegoutvec)
    ix_sus<-which(flag_sus)
    rm(flag_sus,flag_to_check,dev,devperc,devout)
    rm(doit,thrvec,thrposvec,thrnegvec,perc_minvalvec,thrpercvec)
    rm(thrpospercvec,thrnegpercvec,throutvec,thrposoutvec,thrnegoutvec)
    # set dqcflag
    if (length(ix_sus)>0) dqcflag[ix_sus]<-argv$fge.code
    rm(ix_sus)
  }  else {
    print("no valid observations left, no first-guess check")
  }
  if (argv$verbose | argv$debug) {
    print(paste("# observations that fail the first-guess check (ens)=",
                length(which(dqcflag==argv$fge.code))))
    print("+---------------------------------+")
  }
  rm(doit)
  if (argv$debug) 
    save.image(file.path(argv$debug.dir,"dqcres_fge.RData")) 
}
#
#-----------------------------------------------------------------------------
# SCT - Spatial Consistency Test
# NOTE: keep-listed stations are used but they canNOT be flagged here
if (argv$verbose | argv$debug) 
  print(paste0("SCT (",argv$sct.code,")"))
nprev<-0
# set doit vector
doit<-vector(length=ndata,mode="numeric"); doit[]<-NA
eps2.sct<-vector(length=ndata,mode="numeric"); eps2.sct[]<-NA
T2vec<-vector(length=ndata,mode="numeric"); T2vec[]<-NA
T2posvec<-vector(length=ndata,mode="numeric"); T2posvec[]<-NA
T2negvec<-vector(length=ndata,mode="numeric"); T2negvec[]<-NA
for (f in 1:nfin) {
  if (!any(data$prid==argv$prid[f])) next
  aux<-which(data$prid==argv$prid[f])
  doit[aux]<-argv$doit.sct[f]
  eps2.sct[aux]<-argv$eps2.sct[f]
  T2vec[aux]<-argv$thr.sct[f]
  T2posvec[aux]<-argv$thrpos.sct[f]
  T2negvec[aux]<-argv$thrneg.sct[f]
  rm(aux)
}
# test
for (i in 1:argv$i.sct) {
  # use only (probably) good observations with doit!=0
  ix<-which( (is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0 )
  if (length(ix)>0) {
    t0a<-Sys.time()
    # SCT station-by-station
    if (argv$stn_by_stn.sct) {
      if (i>1) break
      sct_loop<-T
      cont_sct_loop<-0
      obs_to_check<-rep(T,length(dqcflag))
      obs_to_check[!is.na(doit) & doit!=1]<-F
      obs_to_use<-rep(T,length(dqcflag))
      if (argv$transf.sct) {
        yo_sct<-boxcox(x=data$value,lambda=argv$boxcox.lambda)
      } else {
        yo_sct<-data$value
      }
      fg_min<-ifelse(argv$transf.sct,
                     boxcox(x=argv$vmin,lambda=argv$boxcox.lambda),
                     argv$vmin)
      fg_max<-ifelse(argv$transf.sct,
                     boxcox(x=argv$vmax,lambda=argv$boxcox.lambda),
                     argv$vmax)
      if (argv$usefge.sct) {
        if (argv$transf.sct) {
          b_sct<-boxcox(x=fge.mu,lambda=argv$boxcox.lambda)
        } else {
          b_sct<-fge.mu
        }
        argv$fglab.sct<-NA
        obs_to_use<-!is.na(b_sct) & !is.nan(b_sct) & is.finite(b_sct)
      } else if (argv$usefg.sct) {
        if (argv$transf.sct) {
          b_sct<-boxcox(x=fg,lambda=argv$boxcox.lambda)
        } else {
          b_sct<-fg
        }
        argv$fglab.sct<-NA
        obs_to_use<-!is.na(b_sct) & !is.nan(b_sct) & is.finite(b_sct)
      }
      while (sct_loop) {
        cont_sct_loop<-cont_sct_loop+1
        t00a<-Sys.time()
        if (cont_sct_loop>(length(ix)/2)) break
        ixg<-which( (is.na(dqcflag) | dqcflag==argv$keep.code) & 
                    obs_to_check & obs_to_use )
        nobs_to_check<-length(ixg)
        xgrid_spint<-x[ixg]
        ygrid_spint<-y[ixg]
        zgrid_spint<-z[ixg]
        lafgrid_spint<-laf[ixg]
        yo_to_check<-yo_sct[ixg]
        if (argv$usefge.sct | argv$usefg.sct) xb_spint<-b_sct[ixg]
        ixo<-which( (is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0 &
                    obs_to_use )
        xobs_spint<-x[ixo]
        yobs_spint<-y[ixo]
        zobs_spint<-z[ixo]
        lafobs_spint<-laf[ixo]
        eps2_spint<-eps2.sct[ixo]
        nobs<-length(ixo)
        yo_spint<-yo_sct[ixo]
        if (argv$usefge.sct | argv$usefg.sct) yb_spint<-b_sct[ixo]
        # CV-analysis and variances at station points
        if (!is.na(argv$cores)) {
          arr<-t(mcmapply(oi_var_gridpoint_by_gridpoint,
                          1:nobs_to_check,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          dh=argv$DhorMin.sct,
                          box_o_nearest_halfwidth=argv$box_o_nearest_halfwidth.sct,
                          dz=argv$Dver.sct,
                          lafmin=argv$lafmin.sct,
                          dh_adaptive=T,
                          corr=argv$corr.sct,
                          pmax=argv$pmax.sct,
                          fg=argv$fglab.sct,
                          fg_gamma=argv$fg_gamma.sct,
                          fg_min=fg_min,
                          fg_max=fg_max,
                          succ_corr=argv$succ_corr.sct,
                          y_elab=F,
                          loocv=T,
                          o_errvar_min=argv$o_errvar_min.sct,
                          o_errvar_max=argv$o_errvar_max.sct,
                          xa_errvar_max=argv$xa_errvar_max.sct,
                          xa_errvar_min=argv$xa_errvar_min))
        # no-multicores
        } else {
          arr<-t(mapply(oi_var_gridpoint_by_gridpoint,
                        1:nobs_to_check,
                        SIMPLIFY=T,
                        dh=argv$DhorMin.sct,
                        box_o_nearest_halfwidth=argv$box_o_nearest_halfwidth.sct,
                        dz=argv$Dver.sct,
                        lafmin=argv$lafmin.sct,
                        dh_adaptive=T,
                        corr=argv$corr.sct,
                        pmax=argv$pmax.sct,
                        fg=argv$fglab.sct,
                        fg_gamma=argv$fg_gamma.sct,
                        fg_min=fg_min,
                        fg_max=fg_max,
                        succ_corr=argv$succ_corr.sct,
                        y_elab=F,
                        loocv=T,
                        o_errvar_min=argv$o_errvar_min.sct,
                        o_errvar_max=argv$o_errvar_max.sct,
                        xa_errvar_max=argv$xa_errvar_max.sct,
                        xa_errvar_min=argv$xa_errvar_min))
        }
        yav<-arr[,1]
        yav_errvar<-arr[,2]
        yo_errvar<-arr[,3]
        dh_ref<-arr[,7]
        rm(arr)
        pog<-(yav-yo_to_check)**2/(yav_errvar+yo_errvar)
        flag_sus<-rep(F,nobs_to_check)
        if (any(!is.na(T2posvec[ixg]))) 
          flag_sus<-flag_sus | 
                    (!is.na(T2posvec[ixg]) & !is.na(pog) & 
                     pog>T2posvec[ixg] & (yo_to_check-yav)>0)
        if (any(!is.na(T2negvec[ixg]))) 
          flag_sus<-flag_sus | 
                    (!is.na(T2negvec[ixg]) & !is.na(pog) & 
                     pog>T2negvec[ixg] & (yo_to_check-yav)<0)
        if (any(!is.na(T2vec[ixg]))) 
          flag_sus<-flag_sus | 
                    (!is.na(T2vec[ixg]) & !is.na(pog) & 
                     pog>T2vec[ixg])
        ix_sus<-which(flag_sus)
        rm(flag_sus)
        t01a<-Sys.time()
        print(paste(" stn-by-stn iteration=",cont_sct_loop,
                    "/#obs to check=",length(ixg),
                    "/time",round(t01a-t00a,1),attr(t01a-t00a,"unit")))
        if (length(ix_sus)==0) {
          sct_loop<-F
        } else {
          find_the_largeErr<-function(i,dist_max){
            dist<-sqrt((xgrid_largeErr[i]-xgrid_largeErr)**2+
                       (ygrid_largeErr[i]-ygrid_largeErr)**2)
            ix_near<-which(dist<dist_max)
            ifelse(any(pog_largeErr[ix_near]>pog_largeErr[i]),F,T)
          }
          xgrid_largeErr<-xgrid_spint[ix_sus]
          ygrid_largeErr<-ygrid_spint[ix_sus]
          pog_largeErr<-pog[ix_sus]
          largeErr<-mapply(find_the_largeErr,
                           1:length(pog_largeErr),
                           SIMPLIFY=T,
                           dist_max=dh_ref[ix_sus])
          dqcflag[ixg[ix_sus[which(largeErr)]]]<-argv$sct.code
          if (!any(largeErr)) sct_loop<-F 
          if (any(!is.na(T2vec[ixg]))) { 
            ix_not_to_check<-which(pog<(T2vec[ixg]/4) & 
                                   !is.na(T2vec[ixg]))
            if (length(ix_not_to_check)>0) 
              obs_to_check[ixg[ix_not_to_check]]<-F
            ix_superBad<-which(pog>(4*T2vec[ixg]) & 
                               !is.na(T2vec[ixg]))
            if (length(ix_superBad)>0) 
              dqcflag[ixg[ix_superBad]]<-argv$sct.code
            rm(ix_not_to_check,ix_superBad)
          }
          if (any(!is.na(T2posvec[ixg]))) {
            ix_not_to_check<-which(pog<(T2posvec[ixg]/4) & 
                                   !is.na(T2posvec[ixg]) & 
                                   (yo_to_check-yav)>0)
            if (length(ix_not_to_check)>0) 
              obs_to_check[ixg[ix_not_to_check]]<-F
            ix_superBad<-which(pog>(4*T2posvec[ixg]) & 
                               !is.na(T2posvec[ixg]) &
                               (yo_to_check-yav)>0)
            if (length(ix_superBad)>0) 
              dqcflag[ixg[ix_superBad]]<-argv$sct.code
            rm(ix_not_to_check,ix_superBad)
          } 
          if (any(!is.na(T2negvec[ixg]))) {
            ix_not_to_check<-which(pog<(T2negvec[ixg]/4) & 
                                   !is.na(T2negvec[ixg]) & 
                                   (yo_to_check-yav)<0)
            if (length(ix_not_to_check)>0) 
              obs_to_check[ixg[ix_not_to_check]]<-F
            ix_superBad<-which(pog>(4*T2negvec[ixg]) & 
                               !is.na(T2negvec[ixg]) &
                               (yo_to_check-yav)<0)
            if (length(ix_superBad)>0) 
              dqcflag[ixg[ix_superBad]]<-argv$sct.code
            rm(ix_not_to_check,ix_superBad)
          }
          if (!exists("yav_prev")) {
            yav_prev<-dqcflag; yav_prev[]<-NA
            yav_prev[ixg]<-yav 
          } else {
            if (any(yav_prev[ixg]==yav)) 
              obs_to_check[ixg[which(yav_prev[ixg]==yav)]]<-F
            yav_prev[ixg]<-yav 
          }
        }
      } # end sct_loop
      if (cont_sct_loop>(length(ix)/2)) {
        print("Warning: SCT loop stopped. Too many iterations. Better check this out.")
      } 
      # compute coefficients of representativeness (if needed)
      if (any(is.na(argv$const.corep))) {
        t00a<-Sys.time()
        ixg<-which(is.na(dqcflag))
        nobs_to_check<-length(ixg)
        xgrid_spint<-x[ixg]
        ygrid_spint<-y[ixg]
        zgrid_spint<-z[ixg]
        lafgrid_spint<-laf[ixg]
        yo_to_check<-yo_sct[ixg]
        if (argv$usefge.sct | argv$usefg.sct) xb_spint<-b_sct[ixg]
        ixo<-ixg
        xobs_spint<-x[ixo]
        yobs_spint<-y[ixo]
        zobs_spint<-z[ixo]
        lafobs_spint<-laf[ixo]
        eps2_spint<-eps2.sct[ixo]
        nobs<-length(ixo)
        yo_spint<-yo_sct[ixo]
        if (argv$usefge.sct | argv$usefg.sct) yb_spint<-b_sct[ixo]
          # CV-analysis and variances at station points
        if (!is.na(argv$cores)) {
          arr<-t(mcmapply(oi_var_gridpoint_by_gridpoint,
                          1:nobs_to_check,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          dh=argv$DhorMin.sct,
                          box_o_nearest_halfwidth=argv$box_o_nearest_halfwidth.sct,
                          dz=argv$Dver.sct,
                          lafmin=argv$lafmin.sct,
                          dh_adaptive=T,
                          corr=argv$corr.sct,
                          pmax=argv$pmax.sct,
                          fg=argv$fglab.sct,
                          fg_gamma=argv$fg_gamma.sct,
                          fg_min=fg_min,
                          fg_max=fg_max,
                          succ_corr=argv$succ_corr.sct,
                          y_elab=F,
                          loocv=T,
                          o_errvar_min=argv$o_errvar_min.sct,
                          o_errvar_max=argv$o_errvar_max.sct,
                          xa_errvar_max=argv$xa_errvar_max.sct,
                          xa_errvar_min=argv$xa_errvar_min))
        # no-multicores
        } else {
          arr<-t(mapply(oi_var_gridpoint_by_gridpoint,
                        1:nobs_to_check,
                        SIMPLIFY=T,
                        dh=argv$DhorMin.sct,
                        box_o_nearest_halfwidth=argv$box_o_nearest_halfwidth.sct,
                        dz=argv$Dver.sct,
                        lafmin=argv$lafmin.sct,
                        dh_adaptive=T,
                        corr=argv$corr.sct,
                        pmax=argv$pmax.sct,
                        fg=argv$fglab.sct,
                        fg_gamma=argv$fg_gamma.sct,
                        fg_min=fg_min,
                        fg_max=fg_max,
                        succ_corr=argv$succ_corr.sct,
                        y_elab=F,
                        loocv=T,
                        o_errvar_min=argv$o_errvar_min.sct,
                        o_errvar_max=argv$o_errvar_max.sct,
                        xa_errvar_max=argv$xa_errvar_max.sct,
                        xa_errvar_min=argv$xa_errvar_min))
        }
        yo_errvar<-arr[,3]
        rm(arr)
        corep[ixo]<-yo_errvar/mean(yo_errvar)
        t01a<-Sys.time()
        print(paste(" repres coeff=",cont_sct_loop,
                    "/#obs=",length(ixg),
                    "/time",round(t01a-t00a,1),attr(t01a-t00a,"unit")))
      } # end calculation of representativeness coefficients
      if (exists("arr")) rm(arr)
      if (exists("ixo")) rm(ixo,xobs_spint,yobs_spint,zobs_spint,lafobs_spint,
                            eps2_spint,nobs,yo_spint)
      if (exists("yb_spint")) rm(yb_spint)
      if (exists("ixg")) rm(ixg,nobs_to_check,xgrid_spint,ygrid_spint,
                           zgrid_spint,lafgrid_spint,yo_to_check)
      if (exists("xb_spint")) rm(xb_spint)
      if (exists("yo_sct")) rm(yo_sct)
      if (exists("b_sct")) rm(b_sct)
      if (exists("yav_prev")) rm(yav_prev)
      if (exists("yav")) rm(yav)
      if (exists("yav_errvar")) rm(yav_errvar)
      if (exists("yo_errvar")) rm(yo_errvar)
      if (exists("dh_ref")) rm(dh_ref)
      # 
    # SCT split grid into boxes
    } else {
      # set min and max for the background values
      sctvmin<-ifelse(argv$variable=="RR",-1./argv$boxcox.lambda,
                                          argv$vmin)
      sctvmax<-ifelse(argv$variable=="RR",boxcox(argv$vmax,argv$boxcox.lambda),
                                          argv$vmax)
      # create the grid for SCT. SCT is done independently in each box
      # NOTE: box size around 100Km should be ok
      if (i==1) {
        r<-raster(e,ncol=argv$grid.sct[2],nrow=argv$grid.sct[1])
        xy<-xyFromCell(r,1:ncell(r))
        xr<-xy[,1]
        yr<-xy[,2]
        ir<-1:ncell(r)
        r[]<-1:ncell(r)
      }
      # define global 1D vector used in sct (1D for fast access)
      xtot<-x[ix]
      ytot<-y[ix]
      ztot<-z[ix]
      if (argv$variable=="RR") {
        ttot<-boxcox(x=data$value[ix],lambda=argv$boxcox.lambda)
      } else {
        ttot<-data$value[ix]
      }
      pridtot<-data$prid[ix]
      laftot<-laf[ix]
      # assign each station to the corresponding box
      if (argv$debug) itotdeb<-extract(r,cbind(x,y))
      itot<-extract(r,cbind(xtot,ytot))
      # count the number of observations in each box
      rnobs<-rasterize(cbind(xtot,ytot),r,ttot,fun=function(x,...)length(x))
      nr<-getValues(rnobs)
      # create the 4D array for the function call via apply
      ixyn<-cbind(ir,xr,yr,nr)
      # SCT within each (grid) box
      # NOTE: dqcflag is updated in "sct" function
      out<-apply(ixyn,
                 FUN=sct,MARGIN=1,
                 nmin=argv$n.sct,
                 dzmin=argv$dz.sct,
                 Dhmin=argv$DhorMin.sct,
                 Dz=argv$Dver.sct,
                 eps2=argv$eps2.sct,
                 T2=argv$thr.sct,
                 T2pos=argv$thrpos.sct,
                 T2neg=argv$thrneg.sct,
                 sus.code=argv$sct.code,
                 faster=argv$fast.sct)
    } # endif SCT stn-by-stn or boxes
  } else {
    print("no valid observations left, no SCT")
  }
  ncur<-length(which(dqcflag==argv$sct.code))
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    str<-"("
    for (f in 1:nfin) {
      if (f>1) str<-paste0(str,"; ") 
      aux<-which(dqcflag==argv$sct.code & data$prid==argv$prid[f])
      str<-paste0(str,"prid",argv$prid[f],"=",length(aux))
      rm(aux)
    }
    str<-paste0(str,")")
    print(paste("SCT, iteration=",i,
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    print(paste("# suspect observations=",ncur-nprev,str))
    rm(str)
  }
  if ((ncur-nprev)==0) break
  nprev<-ncur
}
rm(doit)
if (argv$verbose | argv$debug) 
  print("+---------------------------------+")
if (argv$debug) 
  save.image(file.path(argv$debug.dir,"dqcres_sct.RData")) 
#
# coefficient of observation representativeness
#-----------------------------------------------------------------------------
# corep has been set by function sct to the observation error variance
if (!any(is.na(argv$const.corep))) {
  for (f in 1:nfin) {
    if (any(data$prid[ix]==argv$prid[f])) {
      ip<-which(data$prid[ix]==argv$prid[f])
      if (length(ip)>0) corep[ix[ip]]<-argv$const.corep[f]
    } else {
      print(paste("provider ",argv$prid[f],
      ": no valid data found to compute the coefficient of representativeness",sep=""))
    }
  }
} else {
  ix<-which(!is.na(corep) & (is.na(dqcflag) | dqcflag==argv$keep.code)) 
  if (length(ix)>0) {
    qmn<-0.25
    qmx<-0.75
    qav<-0.5
    acorep<-abs(corep[ix])
    # ecdf(x)(x) here should give us something similar to rank(x)/length(x)
    qcorep<-ecdf(acorep)(acorep)
    if (any(qcorep<qmn)) qcorep[which(qcorep<qmn)]<-qmn
    if (any(qcorep>qmx)) qcorep[which(qcorep>qmx)]<-qmx
    for (f in 1:nfin) {
      if (any(data$prid[ix]==argv$prid[f])) {
        ip<-which(data$prid[ix]==argv$prid[f] & qcorep<=qav)
        if (length(ip)>0)      
          corep[ix[ip]]<-argv$min.corep[f]+
            (argv$mean.corep[f]-argv$min.corep[f])*(qcorep[ip]-qmn)/(qav-qmn)
        ip<-which(data$prid[ix]==argv$prid[f] & qcorep>qav)
        if (length(ip)>0)      
          corep[ix[ip]]<-argv$mean.corep[f]+
            (argv$max.corep[f]-argv$mean.corep[f])*(qcorep[ip]-qav)/(qmx-qav)
      } else {
        print(paste("provider ",argv$prid[f],
        ": no valid data found to compute the coefficient of representativeness",sep=""))
      }
    }
  } else {
    print("no valid first guess for the observation error variance found")
  } 
}
if (argv$debug) 
  save.image(file.path(argv$debug.dir,"dqcres_sct.RData")) 
#-----------------------------------------------------------------------------
# cool test (Check fOr hOLes in the field)
if (argv$cool) {
  nprev<-0
  if (argv$verbose | argv$debug) {
    print(paste0("cool test (",argv$cool.code,")"))
  }
  # set doit vector
  doit<-vector(length=ndata,mode="numeric")
  doit[]<-NA
  for (f in 1:nfin) doit[data$prid==argv$prid[f]]<-argv$doit.cool[f]
  ix<-which((is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0)
  ptmp<-length(ix)
  if (ptmp<1) {
    print("cool test  no valid observations left, no test")
  } else {
    rgrid_cool<-raster(ext=e,resolution=argv$grid_res.cool)
    rgrid_cool[]<-NA
    xygrid_cool<-xyFromCell(rgrid_cool,1:ncell(rgrid_cool))
    xgrid_cool<-xygrid_cool[,1]
    ygrid_cool<-xygrid_cool[,2]
    rm(xygrid_cool)
    ngrid_cool<-length(ygrid_cool)
    if (!exists("xobs_cool_aux")) xobs_cool_aux<-integer(0)
    if (!exists("yobs_cool_aux")) yobs_cool_aux<-integer(0)
    if (!exists("pridobs_cool_aux")) pridobs_cool_aux<-integer(0)
    if (!exists("yo_cool_aux")) yo_cool_aux<-integer(0)
    # test
    for (i in 1:argv$i.cool) {
      t0a<-Sys.time()
      for (j in 1:length(argv$thres.cool)) {
        # use only (probably) good observations
        ix<-which((is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0)
        if (length(ix)>0) {
          xobs_to_check_cool<-x[ix]
          yobs_to_check_cool<-y[ix]
          pridobs_to_check_cool<-data$prid[ix]
          xobs_cool<-c(x[ix],xobs_cool_aux)
          yobs_cool<-c(y[ix],yobs_cool_aux)
          yo_cool<-c(data$value[ix],yo_cool_aux)
          rgrid1<-rasterize(x=cbind(xobs_cool,yobs_cool),
                            y=rgrid_cool,
                            field=yo_cool,fun=mean,na.rm=T)
          ix1<-which(!is.na(getValues(rgrid1)))
          xy1<-xyFromCell(rgrid1,ix1)
          xobs_cool<-xy1[,1]
          yobs_cool<-xy1[,2]
          yo_cool<-getValues(rgrid1)[ix1]
          if (!is.na(argv$cores)) {
            arr<-t(mcmapply(spint_cool,
                            1:ngrid_cool,
                            mc.cores=argv$cores,
                            SIMPLIFY=T,
                            thres=argv$thres.cool[j],
                            condition=argv$condition.cool[j],
                            dh_max=argv$dh_max.cool))
          # no-multicores
          } else {
            arr<-t(mapply(spint_cool,
                          1:ngrid_cool,
                          SIMPLIFY=T,
                          thres=argv$thres.cool[j],
                          condition=argv$condition.cool[j],
                          dh_max=argv$dh_max.cool))
          }
          rgrid_cool[]<-arr
          rc<-clump(rgrid_cool)
          dc<-getValues(rc)
          freq_rc<-freq(rc)
          ytmp<-extract(rc,cbind(xobs_to_check_cool,yobs_to_check_cool),na.rm=T)
          for (k in 1:length(freq_rc[,1])) {
            if (is.na(freq_rc[k,1])) next
            flag_k<-ytmp==freq_rc[k,1] & !is.na(ytmp)
            if (length(which(flag_k))<n.cool[j,1]) {
              for (f in 1:nfin) {
                ixijkf<-which(flag_k & pridobs_to_check_cool==argv$prid[f])
                if (length(ixijkf)<n.cool[j,(f+1)]) {
                  dqcflag[ix[ixijkf]]<-argv$cool.code
                }
              } 
            }
          }
        } else {
          print("no valid observations left, no cool test" )
        }
      } # end of loop over threshold that define events
      ncur<-length(which(dqcflag==argv$cool.code & !is.na(dqcflag)))
      if (argv$verbose | argv$debug) {
        t1a<-Sys.time()
        print(paste("cool test. iteration=",i,
                    "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
        print(paste("# suspect observations=",ncur-nprev))
      } 
      if ((ncur-nprev)==0) break
      nprev<-ncur
    } # end of loop over test iterations
    rm(doit)
    if (argv$verbose | argv$debug) 
      print("+---------------------------------+")
    if (argv$debug) 
      save.image(file.path(argv$debug.dir,"dqcres_cool.RData")) 
  } # end of "if (length(mask_cool)==0)"
} # end of if (argv$cool)
#
#-----------------------------------------------------------------------------
# check for isolated stations
# use only (probably) good observations
# NOTE: keep-listed stations canNOT be flagged here
# set doit vector
doit<-vector(length=ndata,mode="numeric")
doit[]<-NA
for (f in 1:nfin) doit[data$prid==argv$prid[f]]<-argv$doit.iso[f]
#
ix<-which(is.na(dqcflag) & doit!=0)
if (length(ix)>0) {
  # define global 1D vector used in nstat (1D for fast access)
  xtot<-x[ix]
  ytot<-y[ix]
  xy<-cbind(xtot,ytot)
  ns<-apply(xy,FUN=nstat,MARGIN=1,drmin=argv$dr.isol)
  sus<-which(ns<argv$n.isol & doit[ix]==1)
  # set dqcflag
  if (length(sus)>0) dqcflag[ix[sus]]<-argv$isol.code
} else {
  print("no valid observations left, no check for isolated observations")
}
rm(doit)
if (argv$verbose | argv$debug) {
  print(paste("# isolated observations=",length(which(dqcflag==argv$isol.code))))
  print("+---------------------------------+")
}
if (argv$debug) 
  save.image(file.path(argv$debug.dir,"dqcres_iso.RData")) 
#
#-----------------------------------------------------------------------------
# observations not flagged are assumed to be good observations 
dqcflag[is.na(dqcflag)]<-0
if (argv$verbose | argv$debug) {
  nsus<-length(which(dqcflag!=0))
  nok<-length(which(dqcflag==0))
  nsus_notna<-length(which(dqcflag!=0 & !is.na(data$value)))
  nok_notna<-length(which(dqcflag==0 & !is.na(data$value)))
  nnotna<-length(which(!is.na(data$value)))
  nna<-length(which(is.na(data$value)))
  print("summary:")
  print(" #  NAs, number of observations with no value (NAs)")
  print(" #  sus, number of suspicious observations or no-metadata")
  print(" # good, number of good observations")
  print(" NOTE for sus and good, the statistics consider only observations not NAs")
  print("summary:")
  print(paste0(" #  NAs= ",nna))
  print(paste0(" #  sus= ",nsus_notna," [",round(100*nsus_notna/nnotna,0),"%]"))
  print(paste0(" # good= ",nok," [", round(100*nok_notna/nnotna,0), "%]"))
   if (nfin>1) {
    for (f in 1:nfin) {
      nsus<-length(which(dqcflag!=0 & data$prid==argv$prid[f]))
      nok<-length(which(dqcflag==0 & data$prid==argv$prid[f]))
      nsus_notna<-length(which(dqcflag!=0 & 
                               !is.na(data$value) & 
                               data$prid==argv$prid[f]))
      nok_notna<-length(which(dqcflag==0 & 
                              !is.na(data$value) & 
                              data$prid==argv$prid[f]))
      nnotna<-length(which(!is.na(data$value) & 
                           data$prid==argv$prid[f]))
      nna<-length(which(is.na(data$value) & data$prid==argv$prid[f]))
      print(paste("--> summary provider",argv$prid[f]))
      print(paste0("  #  NAs= ",nna))
      print(paste0("  #  sus= ",nsus_notna," [",round(100*nsus_notna/nnotna,0),"%]"))
      print(paste0("  # good= ",nok," [", round(100*nok_notna/nnotna,0), "%]"))

    }
  }
  print("+---------------------------------+")
}
#
#-----------------------------------------------------------------------------
# Include radar-derived precipitation in the output file
if (argv$radarout) {
  if (argv$verbose | argv$debug) 
    print("include radar-derived precipitation in the output file")
  # (optional) aggregate radar data onto a coarser grid
  if (!is.na(argv$radarout.aggfact) & argv$radarout.aggfact>1) {
    raux<-aggregate(rrad,
                    fact=argv$radarout.aggfact,
                    na.rm=T,
                    expand=T, 
                    fun=mean)
    rrad<-raux
    rm(raux)
  }
  drad<-getValues(rrad)
  # get radar-point coordinates into output CRS 
  ix<-which(!is.na(drad)) 
  if (length(ix)) {
    radxy<-as.data.frame(xyFromCell(rrad,ix))
    names(radxy)<-c("x","y")
    coordinates(radxy)<-c("x","y")
    proj4string(radxy)<-CRS(argv$proj4fg)
    radxy.from<-as.data.frame(spTransform(radxy,CRS=argv$proj4_output_files))  
    radx.from<-radxy.from[,1]
    rady.from<-radxy.from[,2]
    radrr<-drad[ix]
  } else {
    radx.from<-integer(0)
    rady.from<-integer(0)
    radrr<-integer(0)
  }
}
#
#-----------------------------------------------------------------------------
# write the output file
if (argv$verbose | argv$debug) {
  print("write the output file")
  print(argv$output)
}
if (argv$proj4_output_files!=argv$proj4_input_obsfiles) {
    xy<-as.data.frame(x=data$lon,y=data$lat)
    coordinates(xy)<-c("x","y")
    proj4string(xy)<-CRS(argv$proj4_input_obsfiles)
    xyt<-as.data.frame(spTransform(xy,CRS=argv$proj4_output_files))  
    xout<-xyt[,1]
    yout<-xyt[,2]
    rm(xy,xyt)
} else {
  xout<-data$lon
  yout<-data$lat
}
varidx.out<-varidx
if (any(!is.na(argv$varname.opt))) 
  varidx.out<-c(varidx,varidx.opt[which(!is.na(varidx.opt))]) 
dataout<-array(data=NA,
               dim=c(length(yout),(length(varidx.out)+4)))
ord.varidx.out<-order(varidx.out)
str<-vector()
for (s in 1:length(ord.varidx.out)) {
  varidx.out.s<-varidx.out[ord.varidx.out[s]]
  pos.s<-which(varidx.out==varidx.out.s)
  if (pos.s>4) {
    posopt.s<-which(varidx.opt==varidx.out.s & !is.na(varidx.opt))
    posopt.nona.s<-which(varidx.opt[which(!is.na(varidx.opt))]==varidx.out.s)
    str[s]<-argv$varname.opt[posopt.s]
    dataout[,s]<-dataopt[,posopt.nona.s]
  } else if (pos.s==1) {
    str[s]<-argv$varname.x.out
    dataout[,s]<-round(yout,argv$xy.dig.out)
  } else if (pos.s==2) {
    str[s]<-argv$varname.y.out
    dataout[,s]<-round(xout,argv$xy.dig.out)
  } else if (pos.s==3) {
    str[s]<-argv$varname.elev.out
    dataout[,s]<-round(z,argv$elev.dig.out)
  } else if (pos.s==4) {
    str[s]<-argv$varname.value.out
    dataout[,s]<-round(data$value,argv$value.dig.out)
  }
}
str[s+1]<-argv$varname.prid
dataout[,(s+1)]<-data$prid
str[s+2]<-argv$varname.dqc
dataout[,(s+2)]<-dqcflag
str[s+3]<-argv$varname.sct
dataout[,(s+3)]<-round(sctpog,2)
str[s+4]<-argv$varname.rep
dataout[,(s+4)]<-round(corep,5)
if (argv$radarout) {
  if (length(radrr)>0) {
    datarad<-array(data=NA,dim=c(length(radrr),(length(varidx.out)+4)))
    datarad[,which(str==argv$varname.x.out)]<-round(rady.from,argv$xy.dig.out)
    datarad[,which(str==argv$varname.y.out)]<-round(radx.from,argv$xy.dig.out)
    datarad[,which(str==argv$varname.value.out)]<-round(radrr,argv$value.dig.out)
    datarad[,which(str==argv$varname.prid)]<-rep(argv$radarout.prid,length(radrr))
    datarad[,which(str==argv$varname.dqc)]<-rep(0,length(radrr))
    dataout<-rbind(dataout,datarad)
  }
}
dataout<-as.data.frame(dataout,stringsAsFactors=F)
names(dataout)<-str
write.table(file=argv$output,
            dataout,
            quote=F,
            col.names=str,
            row.names=F,
            dec=".",
            sep=argv$separator.out)
if (argv$verbose | argv$debug) 
  print("+---------------------------------+")
#
#-----------------------------------------------------------------------------
# Normal exit
t1<-Sys.time()
if (argv$verbose | argv$debug) 
 print(paste("normal exit /time",round(t1-t0,1),attr(t1-t0,"unit")))
q(status=0)
