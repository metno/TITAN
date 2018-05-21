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

#+ call the C-functions
`OI_RR_fast`<-function(yo.sel,
                       yb.sel,
                       xb.sel,
                       xgrid.sel,
                       ygrid.sel,
                       zgrid.sel,
                       VecX.sel,
                       VecY.sel,
                       VecZ.sel,
                       Dh.cur,
                       Dz.cur) {
#------------------------------------------------------------------------------
  no<-length(yo.sel)
  ng<-length(xb.sel)
  xa.sel<-vector(mode="numeric",length=ng)
  vec<-vector(mode="numeric",length=no)
  d<-yo.sel-yb.sel
  out<-.C("oi_rr_first",no=as.integer(no), 
                        innov=as.double(d),
                        SRinv=as.numeric(InvD),
                        vec=as.double(vec) ) 
  vec[1:no]<-out$vec[1:no]
  rm(out)
  out<-.C("oi_rr_fast",ng=as.integer(ng),
                       no=as.integer(no),
                       xg=as.double(xgrid.sel),
                       yg=as.double(ygrid.sel),
                       zg=as.double(zgrid.sel),
                       xo=as.double(VecX.sel),
                       yo=as.double(VecY.sel),
                       zo=as.double(VecZ.sel),
                       Dh=as.double(Dh.cur),
                       Dz=as.double(Dz.cur),
                       xb=as.double(xb.sel),
                       vec=as.double(vec),
                       xa=as.double(xa.sel) )
  out$xa[1:ng]
}


# auxiliary function to keep/blacklist observations
setCode_lonlat<-function(lonlat,code) {
# lonlat. vector. 1=lon; 2=lat
  ix<-datatmp$lon==lonlat[1] & datatmp$lat==lonlat[2]
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
statSpat<-function(ixyzt,drmin,gamma=-0.0065) {
# input
# ixyz=vector. 1=index;2=x;3=y;4=z;5=t
# output
# 1=nobs;2=maxVertDist[m];3=mean(temp);4=sd(temp)
# NOTE: temperature is adjusted for elevation differences (gamma=-0.0065 K/m)
#------------------------------------------------------------------------------
  if (any(is.na(ixyzt))) return(c(NA,NA,NA,NA))
  i<-which(abs(ixyzt[2]-xtot)<=drmin & 
           abs(ixyzt[3]-ytot)<=drmin)
  if (length(i)==1) return(c(1,NA,NA,NA))
  dz<-ztot[i]-ixyzt[4]
  if (argv$variable=="T") {
    tcor<-ttot[i]-gamma*dz
  } else {
    tcor<-ttot[i]
  }
  tmean<-mean(tcor) 
  tsd<-sd(tcor) 
  dz_mx<-max(abs(dz))
  return(c(length(i),round(dz_mx,0),round(tmean,1),round(tsd,3)))
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
#  rr.cor<-rr
#  ix<-which(!is.na(rr) & rr>=0.)
#  if (length(ix)>0) rr.cor[ix]<-rr[ix]*1./f[ix]
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
#                   2/3=easting/northing coord (center of the box);
#                   4=number of stations within the box
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
  # something strange with the number of stations
  if (is.na(ixynp[4]) | is.null(ixynp[4]) | !is.finite(ixynp[4]) ) return(NA)
  # j, index for the stations in the box
  if (argv$usefge.sct) {
    j<-which(itot==ixynp[1] & !is.na(fge.mu[ix]))
  } else if (argv$usefg.sct) {
    j<-which(itot==ixynp[1] & !is.na(fg[ix]))
  } else {
    j<-which(itot==ixynp[1])
  }
  if (length(j)==0) return(-1)
  # case of just one station
  if (ixynp[4]==1) {
    sctpog[ix[j]]<--1
    assign("sctpog",sctpog,envir=.GlobalEnv)
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
#    j<-which(itot==ixynp[1] & !is.na(fge.mu[ix]))
    tb<-fge.mu[ix[j]]
  } else if (argv$usefg.sct) {
#    j<-which(itot==ixynp[1] & !is.na(fg[ix]))
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
#  if (any(is.na(tb))) return(NULL)
  if (any(tb<sctvmin)) tb[which(tb<sctvmin)]<-sctvmin
  if (any(tb>sctvmax)) tb[which(tb>sctvmax)]<-sctvmax
  # OI for SCT (Lussana et al., 2010)
  # initialize variables
  # expand eps2 in eps2.vec
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
  # innvoation
  d<-topt-tb
  # select observations used in the test
  sel<-which(is.na(dqctmp) | dqctmp==argv$keep.code)
  # select observations to test 
  sel2check<-which(is.na(dqctmp) & doittmp==1)
  first<-T
  # loop over SCT iterations 
  # NOTE: SCT flags at most one observation, iterate until no observations fail
  while (length(sel)>1) { 
#    # update selection
#    sel<-which(is.na(dqctmp) | dqctmp==argv$keep.code)
#    sel2check<-which(is.na(dqctmp))
#    if (length(sel2check)==0) break
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
    if (sig2o<0.01) sig2o<-0.01       # safe threshold  
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
  # coefficient of observation representativeness
  # this call to ecdf(x)(x) should be the same as rank(x)/length(x)
  corep[ix[j[sel]]]<-(d[sel]*(-ares))/sig2o
  assign("corep",corep,envir=.GlobalEnv)
  # debug: begin
#  if (argv$debug) {
#  }
  # debug: end
  return(length(which(dqctmp==sus.code)))
}


#+ find the n-th largest element from each matrix row 
findRow <- function (x,n) {   
# order by row number then by value
  y<-t(x)
  array(y[order(col(y), y)], dim(y))[nrow(y) - (n-1), ]
}

#+ mean radial distance between an observation and its k-th closest obs
#dobs_fun<-function(ixynp,k) {
dobs_fun<-function(obs,k) {
# NOTE: k=1 will return 0 everywhere (1st closest obs is the obs itself)
  nobs<-length(obs$x)
  if (nobs<k) return(NA)
  # distance matrices
  disth<-(outer(obs$x,obs$x,FUN="-")**2.+
          outer(obs$y,obs$y,FUN="-")**2.)**0.5
  dobsRow<-findRow(x=disth,n=(nobs-k+1))
  mean(dobsRow)
}

#+ remove clumps determined by less than n observations
remove_ltNobsClumps<-function(r,n,reverse=F,obs) {
  obsflag<-obs$x
  obsflag[]<-0
  d<-getValues(r)
  if (reverse) {
    dd<-d
    d[dd==0]<-1
    d[dd>0]<-0
    r[]<-d
  }
  rc<-clump(r)
  dc<-getValues(rc)
  f<-freq(rc)
  ytmp<-extract(rc,cbind(obs$x,obs$y),na.rm=T)
  sel<-which(!is.na(ytmp))
  for (i in 1:length(f[,1])) {
    if (is.na(f[i,1])) next
    if (any(ytmp[sel]==f[i,1])) {
      ix<-which(ytmp[sel]==f[i,1])
      if (length(ix)<n) {
        dc[dc==f[i,1]]<-NA
        obsflag[sel[ix]]<-1
      }
    }
  }
  dc[is.na(dc)]<-(-1)
  if (any(dc==0)) d[which(dc==0)]<-0
  if (reverse) {
    dd<-d
    d[dd==0]<-1
    d[dd>0]<-0
  }
  r[]<-d
  return(list(r=r,obsflag=obsflag))
}


# + puddle - puddle check
puddle<-function(obs,
                 thres=1,
                 gt_or_lt="gt", #event definition: greater than or less than...
                 Dh=25, # km
                 eps2=0.1,
                 n.eveYES_in_the_clump=2,
                 n.eveNO_in_the_clump=2) {
#------------------------------------------------------------------------------
  p<-length(obs$x)
  ydqc.flag<-vector(mode="numeric",length=p)
  ydqc.flag[]<-0
# eix0: vector of positions to observations that are equivalent to eventNO
# eix1: vector of positions to observations that are equivalent to eventYES
  if (gt_or_lt=="gt") {
    e0<-obs$yo<=thres
    e1<-obs$yo >thres
  } else if (gt_or_lt=="lt") {
    e0<-obs$yo>=thres
    e1<-obs$yo< thres
  } else if (gt_or_lt=="le") {
    e0<-obs$yo >thres
    e1<-obs$yo<=thres
  } else if (gt_or_lt=="ge") {
    e0<-obs$yo <thres
    e1<-obs$yo>=thres
  }
  eix0<-which(e0)
  eix1<-which(e1)
  #
  p0<-length(eix0)
  p1<-length(eix1)
  if ( (p0<=1) | (p1<=1) ) return(NULL)
  # distance matrix
  Disth<-matrix(ncol=p,nrow=p,data=0.)
  Disth<-(outer(obs$y,obs$y,FUN="-")**2.+
          outer(obs$x,obs$x,FUN="-")**2.)**0.5/1000.
  # assign grid points to event yes/no
  t00<-Sys.time()
  #
  D<-exp(-0.5*(Disth[eix0,eix0]/Dh)**2.)
  diag(D)<-diag(D)+eps2
  InvD<-chol2inv(chol(D))
  assign("InvD",InvD,envir=.GlobalEnv)
  rm(D)
  xidi0<-OI_RR_fast(yo.sel=rep(1,p0),
                    yb.sel=rep(0,p0),
                    xb.sel=rep(0,length(xgrid_puddle)),
                    xgrid.sel=xgrid_puddle,
                    ygrid.sel=ygrid_puddle,
                    zgrid.sel=zgrid_puddle,
                    VecX.sel=obs$x[eix0],
                    VecY.sel=obs$y[eix0],
                    VecZ.sel=rep(0,p0),
                    Dh.cur=Dh,
                    Dz.cur=1000000) 
  t11<-Sys.time()
#  if (argv$verbose | argv$debug) 
#    print(paste("xIDI event NO, time=",round(t11-t00,1),attr(t11-t00,"unit")))
  t00<-Sys.time()
  D<-exp(-0.5*(Disth[eix1,eix1]/Dh)**2.)
  diag(D)<-diag(D)+eps2
  InvD<-chol2inv(chol(D))
  assign("InvD",InvD,envir=.GlobalEnv)
  rm(D)
  xidi1<-OI_RR_fast(yo.sel=rep(1,p1),
                    yb.sel=rep(0,p1),
                    xb.sel=rep(0,length(xgrid_puddle)),
                    xgrid.sel=xgrid_puddle,
                    ygrid.sel=ygrid_puddle,
                    zgrid.sel=zgrid_puddle,
                    VecX.sel=obs$x[eix1],
                    VecY.sel=obs$y[eix1],
                    VecZ.sel=rep(0,p1),
                    Dh.cur=Dh,
                    Dz.cur=1000000) 
  t11<-Sys.time()
#  if (argv$verbose | argv$debug) 
#    print(paste("xIDI eventYES, time=",round(t11-t00,1),attr(t11-t00,"unit")))
  rm(InvD,envir=.GlobalEnv)
# x1=1 for grid points where event is YES
  x1<-xidi1
  x1[]<-0
  if (any(xidi1>=xidi0)) {
    x1[which(xidi1>=xidi0)]<-1
  } else {
    return(NULL)
  }
# -- remove YES(NO) observations close to a lot of NO(YES) ones --
# observations might be correct but they will create unrealistic patterns, 
# either a "puudle"/hole or a isolated max, in the analysis field because the local 
# observation density is not sufficient to properly describe such small-scale
# features
  r1<-rgrid_puddle
  r1[]<-NA
  r1[mask_puddle]<-x1
#  ytmp<-extract(r1,cbind(obs$x,obs$y),buffer=(2*mean(res(r1))),fun=mean,na.rm=T)
  ytmp<-extract(r1,cbind(obs$x,obs$y))
  # case 1. observation YES(NO) in an area with significant prevalence of NO(YES)
  cond<-(ytmp==0 & e1) | (ytmp==1 & e0)
  if (any(cond)) ydqc.flag[which(cond)]<-1
  # case 2. puddle of YES surrounded by NO
  out<-remove_ltNobsClumps(r=r1,
                           n=n.eveYES_in_the_clump,
                           reverse=F,
                           obs=data.frame(x=obs$x,
                                          y=obs$y))  
  if (any(out$obsflag==1)) ydqc.flag[which(out$obsflag==1)]<-1
  rm(out)
  # case 3. puddle of NO surrounded by YES
  out<-remove_ltNobsClumps(r=r1,
                           n=n.eveNO_in_the_clump,
                           reverse=T,
                           obs=data.frame(x=obs$x,
                                          y=obs$y))  
  if (any(out$obsflag==1)) ydqc.flag[which(out$obsflag==1)]<-1
  rm(out)
  ydqc.flag
}

# + STEVE - iSolaTed EVent tEst
steve<-function(obs,
                thres=1,
                gt_or_lt="gt", #event definition: greater than or less than...
                pmax=20, 
                dmax=150000,
                n.sector=16,
                n.connected_eveYES=1,
                frac.eveYES_in_the_clump=1,
                dmin.next_eveYES=150000) {
#------------------------------------------------------------------------------
  p<-length(obs$x)
  pmax<-min(p,pmax)
  psub<-vector(mode="numeric",length=p)
  psub[]<-NA
  sector.angle<-360/n.sector
  ydqc.flag<-vector(mode="numeric",length=p)
  ydqc.flag[]<-0
# eix: vector of positions to observations that are equivalent to eventYES
  if (gt_or_lt=="gt") {
    eix<-which(obs$yo>thres)
  } else if (gt_or_lt=="lt") {
    eix<-which(obs$yo<thres)
  } else if (gt_or_lt=="le") {
    eix<-which(obs$yo<=thres)
  } else if (gt_or_lt=="ge") {
    eix<-which(obs$yo>=thres)
  }
  # no events= no isolated events
  if (length(eix)==0) return(ydqc.flag)
  # group eventsYES/NO into clumps of connected eventYES
  # vectors:
  # +neve.vec>number of eventYES clumps
  # +Leve.vec[i]>(i=1,..,neve.vec) number of stns in i-th clump 
  # +eve.vec[i,n]-> n-th stn (i.e. pointers to obs) in i-th clump
  #  matrix(i,j) i=1,neve.vec j=1,Leve.vec[i]
  #  note: a clumps of connected eventYES points may include eventNO points too
  #        provided that an eventNO point is connected to eventYES points
  eve.vec<-matrix(data=NA,ncol=p,nrow=p)
  Leve.vec<-vector(mode="numeric",length=p)
  Leve.vec[]<-NA
  neve.vec<-0
  eve.aux<-matrix(data=NA,ncol=p,nrow=p)
  Leve.aux<-vector(mode="numeric",length=p)
  Leve.aux[]<-NA
  neve.aux<-0
  for (b in eix) {  # START: Cycle over eventYES points 
    # a. identify the closest points to the b-th eventYES point 
    #  outputs: -close2b-> pointer to the closest points
    #            constraints: distance must be less than Lsubsample.DHmax 
    #                         max number of points allowed is Lsubsample.max
    #           -psub[b]-> actual number of points used
    Disth<-sqrt((obs$y[b]-obs$y)**2.+(obs$x[b]-obs$x)**2.)
    if (any(Disth<=dmax)) {
      ix<-which(Disth<=dmax)
      psub[b]<-min(pmax,length(ix))
      c2b<-order(Disth,decreasing=F)[1:psub[b]]
    } else {
      psub[b]<-0
    }
    if (psub[b]<3) next
    close2b<-c2b[1:psub[b]]
    rm(c2b)
    # c. Establish (local-) triangulation (Delauney)
    #    note: all points are used, despite the YES/NO thing
    aux<-cbind(obs$x[close2b],obs$y[close2b])
    nokloop<-1
    # tri.mesh is not happy with duplicates, then we use the jitter option
    jit<-1/mean(c(diff(range(obs$x[close2b])),
                  diff(range(obs$x[close2b]))))
    aux<-cbind(obs$x[close2b],obs$y[close2b])
    nokloop<-1
    while (anyDuplicated(aux)) {
      ixdup<-which(duplicated(aux))
      ndup<-length(ixdup)
      obs$x[close2b][ixdup]<-obs$x[close2b][ixdup]+runif(ndup)
      obs$y[close2b][ixdup]<-obs$y[close2b][ixdup]+runif(ndup)
      aux<-cbind(obs$x[close2b],obs$y[close2b])
      nokloop<-nokloop+1
      if (nokloop>10) {
        print("ERROR: check the input data for too many duplicated stations")
        quit(status=1)
      }
    }
    options(warn = 1, scipen = 999)
    tri.rr<-tri.mesh(obs$x[close2b],obs$y[close2b],
                     jitter=jit,jitter.iter=100,jitter.random=T)
    options(warn = 2, scipen = 999)
    # d. identify all the points (wheter YES or NO) which belongs
    #    to adjacent nodes (respect to the b-th YES point node)
    #  procedure: tri.find-> returns the three vertex indexes of triangle
    #   containing the specified point. To find all the neighbouring 
    #   triangles just span the area surrounding the b-th YES point
    #   (circle of radius=1Km)
    #  output: nodes-> point position respect to obs
    nodes<-vector(mode="numeric")
    tri.loc<-tri.find(tri.rr,obs$x[b],obs$y[b])
    # note: due to the fact that (obs$x[b],obs$y[b]) belongs to the mesh
    # used to establish the triangulation, ...i1,i2,i3 must be >0
    nodes<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
    for (cont.clock in 1:n.sector) {
      x.aux<-obs$x[b]+sin(sector.angle*(cont.clock-1)*pi/180.)
      y.aux<-obs$y[b]+cos(sector.angle*(cont.clock-1)*pi/180.)
      tri.loc<-tri.find(tri.rr,x.aux,y.aux)
      nodes.aux<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
      nodes<-c(nodes,
               nodes.aux[which( (!(nodes.aux %in% nodes)) &
                                  (nodes.aux>0) )])
    }
    rm(x.aux,y.aux,tri.loc,nodes.aux)
    lnodes<-length(nodes)
    # e. update the temporary structure used to identify eve 
    #    merge in a (new) eve all the (old) temporary eve (if any) 
    #    which contain at least one point present in nodes (avoiding
    #    repetition); erase (i.e. set to NAs) the (old) temporary eve merged
    aux<-vector()
    if (neve.aux>0) {
      for (eve in 1:neve.aux) {
        if (Leve.aux[eve]==0) next
        if ( any(nodes %in% eve.aux[eve,1:Leve.aux[eve]]) ) {
          aux<-c(aux[which(!(aux %in% eve.aux[eve,1:Leve.aux[eve]]))],
                 eve.aux[eve,1:Leve.aux[eve]])
          eve.aux[eve,]<-NA
          Leve.aux[eve]<-0
        }
      }
    }
    neve.aux<-neve.aux+1
    Leve.aux[neve.aux]<-length(c(aux,nodes[which(!(nodes%in%aux))]))
    eve.aux[neve.aux,1:Leve.aux[neve.aux]]<-c(aux,nodes[which(!(nodes%in%aux))])
  } # END: Cycle over the eventYES points
# Reorganise eventYES clump labels
  y.eve<-vector(length=p)
  y.eve[1:p]<-NA
  neve.vec<-0
  if (neve.aux>0) {
    for (eve in 1:neve.aux) {
      if (Leve.aux[eve]==0) next
      neve.vec<-neve.vec+1
      Leve.vec[neve.vec]<-Leve.aux[eve]
      eve.vec[neve.vec,1:Leve.vec[neve.vec]]<-eve.aux[eve,1:Leve.aux[eve]]
      y.eve[eve.vec[neve.vec,1:Leve.vec[neve.vec]]]<-neve.vec
    }
  }
  rm(eve.aux,Leve.aux,neve.aux)
# identify eventYES points surrounded only by eventNO points (close enough)
  if (neve.vec>0) {
    for (j in 1:neve.vec) {
      eve.j<-eve.vec[j,1:Leve.vec[j]]
      ix.eY.pts<-which(eve.j %in% eix)
      n.eY<-length(ix.eY.pts)
      frac.eY<-n.eY/Leve.vec[j]
      if (n.eY<n.connected_eveYES & frac.eY<frac.eveYES_in_the_clump) {
        eve.j.eY<-eve.j[ix.eY.pts]
        aux.dist<-(outer(obs$y[eve.j.eY],obs$y[eve.j],FUN="-")**2.+
                   outer(obs$x[eve.j.eY],obs$x[eve.j],FUN="-")**2.)**0.5
        min.dist.from.eY.stn<-min(aux.dist[aux.dist>0])
        if (min.dist.from.eY.stn<dmin.next_eveYES) 
          ydqc.flag[eve.j.eY]<-200
      }
    }
    rm(eve.j,ix.eY.pts)
  }
  if (exists("aux.dist")) rm(aux.dist,min.dist.from.eY.stn,eve.j.eY)
  ydqc.flag
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
    vals<-ncvar_get(defs, varid=varid)
    # close file
    nc_close(defs)
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
  if (is.null(proj4)) {
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

# debug
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

#+
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

#==============================================================================
#  MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN
#==============================================================================
t0<-Sys.time()
#
proj4default<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
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
# titan path (use for the puddle check)
p <- add_argument(p, "--titan_path",
                  help="path to the directory where the TITAN code is",
                  type="character",
                  default=NULL,
                  short="-tip")
#.............................................................................. 
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
   help="offset applied to the input files (one for each provider, default=0)",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-io")
p <- add_argument(p, "--input.negoffset",
                  help="sign for the offsets (1=negative; 0=positive, def=0)",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-ion")
p <- add_argument(p, "--input.cfact",
   help="correction factor applied to the input files (one for each provider, default=1)",
                  type="numeric",
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
p <- add_argument(p, "--steve.code",
                  help="quality code returned in case of STEVE fails",
                  type="numeric",
                  default=9,
                  short="-stevec")
p <- add_argument(p, "--ccrrt.code",
                  help=paste("quality code returned in case of precipitation",
                             "and temperature crosscheck fails"),
                  type="numeric",
                  default=11,
                  short="-ccrrtc")
p <- add_argument(p, "--puddle.code",
             help=paste("quality code returned in case of puddle check fails"),
                  type="numeric",
                  default=12,
                  short="-ccrrtc")
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
# transformation between coordinate reference systems
p <- add_argument(p, "--spatconv",
                  help="flag for conversion of spatial coordinates before running the data quality checks",
                  flag=T,
                  short="-c")
p <- add_argument(p, "--proj4from",
                  help="proj4 string for the original coordinate reference system",
                  type="character",
                  default="+proj=longlat +datum=WGS84",
                  short="-pf")
p <- add_argument(p, "--proj4to",
                  help="proj4 string for the coordinate reference system where the DQC is performed",
                  type="character",
                  default=proj4default,
                  short="-pt")
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
                 help="number of decimal digits for latitude and longitude in the output file",
                 type="numeric",
                 default=5,
                 short="-lldo")
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
                  help="latitude variable name in the output file",
                  type="character",
                  default="lat")
p <- add_argument(p, "--varname.lon",
                  help="character vector, longitude variable name(s) in the input file (default ''lon'')",
                  type="character",
                  nargs=Inf,
                  short="-vlon")
p <- add_argument(p, "--varname.lon.out",
                  help="longitude variable name in the output file",
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
                  help="minimum allowed value [units of the variable specified]",
                  type="numeric",
                  default=-50)
p <- add_argument(p, "--vmax",
                  help="maximum allowed value [units of the variable specified]",
                  type="numeric",
                  default=40)
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
                  help="minimum allowed value [units of the variable specified]",
                  type="numeric",
                  nargs=12,
                  default=c(45,45,40,35,20,15,10,15,15,20,35,45))
p <- add_argument(p, "--vmax.clim",
                  help="maximum allowed value [units of the variable specified]",
                  type="numeric",
                  nargs=12,
                  default=c(20,20,25,25,35,35,40,40,35,30,25,20))
p <- add_argument(p, "--month.clim",
                  help="month (number 1-12)",
                  type="numeric",
                  short="-mC",
                  default=NA)
p <- add_argument(p, "--vminsign.clim",
                  help="minimum allowed value, sign [1=neg, 0=pos]",
                  type="numeric",
                  nargs=12,
                  default=c(1,1,1,1,1,1,1,1,1,1,1,1))
p <- add_argument(p, "--vmaxsign.clim",
                  help="maximum allowed value, sign [1=neg, 0=pos]",
                  type="numeric",
                  nargs=12,
                  default=c(0,0,0,0,0,0,0,0,0,0,0,0))
#.............................................................................. 
# Buddy-check
p <- add_argument(p, "--dr.buddy",
                  help="perform the buddy-check in a dr-by-dr square-box around each observation [m]",
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
                  help=paste("nrow ncol (i.e. number_of_rows number_of_columns).",
                  "SCT is performed independently over several boxes.",
                  "The regular nrow-by-ncol grid is used to define those",
                  "rectangular boxes where the SCT is performed."),
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
p <- add_argument(p, "--laf.file",
                  help="land area fraction file (netCDF in kilometric coordinates)",
                  type="character",
                  default=NULL,
                  short="-lfS")
p <- add_argument(p, "--proj4laf",
                  help="proj4 string for the laf",
                  type="character",
                  default=proj4default,
                  short="-pl")
p <- add_argument(p, "--laf.varname",
                  help="variable name in the netCDF file",
                  type="character",
                  default="land_area_fraction",
                  short="-lfv")
p <- add_argument(p, "--laf.topdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the laf upside down",
                  type="logical",
                  default=FALSE,
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
p <- add_argument(p, "--fast.sct",
                  help=paste("faster spatial consistency test. Allow for",
                      "flagging more than one observation simulataneously"),
                  flag=T,
                  short="-fS")
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
                  default=proj4default,
                  short="-pd")
p <- add_argument(p, "--dem.varname",
                  help="variable name in the netCDF file",
                  type="character",
                  default="altitude",
                  short="-dmv")
p <- add_argument(p, "--dem.topdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the dem upside down",
                  type="logical",
                  default=FALSE,
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
# STEVE - isolated event test
p <- add_argument(p, "--steve",
                  help="do the isolated event test",
                  flag=T,
                  short="-dST")
p <- add_argument(p, "--i.steve",
                  help="number of STEVE iterations",
                  type="integer",
                  default=1,
                  short="-iST")
p <- add_argument(p, "--thres.steve",
                  help=paste("numeric vector with the thresholds used to",
                  " define events (same units of the specified variable).",
                  "Each threshold defines two events: (i) less than the",
                  "threshold and (ii) greater or equal to it."),
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--pmax_lt.steve",
                  help=paste("numeric vector, used for events of the type",
                  "''less than the threshold'' with the maximum number of",
                  "stations defining (in combination with dmax) the",
                  "neighbourhood where to search for",
                  "clumps of connected events (default is 20)"),
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--dmax_lt.steve",
                  help=paste("numeric vector, used for events of the type",
                  "''less than the threshold'' with the maximum distance",
                  "between stations defining (in combination with pmax) the",
                  "neighbourhood where to search for",
                  "clumps of connected events (default is 150000)"),
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--n_lt.steve",
                  help=paste("numeric vector, used for events of the type",
                  "''less than the threshold'' with the minimum expected",
                  "number of connected events. If an event has less than",
                  "this number of connected events, then flag it is as a",
                  "suspect isolated event. The connection between events is",
                  "estimated through a Delaunay triangulation.",
                  "(default is 4)"),
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--frac_lt.steve",
                  help=paste("numeric vector, used for events of the type",
                  "''less than the threshold'' with the maximum expected",
                  "fraction of events in a clump of connected events",
                  "that could include an isolated event.",
                  "If a clump of connected events has more than",
                  "this fraction of events, then it is assumed that",
                  "the situation is too ''patchy'' to flag an isolated event",
                  "with enough confidence (deafult is 0.2)."),
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--dmin_next_lt.steve",
                  help=paste("numeric vector, used for events of the type",
                  "''less than the threshold'' with the maximum expected",
                  "distance between an isolated event and the closest event.",
                  "If a possible isolated event has an other event at a",
                  "distance less than the one defined here then it is not",
                  "considered as an isolated event",
                  "(deafult is 150000)."),
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--pmax_ge.steve",
                  help=paste("numeric vector, used for events of the type",
      "''greater or equal than the threshold'' with the maximum number of",
                  "stations defining (in combination with dmax) the",
                  "neighbourhood where to search for",
                  "clumps of connected events (default is 20)"),
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--dmax_ge.steve",
                  help=paste("numeric vector, used for events of the type",
      "''greater or equal than the threshold'' with the maximum distance",
                  "between stations defining (in combination with pmax) the",
                  "neighbourhood where to search for",
                  "clumps of connected events (default is 150000)"),
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--n_ge.steve",
                  help=paste("numeric vector, used for events of the type",
      "''greater or equal than the threshold'' with the minimum expected",
                  "number of connected events. If an event has less than",
                  "this number of connected events, then flag it is as a",
                  "suspect isolated event. The connection between events is",
                  "estimated through the Delaunay triangulation",
                  "(default is 4)"),
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--frac_ge.steve",
                  help=paste("numeric vector, used for events of the type",
      "''greater or equal than the threshold'' with the maximum expected",
                  "fraction of events in a clump of connected events",
                  "that could include an isolated event.",
                  "If a clump of connected events has more than",
                  "this fraction of events, then it is assumed that",
                  "the situation is too ''patchy'' to flag an isolated event",
                  "with enough confidence (deafult is 0.2)."),
                  type="numeric",
                  default=NA,
                  nargs=Inf)
p <- add_argument(p, "--dmin_next_ge.steve",
                  help=paste("numeric vector, used for events of the type",
      "''greater or equal than the threshold'' with the maximum expected",
                  "distance between an isolated event and the closest event.",
                  "If a possible isolated event has an other event at a",
                  "distance less than the one defined here then it is not",
                  "considered as an isolated event",
                  "(deafult is 100000)."),
                  type="numeric",
                  default=NA,
                  nargs=Inf)
#.............................................................................. 
# precipitation and temperature cross-check
p <- add_argument(p, "--ccrrt",
           help="precipitation (in-situ) and temperature (field) cross-check",
                  flag=T,
                  short="-ccrrtf")
p <- add_argument(p, "--ccrrt.tmin",
                  help="temperature thresholds (vector)",
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
# precipitation and temperature cross-check
p <- add_argument(p, "--puddle",
                  help="puddle check",
                  flag=T,
                  short="-pudf")
p <- add_argument(p, "--dobs_k.puddle",
                  help="in the calculations of the de-correlation lenght scale, consider the observations up to the \"dobs_k\" closest one. Furthermore, the grid resolution of the grid depends on \"dobs_k\"",
                  type="numeric",
                  default=5,
                  short="-putk")
p <- add_argument(p, "--Dh_min.puddle",
                  help="minimum allowed value for the de-correlation lenght scale (in km)",
                  type="numeric",
                  default=5,
                  short="-putd")
p <- add_argument(p, "--i.puddle",
                  help="number of puddle-check iterations",
                  type="integer",
                  default=1,
                  short="-ipudT")
p <- add_argument(p, "--thres.puddle",
                  help="numeric vector with the thresholds used to define events (same units of the specified variable). Each threshold defines two events: (i) less than the threshold and (ii) greater or equal to it.",
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-tpud")
p <- add_argument(p, "--eps2.puddle",
                  help="parameter related to the IDI elaboration (internal OI call). range of values=0.01,...,1. in case of value close to 1, the OI produces a smoothed field. if tha value is close to 0 the OI field tends to fit exactly the observations.",
                  type="numeric",
                  default=0.1,
                  short="-pute2")
p <- add_argument(p, "--n_lt.puddle",
                  help="used for events of the type ''less than the threshold''. minimum acceptable number of events in a clump of connected events. If a clump of connected events has less than this number of events (i.e., it is just a ''puddle'' of cells), then it cannot be properly resolved by the observational network. As a result, the observations causing those small-scale events are assumed to be affected by large representativeness errors and flagged as ''suspect''.",
                  type="integer",
                  default=2,
                  short="-pudnlt")
p <- add_argument(p, "--n_ge.puddle",
                  help="used for events of the type ''greater or equal than the threshold''. minimum acceptable number of events in a clump of connected events. If a clump of connected events has less than this number of events (i.e., it is just a ''puddle'' of cells), then it cannot be properly resolved by the observational network. As a result, the observations causing those small-scale events are assumed to be affected by large representativeness errors and flagged as ''suspect''.",
                  type="integer",
                  default=2,
                  short="-pudnge")
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
p <- add_argument(p, "--doit.steve",
                  help=paste("customize the STEVE application.",comstr),
                  type="numeric",
                  default=NA,
                  nargs=Inf,
                  short="-dost")
p <- add_argument(p, "--doit.puddle",
                  help=paste("customize the puddle check application.",comstr),
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
                  default=proj4default,
                  short="-rrwtp")
p <- add_argument(p, "--t2m.varname",
                  help="air temperature variable name in the netCDF file",
                  type="character",
                  default=NULL,
                  short="-rrwtv")
p <- add_argument(p, "--t2m.topdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the file upside down",
                  type="logical",
                  default=FALSE,
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
                  type="logical",
                  default=FALSE,
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
                  default=proj4default,
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
                  type="logical",
                  default=FALSE,
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
     help="maximum allowed deviation between observation and fg",
                  type="numeric",
                  default=50)
p <- add_argument(p, "--thrpos.fg",
     help="maximum allowed deviation between observation and fg (if obs>fg)",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--thrneg.fg",
     help="maximum allowed deviation between observation and fg (if obs<fg)",
                  type="numeric",
                  default=NA)
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
                  default=proj4default,
                  short="-pfg")
p <- add_argument(p, "--usefg.sct",
         help="use the first-guess field provided as the SCT-background",
                  flag=T)
p <- add_argument(p, "--fg.varname",
                  help="variable name in the netCDF file",
                  type="character",
                  default="land_area_fraction",
                  short="-fgv")
p <- add_argument(p, "--fg.topdown",
                  help="logical, netCDF topdown parameter. If TRUE then turn the fg upside down",
                  type="logical",
                  default=FALSE,
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
                  type="logical",
                  default=FALSE,
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
     help="maximum allowed deviation between observation and fg",
                  type="numeric",
                  default=50)
p <- add_argument(p, "--thrpos.fge",
     help="maximum allowed deviation between observation and fg (if obs>fg)",
                  type="numeric",
                  default=NA)
p <- add_argument(p, "--thrneg.fge",
     help="maximum allowed deviation between observation and fg (if obs<fg)",
                  type="numeric",
                  default=NA)
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
                  default=proj4default,
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
                  type="logical",
                  default=FALSE,
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
                  type="logical",
                  default=FALSE,
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
p <- add_argument(p, "--csd.fge",
  help="coefficient, outlier detection based on ensemble standard deviation",
                  type="numeric",
                  default=5,
                  short="-fgecsd")
p <- add_argument(p, "--infsd.fge",
  help="inflaction coefficient, outlier detection based on ensemble standard deviation",
                  type="numeric",
                  default=3,
                  short="-fgeisd")
p <- add_argument(p, "--ciqr.fge",
  help="coefficient, outlier detection based on ensemble standard deviation",
                  type="numeric",
                  default=1.5,
                  short="-fgeciqr")
p <- add_argument(p, "--infiqr.fge",
  help="inflaction coefficient, outlier detection based on ensemble standard deviation",
                  type="numeric",
                  default=3,
                  short="-fgeiiqr")
p <- add_argument(p, "--sdmin.fge",
                  help="ensemble standard deviation, minimum value",
                  type="numeric",
                  default=NULL,
                  short="-fgesdmn")
p <- add_argument(p, "--iqrmin.fge",
                  help="ensemble inter-quartile range, minimum value",
                  type="numeric",
                  default=NULL,
                  short="-fgesdmn")
#.............................................................................. 
# Timestamp valid for all the netcdf files
p <- add_argument(p, "--timestamp",
                  help="timestamp, valid for all the netCDF file (YYYYMMDDHH00)",
                  type="character",
                  default=NA,
                  short="-t")
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
# set offsets and correction factors
if (is.na(argv$input.offset)) argv$input.offset<-rep(0,nfin)
if (is.na(argv$input.negoffset)) argv$input.negoffset<-rep(0,nfin)
if (is.na(argv$input.cfact)) argv$input.cfact<-rep(1,nfin)
if (is.na(argv$input.negcfact)) argv$input.negcfact<-rep(0,nfin)
argv$input.offset<-argv$input.offset*(-1)**argv$input.negoffset
argv$input.cfact<-argv$input.cfact*(-1)**argv$input.negcfact
# check variable
if (!(argv$variable %in% c("T","RH","RR","SD"))) {
  print("variable must be one of T, RH, RR, SD")
  quit(status=1)
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
      if (is.na(argv$fg.dimnames)) {
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
      if (is.na(argv$fg.dimnames)) {
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
      if (is.na(argv$fg.dimnames)) {
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
  print("ERROR: \"--spatconv\" (-c) option must be used onthe command line")
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
# no more needed?
#  if ((argv$dem | argv$dem.fill) &
#      (argv$proj4fg!=argv$proj4dem |
#       !is.character(argv$proj4fg))) {
#    print("ERROR: anomalies found in the proj4 strings:")
#    print(paste("proj4 fg=",argv$proj4fg))
#    print(paste("proj4 dem=",argv$proj4dem))
#    quit(status=1)
#  }
  if ( argv$variable=="T" &
       !file.exists(argv$fg.demfile)) {
    print("ERROR: for temperature, a digital elevation model must be specified together with a first-guess file")
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
# observation representativeness
if (any(is.na(argv$mean.corep)) | 
    any(is.na(argv$min.corep))  | 
    any(is.na(argv$max.corep)) ) {
  print("++WARNING")
  print("parameters related to the coefficient of observation representativeness are not properly specified")
  print("--mean.corep --min.corep and --max.corep should be specified")
  print("Because they are not specified then it is assumed that the coefficient of observation representativeness is not considered an interesting output anf the corep parameters will be set to default values (min.corep=0.9 mean.corep=1 max.corep=1.1")
  argv$min.corep<-0.9
  argv$mean.corep<-1
  argv$max.corep<-1.1
}
if (length(argv$min.corep)!=nfin) 
  argv$eps2.sct<-rep(argv$min.corep[1],length=nfin)
if (length(argv$mean.corep)!=nfin) 
  argv$eps2.sct<-rep(argv$mean.corep[1],length=nfin)
if (length(argv$max.corep)!=nfin) 
  argv$eps2.sct<-rep(argv$max.corep[1],length=nfin)
# precip and temperature crosscheck
if (length(argv$ccrrt.tmin)!=nfin) 
  argv$ccrrt.tmin<-rep(argv$ccrrt.tmin[1],length=nfin)
aux<-vector(length=nfin,mode="numeric")
for (i in 1:nfin) aux[i]<-as.numeric(gsub("_","-",argv$ccrrt.tmin[i]))
argv$ccrrt.tmin<-aux
rm(aux)
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
  argv$thr.sct<-16
}
if ( length(argv$thr.sct)!=nfin & 
     is.na(argv$thrneg.sct[1])) 
  argv$thr.sct<-rep(argv$thr.sct[1],length=nfin)
# eps2
if (any(is.na(argv$eps2.sct))) {
  print("++ WARNING")
  print("eps2.sct should be specified and it must not contain NAs")
  print(" because either it has not been specified or it has been set to NA, ")
  print(" then TITAN will use the default value of 0.5")
  argv$eps2.sct<-0.5
}
if (length(argv$eps2.sct)!=nfin) 
  argv$eps2.sct<-rep(argv$eps2.sct[1],length=nfin)
# doit flags
if (any(is.na(argv$doit.buddy))) argv$doit.buddy<-rep(1,length=nfin)
if (any(is.na(argv$doit.sct))) argv$doit.sct<-rep(1,length=nfin)
if (any(is.na(argv$doit.clim))) argv$doit.clim<-rep(1,length=nfin)
if (any(is.na(argv$doit.dem))) argv$doit.dem<-rep(1,length=nfin)
if (any(is.na(argv$doit.isol))) argv$doit.isol<-rep(1,length=nfin)
if (any(is.na(argv$doit.steve))) argv$doit.steve<-rep(1,length=nfin)
if (any(is.na(argv$doit.fg))) argv$doit.fg<-rep(1,length=nfin)
if (any(is.na(argv$doit.fge))) argv$doit.fge<-rep(1,length=nfin)
if (any(is.na(argv$doit.puddle))) argv$doit.puddle<-rep(1,length=nfin)
if (any(!(argv$doit.buddy %in% c(0,1,2)))) {
  print("doit.buddy must contain only 0,1,2")
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
if (any(!(argv$doit.steve %in% c(0,1,2)))) {
  print("doit.steve must contain only 0,1,2")
  quit(status=1)
}
if (any(!(argv$doit.fg %in% c(0,1,2)))) {
  print("doit.fg must contain only 0,1,2")
  quit(status=1)
}
if (any(!(argv$doit.puddle %in% c(0,1,2)))) {
  print("doit.puddle must contain only 0,1,2")
  quit(status=1)
}
# set the thresholds for the plausibility check
if (!is.na(argv$tmin) & is.na(argv$vmin)) argv$vmin<-argv$tmin
if (!is.na(argv$tmax) & is.na(argv$vmax)) argv$vmax<-argv$tmax
if (any(is.na(argv$vmin.clim)) & !any(is.na(argv$tmin.clim))) 
  argv$vmin.clim<-argv$tmin.clim
if (any(is.na(argv$vmax.clim)) & !any(is.na(argv$tmax.clim))) 
  argv$vmax.clim<-argv$tmax.clim
argv$vmin<-argv$vmin*(-1)**argv$vminsign
argv$vmax<-argv$vmax*(-1)**argv$vmaxsign
argv$vmin.clim<-argv$vmin.clim*(-1)**argv$vminsign.clim
argv$vmax.clim<-argv$vmax.clim*(-1)**argv$vmaxsign.clim
# puddle checks
if (argv$puddle) {
  if (!file.exists(file.path(argv$titan_path,"oi_rr","oi_rr_first.so"))) {
    print("ERROR: puddle check, file not found")
    print(file.path(argv$titan_path,"oi_rr","oi_rr_first.so"))
    quit(status=1)
  }
  if (!file.exists(file.path(argv$titan_path,"oi_rr","oi_rr_fast.so"))) {
    print("ERROR: puddle check, file not found")
    print(file.path(argv$titan_path,"oi_rr","oi_rr_fast.so"))
    quit(status=1)
  }
  # load external C functions
  dyn.load(file.path(argv$titan_path,"oi_rr","oi_rr_first.so"))
  dyn.load(file.path(argv$titan_path,"oi_rr","oi_rr_fast.so"))
  if (any(is.na(argv$thres.puddle))) {
    print("ERROR: puddle check, thres.puddle must be specified")
    quit(status=1)
  }
}
# STEVE checks
if (argv$steve) {
  suppressPackageStartupMessages(library("tripack"))
  if (any(is.na(argv$thres.steve))) {
    print("ERROR: STEVE, thres.steve must be specified")
    quit(status=1)
  }
  if (any(is.na(argv$pmax_lt.steve)))
    argv$pmax_lt.steve<-rep(20,length(argv$thres.steve))
  if (any(is.na(argv$dmax_lt.steve))) 
    argv$dmax_lt.steve<-rep(150000,length(argv$thres.steve))
  if (any(is.na(argv$n_lt.steve))) 
    argv$n_lt.steve<-rep(4,length(argv$thres.steve))
  if (any(is.na(argv$frac_lt.steve))) 
    argv$frac_lt.steve<-rep(0.2,length(argv$thres.steve))
  if (any(is.na(argv$dmin_next_lt.steve))) 
    argv$dmin_next_lt.steve<-rep(150000,length(argv$thres.steve))
  if (any(is.na(argv$pmax_ge.steve))) 
    argv$pmax_ge.steve<-rep(20,length(argv$thres.steve))
  if (any(is.na(argv$dmax_ge.steve))) 
    argv$dmax_ge.steve<-rep(150000,length(argv$thres.steve))
  if (any(is.na(argv$n_ge.steve))) 
    argv$n_ge.steve<-rep(4,length(argv$thres.steve))
  if (any(is.na(argv$frac_ge.steve))) 
    argv$frac_ge.steve<-rep(0.2,length(argv$thres.steve))
  if (any(is.na(argv$dmin_next_ge.steve))) 
    argv$dmin_next_ge.steve<-rep(100000,length(argv$thres.steve))
  if (length(argv$thres.steve)!=length(argv$pmax_lt.steve)) {
    print("ERROR: STEVE, the number of input parameters differ (pmax_lt.steve)")
    quit(status=1)
  }
  if (length(argv$thres.steve)!=length(argv$dmax_lt.steve)) {
    print("ERROR: STEVE, the number of input parameters differ (dmax_lt.steve)")
    quit(status=1)
  }
  if (length(argv$thres.steve)!=length(argv$n_lt.steve)) {
    print("ERROR: STEVE, the number of input parameters differ (n_lt.steve)")
    quit(status=1)
  }
  if (length(argv$thres.steve)!=length(argv$frac_lt.steve)) {
    print("ERROR: STEVE, the number of input parameters differ (frac_lt.steve)")
    quit(status=1)
  }
  if (length(argv$thres.steve)!=length(argv$dmin_next_lt.steve)) {
    print("ERROR: STEVE, the number of input parameters differ (dmin_next_lt.steve)")
    quit(status=1)
  }
  if (length(argv$thres.steve)!=length(argv$pmax_ge.steve)) {
    print("ERROR: STEVE, the number of input parameters differ (pmax_ge.steve)")
    quit(status=1)
  }
  if (length(argv$thres.steve)!=length(argv$dmax_ge.steve)) {
    print("ERROR: STEVE, the number of input parameters differ (dmax_ge.steve)")
    quit(status=1)
  }
  if (length(argv$thres.steve)!=length(argv$n_ge.steve)) {
    print("ERROR: STEVE, the number of input parameters differ (n_ge.steve)")
    quit(status=1)
  }
  if (length(argv$thres.steve)!=length(argv$frac_ge.steve)) {
    print("ERROR: STEVE, the number of input parameters differ (frac_ge.steve)")
    quit(status=1)
  }
  if (length(argv$thres.steve)!=length(argv$dmin_next_ge.steve)) {
    print("ERROR: STEVE, the number of input parameters differ (dmin_next_ge.steve)")
    quit(status=1)
  }
}
# wind-induced undercatch of precipitation, check consistency of inputs
if (argv$rr.wcor & argv$variable!="RR") {
  print(paste("ERROR: wind-induced correction for undercatch is implemented",
              "for precipitation only"))
  quit(status=1)
}
# set the timestamp
if (!is.na(argv$timestamp)) {
  argv$fg.t<-argv$timestamp
  argv$fge.t<-argv$timestamp
  argv$wind.t<-argv$timestamp
  argv$t2m.t<-argv$timestamp
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
  # varidx is used also in the output session
  varidx<-match(c(argv$varname.lat[f],
                  argv$varname.lon[f],
                  argv$varname.elev[f],
                  argv$varname.value[f]),
                names(datain))
  if (any(is.na(varidx))) {
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
  datatmp<-data.frame(datain[,varidx])
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
  if (any(!is.na(argv$blacklist.idx)) & any(argv$blacklist.fidx==f)) {
    aux[argv$blacklist.idx[which(argv$blacklist.fidx==f)]]<-argv$black.code  
  }
  if (any(!is.na(argv$blacklist.lat)) & any(argv$blacklist.fll==f)) {
    out<-apply(cbind(argv$blacklist.lon[argv$blacklist.fll==f],
                     argv$blacklist.lat[argv$blacklist.fll==f])
               ,FUN=setCode_lonlat,MARGIN=1,code=argv$black.code)
    rm(out)
  }
  if (any(!is.na(argv$keeplist.idx)) & any(argv$keeplist.fidx==f)) {
    aux[argv$keeplist.idx[which(argv$keeplist.fidx==f)]]<-argv$keep.code  
  }
  if (any(!is.na(argv$keeplist.lat)) & any(argv$keeplist.fll==f)) {
    out<-apply(cbind(argv$keeplist.lon[argv$keeplist.fll==f],
                     argv$keeplist.lat[argv$keeplist.fll==f])
               ,FUN=setCode_lonlat,MARGIN=1,code=argv$keep.code)
    rm(out)
  }
  if (first) {
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
      dataopt<-array(data=NA,dim=c(ndatatmp,
                                   length(argv$varname.opt)))
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
      if (any(varidx.opt.check!=varidx.opt | 
              (is.na(varidx.opt.check) & !is.na(varidx.opt)) |
              (is.na(varidx.opt) & !is.na(varidx.opt.check)) )) {
        print("ERROR the header of file")
        print(argv$input.files[f])
        print("is different from the header of the first file")
        print(argv$input.files[1])
        quit(status=1)
      }
      dataopt<-rbind(dataopt,
        datain[,varidx.opt[which(!is.na(varidx.opt))],drop=F])
    }
  }
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
  meta<-!is.na(data$lat[ix]) & 
        !is.na(data$lon[ix]) &
        !is.na(z[ix]) & 
        z[ix]>=argv$zmin & 
        z[ix]<=argv$zmax &
        !is.na(data$value[ix]) 
  if (any(!meta)) dqcflag[ix[which(!meta)]]<-argv$nometa.code
} else {
  print("no valid observations left, no metadata check")
}
if (argv$verbose | argv$debug) {
  print("test for no metdata")
  print(paste("# observations lacking metadata and/or NAs=",length(which(dqcflag==argv$nometa.code & !is.na(argv$nometa.code)))))
  print(paste("  # NAs             =",length(which(is.na(data$value[ix])))))
  print(paste("  # lon-lat missing =",length(which(is.na(data$lat[ix]) | 
                                                  is.na(data$lon[ix])))))
  print(paste("  # z missing       =",length(which(is.na(z[ix]))))) 
  print(paste("  # z out of range  =",length(which(!is.na(z[ix]) & 
                                    (z[ix]<argv$zmin | z[ix]>argv$zmax))))) 
  print("+---------------------------------+")
}
rm(ix)
if (exists("meta")) rm(meta)
#
#-----------------------------------------------------------------------------
# coordinate transformation
if (argv$spatconv) {
  if (argv$debug) print("conversion of spatial coordinates")
  coord<-SpatialPoints(cbind(data$lon,data$lat),
                       proj4string=CRS(argv$proj4from))
  coord.new<-spTransform(coord,CRS(argv$proj4to))
  xy.new<-coordinates(coord.new)
  x<-round(xy.new[,1],0)
  y<-round(xy.new[,2],0)
  xp<-expand.grid(c(argv$lonmin,argv$lonmax),c(argv$latmin,argv$latmax))
  coord<-SpatialPoints(xp,
                       proj4string=CRS(argv$proj4from))
  coord.new<-spTransform(coord,CRS(argv$proj4to))
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
  if (argv$verbose | argv$debug)
    print("read digital elevation model")
  ti<-nc4.getTime(argv$dem.file)
  raux<-try(nc4in(nc.file=argv$dem.file,
                  nc.varname=argv$dem.varname,
                  topdown=argv$dem.topdown,
                  out.dim=list(ndim=argv$dem.ndim,
                               tpos=argv$dem.tpos,
                               epos=NULL,
                               names=argv$dem.dimnames),
                  proj4=argv$proj4dem,
                  nc.proj4=list(var=argv$dem.proj4_var,
                                att=argv$dem.proj4_att),
                  selection=list(t=ti[1],e=NULL)))
  if (is.null(raux)) {
    print("ERROR while reading file:")
    print(argv$dem.file)
    quit(status=1)
  }
  rdem<-raux$stack
  rm(raux,ti)
  if (argv$proj4dem!=argv$proj4to) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4from))
    coord.new<-spTransform(coord,CRS(argv$proj4dem))
    xy.tmp<-coordinates(coord.new)
    zdem<-extract(rdem,xy.tmp)
    rm(coord,coord.new,xy.tmp)
  } else {
    zdem<-extract(rdem,cbind(x,y))
  }
  # fill missing elevation with dem
  if (argv$dem.fill) {
    iz<-which( (is.na(z) & !is.na(zdem)) | 
               (!is.na(zdem) & (z<argv$zmin | z>argv$zmax)) )
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
  if (argv$proj4laf!=argv$proj4to) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4from))
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
if (argv$debug) {
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
               r=rt2m,x=x,y=y,proj4=argv$proj4to,proj4plot=argv$proj4t2m)
  coord<-SpatialPoints(cbind(data$lon,data$lat),
                       proj4string=CRS(argv$proj4from))
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
  if (argv$proj4t2m!=argv$proj4to) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4from))
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
               r=rt2mdem,x=x,y=y,proj4=argv$proj4to,proj4plot=argv$proj4t2m)
  t2m<-t2m+argv$gamma.standard*(z-zt2mdem)
  # cross-check
  for (f in 1:nfin) 
    dqcflag[which(data$prid==argv$prid[f] & 
                  t2m<argv$ccrrt.tmin[f])]<-argv$ccrrt.code
  #
  if (argv$verbose | argv$debug) {
    print("precipitaton and temperature  crosscheck")
    print(paste("temp thresholds =",toString(argv$ccrrt.tmin)))
    print(paste("# suspect observations=",
          length(which(dqcflag==argv$ccrrt.code))))
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
                       proj4string=CRS(argv$proj4from))
  coord.new<-spTransform(coord,CRS(argv$proj4t2m))
  xy.tmp<-coordinates(coord.new)
  t2m<-extract(rt2m,xy.tmp,method="bilinear")
  t2m<-argv$t2m.offset*(-1)**(argv$t2m.negoffset)+
       t2m*argv$t2m.cfact*(-1)**(argv$t2m.negcfact)
  if (argv$debug)
    plot_debug(ff=file.path(argv$debug.dir,"rrwcor_t2m.png"),
               r=rt2m,x=x,y=y,proj4=argv$proj4to,proj4plot=argv$proj4t2m)
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
  if (argv$proj4t2m!=argv$proj4to) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4from))
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
               r=rt2mdem,x=x,y=y,proj4=argv$proj4to,proj4plot=argv$proj4t2m)
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
                 r=ru,x=x,y=y,proj4=argv$proj4to,proj4plot=argv$proj4wind)
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
                 r=rv,x=x,y=y,proj4=argv$proj4to,proj4plot=argv$proj4wind)
    rwind<-rv
    rwind[]<-sqrt(getValues(ru)**2+getValues(rv)**2)
    rm(ru,rv)
  #  case of wrong/no wind varname in the file
  } else {
    print(paste("ERROR precipitation correction for the wind undercatch:",
                " wind varname has not been specified"))
    quit(status=1)
  }
  if (argv$proj4wind!=argv$proj4to) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4from))
    coord.new<-spTransform(coord,CRS(argv$proj4wind))
    ws10m<-extract(rwind,cbind(x,y),method="bilinear")
    xy.tmp<-coordinates(coord.new)
    rm(coord,coord.new,xy.tmp)
  } else {
    ws10m<-extract(rwind,cbind(x,y),method="bilinear")
  }
  if (argv$debug)
    plot_debug(ff=file.path(argv$debug.dir,"rrwcor_ws.png"),
               r=rwind,x=x,y=y,proj4=argv$proj4to,proj4plot=argv$proj4wind)
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
    ix<-which(data$value>=0 & !is.na(data$value) &
              is.na(dqcflag))
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
      print("ERROR while reading file:")
      print(argv$fg.file)
      quit(status=1)
    }
    rfg<-raux$stack
  }
  rm(raux,ti,fg.e,fg.epos)
  # radar fg, data quality control 
  if (argv$fg.type=="radar") {
    if (argv$verbose | argv$debug)
      print("Read radar and do the radar-DQC")
    t0a<-Sys.time()
    suppressPackageStartupMessages(library("igraph"))
    dfg<-getValues(rfg)
    # a. remove not plausible values
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
    raux<-rfg
    daux<-boxcox(x=dfg,lambda=0.5)
    raux[]<-daux
    radardqc.outl.fact<-51
    # aggregate over boxes of fact x fact cells, take the mean
    avg<-getValues(
          crop(
           focal(
            extend(raux,c(radardqc.outl.fact,radardqc.outl.fact)),
            w=matrix(1,radardqc.outl.fact,radardqc.outl.fact),fun=mean,na.rm=T),
          raux) )
    stdev<-getValues(
            crop(
             focal(
              extend(raux,c(radardqc.outl.fact,radardqc.outl.fact)),
              w=matrix(1,radardqc.outl.fact,radardqc.outl.fact),fun=sd,na.rm=T),
            raux) )
    ix<-which(stdev>0)
    # outliers are defined as in Lanzante,1997: abs(value-mean)/st.dev > 5
    suspect<-which((abs(daux[ix]-avg[ix])/stdev[ix])>5) 
    if (length(suspect)>0) dfg[ix[suspect]]<-NA
    rfg[]<-dfg
    rm(raux,daux,avg,stdev,ix,suspect,dfg)
    if (argv$radarout) rrad<-rfg
    t1a<-Sys.time()
    print(paste("time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
  }
  if (argv$proj4fg!=argv$proj4to) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4from))
    coord.new<-spTransform(coord,CRS(argv$proj4fg))
    xy.tmp<-coordinates(coord.new)
    fg<-extract(rfg,xy.tmp,method="bilinear")
    rm(coord,coord.new,xy.tmp)
  } else {
    fg<-extract(rfg,cbind(x,y),method="bilinear")
  }
  fg<-argv$fg.offset*(-1)**(argv$fg.negoffset)+
      fg*argv$fg.cfact*(-1)**(argv$fg.negcfact)
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
    if (argv$proj4fg!=argv$proj4to) {
      coord<-SpatialPoints(cbind(data$lon,data$lat),
                           proj4string=CRS(argv$proj4from))
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
    proj4string(xy.tmp)<-CRS(argv$proj4to)
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
    if (argv$proj4fge!=argv$proj4to) {
      coord<-SpatialPoints(cbind(data$lon,data$lat),
                           proj4string=CRS(argv$proj4from))
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
  ei<-nc4.getDim(argv$fge.file,varid=argv$fge.dimnames[argv$fge.epos])
  if (is.na(argv$fge.epos)) {
    print("ERROR fge.epos must have a valid value")
    quit(status=1)
  }
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
                      selection=list(t=c(tminus1h,argv$fge.t),e=ei[ens])))
      if (is.null(raux)) {
        print("ERROR while reading file:")
        print(argv$fge.file)
        quit(status=1)
      }
      rfge<-raster(raux$stack,"layer.2")-raster(raux$stack,"layer.1")
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
    # for ''output radar data'' usage: save the ensemble median
    if (argv$radarout) {
      dfge<-getValues(rfge)
      if (any(!is.na(dfge))) {
        if (first) {
          gdata<-array(data=NA,dim=c(length(dfge),length(ei)))
          rfge.grid<-rfge
          first<-F
        }
        gdata[,ens]<-getValues(rfge)
      }
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
    if (argv$proj4fge!=argv$proj4to) {
      coord<-SpatialPoints(cbind(data$lon,data$lat),
                           proj4string=CRS(argv$proj4from))
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
    } else if (argv$variable=="RR") {
      edata[,ens]<-boxcox(x=edata[,ens],lambda=argv$boxcox.lambda)
    }
    rm(rfge)
  } # end of cycle over ensemble members
  # compute aggregate quantities
  fge.mu<-rowMeans(edata,na.rm=T)
  fge.sd<-apply(edata,MARGIN=1,FUN=function(x){sd(x,na.rm=T)})
  fge.q50<-apply(edata,MARGIN=1,FUN=function(x){as.numeric(quantile(x,probs=0.5,na.rm=T))})
  fge.q25<-apply(edata,MARGIN=1,FUN=function(x){as.numeric(quantile(x,probs=0.25,na.rm=T))})
  fge.q75<-apply(edata,MARGIN=1,FUN=function(x){as.numeric(quantile(x,probs=0.75,na.rm=T))})
  if (argv$debug)
    save.image(file.path(argv$debug.dir,"input_data_fge.RData")) 
  rm(edata)
  #
  if (argv$radarout) {
    fge.gq50<-apply(gdata,MARGIN=1,FUN=function(x){as.numeric(quantile(x,probs=0.5,na.rm=T))})
    rm(gdata)
  }
  # debug
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
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
  meta<-!is.na(data$lat[ix]) & 
        !is.na(data$lon[ix]) &
        !is.na(z[ix]) & 
        z[ix]>=argv$zmin & 
        z[ix]<=argv$zmax &
        !is.na(data$value[ix]) 
  if (any(!meta)) dqcflag[ix[which(!meta)]]<-argv$nometa.code
} else {
  print("no valid observations left, no metadata check")
}


if (argv$verbose | argv$debug) {
  print("test for no metdata (2nd round)")
  print(paste("# observations lacking metadata and/or NAs=",length(which(dqcflag==argv$nometa.code & !is.na(argv$nometa.code)))))
  print(paste("  # NAs             =",length(which(is.na(data$value[ix])))))
  print(paste("  # lon-lat missing =",length(which(is.na(data$lat[ix]) | 
                                                  is.na(data$lon[ix])))))
  print(paste("  # z missing       =",length(which(is.na(z[ix]))))) 
  print(paste("  # z out of range  =",length(which(!is.na(z[ix]) & 
                                    (z[ix]<argv$zmin | z[ix]>argv$zmax))))) 
  print("+---------------------------------+")
}
#if (argv$verbose | argv$debug) {
#  print("test for no metdata (2nd round)")
##  print(paste(data$lat[which(!meta)],data$lon[which(!meta)],
##              z[which(!meta)],data$value[which(!meta)]))
#  print(paste("# observations lacking metadata (tot, 1st+2nd round)=",
#         length(which(dqcflag==argv$nometa.code))))
#  print("+---------------------------------+")
#}
rm(ix)
if (exists("meta")) rm(meta)
if (argv$debug) 
  save.image(file.path(argv$debug.dir,"dqcres_meta.RData")) 
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
  print(paste("# <min=",length(which(dqcflag==argv$p.code & data$value<argv$vmin))))
  print(paste("# >max=",length(which(dqcflag==argv$p.code & data$value>argv$vmax))))
  print(paste("# suspect observations=",length(which(dqcflag==argv$p.code))))
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
# buddy check 
#  compare each observation against the average of neighbouring observations 
# NOTE: keep-listed stations are used but they canNOT be flagged here
if (argv$verbose | argv$debug) nprev<-0
# set doit vector
doit<-vector(length=ndata,mode="numeric")
doit[]<-NA
for (f in 1:nfin) doit[data$prid==argv$prid[f]]<-argv$doit.buddy[f]
# test
print(paste0("buddy-check (",argv$buddy.code,")"))
for (i in 1:argv$i.buddy) {
  # use only (probably) good observations with doit!=0
  ix<-which( (is.na(dqcflag) | dqcflag==argv$keep.code) &
             doit!=0 )
  t0a<-Sys.time()
  if (length(ix)>0) {
    # define global 1D vector used in statSpat (1D for fast access)
    xtot<-x[ix]
    ytot<-y[ix]
    ztot<-as.numeric(z[ix])
    if (argv$variable=="RR") {
      ttot<-boxcox(x=data$value[ix],lambda=argv$boxcox.lambda)
    } else {
      ttot<-data$value[ix]
    }
    # apply will loop over this 4D array
    ixyzt_tot<-cbind(1:length(xtot),xtot,ytot,ztot,ttot)
    stSp_3km<-apply(ixyzt_tot,FUN=statSpat,MARGIN=1,drmin=argv$dr.buddy)
    # probability of gross error
    pog<-abs(ttot-stSp_3km[3,])/stSp_3km[4,]
    # suspect if: 
    sus<-which( (pog>argv$thr.buddy & 
                 stSp_3km[1,]>argv$n.buddy & 
                 stSp_3km[2,]<argv$dz.buddy &
                 is.na(dqcflag[ix])) &
                 doit[ix]==1 )
    # set dqcflag
    if (length(sus)>0) dqcflag[ix[sus]]<-argv$buddy.code
  } else {
    print("no valid observations left, no buddy check")
  }
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    print(paste("buddy-check, iteration=",i,
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    ncur<-length(which(dqcflag==argv$buddy.code))
    print(paste("# suspect observations=",ncur-nprev))
    nprev<-length(which(dqcflag==argv$buddy.code))
  }
}
rm(doit)
if (argv$debug) 
  save.image(file.path(argv$debug.dir,"dqcres_buddy.RData")) 
if (argv$verbose | argv$debug) 
  print("+---------------------------------+")
if (exists("stSp_3km")) rm(stSp_3km)
if (exists("sus")) rm(sus)
if (exists("pog")) rm(pog)
if (exists("xtot")) rm(xtot)
if (exists("ytot")) rm(ytot)
if (exists("ztot")) rm(ztot)
if (exists("ixyzt_tot")) rm(ixyzt_tot)
#
#-----------------------------------------------------------------------------
# STEVE - isolated event test (YES/NO)
if (argv$steve) {
  if (argv$verbose | argv$debug) {
    nprev<-0
    print(paste0("STEVE (",argv$steve.code,")"))
  }
  # set doit vector
  doit<-vector(length=ndata,mode="numeric")
  doit[]<-NA
  for (f in 1:nfin) doit[data$prid==argv$prid[f]]<-argv$doit.steve[f]
  # test
  for (i in 1:argv$i.steve) {
    t0a<-Sys.time()
    for (j in 1:length(argv$thres.steve)) {
      # use only (probably) good observations
      ix<-which((is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0)
      if (length(ix)>0) {
        obs<-data.frame(x[ix],y[ix],data$value[ix])
        aux<-steve(obs=data.frame(x=x[ix],y=y[ix],yo=data$value[ix]),
                   thres=argv$thres.steve[j],
                   gt_or_lt="lt",
                   pmax=argv$pmax_lt.steve[j],
                   dmax=argv$dmax_lt.steve[j],
                   n.sector=16,
                   n.connected_eveYES=argv$n_lt.steve[j],
                   frac.eveYES_in_the_clump=argv$frac_lt.steve[j],
                   dmin.next_eveYES=argv$dmin_next_lt.steve[j])
        sus<-which(aux!=0 & is.na(dqcflag[ix]) & doit[ix]==1)
        # set dqcflag
        if (length(sus)>0) dqcflag[ix[sus]]<-argv$steve.code
      } else {
        print("no valid observations left, no STEVE")
      }
      ix<-which((is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0)
      if (length(ix)>0) {
        obs<-data.frame(x[ix],y[ix],data$value[ix])
        aux<-steve(obs=data.frame(x=x[ix],y=y[ix],yo=data$value[ix]),
                   thres=argv$thres.steve[j],
                   gt_or_lt="ge",
                   pmax=argv$pmax_ge.steve[j],
                   dmax=argv$dmax_ge.steve[j],
                   n.sector=16,
                   n.connected_eveYES=argv$n_ge.steve[j],
                   frac.eveYES_in_the_clump=argv$frac_ge.steve[j],
                   dmin.next_eveYES=argv$dmin_next_ge.steve[j])
        sus<-which(aux!=0 & is.na(dqcflag[ix]) & doit[ix]==1)
        # set dqcflag
        if (length(sus)>0) dqcflag[ix[sus]]<-argv$steve.code
      } else {
        print("no valid observations left, no STEVE")
      }
    } # end of loop over threshold that define events
    if (argv$verbose | argv$debug) {
      t1a<-Sys.time()
      print(paste("STEVE, iteration=",i,
                  "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
      ncur<-length(which(dqcflag==argv$steve.code))
      print(paste("# suspect observations=",ncur-nprev))
      nprev<-ncur
    }
  }
  rm(doit)
  if (argv$verbose | argv$debug) 
    print("+---------------------------------+")
  if (argv$debug) 
    save.image(file.path(argv$debug.dir,"dqcres_steve.RData")) 
}
#
#-----------------------------------------------------------------------------
# puddle check (first round)
if (argv$puddle) {
  if (argv$verbose | argv$debug) {
    nprev<-0
    print(paste0("puddle-check (",argv$puddle.code,")"))
  }
  # set doit vector
  doit<-vector(length=ndata,mode="numeric")
  doit[]<-NA
  for (f in 1:nfin) doit[data$prid==argv$prid[f]]<-argv$doit.puddle[f]
  ix<-which((is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0)
  ptmp<-length(ix)
  if (ptmp<1) {
    print("puddle check: no valid observations left, no test")
    mask_puddle<-integer(0)
  } else {
    if (!exists("mask_puddle")) {
      # prepare mask with gridpoints where the calculation can be done
      t00<-Sys.time()
      xytmp<-cbind(x[ix],y[ix])
      Dh_puddle<-dobs_fun(obs=data.frame(x=x[ix],y=y[ix]),
                          k=argv$dobs_k.puddle)
      Dh_puddle<-max(argv$Dh_min.puddle,(Dh_puddle/1000))
      rm(xytmp)
      # define raster grid
      grid_res<-as.integer((Dh_puddle/argv$dobs_k.puddle)*1000)
      if (argv$verbose | argv$debug) {
        print(paste("Dh (to the closest",argv$dobs_k.puddle,"neighbour)=",
              round(Dh_puddle,1),"km (grid res=",round(grid_res,0),"m)"))
      }
      rgrid_puddle<-raster(ext=e,resolution=grid_res)
      rgrid_puddle[]<-NA
      xygrid_puddle<-xyFromCell(rgrid_puddle,1:ncell(rgrid_puddle))
      xgrid_puddle<-xygrid_puddle[,1]
      ygrid_puddle<-xygrid_puddle[,2]
      rm(xygrid_puddle)
      # NOTE: z is not used here, but we are ready to introduce it for future versions
      zgrid_puddle<-rep(0,length(ygrid_puddle))
      D<-exp(-0.5*(((outer(x[ix],x[ix],FUN="-")**2.+
                     outer(y[ix],y[ix],FUN="-")**2.)**0.5/1000.)/Dh_puddle)**2.)
      diag(D)<-diag(D)+argv$eps2.puddle
#      t00b<-Sys.time()
      InvD<-chol2inv(chol(D))
#      t01b<-Sys.time()
#      print(paste("InvD time",round(t01b-t00b,1),attr(t01b-t00b,"unit")))
      rm(D)
#      t00b<-Sys.time()
      xidi<-OI_RR_fast(yo.sel=rep(1,ptmp),
                       yb.sel=rep(0,ptmp),
                       xb.sel=rep(0,length(xgrid_puddle)),
                       xgrid.sel=xgrid_puddle,
                       ygrid.sel=ygrid_puddle,
                       zgrid.sel=zgrid_puddle,
                       VecX.sel=x[ix],
                       VecY.sel=y[ix],
                       VecZ.sel=rep(0,ptmp),
                       Dh.cur=Dh_puddle,
                       Dz.cur=1000000)
#      t01b<-Sys.time()
#      print(paste("OI_RR_fast",round(t01b-t00b,1),attr(t01b-t00b,"unit")))
      rm(InvD)
      mask_puddle<-which(xidi>=0.00001)
      t11<-Sys.time()
      if (argv$verbose | argv$debug) {
        print(paste("prepare mask, time=",round(t11-t00,1),attr(t11-t00,"unit")))
        print(paste("tot points=",length(xgrid_puddle)))
        print(paste("masked points=",length(mask_puddle)))
      }
    } # end "if (!exists("mask_puddle"))"
  }
  if (length(mask_puddle)==0) {
    print("puddle check: no test")
  } else {
    xgrid_puddle<-xgrid_puddle[mask_puddle]
    ygrid_puddle<-ygrid_puddle[mask_puddle]
    zgrid_puddle<-zgrid_puddle[mask_puddle]
    # test
    for (i in 1:argv$i.puddle) {
      t0a<-Sys.time()
      for (j in 1:length(argv$thres.puddle)) {
        # use only (probably) good observations
        ix<-which((is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0)
        if (length(ix)>0) {
          aux<-puddle(obs=data.frame(x=x[ix],y=y[ix],yo=data$value[ix]),
                      thres=argv$thres.puddle[j],
                      gt_or_lt="lt",
                      Dh=Dh_puddle,
                      eps2=argv$eps2.puddle,
                      n.eveYES_in_the_clump=argv$n_lt.puddle[j],
                      n.eveNO_in_the_clump=argv$n_ge.puddle[j])
          sus<-which(aux!=0 & is.na(dqcflag[ix]) & doit[ix]==1)
          # set dqcflag
          if (length(sus)>0) dqcflag[ix[sus]]<-argv$puddle.code
        } else {
          print("no valid observations left, no puddle check")
        }
      } # end of loop over threshold that define events
      if (argv$verbose | argv$debug) {
        t1a<-Sys.time()
        print(paste("puddle plot, iteration=",i,
                    "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
        ncur<-length(which(dqcflag==argv$puddle.code))
        print(paste("# suspect observations=",ncur-nprev))
        nprev<-ncur
      } 
    } # end of loop over test iterations
    rm(doit)
    if (argv$verbose | argv$debug) 
      print("+---------------------------------+")
    if (argv$debug) 
      save.image(file.path(argv$debug.dir,"dqcres_puddle.RData")) 
  } # end of "if (length(mask_puddle)==0)"
} # end of if (argv$puddle)
#
#-----------------------------------------------------------------------------
# check against a first-guess (deterministic)
if (argv$fg) {
  if (argv$verbose | argv$debug) {
    nprev<-0
    print(paste0("first-guess check det (",argv$fg.code,")"))
  }
  # set doit vector
  doit<-vector(length=ndata,mode="numeric")
  doit[]<-NA
  for (f in 1:nfin) doit[data$prid==argv$prid[f]]<-argv$doit.fg[f]
  # use only (probably) good observations
  ix<-which(is.na(dqcflag) & doit!=0)
  if (length(ix)>0) {
    if (is.na(argv$thrpos.fg) | is.na(argv$thrneg.fg)) {
      sus<-which(!is.na(fg) & 
                 is.na(dqcflag) & 
                 (abs(data$value-fg)>argv$thr.fg) &
                 doit==1)
    } else {
      sus<-which(!is.na(fg) & 
                 is.na(dqcflag) & 
                 ( ((data$value-fg)>argv$thrpos.fg) | 
                   ((fg-data$value)>argv$thrneg.fg) ) &
                 doit==1)
    }
    # set dqcflag
    if (length(sus)>0) dqcflag[sus]<-argv$fg.code
  }  else {
    print("no valid observations left, no first-guess check")
  }
  if (argv$verbose | argv$debug) {
    print(paste("# observations that fail the first-guess check (det)=",
                length(which(dqcflag==argv$fg.code))))
    print("+---------------------------------+")
  }
  rm(doit)
  if (argv$debug) 
    save.image(file.path(argv$debug.dir,"dqcres_fg.RData")) 
}
#
#-----------------------------------------------------------------------------
# check against a first-guess (ensemble)
if (argv$fge) {
  # set doit vector
  doit<-vector(length=ndata,mode="numeric")
  doit[]<-NA
  for (f in 1:nfin) doit[data$prid==argv$prid[f]]<-argv$doit.fge[f]
  # use only (probably) good observations
  ix<-which(is.na(dqcflag) & doit!=0)
  if (length(ix)>0) {
    if (argv$variable=="RR") {
      ttot<-boxcox(x=data$value,lambda=argv$boxcox.lambda)
    } else {
      ttot<-data$value
    }
    fge.sd<-pmax(fge.sd,argv$sdmin.fge,na.rm=T)
    fge.iqr<-pmax((fge.q75-fge.q25),argv$iqrmin.fge,na.rm=T)
    sus<-which( ( 
    (abs(fge.mu-ttot)>(argv$csd.fge*argv$infsd.fge*fge.sd))    | 
    ((ttot-fge.q75)>(argv$ciqr.fge*argv$infiqr.fge*fge.iqr))   |
    ((fge.q25-ttot)>(argv$ciqr.fge*argv$infiqr.fge*fge.iqr)) ) &
    !is.na(fge.mu) & 
    !is.na(fge.sd) & 
    !is.na(fge.iqr) & 
    is.na(dqcflag) &
    doit==1)
    # set dqcflag
    if (length(sus)>0) dqcflag[sus]<-argv$fge.code
  }  else {
    print("no valid observations left, no ensemble first-guess check")
  }
  if (argv$verbose | argv$debug) {
    print(paste("# observations that fail the first-guess check (ens)=",
                length(which(dqcflag==argv$fge.code))))
    print("+---------------------------------+")
  }
  rm(doit)
  if (argv$debug) 
    save.image(file.path(argv$debug.dir,"dqcres_fge.RData")) 
  if (exists("ttot")) rm(ttot)
  if (exists("ix")) rm(ix)
}
#
#-----------------------------------------------------------------------------
# SCT - Spatial Consistency Test
# NOTE: keep-listed stations are used but they canNOT be flagged here
if (argv$verbose | argv$debug) nprev<-0
# set doit vector
doit<-vector(length=ndata,mode="numeric")
doit[]<-NA
for (f in 1:nfin) doit[data$prid==argv$prid[f]]<-argv$doit.sct[f]
# set min and max for the background values
sctvmin<-ifelse(argv$variable=="RR",-1./argv$boxcox.lambda,
                                    argv$vmin)
sctvmax<-ifelse(argv$variable=="RR",boxcox(argv$vmax,argv$boxcox.lambda),
                                    argv$vmax)
# test
for (i in 1:argv$i.sct) {
  # use only (probably) good observations with doit!=0
  ix<-which( (is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0 )
  if (length(ix)>0) {
    t0a<-Sys.time()
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
    # define global 1D vector used in statSpat (1D for fast access)
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
  } else {
    print("no valid observations left, no SCT")
  }
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    print(paste("SCT, iteration=",i,
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    ncur<-length(which(dqcflag==argv$sct.code))
    print(paste("# suspect observations=",ncur-nprev))
    nprev<-length(which(dqcflag==argv$sct.code))
  }
}
rm(doit)
if (argv$verbose | argv$debug) 
  print("+---------------------------------+")
#
# coefficient of observation representativeness
#-----------------------------------------------------------------------------
# corep has been set by function sct to the observation error variance
qmn<-0.25
qmx<-0.75
qav<-0.5
ix<-which(!is.na(corep) & (is.na(dqcflag) | dqcflag==argv$keep.code)) 
if (length(ix)>0) {
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
if (argv$debug) 
  save.image(file.path(argv$debug.dir,"dqcres_sct.RData")) 
#
#-----------------------------------------------------------------------------
# check elevation against dem 
# NOTE: keep-listed stations canNOT be flagged here
if (argv$dem) {
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
    print(paste("# observations too far from dem=",
                length(which(dqcflag==argv$dem.code))))
    print("+---------------------------------+")
  }
  rm(doit)
}
if (argv$debug) 
  save.image(file.path(argv$debug.dir,"dqcres_demcheck.RData")) 
#
#-----------------------------------------------------------------------------
# STEVE - isolated event test (YES/NO)
if (argv$steve) {
  t0a<-Sys.time()
  if (argv$verbose | argv$debug) nprev<-0
  # set doit vector
  doit<-vector(length=ndata,mode="numeric")
  doit[]<-NA
  for (f in 1:nfin) doit[data$prid==argv$prid[f]]<-argv$doit.steve[f]
  # test
  for (i in 1:argv$i.steve) {
    t0a<-Sys.time()
    for (j in 1:length(argv$thres.steve)) {
      # use only (probably) good observations
      ix<-which((is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0)
      if (length(ix)>0) {
        obs<-data.frame(x[ix],y[ix],data$value[ix])
        aux<-steve(obs=data.frame(x=x[ix],y=y[ix],yo=data$value[ix]),
                   thres=argv$thres.steve[j],
                   gt_or_lt="lt",
                   pmax=argv$pmax_lt.steve[j],
                   dmax=argv$dmax_lt.steve[j],
                   n.sector=16,
                   n.connected_eveYES=argv$n_lt.steve[j],
                   frac.eveYES_in_the_clump=argv$frac_lt.steve[j],
                   dmin.next_eveYES=argv$dmin_next_lt.steve[j])
        sus<-which(aux!=0 & is.na(dqcflag[ix]) & doit[ix]==1)
        # set dqcflag
        if (length(sus)>0) dqcflag[ix[sus]]<-argv$steve.code
      } else {
        print("no valid observations left, no STEVE")
      }
      ix<-which((is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0)
      if (length(ix)>0) {
        obs<-data.frame(x[ix],y[ix],data$value[ix])
        aux<-steve(obs=data.frame(x=x[ix],y=y[ix],yo=data$value[ix]),
                   thres=argv$thres.steve[j],
                   gt_or_lt="ge",
                   pmax=argv$pmax_ge.steve[j],
                   dmax=argv$dmax_ge.steve[j],
                   n.sector=16,
                   n.connected_eveYES=argv$n_ge.steve[j],
                   frac.eveYES_in_the_clump=argv$frac_ge.steve[j],
                   dmin.next_eveYES=argv$dmin_next_ge.steve[j])
        sus<-which(aux!=0 & is.na(dqcflag[ix]) & doit[ix]==1)
        # set dqcflag
        if (length(sus)>0) dqcflag[ix[sus]]<-argv$steve.code
      } else {
        print("no valid observations left, no STEVE")
      }
    } # end of loop over threshold that define events
    if (argv$verbose | argv$debug) {
      t1a<-Sys.time()
      print(paste("STEVE, iteration=",i,
                  "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
      ncur<-length(which(dqcflag==argv$steve.code))
      print(paste("# suspect observations=",ncur-nprev))
      nprev<-length(which(dqcflag==argv$steve.code))
    }
  }
  rm(doit)
  if (argv$verbose | argv$debug) 
    print("+---------------------------------+")
  if (argv$debug) 
    save.image(file.path(argv$debug.dir,"dqcres_steve2.RData")) 
}
#
#-----------------------------------------------------------------------------
# puddle check (second round)
if (argv$puddle) {
  if (argv$verbose | argv$debug) {
    nprev<-0
    print(paste0("puddle-check 2 (",argv$puddle.code,")"))
  }
  # set doit vector
  doit<-vector(length=ndata,mode="numeric")
  doit[]<-NA
  for (f in 1:nfin) doit[data$prid==argv$prid[f]]<-argv$doit.puddle[f]
  ix<-which((is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0)
  ptmp<-length(ix)
  if (ptmp<1) {
    print("puddle check: no valid observations left, no test")
    mask_puddle<-integer(0)
  } else {
#    if (!exists("mask_puddle")) {
      # prepare mask with gridpoints where the calculation can be done
      t00<-Sys.time()
      xytmp<-cbind(x[ix],y[ix])
      Dh_puddle<-dobs_fun(obs=data.frame(x=x[ix],y=y[ix]),
                          k=argv$dobs_k.puddle)
      Dh_puddle<-max(argv$Dh_min.puddle,(Dh_puddle/1000))
      rm(xytmp)
      # define raster grid
      grid_res<-as.integer((Dh_puddle/argv$dobs_k.puddle)*1000)
      if (argv$verbose | argv$debug) {
        print(paste("Dh (to the closest",argv$dobs_k.puddle,"neighbour)=",
              round(Dh_puddle,1),"km (grid res=",round(grid_res,0),"m)"))
      }
      rgrid_puddle<-raster(ext=e,resolution=grid_res)
      rgrid_puddle[]<-NA
      xygrid_puddle<-xyFromCell(rgrid_puddle,1:ncell(rgrid_puddle))
      xgrid_puddle<-xygrid_puddle[,1]
      ygrid_puddle<-xygrid_puddle[,2]
      rm(xygrid_puddle)
      # NOTE: z is not used here, but we are ready to introduce it for future versions
      zgrid_puddle<-rep(0,length(ygrid_puddle))
      D<-exp(-0.5*(((outer(x[ix],x[ix],FUN="-")**2.+
                     outer(y[ix],y[ix],FUN="-")**2.)**0.5/1000.)/Dh_puddle)**2.)
      diag(D)<-diag(D)+argv$eps2.puddle
      InvD<-chol2inv(chol(D))
      rm(D)
      xidi<-OI_RR_fast(yo.sel=rep(1,ptmp),
                       yb.sel=rep(0,ptmp),
                       xb.sel=rep(0,length(xgrid_puddle)),
                       xgrid.sel=xgrid_puddle,
                       ygrid.sel=ygrid_puddle,
                       zgrid.sel=zgrid_puddle,
                       VecX.sel=x[ix],
                       VecY.sel=y[ix],
                       VecZ.sel=rep(0,ptmp),
                       Dh.cur=Dh_puddle,
                       Dz.cur=1000000)
      rm(InvD)
      mask_puddle<-which(xidi>=0.00001)
      t11<-Sys.time()
      if (argv$verbose | argv$debug) {
        print(paste("prepare mask, time=",round(t11-t00,1),attr(t11-t00,"unit")))
        print(paste("tot points=",length(xgrid_puddle)))
        print(paste("masked points=",length(mask_puddle)))
      }
#    } # end "if (!exists("mask_puddle"))"
  }
  if (length(mask_puddle)==0) {
    print("puddle check: no test")
  } else {
    xgrid_puddle<-xgrid_puddle[mask_puddle]
    ygrid_puddle<-ygrid_puddle[mask_puddle]
    zgrid_puddle<-zgrid_puddle[mask_puddle]
    # test
    for (i in 1:argv$i.puddle) {
      t0a<-Sys.time()
      for (j in 1:length(argv$thres.puddle)) {
        # use only (probably) good observations
        ix<-which((is.na(dqcflag) | dqcflag==argv$keep.code) & doit!=0)
        if (length(ix)>0) {
          aux<-puddle(obs=data.frame(x=x[ix],y=y[ix],yo=data$value[ix]),
                      thres=argv$thres.puddle[j],
                      gt_or_lt="lt",
                      Dh=Dh_puddle,
                      eps2=argv$eps2.puddle,
                      n.eveYES_in_the_clump=argv$n_lt.puddle[j],
                      n.eveNO_in_the_clump=argv$n_ge.puddle[j])
          sus<-which(aux!=0 & is.na(dqcflag[ix]) & doit[ix]==1)
          # set dqcflag
          if (length(sus)>0) dqcflag[ix[sus]]<-argv$puddle.code
        } else {
          print("no valid observations left, no puddle check")
        }
      } # end of loop over threshold that define events
      if (argv$verbose | argv$debug) {
        t1a<-Sys.time()
        print(paste("puddle plot, iteration=",i,
                    "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
        ncur<-length(which(dqcflag==argv$puddle.code))
        print(paste("# suspect observations=",ncur-nprev))
        nprev<-ncur
      } 
    } # end of loop over test iterations
    rm(doit)
    if (argv$verbose | argv$debug) 
      print("+---------------------------------+")
    if (argv$debug) 
      save.image(file.path(argv$debug.dir,"dqcres_puddle2.RData")) 
  } # end of "if (length(mask_puddle)==0)"
} # end of if (argv$puddle)
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
  # define global 1D vector used in statSpat (1D for fast access)
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
  print("summary")
  print(paste("# total suspect & no-metadata observations=",
              length(which(dqcflag!=0))," [",
              round(100*length(which(dqcflag!=0))/ndata,0),
              "%]",sep="") )
   print(paste("# total suspect & no-metadata observations (statistics over not NAs only)=",
              length(which(dqcflag!=0 & !is.na(data$value)))," [",
        round(100*length(which(dqcflag!=0 & !is.na(data$value)))/length(which(!is.na(data$value))),0),
              "%]",sep="") )
  print(paste("# total good observations=",
              length(which(dqcflag==0))," [",
              round(100*length(which(dqcflag==0))/ndata,0),
              "%]",sep="") )
  print(paste("# total good observations (statistics over not NAs only)=",
              length(which(dqcflag==0))," [",
   round(100*length(which(dqcflag==0 & !is.na(data$value)))/length(which(!is.na(data$value))),0),
              "%]",sep="") )
  print("+---------------------------------+")
}
#
#-----------------------------------------------------------------------------
# Include radar-derived precipitation in the output file
#  abbreviation fge = first-guess ensemble
if (argv$radarout) {
  if (argv$verbose | argv$debug) 
    print("include radar-derived precipitation in the output file")
  drad<-getValues(rrad) #fg-vec & fg-grid
  # get radar-point coordinates into fge CRS (only for not-NA points)
  ix1<-which(!is.na(drad)) # indx over fg-vec
  radxy<-as.data.frame(xyFromCell(rrad,ix1)) #ix1-vec
  names(radxy)<-c("x","y")
  coordinates(radxy)<-c("x","y")
  proj4string(radxy)<-CRS(argv$proj4fg)
  radxy.fge<-spTransform(radxy,CRS=argv$proj4fge) #ix1-vec
  radxy.fge<-as.data.frame(radxy.fge)
  # aggregate radar values onto the fge grid
  raux<-rasterize(radxy.fge, #fge-grid
                  rfge.grid,
                  drad[ix1],
                  fun=mean)
  raux[raux<=0]<-NA 
  daux<-getValues(raux) #fge-vec
  if (argv$debug) 
    plot_debug(ff=file.path(argv$debug.dir,"radar_1_fgegrid.png"),
               r=raux,x=x,y=y,proj4=argv$proj4to,proj4plot=argv$proj4fge)
  # keep points where radar prec is greater than the ensemble median
  ix2<-which(!is.na(daux) & daux<fge.gq50) # indx over fge-vec
  if (length(ix2)>0) daux[ix2]<-NA
  raux[]<-daux
  if (argv$debug) { 
    plot_debug(ff=file.path(argv$debug.dir,"radar_2_fgegrid.png"),
               r=raux,x=x,y=y,proj4=argv$proj4to,proj4plot=argv$proj4fge)
    ixdeb<-which(!is.na(fge.gq50) & !is.na(daux))
    png(file=file.path(argv$debug.dir,"radar_vs_fgegq50.png"))
    plot(daux[ixdeb],fge.gq50[ixdeb],xlab="radar (mm)",
         ylab="numerical model, ensemble median (mm)")
    dev.off()
  }
  ix1.1<-which(!is.na(extract(raux,radxy.fge))) #indx over ix1-vec
  ix1.2<-which(is.na(extract(raux,radxy.fge))) #indx over ix1-vec
  rm(raux,daux)
  # case of valid radar-data points found
  if (length(ix1.1)>0) {
    # remove patches of connected cells that are too small
    drad[ix1[ix1.2]]<-NA
    rrad[]<-drad
    rclump<-clump(rrad)
    fr<-freq(rclump)
    ix<-which(!is.na(fr[,2]) & fr[,2]<=10)
    if (length(ix)>0) drad[getValues(rclump) %in% fr[ix,1]]<-NA
    rrad[]<-drad
    rm(rclump,fr,ix)
    # (optional) aggregate radar data onto a coarser grid
    if (argv$radarout.aggfact>1) {
      raux<-aggregate(rrad,fact=argv$radarout.aggfact,na.rm=T,expand=T,fun=mean)
      rrad<-raux
      rm(raux)
      drad<-getValues(rrad)
    }
    if (argv$debug) 
      plot_debug(ff=file.path(argv$debug.dir,"radar_3_fggrid.png"),
                 r=rrad,x=x,y=y,proj4=argv$proj4to,proj4plot=argv$proj4fg)
    # prepare list of valid radar-points in output CRS
    ix1<-which(!is.na(drad)) # indx over fg-vec
    radxy<-as.data.frame(xyFromCell(rrad,ix1)) #ix1-vec
    names(radxy)<-c("x","y")
    coordinates(radxy)<-c("x","y")
    proj4string(radxy)<-CRS(argv$proj4fg)
    radxy.from<-as.data.frame(spTransform(radxy,CRS=argv$proj4from))  
    radx.from<-radxy.from[,1]
    rady.from<-radxy.from[,2]
    radrr<-drad[ix1]
    if (argv$debug) {
      ixdeb<-which(dqcflag==0 & !is.na(dqcflag))
      png(file=file.path(argv$debug.dir,"radar_4_ll.png"),
                         width=800,height=800)
      plotp(x=c(data$lon[ixdeb],radx.from),
            y=c(data$lat[ixdeb],rady.from),
            val=c(data$value[ixdeb],radrr),
            br=c(0,0.1,0.5,1,2,3,4,5,7,10,15,20,50,100),
            col=c("gray",rev(rainbow(12))),
            map=NULL,map.br=NULL,map.col=NULL,
            xl=range(radx.from),yl=range(rady.from))
      dev.off() 
    }
  # case of no valid radar-data points found
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
varidx.out<-varidx
if (any(!is.na(argv$varname.opt))) 
  varidx.out<-c(varidx,varidx.opt[which(!is.na(varidx.opt))]) 
dataout<-array(data=NA,
               dim=c(length(data$lat),(length(varidx.out)+4)))
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
    str[s]<-argv$varname.lat.out
    dataout[,s]<-round(data$lat,argv$latlon.dig.out)
  } else if (pos.s==2) {
    str[s]<-argv$varname.lon.out
    dataout[,s]<-round(data$lon,argv$latlon.dig.out)
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
    datarad[,which(str==argv$varname.lat.out)]<-round(rady.from,argv$latlon.dig.out)
    datarad[,which(str==argv$varname.lon.out)]<-round(radx.from,argv$latlon.dig.out)
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
