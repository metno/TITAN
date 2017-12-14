# TITAN - Temperature spatIal daTa quAlity coNtrol
# command line
# >R --vanilla file_input file_output < titan.R
# 
# arguments:
# file_input, text file with colums separated by ";"
# column names = lat,lon,elev,value (any order is accepted)
#
# file_output, text file with colums separated by ";"
#  same number of rows as file_in
#  header:
#   lat;lon;x;y;z;val;dqc;
#  data quality control (dqc) codes:
#  0 = ok
#  1 = missing metadata
#  2 = plausibility test failed
#  3 = buddy check failed
#  4 = SCT failed
#  5 = isolated station 
#-----------------------------------------------------------------------------
library(sp)
library(raster)
library(rgdal)
library(dotnc)
#options(warn = 2, scipen = 999)
options(scipen = 999)
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# FUNCTIONS

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
statSpat<-function(ixyzt,drmin) {
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
  tcor<-ttot[i]+0.0065*dz
  tmean<-mean(tcor) 
  tsd<-sd(tcor) 
  dz_mx<-max(abs(dz))
  return(c(length(i),round(dz_mx,0),round(tmean,1),round(tsd,3)))
}

#+ vertical profile of temperature (Frei, 2014)
tvertprof<-function(z,t0,gamma,a,h0,h1i) {
# ref:
# Frei, C. (2014). Interpolation of temperature in a mountainous region 
#  using nonlinear profiles and non‐Euclidean distances.
#  International Journal of Climatology, 34(5), 1585-1605.
# input
#  z= array. elevations [m]
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
sct<-function(ixyn,
              nmin=50,dzmin=30,
              Dhmin=10,Dz=200,eps2=0.5,
              T2=16,sus.code=4) {
# ref:
#  Lussana, C., Uboldi, F., & Salvati, M. R. (2010). A spatial consistency 
#   test for surface observations from mesoscale meteorological networks.
#   Quarterly Journal of the Royal Meteorological Society, 136(649), 1075-1088.
# input
#  ixyn= vector(4). 1=box identifier on the grid; 
#                   2/3=easting/northing coord (center of the box);
#                   4=number of stations within the box
#  NOTE: stations are associated to a box before running this function
#  nmin= numeric. minimum number of stations to fit a vertical profile
#  dzmin= numeric. minimum elevation range to fit a vertical profile [m]
#  Dhmin= numeric. minimum value for OI horizontal decorellation length [km]
#  Dz= numeric. OI vertical decorellation length [m]
#  eps2= numeric. OI ratio between obs_err_variance/backg_err_variance
#  T2=numeric. SCT threshold. (obs-pred)^2/(varObs+varPred)^2 > T2, suspect!
#  sus.code=numeric. identifier code for suspect observation
# output
#  number of rejected stations. NA if the function has not been applied
#  NOTE: "dqcflag" global variable is also updated
#------------------------------------------------------------------------------
  # few stations in the box = no SCT
  if (ixyn[4]<nmin | is.na(ixyn[4])) return(NA)
  # j, index for the stations in the box
  j<-which(itot==ixyn[1])
  # "zopt" and "topt" are set as global variables, so that they can be used by 
  #   other functions in the optimizaiton procedure
  assign("zopt",ztot[j],envir=.GlobalEnv)
  dz<-as.numeric(quantile(zopt,probs=0.95))-
      as.numeric(quantile(zopt,probs=0.05))
  assign("topt",ttot[j],envir=.GlobalEnv)
  # if region is too flat, then background is the average
  if (dz<dzmin) {
    tb<-topt
    tb[]<-mean(topt)
  # otherwise, fit a vertical profile of temperature
  } else {
    par<-c(mean(topt),-0.0065,5,
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
  # OI for SCT (Lussana et al., 2010)
  # initialize variables
  #  dqc flags
  dqctmp<-dqcflag[ix[j]]
  #  probability of gross error (pog)
  pog<-dqctmp
  pog[]<-NA
  # distance matrices
  disth<-(outer(xtot[j],xtot[j],FUN="-")**2.+
          outer(ytot[j],ytot[j],FUN="-")**2.)**0.5/1000.
  distz<-abs(outer(ztot[j],ztot[j],FUN="-"))
  # set to optimal Dh
  Dh<-max(Dhmin,
          as.numeric(quantile(disth,probs=0.1)))
  # background error correlation matrix
  S<-exp(-0.5*(disth/Dh)**2.-0.5*(distz/Dz)**2.)
  # S+eps2I
  diag(S)<-diag(S)+eps2
  # innvoation
  d<-topt-tb
  # select observations to test 
  sel<-which(is.na(dqctmp))
  first<-T
  # loop over SCT iterations 
  # NOTE: SCT flags at most one observation, iterate until no observations fail
  while (length(sel)>nmin) { 
    # update selection
    sel<-which(is.na(dqctmp))
    # first iteration, inver the matrix
    if (first) {
      SRinv<-chol2inv(chol(S))
      first<-F
    } else {
      # Update inverse matrix (Uboldi et al 2008, Appendix and erratum)
      aux<-SRinv
      SRinv<-aux[-indx,-indx]-
             (tcrossprod(aux[indx,-indx],aux[-indx,indx]))*Zinv[indx]
      rm(aux)
    }
    # next tree lines: compute cvres=(obs - CrossValidation_prediction)
    Zinv<-1/diag(SRinv)
    SRinv.d<-crossprod(SRinv,d[sel])
    cvres<-Zinv*SRinv.d
    # pog=cvres/(sig2obs+sig2CVpred)
    pog[sel]<-cvres**2/Zinv
    # check if any obs fails the test
    if (any(pog[sel]>T2)) {
      # flag as suspect only the observation with the largest cvres 
      indx<-which.max(pog[sel])
      dqctmp[sel[indx]]<-sus.code
      # update global variable with flags
      dqcflag[ix[j[sel[indx]]]]<-sus.code
      assign("dqcflag",dqcflag,envir=.GlobalEnv)
#      print("yo yb cvres cvres/sigma sig2o+sig2cv srinv.d")
#      print(cbind(topt[sel][indx],
#                  tb[sel][indx],
#                  cvres[indx],
#                  pog[sel[indx]],
#                  round(Zinv[indx],4), 
#                  round(SRinv[indx],4) )) 
    } else {
      break
    }
  } # end cycle SCT model
  # debug: begin
  if (debug) {
    susi<-which(!is.na(dqctmp))
    png(file=paste("vert_",formatC(ixyn[1],width=5,flag="0"),".png",sep=""),
        width=800,height=800)
    plot(topt,zopt,xlim=c(max(-30,min(topt)),min(30,max(topt))),
                   ylim=c(0,min(2000,max(zopt))) )
    zz<-seq(0,2000,by=0.1)
    if (!(dz<dzmin)) {
      points(tvertprof(zz,t0=opt$par[1],gamma=opt$par[2],
                       a=opt$par[3],h0=opt$par[4],h1i=opt$par[5]),
                       zz,col="blue",pch=19,cex=0.5)
    } else {
      points(rep(tb[1],length(zz)),zz,col="gold",pch=19,cex=0.5)
    }
    points(topt[susi],zopt[susi],pch=19,col="red")
    abline(h=seq(-1000,10000,by=100),col="gray",lty=2)
    abline(h=0,lwd=2,col="black")
    dev.off()
    png(file=paste("horz_",formatC(ixyn[1],width=5,flag="0"),".png",sep=""),
        width=800,height=800)
    plot(xtot[j],ytot[j])
    dem1<-crop(dem$raster,
               extent(c(ixyn[2]-100000,
                        ixyn[2]+100000,
                        ixyn[3]-100000,
                        ixyn[3]+100000
                        )))
    image(dem1,add=T,
          breaks=c(0,10,25,50,100,250,500,750,1000,1250,1500,1750,2000,2500,3000),
          col=gray.colors(14))
    points(xtot[j],ytot[j],pch=19,col="blue")
    points(xtot[j[susi]],ytot[j[susi]],pch=19,col="red")
    dev.off()
  }
  # debug: end
  return(length(which(dqctmp==sus.code)))
}
#==============================================================================
#  MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN
#==============================================================================
t0<-Sys.time()
# constants
# define domains corners
#  NOTE: setup to have Oslo in a single box
lonl<-c(5,28)
latl<-c(53.25,71.8)
# proj4 strings
proj4.wgs84<-"+proj=longlat +datum=WGS84"
proj4.lcc<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
#
# dqc constants
min.elev<-0 # m amsl
max.elev<-2500 # m amsl
# plausibility test
# TODO: monthly thresholds
tmin<--50
tmax<-40
# buddy check
drmin.buddy<-3000 # m
pog.max.buddy<-5 # probability of gross error, theshold
nstat.min.buddy<-5 # minimum number of stations
dz.max.buddy<-30 # m, elevation range (no check if elevation > dz.max.buddy) 
# SCT (more parameter definitions, see function "sct")
ncol.sct<-20 # define grid for SCT
nrow.sct<-20 # define grid for SCT
# isolated stations
drmin.iso<-25000 # m
ns.min.iso<-10 # number of stations, less than that and you are alone 
#
debug<-F
#
#-----------------------------------------------------------------------------
# Read command line arguments
arguments <- commandArgs()
print("Arguments")
print(arguments)
fin<-arguments[3]
fout<-arguments[4]
#
#-----------------------------------------------------------------------------
# read the dem (used for debugging/plots)
if (debug) {
  fdem<-"/lustre/storeB/project/metkl/klinogrid/geoinfo/meps_gmted2010_1km_topo_topdown.nc"
  dem<-nc4in_easy(file=fdem,
                  var.name="altitude",
                  proj.name="projection_lambert",
                  proj.att="proj4",topdown=F)
}
#
#-----------------------------------------------------------------------------
# read data
data<-read.table(file=fin,header=T,sep=";",stringsAsFactors=F,strip.white=T)
data$lat<-as.numeric(data$lat)
data$lon<-as.numeric(data$lon)
data$elev<-as.numeric(data$elev)
z<-data$elev
data$value<-as.numeric(data$value)
ndata<-length(data$lat)
#
#-----------------------------------------------------------------------------
# define dqcflag vector
dqcflag<-vector(mode="numeric",length=ndata)
dqcflag[]<-NA
#
#-----------------------------------------------------------------------------
# test for no metadata
meta<-!is.na(data$lat) & 
      !is.na(data$lon) &
      !is.na(data$elev) &
      data$elev>=min.elev & data$elev<=max.elev &
      !is.na(data$value)
dqcflag[which(!meta)]<-1
#
#-----------------------------------------------------------------------------
# plausibility test
ix<-which( is.na(dqcflag)  &
           data$value<tmin &
           data$value>tmax)
dqcflag[ix]<-2
#
#-----------------------------------------------------------------------------
# lat-lon to km tranformation
coord<-SpatialPoints(cbind(data$lon,data$lat),
                     proj4string=CRS(proj4.wgs84))
coord.new<-spTransform(coord,CRS(proj4.lcc))
xy.RR<-coordinates(coord.new)
x<-round(xy.RR[,1],0)
y<-round(xy.RR[,2],0)
xp<-expand.grid(lonl,latl)
coord<-SpatialPoints(xp,
                     proj4string=CRS(proj4.wgs84))
coord.new<-spTransform(coord,CRS(proj4.lcc))
# define the extent for the SCT grid
e<-extent(coord.new)
xl<-e[1:2]
yl<-e[3:4]
#
#-----------------------------------------------------------------------------
# buddy check 
#  compare each observation against the average of observations 
#  within a box of /pm 1.5Km
ix<-which(is.na(dqcflag))
# define global 1D vector used in statSpat (1D for fast access)
xtot<-x[ix]
ytot<-y[ix]
ztot<-as.numeric(data$elev[ix])
ttot<-as.numeric(data$value[ix])
# apply will loop over this 4D array
ixyzt_tot<-cbind(1:length(xtot),
                 xtot,ytot,
                 ztot,ttot)
stSp_3km<-apply(ixyzt_tot,FUN=statSpat,MARGIN=1,drmin=drmin.buddy)
# probability of gross error
pog<-(ttot-stSp_3km[3,])**2/stSp_3km[4,]**2
# suspec: more than 5 obs in the box; all obs within 30m in the vertical; pog>5 
sus<-which(pog>pog.max.buddy & 
           stSp_3km[1,]>nstat.min.buddy & 
           stSp_3km[2,]<dz.max.buddy)
# set dqcflag
dqcflag[ix[sus]]<-3
#
#-----------------------------------------------------------------------------
# SCT - Spatial Consistency Test
# create the grid for SCT. SCT is done independently in each box
# NOTE: box size around 100Km should be ok
r<-raster(e,ncol=ncol.sct,nrow=nrow.sct)
print(r)
xy<-xyFromCell(r,1:ncell(r))
xr<-xy[,1]
yr<-xy[,2]
ir<-1:ncell(r)
r[]<-1:ncell(r)
# use only probably good observations
ix<-which(is.na(dqcflag))
# define global 1D vector used in statSpat (1D for fast access)
xtot<-x[ix]
ytot<-y[ix]
ztot<-z[ix]
ttot<-as.numeric(data$value[ix])
dqcflagtot<-dqcflag[ix]
# assign each station to the corresponding box
itot<-extract(r,cbind(xtot,ytot))
# count the number of observations in each box
rnobs<-rasterize(cbind(xtot,ytot),r,ttot,fun=function(x,...)length(x))
nr<-getValues(rnobs)
# create the 4D array for the function call via apply
ixyn<-cbind(ir,xr,yr,nr)
# SCT within each (grid) box
# NOTE: dqcflag is update in "sct" function
out<-apply(ixyn,FUN=sct,MARGIN=1)
#
#-----------------------------------------------------------------------------
# check for isolated stations
#  less than 10 stations in 25km
# use only probably good observations
ix<-which(is.na(dqcflag))
# define global 1D vector used in statSpat (1D for fast access)
xtot<-x[ix]
ytot<-y[ix]
xy<-cbind(xtot,ytot)
ns<-apply(xy,FUN=nstat,MARGIN=1,drmin=drmin.iso)
sus<-which(ns<ns.min.iso)
# set dqcflag
dqcflag[ix[sus]]<-5
#
#-----------------------------------------------------------------------------
# write the output file
cat(file=fout,"lat;lon;x;y;z;val;dqc;\n",append=F)
#
dqcflag[is.na(dqcflag)]<-0
cat(file=fout,
    paste(data$lat,data$lon,
          x,y,z,
          data$value,
          dqcflag,"\n",sep=";"),
    append=T)
#
#-----------------------------------------------------------------------------
# Normal exit
print("normal exit")
t1<-Sys.time()
print(t1-t0)
q(status=0)
