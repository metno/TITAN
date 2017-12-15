#!/usr/bin/env Rscript
# + TITAN - Temperature spatIal daTa quAlity coNtrol
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
              Dhmin=10000,Dz=200,eps2=0.5,
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
#  Dhmin= numeric. minimum value for OI horizontal decorellation length [m]
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
          outer(ytot[j],ytot[j],FUN="-")**2.)**0.5
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
    sctpog[ix[j[sel]]]<-pog[sel]
    assign("sctpog",sctpog,envir=.GlobalEnv)
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
  if (argv$debug) {
    susi<-which(!is.na(dqctmp))
    if (!dir.exists(argv$debug.dir)) 
      dir.create(argv$debug.dir,showWarnings=F,recursive=T)
    f<-file.path(argv$debug.dir,
         paste("vert_",formatC(ixyn[1],width=5,flag="0"),".png",sep=""))
    png(file=f,width=800,height=800)
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
    f<-file.path(argv$debug.dir,
         paste("horz_",formatC(ixyn[1],width=5,flag="0"),".png",sep=""))
    png(file=f,width=800,height=800)
    plot(xtot[j],ytot[j])
#    dem1<-crop(dem$raster,
#               extent(c(ixyn[2]-100000,
#                        ixyn[2]+100000,
#                        ixyn[3]-100000,
#                        ixyn[3]+100000
#                        )))
#    image(dem1,add=T,
#          breaks=c(0,10,25,50,100,250,500,750,1000,1250,1500,1750,2000,2500,3000),
#          col=gray.colors(14))
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
#
# create parser object
p <- arg_parser("titan")
# specify our desired options 
# by default ArgumentParser will add an help option 
p <- add_argument(p, "input",help="input file",type="character")
p <- add_argument(p, "output",help="output file",type="character",
                  default="output.txt")
#
p <- add_argument(p, "--debug",help="debug mode",flag=T,short="-dbg")
p <- add_argument(p, "--debug.dir",help="directory for debug output",
                  type="character",default=".",short="-dbgd")
p <- add_argument(p, "--verbose",help="debug mode",flag=T,short="-v")
# NOTE: lat-lon setup to have Oslo in a single box
p <- add_argument(p, "--lonmin",help="longitude of south-eastern domain corner",
                  type="numeric",default=5,short="-x")
p <- add_argument(p, "--lonmax",help="longitude of south-western domain corner",
                  type="numeric",default=28,short="-X")
p <- add_argument(p, "--latmin",help="latitude of south-eastern domain corner",
                  type="numeric",default=53.25,short="-y")
p <- add_argument(p, "--latmax",help="latitude of north-western domain corner",
                  type="numeric",default=71.8,short="-Y")
#
p <- add_argument(p, "--spatconv",help="flag for conversion of spatial coordinates before running the data quality checks",
                  flag=T,short="-c")
p <- add_argument(p, "--proj4from",help="proj4 string for the original coordinate reference system",
                  type="character",default="+proj=longlat +datum=WGS84",short="-pf")
p <- add_argument(p, "--proj4to",help="proj4 string for the coordinate reference system where the DQC is performed",
                  type="character",default="+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06",short="-pt")
# metadata check
p <- add_argument(p, "--zmin",help="minimum allowed elevation in the domain [m amsl]",
                  type="numeric",default=0,short="-z")
p <- add_argument(p, "--zmax",help="maximum allowed elevation in the domain [m amsl]",
                  type="numeric",default=2500,short="-Z")
# Plausibility check
p <- add_argument(p, "--tmin",help="minimum allowed temperature [K or degC]",
                  type="numeric",default=-50,short="-tP")
p <- add_argument(p, "--tmax",help="maximum allowed temperature [K or degC]",
                  type="numeric",default=40,short="-TP")
# Buddy-check
p <- add_argument(p, "--dr.buddy",help="perform the buddy-check in a dr-by-dr square-box around each observation [m]",
                  type="numeric",default=3000,short="-dB")
p <- add_argument(p, "--i.buddy",help="number of buddy-check iterations",
                  type="integer",default=1,short="-iB")
p <- add_argument(p, "--thr.buddy",help="buddy-check threshold. flag observation if: (obs-pred)^2/var > thr.buddy",
                  type="numeric",default=5,short="-thB")
p <- add_argument(p, "--n.buddy",help="minimum number of neighbouring observations to perform the buddy-check",
                  type="integer",default=5,short="-nB")
p <- add_argument(p, "--dz.buddy",help="maximum allowed range of elevation in a square-box to perform the buddy-check (i.e. no check if elevation > dz.buddy)",
                  type="numeric",default=30,short="-zB")
# isolated stations
p <- add_argument(p, "--dr.isol",help="check for the number of observation in a dr-by-dr square-box around each observation [m]",
                  type="numeric",default=25000,short="-dI")
p <- add_argument(p, "--n.isol",help="threshold (number of neighbouring observations) for the identification of isolated observations.",
                  type="integer",default=10,short="-nI")
# spatial consistency test
p <- add_argument(p, "--grid.sct",help="nrow ncol (i.e. number_of_rows number_of_columns). used to define grid of boxes where the SCT is performed. SCT in each box is independent from the others",
                  type="integer",nargs=2,default=c(20,20),short="-gS")
p <- add_argument(p, "--i.sct",help="number of SCT iterations",
                  type="integer",default=1,short="-iS")
p <- add_argument(p, "--n.sct",help="minimum number of stations in a box to run SCT",
                  type="integer",default=50,short="-nS")
p <- add_argument(p, "--dz.sct",help="minimum range of elevation in a box to run SCT [m]",
                  type="numeric",default=30,short="-zS")
p <- add_argument(p, "--DhorMin.sct",help="OI, minimum allowed value for the horizontal de-correlation lenght (of the background error correlation) [m]",
                  type="numeric",default=10000,short="-hS")
p <- add_argument(p, "--Dver.sct",help="OI, vertical de-correlation lenght  (of the background error correlation) [m]",
                  type="numeric",default=200,short="-vS")
p <- add_argument(p, "--eps2.sct",help="OI, ratio between observation error variance and background error variance",
                  type="numeric",default=0.5,short="-eS")
p <- add_argument(p, "--thr.sct",help="SCT threshold. flag observation if: (obs-Cross_Validation_pred)^2/(varObs+varCVpred) > thr.sct",
                  type="numeric",default=16,short="-tS")
argv <- parse_args(p)
#
# check if input exists
if (!file.exists(argv$input)) {
  print("ERROR: input file not found")
  print(argv$input)
  quit(status=1)
}
if (argv$verbose | argv$debug) print(">> TITAN <<")
#
#-----------------------------------------------------------------------------
# constants
nometa.code<-1
p.code<-2
buddy.code<-3
sct.code<-4
isol.code<-5
#
#-----------------------------------------------------------------------------
# read data
data<-read.table(file=argv$input,header=T,sep=";",
                 stringsAsFactors=F,strip.white=T)
data$lat<-suppressWarnings(as.numeric(data$lat))
data$lon<-suppressWarnings(as.numeric(data$lon))
data$elev<-suppressWarnings(as.numeric(data$elev))
z<-data$elev
data$value<-suppressWarnings(as.numeric(data$value))
ndata<-length(data$lat)
if (ndata==0) {
  print("input file is empty")
  quit(status=0)
}
if (argv$verbose | argv$debug) {
  print(paste("number of observations=",ndata))
}
#
#-----------------------------------------------------------------------------
# define dqcflag vector
dqcflag<-vector(mode="numeric",length=ndata)
dqcflag[]<-NA
sctpog<-vector(mode="numeric",length=ndata)
sctpog[]<-NA
#
#-----------------------------------------------------------------------------
# test for no metadata 
meta<-!is.na(data$lat) & 
      !is.na(data$lon) &
      !is.na(z) & z>=argv$zmin & z<=argv$zmax &
      !is.na(data$value)
if (any(!meta)) dqcflag[which(!meta)]<-nometa.code
if (argv$verbose | argv$debug) {
  print("test for no metdata")
#  print(paste(data$lat[which(!meta)],data$lon[which(!meta)],
#              z[which(!meta)],data$value[which(!meta)]))
  print(paste("# observations lacking metadata=",length(which(!meta))))
}
#
#-----------------------------------------------------------------------------
# plausibility test
ix<-which( is.na(dqcflag)  &
           data$value<argv$tmin &
           data$value>argv$tmax)
if (length(ix)>0) dqcflag[ix]<-p.code
if (argv$verbose | argv$debug) {
  print("plausibility test")
  print(paste("# suspect observations=",length(ix)))
}
#
#-----------------------------------------------------------------------------
# coordinate transformation
if (argv$spatconv) {
  if (argv$debug) print("spatial conversion required")
  coord<-SpatialPoints(cbind(data$lon,data$lat),
                       proj4string=CRS(argv$proj4from))
  coord.new<-spTransform(coord,CRS(argv$proj4to))
  xy.new<-coordinates(coord.new)
  x<-round(xy.new[,1],0)
  y<-round(xy.new[,2],0)
  rm(xy.new)
  xp<-expand.grid(c(argv$lonmin,argv$lonmax),c(argv$latmin,argv$latmax))
  coord<-SpatialPoints(xp,
                       proj4string=CRS(argv$proj4from))
  coord.new<-spTransform(coord,CRS(argv$proj4to))
  # define the extent for the SCT grid
  e<-extent(coord.new)
  xl<-e[1:2]
  yl<-e[3:4]
} else {
  x<-data$lon
  y<-data$lat
  xl<-c(argv$lonmin,argv$lonmax)
  yl<-c(argv$latmin,argv$latmax)
  e<-extent(c(xl,yl))
}
#
#-----------------------------------------------------------------------------
# buddy check 
#  compare each observation against the average of neighbouring observations 
if (argv$verbose | argv$debug) nprev<-0
for (i in 1:argv$i.buddy) {
  # use only (probably) good observations
  ix<-which(is.na(dqcflag))
  if (length(ix)>0) {
    t0a<-Sys.time()
    # define global 1D vector used in statSpat (1D for fast access)
    xtot<-x[ix]
    ytot<-y[ix]
    ztot<-as.numeric(z[ix])
    ttot<-as.numeric(data$value[ix])
    # apply will loop over this 4D array
    ixyzt_tot<-cbind(1:length(xtot),xtot,ytot,ztot,ttot)
    stSp_3km<-apply(ixyzt_tot,FUN=statSpat,MARGIN=1,drmin=argv$dr.buddy)
    # probability of gross error
    pog<-(ttot-stSp_3km[3,])**2/stSp_3km[4,]**2
    # suspect if: 
    sus<-which(pog>argv$thr.buddy & 
               stSp_3km[1,]>argv$n.buddy & 
               stSp_3km[2,]<argv$dz.buddy)
    # set dqcflag
    if (length(sus)>0) dqcflag[ix[sus]]<-buddy.code
  } else {
    print("no valid observations left, no buddy check")
  }
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    print(paste("buddy-check, iteration=",i,
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    ncur<-length(which(dqcflag==buddy.code))
    print(paste("# suspect observations=",ncur-nprev))
    nprev<-length(which(dqcflag==buddy.code))
  }
}
#
#-----------------------------------------------------------------------------
# SCT - Spatial Consistency Test
if (argv$verbose | argv$debug) nprev<-0
for (i in 1:argv$i.sct) {
  # use only (probably) good observations
  ix<-which(is.na(dqcflag))
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
    # NOTE: dqcflag is updated in "sct" function
    out<-apply(ixyn,FUN=sct,MARGIN=1,nmin=argv$n.sct,dzmin=argv$dz.sct,
               Dhmin=argv$DhorMin.sct,Dz=argv$Dver.sct,eps2=argv$eps2.sct,
               T2=argv$thr.sct,sus.code=sct.code)
  } else {
    print("no valid observations left, no SCT")
  }
  if (argv$verbose | argv$debug) {
    t1a<-Sys.time()
    print(paste("SCT, iteration=",i,
                "/time",round(t1a-t0a,1),attr(t1a-t0a,"unit")))
    ncur<-length(which(dqcflag==sct.code))
    print(paste("# suspect observations=",ncur-nprev))
    nprev<-length(which(dqcflag==sct.code))
  }
}
#
#-----------------------------------------------------------------------------
# check for isolated stations
# use only (probably) good observations
ix<-which(is.na(dqcflag))
if (length(ix)>0) {
  # define global 1D vector used in statSpat (1D for fast access)
  xtot<-x[ix]
  ytot<-y[ix]
  xy<-cbind(xtot,ytot)
  ns<-apply(xy,FUN=nstat,MARGIN=1,drmin=argv$dr.isol)
  sus<-which(ns<argv$n.isol)
  # set dqcflag
  if (length(sus)>0)dqcflag[ix[sus]]<-isol.code
} else {
  print("no valid observations left, no check for isolated observations")
}
if (argv$verbose | argv$debug) {
  print(paste("# isolated observations=",length(which(dqcflag==isol.code))))
}
#
#-----------------------------------------------------------------------------
# observations not flagged are assumed to be good observaitons 
dqcflag[is.na(dqcflag)]<-0
if (argv$verbose | argv$debug) {
  print("summary")
  print(paste("# total suspect observations=",
              length(which(dqcflag!=0))," [",
              round(100*length(which(dqcflag!=0))/ndata,0),
              "%]",sep="") )
}
#
#-----------------------------------------------------------------------------
# write the output file
cat(file=argv$output,"prid;lat;lon;elev;value;dqc;sct;rep;\n",append=F)
#
cat(file=argv$output,
    paste(NA,
          round(data$lat,5),
          round(data$lon,5),
          round(z,1),
          round(data$value,1),
          dqcflag,
          round(sctpog,2),
          NA,
          "",sep=";",collapse='\n'),
    append=T)
#
#-----------------------------------------------------------------------------
# Normal exit
t1<-Sys.time()
if (argv$verbose | argv$debug) 
 print(paste("normal exit /time",round(t1-t0,1),attr(t1-t0,"unit")))
q(status=0)
