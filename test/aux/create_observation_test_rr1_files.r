library(ncdf4)
library(rgdal)
library(dotnc)

#
# constants
proj4.wgs84 <- "+proj=longlat +datum=WGS84"
proj4.lcc<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06" 

#
# read the variables from the nc-file
# (note: you may change this part accordingto your nc-tastes)
 
t <- nc4.getTime( "../data/background_test_rr1.nc")
r <- read_dotnc( nc.file="../data/background_test_rr1.nc",
                 nc.varname="precipitation_amount",
                 topdown=T,
                 out.dim=list(ndim=4,tpos=4,epos=3,
                              names=c("x","y","ensemble_member","time")),
                 proj4=proj4.lcc,
                 nc.proj4=NULL,
                 selection=list(t=t,e=0),
                 gridded=T,
                 verbose=F)
r$stack[r$stack<0]<-0

dem <- read_dotnc( nc.file="../data/background_test_rr1.nc",
                 nc.varname="altitude",
                 topdown=T,
                 out.dim=list(ndim=2,
                              names=c("x","y")),
                 proj4=proj4.lcc,
                 nc.proj4=NULL,
                 selection=list(t=t),
                 gridded=T,
                 verbose=F)

#
# generate observations
p     <- c(2000,100)
pGE   <- c(20,3)
nprid <- length(p)
for (i in 1:nprid) {
  x <- runif( p[i], min=extent(r$stack)[1], max=extent(r$stack)[2])
  y <- runif( p[i], min=extent(r$stack)[3], max=extent(r$stack)[4])
  z <- extract( dem$stack, cbind(x,y), method="bilinear")
  val <- extract( r$stack, cbind(x,y), method="bilinear")
  ge <- val; ge[] <- 0
  ix <- sample( 1:p[i], ceiling( p[i]*pGE[i]/100))
  ge[ix] <- 1
  val[ix] <- runif( length(ix), min=0, max=300)  
  coord <- spTransform( SpatialPoints( cbind( x, y), proj4string=CRS(proj4.lcc )), CRS(proj4.wgs84))
  lon <- attr(coord,"coords")[,1]
  lat <- attr(coord,"coords")[,2]
  ffout <- paste0( "../data/observation_test_rr1_prid",
                   formatC(i,width=2,flag="0"),
                   "_p",formatC(p[i],width=5,flag="0"),
                   "_pGE",formatC(pGE[i],width=3,flag="0"),
                   "percent.txt")
  cat( file=ffout, append=F, "sourceId;lon;lat;x_lcc;y_lcc;val;ge\n")
  cat( file=ffout, append=T, 
       paste0( 1:p[i], ";",
               round(lon,6), ";", 
               round(lat,6), ";",
               round(x), ";",
               round(y), ";",
               round(val,1), ";",
               round(ge), "\n"))
}

p     <- c(2000,100)
pGE   <- c(3,1)
nprid <- length(p)
for (i in 1:nprid) {
  x <- runif( p[i], min=extent(r$stack)[1], max=extent(r$stack)[2])
  y <- runif( p[i], min=extent(r$stack)[3], max=extent(r$stack)[4])
  z <- extract( dem$stack, cbind(x,y), method="bilinear")
  val <- extract( r$stack, cbind(x,y), method="bilinear")
  ge <- val; ge[] <- 0
  ixx <- which( val > 0.1)
  ix <- ixx[sample( 1:length(ixx), ceiling( p[i]*pGE[i]/100))]
  ge[ix] <- 1
  val[ix] <- rep( 0, length(ix))
  coord <- spTransform( SpatialPoints( cbind( x, y), proj4string=CRS(proj4.lcc )), CRS(proj4.wgs84))
  lon <- attr(coord,"coords")[,1]
  lat <- attr(coord,"coords")[,2]
  ffout <- paste0( "../data/observation_test_rr1_prid",
                   formatC(i,width=2,flag="0"),
                   "_p",formatC(p[i],width=5,flag="0"),
                   "_zero",formatC(pGE[i],width=3,flag="0"),
                   "percent.txt")
  cat( file=ffout, append=F, "sourceId;lon;lat;x_lcc;y_lcc;val;ge\n")
  cat( file=ffout, append=T, 
       paste0( 1:p[i], ";",
               round(lon,6), ";", 
               round(lat,6), ";",
               round(x), ";",
               round(y), ";",
               round(val,1), ";",
               round(ge), "\n"))
}


#
#
q()
