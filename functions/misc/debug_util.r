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
  if (!is.null(proj4plot)) xy<-spTransform(xy,crs(proj4plot))
  points(xy,cex=0.8,pch=19)
  dev.off()
}

