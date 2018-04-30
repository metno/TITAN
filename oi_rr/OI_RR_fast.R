#+
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

