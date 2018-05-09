conf<-list(
           verbose=T,
           debug=T,
           debug.dir="/home/cristianl/scratch",
           #
           spatconv=T,
           variable="RH",
           # use of geographical info
           dem.fill=T,
           dem.file="/home/cristianl/data/titan/meps_gmted2010_1km_topo_topdown.nc",
           dem.varname="altitude",
           proj4dem="+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06",
           laf.sct=T,
           laf.file="/home/cristianl/data/titan/meps_gmted2010_1km_laf_topdown.nc",
           # fg (MEPS)
           fg=T,
           fg.file="/home/cristianl/data/titan/meps_2_5km_20180403T12Z_test.nc",
           fg.type="meps",
           fg.cfact=100, #MEPS RH 0-1 -> 0-100%
           thr.fg=70,
           # fge (MEPS)
           fge=T,
           fge.file="/home/cristianl/data/titan/meps_2_5km_20180403T12Z_test.nc",
           usefge.sct=T,
           fge.cfact=100,
           fge.type="meps",
           sdmin.fge=10,
           iqrmin.fge=10,
           csd.fge=3,
           infsd.fge=3,
           ciqr.fge=3,
           infiqr.fge=3,
           # buddy-check
           i.buddy=4,
           dr.buddy=3000, #m
           thr.buddy=3,
           n.buddy=5,
           dz.buddy=500, #m
           # steve
#           steve=T,
#           i.steve=3, 
#           thres.steve=c(10,90),
#           pmax_lt.steve=c(20,20),
#           dmax_lt.steve=c(150000,150000),
#           n_lt.steve=c(4,4),
#           frac_lt.steve=c(0.2,0.2),
#           dmin_next_lt.steve=c(150000,150000),
#           pmax_ge.steve=c(20,20),
#           dmax_ge.steve=c(150000,150000),
#           n_ge.steve=c(4,4),
#           frac_ge.steve=c(0.2,0.2),
#           dmin_next_ge.steve=c(100000,100000),
           # SCT
           i.sct=3,
           n.sct=50,
           dz.sct=1,
           DhorMin.sct=5000, #m
           Dver.sct=2000, #m
           eps2.sct=0.1,
           thr.sct=9,
           fast.sct=T,
           min.corep=0.5,
           mean.corep=1,
           max.corep=2,
           # dem
           dem=T,
           dz.dem=1500, #m
           # plausibility
           vmin=0,  #%
           vmax=100 #%
)
