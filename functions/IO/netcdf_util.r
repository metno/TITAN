
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

#+ wraper, read data from nc-file and format it
get_data_from_ncfile<-function(nc.file,
                               nc.varname,
                               nc.t=NA,
                               nc.e=NA,
                               topdown,
                               var.dim,
                               var.de_acc=F,
                               var.de_acc_by="-1 hour",
                               proj4,
                               proj4_from_nc,
                               xy_as_vars=F,
                               x_as_var.varname=NA,
                               y_as_var.varname=NA,
                               xy_as_var.dim=NA,
                               xy_as_var.dh_max=NA,
                               return_raster=F,
                               return_obsloc=T,
                               extent=NA,
                               debug.file=NA,
                               nc.dqc_mode="none") { #"radar_hourly"
#------------------------------------------------------------------------------
  val <- NA
  if (var.dim$ndim==2) { var.dim$tpos<-NULL; var.dim$epos<-NULL }
  if (is.null(var.dim$tpos)) {
    nc.t<-NULL
  } else {
    if (is.na(var.dim$tpos)) {
      nc.t<-NULL
      var.dim$tpos<-NULL
    } else {
      ti<-nc4.getTime(nc.file)
      if (is.na(nc.t)) nc.t<-ti[1]
      if (!(nc.t %in% ti)) 
        boom(paste("ERROR, file",nc.file,",timestamp",nc.t,"not present"))
    }
  }
  if (is.null(var.dim$epos)) {
    nc.e<-NULL
  } else {
    if (is.na(var.dim$epos)) {
      nc.e<-NULL
      var.dim$epos<-NULL
    } else {
      if (is.na(nc.e)) nc.e<-0
    }
  }
  if (var.de_acc) {
    tminus1h<-format(as.POSIXlt(
                     seq(as.POSIXlt(strptime(nc.t,"%Y%m%d%H%M",tz="UTC")),
                         length=2,by=var.de_acc_by),"UTC")[2],
                     "%Y%m%d%H%M",tz="UTC")
    if (!(tminus1h %in% ti)) 
      boom(paste("ERROR, file",nc.file,",timestamp",tminus1h,"not present"))
    raux<-try(nc4in(nc.file=nc.file,
                    nc.varname=nc.varname,
                    topdown=topdown,
                    out.dim=var.dim,
                    proj4=proj4,
                    nc.proj4=proj4_from_nc,
                    selection=list(t=c(tminus1h,nc.t),e=nc.e)))
    if (is.null(raux)) boom(paste("ERROR while reading \"var\" from file:",argv$nc.file))
    r<-raster(raux$stack,"layer.2")-raster(raux$stack,"layer.1")
  } else {
    raux <- try( nc4in( nc.file=nc.file,
                        nc.varname=nc.varname,
                        topdown=topdown,
                        out.dim=var.dim,
                        proj4=proj4,
                        nc.proj4=proj4_from_nc,
                        selection=list(t=nc.t,e=nc.e)))
    r<-raux$stack
  }
  rm(raux)
  proj4_nc<-as.character(crs(r))
  if (nc.dqc_mode=="radar_hourly") {
    if (argv$verbose) print("Read radar and do the radar-DQC")
    t0a<-Sys.time()
    suppressPackageStartupMessages(library("igraph"))
    rval<-getValues(r)
    # a. remove not plausible values
    if (argv$verbose) print(" remove not plausible values")
    radardqc.min<-0
    radardqc.max<-300
    ix<-which( !is.na(rval) & (rval<radardqc.min | rval>radardqc.max) )
    if (length(ix)>0) {
      rval[ix]<-NA
      r[]<-rval
    }
    rm(ix)
    # b. remove patches of connected cells that are too small
    #  check for small and isolated clumps (patches) of connected cells with 
    #  precipitation greater than a predefined threshold
    #   threshold 0 mm/h. remove all the clumps made of less than 100 cells
    #   threshold 1 mm/h. remove all the clumps made of less than 50 cells
    if (argv$verbose) print(" remove small clumps")
    radardqc.clump.thr<-c(0,1)
    radardqc.clump.n<-c(100,50)
    for (i in 1:length(radardqc.clump.thr)) {
      raux<-r
      if (any(rval<=radardqc.clump.thr[i])) 
        raux[which(rval<=radardqc.clump.thr[i])]<-NA
      rclump<-clump(raux)
      fr<-freq(rclump)
      ix<-which(!is.na(fr[,2]) & fr[,2]<=radardqc.clump.n[i])
      if (length(ix)>0) {
        rval[getValues(rclump) %in% fr[ix,1]]<-NA
        r[]<-rval
      }
      rm(raux,fr,ix,rclump)
    }
    # c. remove outliers. Check for outliers in square boxes of 51km by 51km
    raux<-r
    daux<-boxcox(x=rval,lambda=0.5)
    raux[]<-daux
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
    if (ncell(raux)>ncell(r)) {
      raux<-crop(raux,r)
    } else if (ncell(raux)<ncell(r)) {
      raux<-extend(raux,r)
    }
    avg<-getValues(raux)
    raux_agg[]<-NA; raux_agg[ix_aux]<-arr[2,]
    raux<-disaggregate(raux_agg,fact=5,method="bilinear",na.rm=T)
    if (ncell(raux)>ncell(r)) {
      raux<-crop(raux,r)
    } else if (ncell(raux)<ncell(r)) {
      raux<-extend(raux,r)
    }
    stdev<-getValues(raux)
    ix<-which(stdev>0 & !is.na(daux) & !is.na(avg) & !is.na(stdev))
    rm(arr,raux_agg,ix_aux,xrad_aux,yrad_aux,vrad_aux,daux_agg,xyaux)
    # outliers are defined as in Lanzante,1997: abs(value-mean)/st.dev > 5
    suspect<-which((abs(daux[ix]-avg[ix])/stdev[ix])>5) 
    if (length(suspect)>0) rval[ix[suspect]]<-NA
    r[]<-rval
    rm(raux,daux,avg,stdev,ix,suspect,rval)
    t1a<-Sys.time()
    if (argv$verbose) print(paste(" remove outliers - time",round(t1a-t0a,1),
                                                       attr(t1a-t0a,"unit")))
  }
  #
  if (xy_as_vars) {
    ix<-which(!is.na(getValues(r)))
    if (length(ix)==0) {
      val<-rep(NA,ndata)
    } else {
      if (!is.null(xy_as_var.dim$tpos)) {
        if (is.na(xy_as_var.dim$tpos)) {
          nc.t<-NULL
          xy_as_var.dim$tpos<-NULL
        } else {
          ti<-nc4.getTime(nc.file)
          if (is.na(nc.t)) nc.t<-ti[1]
          if (!(nc.t %in% ti)) 
            boom(paste("ERROR, file",nc.file,",timestamp",nc.t,"not present"))
        }
      }
      if (!is.null(xy_as_var.dim$epos)) {
        if (is.na(xy_as_var.dim$epos)) {
          nc.e<-NULL
          xy_as_var.dim$epos<-NULL
        } else {
          if (is.na(nc.e)) nc.e<-0
        }
      }
#      rval<-getValues(r)[ix]
      rval<-getValues(r)
      raux<-try(nc4in(nc.file=nc.file,
                      nc.varname=x_as_var.varname,
                      topdown=topdown,
                      out.dim=xy_as_var.dim,
                      proj4=proj4,
                      nc.proj4=proj4_from_nc,
                      selection=list(t=NULL,e=NULL)))
      if (is.null(raux)) boom(paste("ERROR while reading \"x_as_var\" from file:",argv$nc.file))
#      rx<-getValues(raux$stack)[ix]; rm(raux)
      rx<-getValues(raux$stack); rm(raux)
      raux<-try(nc4in(nc.file=nc.file,
                      nc.varname=y_as_var.varname,
                      topdown=topdown,
                      out.dim=xy_as_var.dim,
                      proj4=proj4,
                      nc.proj4=proj4_from_nc,
                      selection=list(t=NULL,e=NULL)))
      if (is.null(raux)) boom(paste("ERROR while reading \"y_as_var\" from file:",argv$nc.file))
#      ry<-getValues(raux$stack)[ix]; rm(raux)
      ry<-getValues(raux$stack); rm(raux)
      if (is.na(xy_as_var.dh_max)) 
        xy_as_var.dh_max<-10*mean( median(abs(diff(rx)),na.rm=T), 
                                   median(abs(diff(ry)),na.rm=T) )  
      val<-apply(cbind(data$lon,data$lat),FUN=spint_nn,MARGIN=1,
                 rx=rx,ry=ry,rval=rval,dh_max=xy_as_var.dh_max)
    }
  } else if (proj4_nc!=argv$proj4_where_dqc_is_done) {
    coord<-SpatialPoints(cbind(data$lon,data$lat),
                         proj4string=CRS(argv$proj4_input_obsfiles))
    coord.new<-spTransform(coord,CRS(proj4_nc))
    xy.tmp<-coordinates(coord.new)
    if ( return_obsloc) val<-extract(r,xy.tmp)
    rm(coord,coord.new,xy.tmp)
  } else {
    if ( return_obsloc) val<-extract(r,cbind(x,y))
  }
  if (!is.na(debug.file)) save(r,val,file=debug.file)
  if (return_raster) {
    if ( !any( is.na(extent))) {
      coord <- SpatialPoints( cbind( c( extent[1], extent[1], extent[2], extent[2]),
                                     c( extent[3], extent[4], extent[3], extent[4])),
                           proj4string=CRS("+proj=longlat +datum=WGS84"))
      coord.new<-spTransform(coord,CRS(proj4_nc))
      xy.tmp<-coordinates(coord.new)
      extentxy <- as( extent( min(xy.tmp[,1]), max(xy.tmp[,1]), 
                              min(xy.tmp[,2]), max(xy.tmp[,2])), 'SpatialPolygons')
      crs( extentxy) <- crs(r) 
      r <- crop( r, extentxy)
    }
    return(list(raster=r,values=val))
  } else {
    return(val)
  }
}

