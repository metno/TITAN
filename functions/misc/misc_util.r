
# + manage fatal error
boom <- function( str=NA, code=NA) {
  cat("Fatal Error ")
  if ( !is.na(code)) {
    if ( code == 1) cat("file not found ")
  }
  if ( !is.na(str)) cat( str)
  cat("\n")
  quit( status= 1)
}

#+ the end 
rip <- function( str=NA, code=NA, t0=NA) {
  cat( "the End : ")
  if ( !is.na(code) ) {
    if ( code == 1 ) cat( "normal exit : ")
  }
  if ( !is.na(t0)) {
    t1 <- Sys.time()
    cat( paste( "total time=", round(t1-t0,1), attr(t1-t0,"unit")))
  }
  if ( !is.na(str)) cat( str)
  cat("\n")
  quit( status= 0 )
}

#+ convert (input) strings to numeric values
strings_to_numbers <- function( strings,
                                default=NA,
                                strings_dim=1,
                                neg=NA) {
# strings is a vector of characters
#------------------------------------------------------------------------------
  options( warn = 1)
  # string is a one-dimensional
  if ( strings_dim == 1) {
    if ( is.na( strings)) strings <- default
    numbers <- as.numeric( gsub( "_", "-", strings))
    if ( !is.na( neg)) numbers <- numbers * (-1)**as.numeric(neg)
  # string is a vector
  } else {
    if ( any( is.na( strings))) {
      numbers <- rep( default, length=strings_dim)
    } else {
      if ( length( strings) != strings_dim) 
        strings <- rep( strings[1], length=strings_dim)
      aux <- vector( length=strings_dim, mode="numeric")
      for (i in 1:strings_dim) aux[i] <- as.numeric( gsub( "_", "-", strings[i]))
      numbers <- aux
      rm( aux)
      if ( !any( is.na( neg))) {
        if ( length( neg) != strings_dim) neg <- rep( neg[1], length=strings_dim)
        for (i in 1:strings_dim) numbers[i] <- numbers[i] * (-1)**as.numeric(neg[i])
      }
    }
  }
  options( warn = 2)
  numbers
}

# auxiliary function to keep/blacklist observations
setCode_lonlat<-function(lonlat,code,datatmp) {
# lonlat. vector. 1=lon; 2=lat
  ix<-which(datatmp$lon==lonlat[1] & datatmp$lat==lonlat[2])
  if (length(ix)>0)  {
    aux[ix]<-code
    assign("aux",aux,envir=.GlobalEnv)
  }
  return(length(ix))
}

