#!/usr/bin/env Rscript
# + TITAN - spatial data quality control for in-situ observations
# mailto: cristianl@met.no
# https://github.com/metno/TITAN
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
suppressPackageStartupMessages( library( "argparser"))
suppressPackageStartupMessages( library( "sp"))
suppressPackageStartupMessages( library( "raster"))
suppressPackageStartupMessages( library( "rgdal"))

options( warn = 2, scipen = 999)

#
#==============================================================================
#  MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN*MAIN
#==============================================================================

t0 <- Sys.time() # game on

#
#-----------------------------------------------------------------------------
# path to the titan functions is stored in the enviroment var TITANR_FUN

titan_fun_path <- Sys.getenv( "TITANR_PATH")

#
#-----------------------------------------------------------------------------
# load functions

fun_file <- file.path( titan_fun_path, "fun_list.r")
if ( !file.exists( fun_file)) boom( fun_file, code=1)

source( fun_file)

for (fun in fun_list) {
  if ( !file.exists(file.path( titan_fun_path, fun)))
    boom( file.path( titan_fun_path, fun), code=1)
  source( file.path( titan_fun_path, fun))
}

rm( fun_file, fun_list, fun)               

# Alternatively, more elegant but difficult to find missing functions
#for (file in list.files(path = titan_fun_path, pattern = ".r", full.names=T) ) 
#  source(file)
#rm(file)

#
#-----------------------------------------------------------------------------
# define constants

const_file <- file.path( titan_fun_path, "constants.r")

if ( !file.exists( const_file)) boom( const_file, code=1)

source( const_file)

#
#-----------------------------------------------------------------------------
# read command line arguments and/or configuration file

# environment where to store background fields
fg_env    <- new.env( parent=emptyenv())
fg_env$fg <- list()

argv <- argparser()

#
#-----------------------------------------------------------------------------
# Load titanlib

dyn.load( file.path( argv$titanlib_path, paste0( "titanlib", .Platform$dynlib.ext)))

source( file.path( argv$titanlib_path, "titanlib.R"))

#
#-----------------------------------------------------------------------------
# Multi-cores run

if ( !is.na( argv$cores)) {
  suppressPackageStartupMessages( library( "parallel"))
  if ( argv$cores==0) argv$cores <- detectCores()
  cat( paste( "--> multi-core run, cores=", argv$cores, "\n"))
}

#
#-----------------------------------------------------------------------------
# read data
#
res         <- read_data_to_check( argv)
extent      <- res$extent
data        <- res$data
dqcflag     <- res$dqcflag
z           <- res$z
dataopt     <- res$dataopt
varidx      <- res$varidx
varidx.opt  <- res$varidx.opt
rm(res)

#
#-----------------------------------------------------------------------------
# test for no metadata (1st round) 
#  1st round, check for missing metadata in the original data
#  2nd round, check for missing metadata after dem.fill 

dqcflag.bak <- dqcflag # bakup, used in main_read_dem.r (if dem.fill=T)

dqcflag <- metadata_check_r ( argv, data, z, extent, dqcflag)

#
#-----------------------------------------------------------------------------
# coordinate transformation

res <- spatconv( argv, data, extent)
x  <- res$x
y  <- res$y
xl <- res$xl
yl <- res$yl
e  <- res$e
rm(res)

#
#-----------------------------------------------------------------------------
# Read geographical information (optional) 
# digital elevation model

if (argv$dem | argv$dem.fill) { 
  res     <- read_dem( argv, data, z, dqcflag, dqcflag.bak)
  z       <- res$z    # station elevations (either from input or dem) 
  zdem    <- res$zdem # dem elevation
  dqcflag <- res$dqcflag
  rm(res, dqcflag.bak)
}

#
#-----------------------------------------------------------------------------
# precipitation (in-situ) and temperature (field) cross-check

if (argv$ccrrt) {
  res     <- ccrrt( argv, data, z, dqcflag)
  dqcflag <- res$dqcflag
  if (argv$rr.wcor) { t2m <- res$t2m } else { t2m <- NULL }
  rm(res)
}

#
#-----------------------------------------------------------------------------
# Correction for the wind-undercatch of precipitation 

if (argv$rr.wcor) 
  data  <- rr_windcorr( argv, data, z, dqcflag, t2m) 

#
#-----------------------------------------------------------------------------
# Read first guess fields

if ( length( fg_env$fg) > 0) 
  res <- read_fgs( argv, extent)

#
#-----------------------------------------------------------------------------
# test for no metadata (2nd and final) 

dqcflag <- metadata_check_r ( argv, data, z, extent, dqcflag)

#
#-----------------------------------------------------------------------------
# check elevation against dem 
# NOTE: keep-listed stations canNOT be flagged here

if (argv$dem) 
  dqcflag <- check_z_against_dem( argv, data, z, zdem, dqcflag )

#
#-----------------------------------------------------------------------------
# plausibility test
# NOTE: keep-listed stations could be flagged here

dqcflag <- plausibility_test( argv, data, dqcflag)

#
#-----------------------------------------------------------------------------
# climatological check 
# NOTE: keep-listed stations canNOT be flagged here
# use only (probably) good observations

if ( !is.na( argv$month.clim))
  dqcflag <- climatological_check( argv, data, dqcflag)

#
#-----------------------------------------------------------------------------
# SCT for dichotomous (yes/no) variables with the background

if (argv$sct_fg_dual) 
  dqcflag <- sct_fg_dual_r( argv, data, x, y, z, dqcflag)

#
#-----------------------------------------------------------------------------
# SCT for dichotomous (yes/no) variables

if (argv$sct_dual) 
  dqcflag <- sct_dual_r( argv, data, z, dqcflag)

#
#-----------------------------------------------------------------------------
# check against first-guess fields

if (argv$fgt) 
  dqcflag <- fgt_r( argv, data, x, y, z, dqcflag)

#
#-----------------------------------------------------------------------------
# buddy check 
#  compare each observation against the neighbouring observations 
# NOTE: keep-listed stations are used but they canNOT be flagged here

if (argv$buddy)
  dqcflag <- buddy( argv, data, z, dqcflag)

#
#-----------------------------------------------------------------------------
# SCT - Spatial Consistency Test, using background fields
# NOTE: keep-listed stations are used but they canNOT be flagged here

if (argv$sct_fg)
  dqcflag <- sct_fg_resistant( argv, data,  x, y, z, dqcflag)

#
#-----------------------------------------------------------------------------
# SCT - Spatial Consistency Test
# NOTE: keep-listed stations are used but they canNOT be flagged here

if ( argv$sct)
  dqcflag <- sct_resistant_r( argv, data, z, dqcflag)

#
#-----------------------------------------------------------------------------
# check for isolated stations
# use only (probably) good observations

if ( argv$iso)
  dqcflag <- isolation_test( argv, data, dqcflag)

#
#-----------------------------------------------------------------------------
# observations not flagged are assumed to be good observations 

dqcflag <- final_decision( data, dqcflag)

#
#-----------------------------------------------------------------------------
# write the output file

write_output( argv, data, dqcflag, dataopt, varidx, varidx.opt)
 
#
#-----------------------------------------------------------------------------
# Normal exit

rip( code=1, t0=t0) # exit status is 0
