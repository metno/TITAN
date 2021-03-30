#!/usr/bin/env Rscript

# Test TITAN 

system("export TITANR_PATH=$HOME/projects/titanlab/R/functions; ../titan.r --input.files data/observation_test_ta_prid01_p02000_pGE020percent.txt data/observation_test_ta_prid02_p00100_pGE003percent.txt --output.file data/out.txt --config.files ini/ta_test_titan.ini ini/ta_sct_fg.ini --fg.files ini/background_test_ta_det.ini ini/background_test_ta_ens.ini")

# Read DQC output
dat <- read.table(file="data/out.txt",header=T,sep=";",stringsAsFactors=T,strip.white=T)

a <- length( which( dat$ge == 1 & dat$dqc == 2))
c <- length( which( dat$ge == 1 & dat$dqc == 0))
b <- length( which( dat$ge == 0 & dat$dqc == 2))
d <- length( which( dat$ge == 0 & dat$dqc == 0))
rand <- (a+c) * (a+b) / (a+b+c+d)
ets  <- (a-rand) / (a+b+c-rand)
acc  <- (a+d)/(a+b+c+d)
pod  <- a/(a+c)
pofa <- b/(b+d)

print( paste("a(bad) b c d", a,"(",length(which( dat$ge==1)),")", b, c, d, a+b+c+d))
print( paste("acc pod pofa ets", round(acc,2), round(pod,2), round(pofa,2), round(ets,2)))

# plots
par(mar=c(3,3,1,1))
plot( dat$lon, dat$lat, xlab="", ylab="", col="white")
br<-c(-60,seq(-15,3,length=12),60)
col<-rev(rainbow(length(br)-1))
for (i in length(col):1) {
  if ( length( ix<- which( dat$value>br[i] & dat$value<=br[i+1])) > 0) 
    points( dat$lon[ix], dat$lat[ix], pch=21, bg=col[i],col=col[i], cex=2)
  if ( length( ix<- which( dat$value>br[i] & dat$value<=br[i+1] & dat$ge == 1)) > 0) 
    points( dat$lon[ix], dat$lat[ix], pch=21, bg=col[i],col="red", cex=2,lwd=3)
  if ( length( ix<- which( dat$value>br[i] & dat$value<=br[i+1] & dat$dqc == 2)) > 0) 
    points( dat$lon[ix], dat$lat[ix], pch=4, col="red", cex=2)
}

q()
