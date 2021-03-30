#+ Correction for wind-induced undercatch of precipitation (Wolff et al., 2015)
wolff_correction<-function(par,rr,t2m,ws10m) {
# Wolff MA, Isaksen K, Petersen-Øverleir A, Ødemark K, Reitan T, Brækkan R. 
# Derivation of a new continuous adjustment function for correcting 
# wind-induced loss of solid precipitation: results of a Norwegian field study.
# Hydrology and Earth System Sciences. 2015 Feb 1;19(2):951.
#------------------------------------------------------------------------------
  theta<-par[1]
  beta<-par[2]
  tau1<-par[3]
  tau2<-par[4]
  Ttau<-par[5]
  stau<-par[6]
  sig1<-par[7]
  sig2<-par[8]
  Tsig<-par[9]
  ssig<-par[10]
  n<-length(t2m)
  aux0<-exp((t2m-Ttau)/stau)
  aux1<-aux0/(1+aux0)
  f<- ( 1 - tau1 - (tau2-tau1) * aux1 ) * exp(-(ws10m/theta)**beta) + 
     tau1 + (tau2-tau1) * aux1
  aux2<-exp((t2m-Tsig)/ssig)
  s<-sig1+(sig2-sig1)*(aux2)/(1+aux2)
  rr.cor<-rr*1./f
  return(list(f=f,sigma=s,rr.cor=rr.cor))
}

