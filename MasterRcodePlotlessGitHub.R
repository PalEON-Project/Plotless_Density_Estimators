
library(spatstat)

#setwd("C:Users/Cogbill/Documents/Charlie/R in PalEON/RSimMap)
getwd()

P='Title '
# n = number in aggregation pool
n<-100
# S = number of replications
S<-2000
#Cy<-10
# bord = proportion of unit square buffered border
bord<-.10
# q = basic density parameter per unit square
q<-300
# samspa = spacing parameter for sampling in proportion of unit (100 m)
samspa<-.01
# minden = minimum density parameter gradient
minden<-20
# maxden = maximum density parameter in gradient
maxden<-1000
# expn = exponent parameter in gradient 
expn<-2
# mat = MaternII parameter inhibition distance parameter in proportion of of unit (100 m)
mat<-0.015
#  mu = rThomas parameter of number of events within clumps
mu<-2
#  gamma = rStrauss parameter proportion of pairs pruned 
gamma<-0.0
# sigma = rStrauss parameter of inhibition distance proportion of unit (100 m)
sigma<-0.015
#sigma<-mat
# gt - rank order for Morisita g-tree, Prodan, Kleinn-Vilcko estimators
gt<-6
scale<-1
bias<-.015
rancut<-1



rows<-as.integer(1/samspa)
colm<-as.integer(1/samspa)
GX<-rows-4
GY<-colm-4
GX<-10
GY<-10
#n<-GX*GY


distveil<-.00/sqrt(scale)
veilcut<-1

bull<-.005
mind<-minden
kick<-.0
#geocut2<-bias
#geocut3p<-.0



#geocut2<-bias*(1-bias)
#geocut3p<-(bias-geocut2)


nearcut1<-bias
nearcut2<-bias
nearcut3<-bias
nearcut4<-bias
nearPcut1<-bias
nearPcut2<-bias
addrank<-2
baserank<-1




ciel<-maxden
r<-50
# a = exponent parameter in density gradient on x-axis
a<-2
# b =  exponent parameter in density gradient on y-axis
b<-2

#mat<-0.00
#bord<-.10
# qNE = basic density parameter for gradient in 1st quadrant 
qNE<-35
# qSE = basic density parameter for gradient in 2nd quadrant 
qSE<-70
# qSW = basic density parameter for gradient in 3rd quadrant 
qSW<-35
# qNW = basic density parameter for gradient in 4th quadrant 
qNW<-70

#mat<-.015/sqrt(scale)


# matNE ...  = inhibition parameter in quadrant
matNE<-mat
matSE<-mat
matSW<-mat
matNW<-mat

# addNE ... = gradient parameter in quadrant
addNE<-1
addSE<-1
addSW<-1
addNW<-1

#  NEgam ... = gradient parameter in quadrant
NEgam<-1500
SWgam<-1500
SEgam<-1500
NWgam<-1500

#NElam<-500
#SElam<-500
#SWlam<-500
#NWlam<-500



#a<-2
#b<-2
#cc<-1000*a*b/(exp(a)-1)/(exp(b)-1)
#m<-cc/(a*b)*(exp(a)-1)*(exp(b)-1)*pi
#cc<-q*a*b/exp(a)-1/exp(b)-1)
#cc<-q
#range.IPlow<-cc*(exp(a+b)-1)
#lambda=function(x,y){minden+q*exp(a*x)*exp(b*y)}
#lambda<-600
#InH<-rpoispp(lambda,win=owin(c(.75,1),c(0,1)))

# NElam ... = density parameter in gradient quadrant  
NElam=function(x,y){minden+qNE*exp(a*(2*x-1))*exp(a*(2*y-1))}
SElam=function(x,y){minden+qSE*exp(a*(2*x-1))*exp(a*(1-2*y))}
SWlam=function(x,y){minden+qSW*exp(b*(1-2*x))*exp(b*(1-2*y))}
NWlam=function(x,y){minden+qNW*exp(b*(1-2*x))*exp(b*(2*y-1))}



  par(mfrow=c(1,1))
# InNE ...  = MaternII pruned simulated pattern in quadrant
InNE<-rMaternII(NElam, matNE, win=owin(c(.5,1),c(.5,1)))
InSE<-rMaternII(SElam, matSE, win=owin(c(.5,1),c(0,.5)))
InSW<-rMaternII(SWlam, matSW, win=owin(c(0,.5),c(0,.5)))
InNW<-rMaternII(NWlam, matNW, win=owin(c(0,.5),c(.5,1)))
#  = simulated gradient pattern combining 4 quadrants
eg.IP<-superimpose(InNE,InSE,InSW,InNW, W=owin(c(0,1), c(0,1)))

# = simulated pattern of Strauss inhibition parameters density, proportion, inhibition distance
#eg.IP<-rStrauss(q, gamma, mat)

#SSNE<-rSSI(addNE, n=NEgam, win=owin(c(.5,1),c(.5,1)), giveup=1000, x.init=InNE)
#SSSE<-rSSI(addSE, n=SEgam, win=owin(c(.5,1),c(0,.5)),giveup=1000,  x.init=InSE)
#SSSW<-rSSI(addSW, n=SWgam, win=owin(c(0,.5),c(0,.5)), giveup=1000, x.init=InSW)
#SSNW<-rSSI(addNW, n=NWgam, win=owin(c(0,.5),c(.5,1)),giveup=1000,  x.init=InNW)
#  = simulated SSI gradient pattern combining 4 quadrants
#eg.IP<-superimpose(SSNE,SSSE,SSSW,SSNW, W=owin(c(0,1), c(0,1))) 

#  = simulated pattern of homogenous Poisson process parameters density
#eg.IP<-rpoispp(q)

# = simulated pattern of Matern II inhibition parameters density, inhibition distance
#eg.IP<-rMaternII(q,mat)

# = simulated pattern of rThomas clumping parameters density, sd of distance, number in clump
#eg.IP<-rThomas(q,mat,mu)



##plot(SSNE)
##plot(SSSE)
##plot(SSSW)
##plot(SSNW)
par(mfrow=c(1,1),mar=c(2,2,1,1))

plot(eg.IP)
win<-owin(c(bord,1-bord),c(bord,1-bord))
#w<-owin(c(.1,.9),c(.1,.9))
#clarkevans.test(InH, clipregion=win    , correction="guard")
#dual_InH<-InH
#plot(dual_InH)

rect(bord,bord,bord,1-bord,col='red',lwd=2)
rect(bord,1-bord,1-bord,1-bord,col='red',lwd=2)
rect(1-bord,1-bord,1-bord,bord,col='red',lwd=2)
rect(1-bord,bord,bord,bord,col='red',lwd=2)

#clarkevans.test(InH,clipregion=owin(c(bord,1-bord),c(bord,1-bord)))


#eg.IP<-rStrauss(300, 0, R=.03)

#BNE<-rSSI(matNE, n=Inf , win=owin(c(.5,1),c(.5,1)), giveup=1000, x.init=InNE)

#BSE<-rSSI(matSE, n=Inf , win=owin(c(.5,1),c(0,.5)), giveup=1000, x.init=InSE)

#BSW<-rSSI(matSW, n=Inf , win=owin(c(0,.5),c(0,.5)), giveup=1000, x.init=InSW)

#BNW<-rSSI(matNW, n=Inf , win=owin(c(0,.5),c(.5,1)), giveup=1000, x.init=InNW)




#SSI<-rSSI(r=sigma,n=Inf,giveup=1000, win=owin(c(0,.5), c(0,1)))
#
#InNE<-rpoispp(NElam, win=owin(c(.5,1), c(0,1)))
#InSE<-rpoispp(SElam, win=owin(c(.5,1), c(0,1)))
#InSW<-rpoispp(SWlam, win=owin(c(.5,1), c(0,1)))
#InNW<-rpoispp(NWlam, win=owin(c(.5,1), c(0,1)))


#plot(InSE)
#plot(InSW)
#plot(InNW)


                 
#par(mfrow=c(1,1))


#eg.IP<-rStrauss(q, gamma, sigma)
#eg.IP<-rpoispp(q)
#eg.IP<-rpoispp(lambda)
#eg.IP<-rSSI(r=sigma, n=ciel, giveup=1000, win=owin(c(0,1),c(0,1)), x.init=InH)
 
#plot(eg.IP)

#eg.IP<-rpoispp(lambda=function(x,y){cc*(a*x)})
#zeta<-as.im(function(x,y){2*exp(4*x-1)},owin(c(0,1),c(0,1)))
# 
#   SIMULATION of outcome homogenous Poisson process intensity of (q) per unit sq 
#   (if unit is 100m ==  on 100X100m square per /ha  or q t/ha)
print('t/ha')

# el ...  = east, west, south, & north proportion in border
el<-1-bord
wl<-bord
sl<-bord
nl<-1-bord


#min den plots=g, max den plots=e;  min focus =v; max focus=w
g<-0
e<-1500
u<-0
w<-1500


#zeta<-as.im(function(x,y){2*exp(2*x-1)},owin(c(0,1),c(0,1)))
#eg.IP<-rThomas(q,sigma, zeta)
#eg.IP<-rpoispp(lambda=function(x,y){cc*exp(a*x+b*y)})

# Correction constant for two nearest quadrats


#  std.dat  = empirical pattern imported from external file
#std.dat<-ppp(Raw[ ,2],Raw[ ,3],c(min(Raw[,2]),max(Raw[,2])),
 #            c(min(Raw[,3]),max(Raw[,3])))
#summary(std.dat)

#Stanclean<-unique(Stand)
#any(duplicated(Stanclean))
#sim.dat<-rpoispp(q/pi)

# E.g. unbiased survey

#eg<-dist.pts2(std.dat$x,std.dat$y)

#library(fields)
#library(gtools)



#sim.dat<-rpoispp(q)

# E.g. unbiased survey

# Definition and Initilation of vectors
truen<-rep(0,n)
true<-rep(0,S)
rbias<-function(est,true){mean(est/true)}
diff<-function(est,true){est-true}
me<-function(error){mean(error)}
rme<-function(error,true){mean(error/true)}
mae<-function(error){mean(abs(error))}
rmae<-function(error,true){mean(abs(error)/true)}
rrmse<-function(error,true){sqrt(mean(error^2/true^2))}


#par(mfrow=c(1,1))

# The forest

#plot(sim.dat$y~sim.dat$x,pch=20,cex=.5)

#intensityfn.HP<-matrix(m/pi,length(xvals),length(yvals))

# (1) Homogeneous Poisson point process

#m<-q*pi

#eg.HP<-rpoispp(m/pi)


# The survey point



# The nearest tree (blue) and the biased nearest tree (green) 
# in the opposite half (the cyan tree should have been sampled,
# but the green tree was sampled instead)
#points(eg$surveypt[2]~eg$surveypt[1],col='red',pch=20)
#abline(v=eg$surveypt[1],col='gray', lty=2)
#points(eg$point1[2]~eg$point1[1],col='blue',pch=20)
#points(eg$point2[2]~eg$point2[1],col='green',pch=20)
#points(eg$censoredpoint[2]~eg$censoredpoint[1],col='cyan',pch=20)

############################################################

#-----Point Half-----Point Half-----Point Half-----Point Half-----Point Half-----
#-----Point Half-----Point Half-----Point Half-----Point Half-----Point Half-----


# Simulating different point patterns and 
# estimating stem density in UNBIASED surveys.

# Definition and Initilation of vectors

Poisden<-rep(0,S)
rat2.HP<-rep(0,S)
moratio.HP<-rep(0,S)

rat1.HP<-rep(0,n)
dist1.HP<-rep(0,S)
dist2.HP<-rep(0,S)
rofmean.HP<-rep(0,S)
hom.HP<-rep(0,S)
inhom.HP<-rep(0,S)
cott.HP <- rep(0,S)
shvc.HP <- rep (0,S)
eph.HP <-rep(0,S)
inhom.TQ<-rep(0,S)
hom.TQ<-rep(0,S)
cott.TQ <- rep(0,S)
shvc.TQ <- rep (0,S)
eph.TQ <-rep(0,S)


avrbar.HP <-rep(0,S)
avrsqbar.HP <-rep(0,S)
avrbar.QU <-rep(0,S)
avrsqbar.QU <-rep(0,S)

xvals<-seq(0,1,.01)
yvals<-seq(0,1,.01)
dia.nearest<-rep(0,n)
dist.nearest<-rep(0,n)
dia.dist.second<-rep(0,n)
dist.dist.second<-rep(0,n)
point.dist.mean.HP<-rep(0,n)
ave.dist.point.HP<-rep(0,n)

pp.dat<-matrix(0,n,3)
pp.x<-rep(0,n)
pp.y<-rep(0,n)
ratio.HP<-matrix(0,n,2)

avBA.HP<-rep(0,n)
avdia.HP<-rep(0,n)
avdist.HP<-rep(0,n)
sumofsq.HP<-rep(0,n)
BAmean.HP<-rep(0,S)
Mori.HP<-rep(0,S)
Mori.BA.HP<-rep(0,S)
Cott.BA.HP<-rep(0,S)
Cott.Dia.HP<-rep(0,S)
Cott.Area.HP<-rep(0,S)
Xprod.HP<-rep(0,n)
aver.dist.point.BH2X<-rep(0,n)
aver.dist.point.HP<-rep(0,n)
nearest.Sq.HP<-rep(0,n)
HP.dia<-matrix(0,n,2)
nearest.dist.HP<-rep(0,n)
second.dist.HP<-rep(0,n)
nearest.dia.HP<-rep(0,n)
second.dia.HP<-rep(0,n)
second.HP<-matrix(0,n,4)
point.HP<-matrix(0,n,8)
sum.HP<-matrix(0,S,8)
sec.dia.HP<-rep(0,n)
BH2X.4dist<-matrix(0,n,2)
BX2.2order.BH2X<-rep(0,n)
BH2X<-rep(0,n)
B2X<-rep(0,n)
sumofsq.BH2X<-rep(0,n)
avdist.BH2X<-rep(0,n)
ave.dist.point.BH2X<-rep(0,n)
inhibden<-rep(0,S)
shanks1<-rep(0,S)
trip3<-rep(0,n)
Trip<-rep(0,S)
QU2X<-rep(0,n)
nearest.dist.QU2X<-rep(0,n)
cott.BH2X<-rep(0,n)
hom.BH2X<-rep(0,S)
inhom.BH2X<-rep(0,S)
cott.BH2X<-rep(0,S)
shvc.BH2X<-rep(0,S)
eph.BH2X<-rep(0,S)
BHX<-matrix(0,2,2)
BH2X.2Q2tree<-rep(0,n)   
BH2X.1Q2tree<-rep(0,n)

cott.EQ<- rep(0,S)  
hom.EQ<- rep(0,S)
inhom.EQ <- rep(0,S)
shvc.EQ <- rep(0,S)
EQ.all<-matrix(0,n,3)

avdist.HP2X<-rep(0,n)
avdist.QU<-rep(0,n)
avdist.QU2X<-rep(0,n)
avdist.QUNX<-rep(0,n)
avdist.HPNX<-rep(0,n)

aver.dist.point.HP2X<-rep(0,n)
aver.dist.point.QU<-rep(0,n)
aver.dist.point.QU2X<-rep(0,n)
aver.dist.point.QUNX<-rep(0,n)
aver.dist.point.HPNX<-rep(0,n)

sumofsq.HP2X<-rep(0,n)
sumofsq.QU<-rep(0,n)
sumofsq.QU2X<-rep(0,n)
sumofsq.QUNX<-rep(0,n)
sumofsq.HPNX<-rep(0,n)

rat1.QU2X<-matrix(0,n,3)
rat2X.QU2X<-matrix(0,n,3)
rat.HP2X<-rep(0,n)
rat.HP<-rep(0,n)


dist1.QU2X<-rep(0,S)
dist2.QU2X<-rep(0,S)
dist3.QU2X<-rep(0,S)
dist4.QU2X<-rep(0,S)
eph.HP2X<-rep(0,S)
rat2.QU2X<-rep(0,S)
rat3.QU2X<-rep(0,S)
rat4.QU2X<-rep(0,S)
rat5.QU2X<-rep(0,S)
rat6.QU2X<-rep(0,S)
rat7.QU2X<-rep(0,S)
rat21.QU2X<-rep(0,n)
rat31.QU2X<-rep(0,n)
rat41.QU2X<-rep(0,n)
rat21X.QU2X<-rep(0,n)
rat31X.QU2X<-rep(0,n)
rat41X.QU2X<-rep(0,n)

mor.HP<-rep(0,S)
mor.HP2X<-rep(0,S)
cott.HP2X<-rep(0,S)
hom.HP2X<-rep(0,S)
inhom.HP2X<-rep(0,S)
shvc.HP2X<-rep(0,S)
cott.QU2X<-rep(0,S)
hom.QU2X<-rep(0,S)
inhom.QU2X<-rep(0,S)
shvc.QU2X<-rep(0,S)
eph.QU2X<-rep(0,S)
cott.QUNX<-rep(0,S)
hom.QUNX<-rep(0,S)
inhom.QUNX<-rep(0,S)
shvc.QUNX<-rep(0,S)
eph.QUNX<-rep(0,S)
dist1sq.QU2X<-rep(0,S)
trueden<-rep(0,S)
inhomden<-rep(0,S)
hom.QU<-rep(0,S)
inhom.QU<-rep(0,S)
cott.QU<-rep(0,S)
shvc.QU<-rep(0,S)
eph.QU<-rep(0,S)
cott.QUNX <-rep(0,S) 
hom.QUNX<-rep(0,S)
inhom.QUNX<-rep(0,S)
shvc.QUNX<-rep(0,S)
eph.QUNX<-rep(0,S)
cott.HPNX<-rep(0,S) 
hom.HPNX<-rep(0,S)
inhom.HPNX<-rep(0,S)
shvc.HPNX<-rep(0,S)
eph.HPNX<-rep(0,S)
ClarkN1<-rep(0,S)
MooreN1<-rep(0,S)
Bythmed<-rep(0,S)
Klinemed4<-rep(0,S)
Klinemed2<-rep(0,S)
RoM.HP<-rep(0,S)
RoM.QU<-rep(0,S)
RoM.HP2X<-rep(0,S)
RoM.HP3X<-rep(0,S)
RoM.QU2X<-rep(0,S)
RoM.HPNX<-rep(0,S)
MoR.HP<-rep(0,S)
MoR.HP2X<-rep(0,S)
MoR.HPNX<-rep(0,S)
MoR.HP3X<-rep(0,S)
MoR.HPNE<-rep(0,S)
RoM.HPNE<-rep(0,S)
MoR.QU<-rep(0,S)
MoR.QU2X<-rep(0,S)
min.dist.HP<-rep(0,n)
max.dist.HP<-rep(0,n)
min.dist.HP2X<-rep(0,n)
max.dist.HP2X<-rep(0,n)
  MoR.QU21<-rep(0,S)
  MoR.QU31<-rep(0,S)
  MoR.QU41<-rep(0,S)
  RoM.QU21<-rep(0,S) 
  RoM.QU31<-rep(0,S)
  RoM.QU41<-rep(0,S)
  MoR.QU2X21<-rep(0,S)
  MoR.QU2X31<-rep(0,S)
 MoR.QU2X41<-rep(0,S)
RoM.QU2X41<-rep(0,S) 
  RoM.QU2X21<-rep(0,S)
  RoM.QU2X31<-rep(0,S)
  RoM.QU2X41<-rep(0,S)
  MoR.dist.HP2X<-matrix(0,S,n)                   
MoR.QUN21<-rep(0,S)
MoR.QUN31<-rep(0,S)
MoR.QUN41<-rep(0,S)
RoM.QUN21<-rep(0,S)
RoM.QUN31<-rep(0,S)
RoM.QUN41<-rep(0,S)
cottHP<- rep(0,S)
homHP<-rep(0,S)
cottQU<-rep(0,S)
homQU<-rep(0,S)
cottTQ<-rep(0,S)
homTQ<-rep(09,S)
Clark.N1<-rep(0,S)
Moore.N1<-rep(0,S)
Byth.med<-rep(0,S)
Kline.med4<-rep(0,S)
GTree.PDE<-matrix(0,S,3)

ErrCottHP<-rep(0,Cy)
ErrCottQU<-rep(0,Cy)
ErrPollHP<-rep(0,Cy)
ErrPollQU<-rep(0,Cy)
ErrMoriHP<-rep(0,Cy)
ErrMoriQU<-rep(0,Cy)

#HP.dia<-cbind(x=std.dat$x,y=std.dat$y,dbh=m)
par(mfrow=c(1,1))


decy<-matrix(0,Cy,7)
msecy<-matrix(0,Cy,7)
# Begin cz loop of aggregation 
#for(cz in 1:Cy){
 
 
  
 
#plot(sim.dat$y~sim.dat$x,pch=20)
# Begin S loop for replications
for(s in 1:S){
  # Definition and Initilation of vectors within S loop
  HP1s.dist<-matrix(0,n,2)
  HP2Xs.dist<-matrix(0,n,2)
  HP3Xs.dist<-matrix(0,n,2)
  QU3Xs.dist<-matrix(0,n,4)
  QUs.dist<-matrix(0,n,4)
 QUNX.dist<-matrix(0,n,4)
  QU2Xs.dist<-matrix(0,n,4)
 QUNXs.dist<-matrix(0,n,4)
  QU1.dist<-matrix(0,n,4)
  QU2X.dist<-matrix(0,n,4)
  HP1.dist<-matrix(0,n,2)
  HP2X.dist<-matrix(0,n,2)
 HP3X.dist<-matrix(0,n,2)
 QU3X.dist<-matrix(0,n,4)
 HPNX.dist<-matrix(0,n,2)
 HPNXs.dist<-matrix(0,n,2)
 TQ.all<-matrix(0,n,2)
 minNE<-matrix(0,n,4)
 minSE<-matrix(0,n,4)
 minSW<-matrix(0,n,4)
 minNW<-matrix(0,n,4)
 minEE<-matrix(0,n,2)
 minWW<-matrix(0,n,2)
 HPNE<-matrix(0,n,2)
 HPNE.dist<-matrix(0,n,2)
 m<-q*pi
 #  g-tree pt matrix (0,n,3) dist to g ; g+1 tree; ave g & g+1
 GT.pts<-matrix(0,n,5)
 
 #cc<-q*a*b/(exp(a)-1)/(exp(b)-1)
 #m<-cc/(a*b)*(exp(a)-1)*(exp(b)-1)*pi
 #range.IPlow<-cc*(exp(a+b)-1)
 xvals<-seq(0,1,.01)
 yvals<-seq(0,1,.01)

 
 
 
 # NElam ... = density parameter in gradient quadrant  
 NElam=function(x,y){minden+qNE*exp(a*(2*x-1))*exp(a*(2*y-1))}
 SElam=function(x,y){minden+qSE*exp(a*(2*x-1))*exp(a*(1-2*y))}
 SWlam=function(x,y){minden+qSW*exp(b*(1-2*x))*exp(b*(1-2*y))}
 NWlam=function(x,y){minden+qNW*exp(b*(1-2*x))*exp(b*(2*y-1))}
 
 
 
 par(mfrow=c(1,1))
 # InNE ...  = MaternII pruned simulated pattern in quadrant
 InNE<-rMaternII(NElam, matNE, win=owin(c(.5,1),c(.5,1)))
 InSE<-rMaternII(SElam, matSE, win=owin(c(.5,1),c(0,.5)))
 InSW<-rMaternII(SWlam, matSW, win=owin(c(0,.5),c(0,.5)))
 InNW<-rMaternII(NWlam, matNW, win=owin(c(0,.5),c(.5,1)))
  #  = simulated gradient pattern combining 4 quadrants
eg.IP<-superimpose(InNE, InSE, InSW, InNW)
 #InH<-eg.IP
 
 # = simulated pattern of Strauss inhibition parameters density, proportion, inhibition distance
 #eg.IP<-rStrauss(q, gamma, sigma)
 
 
 #SSNE<-rSSI(addNE, n=NEgam, win=owin(c(.5,1),c(.5,1)), giveup=1000, x.init=InNE)
 #SSSE<-rSSI(addSE, n=SEgam, win=owin(c(.5,1),c(0,.5)),giveup=1000,  x.init=InSE)
 #SSSW<-rSSI(addSW, n=SWgam, win=owin(c(0,.5),c(0,.5)), giveup=1000, x.init=InSW)
 #SSNW<-rSSI(addNW, n=NWgam, win=owin(c(0,.5),c(.5,1)),giveup=1000,  x.init=InNW)
 #  = simulated SSI gradient pattern combining 4 quadrants
 #eg.IP<-superimpose(SSNE,SSSE,SSSW,SSNW, W=owin(c(0,1), c(0,1))) 
 
 #  = simulated pattern of homogenous Poisson process parameters density
 #eg.IP<-rpoispp(q)
 
 # = simulated pattern of Matern II inhibition parameters density, inhibition distance
 #eg.IP<-rMaternII(q,mat)
  
 # = simulated pattern of rThomas clumping parameters density, sd of distance, number in clump
 #eg.IP<-rThomas(q,sigma,mu)
 
 # = simulated pattern Poisson process with homogenous q or inhomogenous vq 
  #vq<-minden+maxden*runif(1)^expn
 #eg.IP<-rpoispp(q)
 #eg.IP<-rpoispp(vq)
 #eg.IP<-rMaternII(q,mat)
 #eg.IP<-rThomas(vq,mat,mu)
 #eg.IP<-rStrauss(vq, gamma, sigma)
 #eg.IP<-rSSI(mat,q)
 #eg.IP<-rSSI(r=sigma, giveup=1000, x.init=InH)
 #eg.IP<-rSSI(r=sigma, n=ciel , giveup=1000, win=owin(c(0,1),c(0,1)), x.init=InH)
 #eg.IP<-rpoispp(lambda=function(x,y){cc*(a*x)})
 #eg.IP<-rpoispp(lambda=function(x,y){minden+q*exp(a*x)})
 
 #AInNE<-rMaternII(qNE,  matNE, win=owin(c(.5,1),c(.5,1)))
 #AInSE<-rMaternII(qSE,  matSE, win=owin(c(.5,1),c(0,.5)))
 #AInSW<-rMaternII(qSW,  matSW, win=owin(c(0,.5),c(0,.5)))
 #AInNW<-rMaternII(qNW,  matNW, win=owin(c(0,.5),c(.5,1)))
 #eg.IP<-superimpose(AInNE,AInSE,AInSW,AInNW, W=owin(c(0,1), c(0,1))) 
 
 #offset<-samspa*runif(1,0,1)  
 
# Begin n loop calculations of n samples of one pattern
 for(k in 1:n){
   g<-as.integer((k/GX)-.00000001)+1
   h<-k-(g-1)*(GX)
   
   #eg.IP<-rpoispp(q)
   #vq<-minden+(maxden*runif(1))^expn
   vq<-q
   #eg.IP<-rpoispp(vq)
  
   #eg.IP<-rThomas(q,mat,mu)
   #eg.IP<-rStrauss(vq, gamma, mat)
   #eg.IP<-rMaternII(vq,mat)
   # inhib = total density
   inhibden[s]<-length(eg.IP$x)
   std.dat<-eg.IP
   IP.dat<-list(X=std.dat$x,Y=std.dat$y)
   pp.dat<-cbind(std.dat$x,std.dat$y)
   
   sample.dat<-std.dat[std.dat$x<=el&std.dat$x>=wl&std.dat$y>=sl&std.dat$y<=nl]
   
   # truen = density in actual sample
   truen[k]<-length(sample.dat$x)/((el-wl) *(nl-sl))
   
    pp.x<-IP.dat$X[!is.na(IP.dat$X)&!is.na(IP.dat$Y)]
    pp.y<-IP.dat$Y[!is.na(IP.dat$X)&!is.na(IP.dat$Y)]
    
   #  survey points either centered, or random  
   # surveypt<-c(.5,.5)
   surveypt<-c(runif(1,min(pp.x)+wl*(max(pp.x)-min(pp.x)),max(pp.x)-(1-el)*(max(pp.x)-min(pp.x))),
           runif(1,min(pp.y)+sl*(max(pp.y)-min(pp.y)),max(pp.y)-(1-nl)*(max(pp.y)-min(pp.y))))           
   #surveypt<-c(bord+offset+(h-1)*samspa, bord+ offset +(g-1)*samspa)
    
   # calculation of distances from sample points to pattern events 
   #nearest<-which.min(as.matrix(dist(rbind(surveypt,cbind(pp.x,pp.y))))[1,-1])
    nearest<-which.min(sqrt(apply((t(cbind(pp.x,pp.y))-surveypt)^2,2,sum)))
    nearest.dist.QU2X[k]<-sqrt(sum((surveypt-cbind(pp.x,pp.y)[nearest,])^2))
    
     
    nearest<-which.min(sqrt(apply((t(cbind(pp.x,pp.y))-surveypt)^2,2,sum)))
    nearest2<- min(sqrt(apply((t(cbind(pp.x,pp.y))-surveypt)^2,2,sum)))
    
    dist.nearest[k]<-sqrt(sum((surveypt-cbind(pp.x,pp.y)[nearest,])^2))
   dist.order<-sort(sqrt(apply((t(cbind(pp.x,pp.y))-surveypt)^2,2,sum)) )   
   dist.one.all<-sqrt(apply((t(cbind(pp.x, pp.y))-surveypt)^2,2,sum))
 
     # load matrix with g & g+1 at point   
       GT.pts[k,1]<-dist.order[gt]
       GT.pts[k,2]<-dist.order[gt+1]
       GT.pts[k,3]<-sum(GT.pts[k,1],GT.pts[k,2])/2
       GT.pts[k,4]<-GT.pts[k,1]^2
       GT.pts[k,5]<-GT.pts[k,3]^2
  
  # calculation of distances from sample points to events in 4 quadrants
    pp.x.NE<-pp.x[pp.x>surveypt[1] & pp.y>surveypt[2]]
    pp.x.SE<-pp.x[pp.x>surveypt[1] & pp.y<surveypt[2]]
    pp.y.NE<-pp.y[pp.x>surveypt[1] & pp.y>surveypt[2]]
    pp.y.SE<-pp.y[pp.x>surveypt[1] & pp.y<surveypt[2]]
     
      
    pp.x.SW<-pp.x[pp.x<surveypt[1] & pp.y<surveypt[2]]
    pp.x.NW<-pp.x[pp.x<surveypt[1] & pp.y>surveypt[2]]
    pp.y.SW<-pp.y[pp.x<surveypt[1] & pp.y<surveypt[2]]
    pp.y.NW<-pp.y[pp.x<surveypt[1] & pp.y>surveypt[2]]
    
    pp.x.EE<-pp.x[pp.x>surveypt[1]]
    pp.x.WW<-pp.x[pp.x<surveypt[1]]
    pp.y.EE<-pp.y[pp.x>surveypt[1]]
    pp.y.WW<-pp.y[pp.x<surveypt[1]]
      
    
    #nearest.NE<-which.min(sqrt(apply((t(cbind(pp.x.NE, pp.y.NE))-surveypt)^2,2,sum)))
    #QU1.dist[k,1]<-sqrt(sum((surveypt-cbind(pp.x.NE,pp.y.NE)[nearest.NE,])^2))
    dist.NE.all<-sqrt(apply((t(cbind(pp.x.NE, pp.y.NE))-surveypt)^2,2,sum))
    dist.order.NE<-order(dist.NE.all)
    minNE[k,1]<-dist.NE.all[dist.order.NE[1]]
    #QU1.dist[k,1]<- dist.NE.all[dist.order.NE[1]]
     
   # distances from events to ordered events in same quadrant
    QU1.dist[k,1]<-ifelse(minNE[k,1]> distveil, dist.NE.all[dist.order.NE[1]],
                          dist.NE.all[dist.order.NE[2]]) 
    QU2X.dist[k,1]<-ifelse(random1[k]<rancut, dist.NE.all[dist.order.NE[2]],
                           dist.NE.all[dist.order.NE[1]])
    QU3X.dist[k,1]<-ifelse(random1[k]<rancut, dist.NE.all[dist.order.NE[3]],
                           dist.NE.all[dist.order.NE[1]])
    
    
   
   
    #nearest.SE<-which.min(sqrt(apply((t(cbind(pp.x.SE, pp.y.SE))-surveypt)^2,2,sum)))
    #QU1.dist[k,2]<-sqrt(sum((surveypt-cbind(pp.x.SE,pp.y.SE)[nearest.SE,])^2))
    
    dist.SE.all<-sqrt(apply((t(cbind(pp.x.SE, pp.y.SE))-surveypt)^2,2,sum))
    dist.order.SE<-order(dist.SE.all)
    minSE[k,1]<-dist.SE.all[dist.order.SE[1]]
    #QU1.dist[k,2]<- dist.SE.all[dist.order.SE[1]]
    QU1.dist[k,2]<-ifelse(minSE[k,1]> distveil, dist.SE.all[dist.order.SE[1]],
                          dist.SE.all[dist.order.SE[2]]) 
    QU2X.dist[k,2]<-ifelse(random2[k]<rancut, dist.SE.all[dist.order.SE[2]],
                           dist.SE.all[dist.order.SE[1]])
    QU3X.dist[k,2]<-ifelse(random2[k]<rancut , dist.SE.all[dist.order.SE[3]],
                           dist.SE.all[dist.order.SE[1]])
    
    
   
        
    #nearest.SW<-which.min(sqrt(apply((t(cbind(pp.x.SW, pp.y.SW))-surveypt)^2,2,sum)))
    #QU1.dist[k,3]<-sqrt(sum((surveypt-cbind(pp.x.SW,pp.y.SW)[nearest.SW,])^2))
    dist.SW.all<-sqrt(apply((t(cbind(pp.x.SW, pp.y.SW))-surveypt)^2,2,sum))
    dist.order.SW<-order(dist.SW.all)
    minSW[k,1]<-dist.SW.all[dist.order.SW[1]]
    #QU1.dist[k,3]<- dist.SW.all[dist.order.SW[1]]
                          
   QU1.dist[k,3]<-ifelse(minSW[k,1]> distveil, dist.SW.all[dist.order.SW[1]],
                          dist.SW.all[dist.order.SW[2]]) 
    QU2X.dist[k,3]<-ifelse(random3[k]<rancut, dist.SW.all[dist.order.SW[2]],
                           dist.SW.all[dist.order.SW[1]])
    QU3X.dist[k,3]<-ifelse(random3[k]<rancut , dist.SW.all[dist.order.SW[3]],
                           dist.SW.all[dist.order.SW[1]])
    
    
     
    #nearest.NW<-which.min(sqrt(apply((t(cbind(pp.x.NW, pp.y.NW))-surveypt)^2,2,sum)))
    #QU1.dist[k,4]<-sqrt(sum((surveypt-cbind(pp.x.NW,pp.y.NW)[nearest.NW,])^2))
    dist.NW.all<-sqrt(apply((t(cbind(pp.x.NW, pp.y.NW))-surveypt)^2,2,sum))
    dist.order.NW<-order(dist.NW.all)
    minNW[k,1]<-dist.NW.all[dist.order.NW[1]]
     #QU1.dist[k,4]<- dist.NW.all[dist.order.NW[1]]
                        
    QU1.dist[k,4]<-ifelse(minNW[k,1]>distveil, dist.NW.all[dist.order.NW[1]],
                          dist.NW.all[dist.order.NW[2]]) 
    QU2X.dist[k,4]<-ifelse(random4[k]<rancut, dist.NW.all[dist.order.NW[2]],
                           dist.NW.all[dist.order.NW[1]])
    QU3X.dist[k,4]<-ifelse(random4[k]<rancut , dist.NW.all[dist.order.NW[3]],
                           dist.NW.all[dist.order.NW[1]])
    
   

    QUs.ord<-order(QU1.dist[k,])
 
    QUs.dist[k,1]<-QU1.dist[k,QUs.ord[1]]
    QUs.dist[k,2]<-QU1.dist[k,QUs.ord[2]]
    QUs.dist[k,3]<-QU1.dist[k,QUs.ord[3]]
    QUs.dist[k,4]<-QU1.dist[k,QUs.ord[4]]
    QU2Xs.ord<-order(QU2X.dist[k,])
  
    QU2Xs.dist[k,1]<-QU2X.dist[k,QU2Xs.ord[1]]
    QU2Xs.dist[k,2]<-QU2X.dist[k,QU2Xs.ord[2]]
    QU2Xs.dist[k,3]<-QU2X.dist[k,QU2Xs.ord[3]]
    QU2Xs.dist[k,4]<-QU2X.dist[k,QU2Xs.ord[4]]
 
   
 
    QUNX.dist[k,1]<-ifelse(random1[k]<bias,QU2X.dist[k,QUs.ord[1]],QU1.dist[k,QUs.ord[1]])
    QUNX.dist[k,2]<-ifelse(random1[k]<bias,QU2X.dist[k,QUs.ord[2]],QU1.dist[k,QUs.ord[2]])
    QUNX.dist[k,3]<-ifelse(random1[k]<bias,QU2X.dist[k,QUs.ord[3]],QU1.dist[k,QUs.ord[3]])
    QUNX.dist[k,4]<-ifelse(random1[k]<bias,QU2X.dist[k,QUs.ord[4]],QU1.dist[k,QUs.ord[4]])
    QUNX.ord<-order(QUNX.dist[k, ])
 
    QUNXs.dist[k,1]<-QUNX.dist[k,QUNX.ord[1]]
    QUNXs.dist[k,2]<-QUNX.dist[k,QUNX.ord[2]]
    QUNXs.dist[k,3]<-QUNX.dist[k,QUNX.ord[3]]
    QUNXs.dist[k,4]<-QUNX.dist[k,QUNX.ord[4]]
   
   
   
   TQ.all[k,1]<-QUs.dist[k,1]  
    #dia.TQ<-c(dia.NE,dia.SE,dia.SW,dia.NW)
   EQ.all[k,1]<-QUs.dist[k,1]
   #TQ.all[k,3]<-dia.TQ[nearest.TQ]   
   
    #TQ.dia<-c(dia.NE,dia.SE,dia.SW,dia.NW)
    TQ.all[k,2]<-QUs.dist[k,2]
   EQ.all[k,2]<-QUs.dist[k,2]
   EQ.all[k,3]<-QUs.dist[k,3]
   
   #second<-which(TQ.all[k,2]==TQ.dist) 
    #TQ.all[k,4]<-dia.TQ[second]
    #avBA.HP[k]<-(pi/4)*mean(nearest.dia.HP[k]^2,sec.dia.HP[k]^2)
    
    
   # calculation of distances from sample points to points in halves
    #nearest.EE<-which.min(sqrt(apply((t(cbind(pp.x.EE, pp.y.EE))-surveypt)^2,2,sum)))
    #HP.dist[k,1]<-sqrt(sum((surveypt-cbind(pp.x.EE,pp.y.EE)[nearest.EE,])^2))
    dist.EE.all<-sqrt(apply((t(cbind(pp.x.EE, pp.y.EE))-surveypt)^2,2,sum))
    dist.order.EE<-order(dist.EE.all)
    minEE[k,1]<-dist.EE.all[dist.order.EE[1]]
    HP1.dist[k,1]<-ifelse(minEE[k,1]<= distveil && random2[k]<veilcut, dist.EE.all[dist.order.EE[2]],
                         dist.EE.all[dist.order.EE[1]]) 
    HP2X.dist[k,1]<-ifelse(random1[k]<veilcut, dist.EE.all[dist.order.EE[2]],
                          dist.NW.all[dist.order.NW[1]])
    HP3X.dist[k,1]<-ifelse(random1[k]<veilcut , dist.EE.all[dist.order.EE[3]],
                          dist.NW.all[dist.order.NW[1]])
   
      
         
    #nearest.WW<-which.min(sqrt(apply((t(cbind(pp.x.WW, pp.y.WW))-surveypt)^2,2,sum)))
    #HP.dist[k,2]<-sqrt(sum((surveypt-cbind(pp.x.WW,pp.y.WW)[nearest.WW,])^2))
    dist.WW.all<-sqrt(apply((t(cbind(pp.x.WW, pp.y.WW))-surveypt)^2,2,sum))
    dist.order.WW<-order(dist.WW.all)
   
    minWW[k,1]<-dist.WW.all[dist.order.WW[1]]
    HP1.dist[k,2]<-ifelse(minWW[k,1]<= distveil && random1[k]<= veilcut, dist.WW.all[dist.order.WW[2]],
                         dist.WW.all[dist.order.WW[1]]) 
    HP2X.dist[k,2]<-ifelse(random2[k]<veilcut, dist.WW.all[dist.order.WW[2]],
                          dist.WW.all[dist.order.WW[1]])
    HP3X.dist[k,2]<-ifelse(random2[k]<veilcut , dist.WW.all[dist.order.WW[3]],
                          dist.WW.all[dist.order.WW[1]])
   
  
   HP1s.dist[k,1]<-min(HP1.dist[k,])
   HP1s.dist[k,2]<-max(HP1.dist[k,])
   HP2Xs.dist[k,1]<-min(HP2X.dist[k,])
   HP2Xs.dist[k,2]<-max(HP2X.dist[k,]) 
   HP3Xs.dist[k,1]<-min(HP3X.dist[k,])
   HP3Xs.dist[k,2]<-max(HP3X.dist[k,]) 
   
  
   
  #print(HP1.dist)
  #print(QU1.dist)
 
 
  
 # calculation of summary statistice in each sample
    avdist.HP[k]<-mean(HP1.dist[k,])
    avdist.QU[k]<-mean(QU1.dist[k,])
    avdist.HP2X[k]<-mean(HP2X.dist[k,])
    avdist.QU2X[k]<-mean(QU2X.dist[k,])
    avdist.QUNX[k]<-mean(QUNXs.dist[k,])
    avdist.HPNX[k]<-mean(HPNXs.dist[k,])
    
    trip3[k]<-sum(QUs.dist[k,1],QUs.dist[k,2]/2)^2
    
    aver.dist.point.HP[k]<-1/((sum(HP1.dist[k,])/2)^2)
    aver.dist.point.QU[k]<-1/((sum(QU1.dist[k,])/2)^2)
    aver.dist.point.HP2X[k]<-1/((sum(HP2X.dist[k,])/2)^2)
    aver.dist.point.QU2X[k]<-1/((sum(QU2X.dist[k,])/2)^2)
    aver.dist.point.QUNX[k]<-1/((sum(QUNXs.dist[k,])/2)^2)
    aver.dist.point.HPNX[k]<-1/((sum(HPNXs.dist[k,])/2)^2)
    
    sumofsq.HP[k]<-sum(HP1.dist[k,]^2)
    sumofsq.QU[k]<-sum(QU1.dist[k,]^2)
    sumofsq.HP2X[k]<-sum(HP2X.dist[k,]^2)
    sumofsq.QU2X[k]<-sum(QU2X.dist[k,]^2)
    sumofsq.QUNX[k]<-sum(QUNXs.dist[k,]^2)
        
  
# end of n loop for k sample points 
  }

# calculation of summary values for each of s replications  
Trip[s]<-0.5/mean(trip3)
avrbar.HP[s]<-mean(avdist.HP)
avrbar.QU[s]<-mean(avdist.QU)
avrsqbar.HP[s]<-mean(sumofsq.HP)
avrsqbar.QU[s]<-mean(sumofsq.QU)


  
# calculation of estimators for each replication
 
  
  cott.HP[s]<-.5/(mean(HP1.dist)^2)
  hom.HP[s]<- (4*n-2)/(pi*sum(HP1.dist^2))
  inhom.HP[s] <- (2/pi)*mean(1/rowSums(HP1.dist^2))
  shvc.HP[s] <- 0.5*mean(1/((rowSums(HP1.dist)/2)^2))
  
  cott.QU[s] <- 1/(mean(QU1.dist))^2    
  hom.QU[s] <- (16*n-4)/(pi*sum(QU1.dist^2))
  inhom.QU[s] <- 12*mean(1/(pi*rowSums(QU1.dist^2)))
  shvc.QU[s] <- mean(1/(rowSums(QU1.dist)/4)^2)

  f<-0.857266
  cott.TQ[s] <- (0.5*f)/(mean(TQ.all)^2)    
  hom.TQ[s]<- (((4*n-2)*.824)/pi)/(sum(TQ.all^2))
  inhom.TQ[s]<- 0.820*(2/pi)*mean(1/rowSums(TQ.all^2))
  shvc.TQ[s] <- f*0.5*mean(1/((rowSums(TQ.all)/2)^2))
 
  ClarkN1[s]<-1/(4*mean(HP1s.dist[,1])^2)
  MooreN1[s]<-(n-1)/(pi*(sum(HP1s.dist[,1]^2)))
  Bythmed[s]<-0.69/(pi*median(HP1s.dist[,1])^2)
  Klinemed4[s]<-1/((median(QU1.dist))^2) 


  GTree.PDE[s,1]<-((gt-1)/pi)*(mean(1/GT.pts[,4]))
  GTree.PDE[s,2]<-((gt-0.5)/pi)*(mean((1/GT.pts[,4])))
  GTree.PDE[s,3]<-((gt-1)/pi)*(mean((1/GT.pts[,5])))
 
  cott.EQ[s] <- 0.5/(mean(EQ.all)^2)    
  hom.EQ[s]<- ((9*n-3)/pi)/sum(EQ.all^2)
  inhom.EQ[s]<- mean(6/(pi*rowSums(EQ.all^2)))
  shvc.EQ[s]<- 0.5*mean(1/((rowSums(EQ.all)/2))^2)

 
  Poisden[s]<-mean(truen)  

  #if(s%%100==0){print(s)}
  #print(s)




# drop through cz loop and bypass of n loop 

#cottQU[s]<- 1/(mean(cott.QU)^2)    
#homQU[s] <- ((4*S-1)/(S*pi))/mean(hom.QU)

#ErrCottHP[cz]<-abs(cottHP-mean(Poisden))
#ErrCottQU[cz]<-abs(cottQU-mean(Poisden))
#ErrPollHP[cz]<-abs(homHP-mean(Poisden))
#ErrPollQU[cz]<-abs(homQU-mean(Poisden))
#ErrMoriHP[cz]<-abs(mean(inhom.HP)-mean(Poisden))
#ErrMoriQU[cz]<-abs(mean(inhom.QU)-mean(Poisden))



#decy[cz,1]<-mean(Poisden)
#decy[cz,2]<-cottHP
#decy[cz,3]<-homHP
#decy[cz,4]<-mean(inhom.HP)
#decy[cz,5]<-cottQU
#decy[cz,6]<-homQU
#decy[cz,7]<-mean(inhom.QU)

#msecy[cz,1]<-mean(truen)
#msecy[cz,2]<-cottHP-mean(Poisden)
#msecy[cz,3]<-homHP-mean(Poisden)
#msecy[cz,4]<-mean(inhom.HP)-mean(Poisden)
#msecy[cz,5]<-cottQU-mean(Poisden)
#msecy[cz,6]<-homQU-mean(Poisden)
#msecy[cz,7]<-mean(inhom.QU)-mean(Poisden)
#end of regular s replications loop
}

#print(mean(true))
if(cz%%100==0){print(cz)}
#Cy loop third summary
#}



Poi<-mean(true)
err<-(Poi-q)
mean(Poisden)
sd(Poisden)
err


cat(P,r,'virtual plotless sampling', density,'t/ha')
print(" density estimate")
cat(mean(Poisden), mean(cott.HP),mean(hom.HP),mean(inhom.HP),mean(shvc.HP), 
    mean(cott.QU),mean(hom.QU),mean(inhom.QU),mean(shvc.QU), 
    mean(cott.TQ),mean(hom.TQ),mean(inhom.TQ),mean(shvc.TQ), 
    mean(ClarkN1),mean(MooreN1),mean(Bythmed),mean(Klinemed4),
    mean(GTree.PDE[,1]),mean(GTree.PDE[,2]),mean(GTree.PDE[,3]))

print ("relative density")

cat(mean(Poisden)/mean(Poisden), mean(cott.HP)/mean(Poisden),mean(hom.HP)/mean(Poisden),mean(inhom.HP)/mean(Poisden),mean(shvc.HP)/mean(Poisden), 
    mean(cott.QU)/mean(Poisden),mean(hom.QU)/mean(Poisden),mean(inhom.QU)/mean(Poisden),mean(shvc.QU)/mean(Poisden), 
    mean(cott.TQ)/mean(Poisden),mean(hom.TQ)/mean(Poisden),mean(inhom.TQ)/mean(Poisden),mean(shvc.TQ)/mean(Poisden), 
    mean(ClarkN1)/mean(Poisden),mean(MooreN1)/mean(Poisden),mean(Bythmed)/mean(Poisden),mean(Klinemed4)/mean(Poisden),
    mean(GTree.PDE[,1])/mean(Poisden), mean(GTree.PDE[,2])/mean(Poisden),mean(GTree.PDE[,3])/mean(Poisden))

print("std dev of aggr density")
cat(sd(Poisden), sd(cott.HP),sd(hom.HP),sd(inhom.HP),sd(shvc.HP), 
    sd(cott.QU),sd(hom.QU),sd(inhom.QU),sd(shvc.QU), 
    sd(cott.TQ),sd(hom.TQ),sd(inhom.TQ),sd(shvc.TQ), 
    sd(ClarkN1),sd(MooreN1),sd(Bythmed),sd(Klinemed4),
    sd(GTree.PDE[,1]),sd(GTree.PDE[,2]),sd(GTree.PDE[,3]))

print ("relative std dev")
cat(sd(Poisden)/mean(Poisden), sd(cott.HP)/mean(cott.HP),sd(hom.HP)/mean(hom.HP),sd(inhom.HP)/mean(inhom.HP),sd(shvc.HP)/mean(shvc.HP), 
    sd(cott.QU)/mean(cott.QU),sd(hom.QU)/mean(hom.QU),sd(inhom.QU)/mean(inhom.QU),sd(shvc.QU)/mean(shvc.QU), 
    sd(cott.TQ)/mean(cott.TQ),sd(hom.TQ)/mean(hom.TQ),sd(inhom.TQ)/mean(inhom.TQ),sd(shvc.TQ)/mean(shvc.TQ), 
    sd(ClarkN1)/mean(ClarkN1),sd(MooreN1)/mean(MooreN1),sd(Bythmed)/mean(Bythmed),sd(Klinemed4)/mean(Klinemed4),
    sd(GTree.PDE[,1])/mean(GTree.PDE[,1]),sd(GTree.PDE[,2])/mean(GTree.PDE[,2]),sd(GTree.PDE[,3])/mean(GTree.PDE[,3]))


