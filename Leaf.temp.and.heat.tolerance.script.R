library(scales)
library(pse)

#----------------------------------------------------------------------------------------------------
# Observed vs. Predicted leaf temperatures & plotting realtionship
#----------------------------------------------------------------------------------------------------

#Retrieve observed leaf temperature and leaf traits data
ltdat=read.csv("Leaf.and.Temp.Data.csv", header=T)

#Retrieve heat tolerance data
HTs=read.csv("/Users/timothyperez/Dropbox/CH1/Leaf.temp-PHT.manuscript/Final.Data/CH1.HTs.csv", header=T)

#Leaf temperature function & using it to predict leaf temperature
Ltemp=function( W.m2, ws, Ta, RH, WE, Alt, g, abs){
  u=ws #windspeedm/s
  #albedo (Table 11.2, Jones) #sand
  A=10 #altitude above sea-level (m)
  d=(WE/100)*.72 #leaf width (m)
  #d=w* #characteristic dimension of leaf(m)
  cp=29.35 #specific heat of air @constant pressure(J*mol^-1*C^-1) -this should change with air temp, and could be better
  y=6.66*10^-4 #psychrometric constant is constants, but does change with temp, so the value could be better
  sb=(5.67*10^-8) #Stefan-Boltzman constant
  es=(0.611)*exp((17.502*Ta)/(Ta+240.97)) #saturation vapor pressure Kpa (eq. 3.8, Jones p. 41)
  ea=es*(RH/100) #vapor pressure of air (eq. 3.11, Jones p. 41)
  D= es-ea #vapor presure defecit (eq. 3.12, Jones p. 41)
  DELT=(17.502*240.97*es)/((240.97+Ta)^2) #slope of the saturation vapor pressure function (eq 3.9, Jones p. 41)
  pa= 101.3*exp(-A/8200) #atmospheric pressure (eq. 3.7, Jones p. 41)
  s=DELT/pa #slope of saturation mole fraction (eq. 3.10, Jones p. 41)
  gHa=1.4*0.135*sqrt(u/d) #boundary layer conductance for heat
  gr= (4*0.98*(5.67*10^-8)*((273.15+Ta)^3))/cp #radiative conductance can't find an equation for specific heat capacity of air in diferent temps or air pressures
  gva=1.4*0.147*sqrt(u/d) #boundary layer conductance for vapor
  gHr=gHa+gr #sum of boundary layer and radiative conductance (mol*m^-2*s^-1)
  gv=(0.5*g*gva)/(g+gva) #conductance to water vapor (eq. 14.2, Jones)
  ym= y*gHr/gv#modified psychrometric constant
  Rld.sky.lawn=0.98*(5.67*10^-8) * (((273.3+Ta)^4)+((273.3-20)^4)) ##(Jones pg. 27)
  #deg2rad <- function(deg) {(deg * pi) / (180)}
  #view=cos(deg2rad(angle))
  Rsw=W.m2*abs*(1+.26) #(Jones pg. 27)
  Rn=Rsw+Rld.sky.lawn-(2*(0.98*((5.67*10^-8)*((273.3+Ta)^4))))
  #
  TL=Ta+(ym/(s+ym))*((Rn/(gHr*cp))-(D/(pa*ym)))
  return(TL)
}
ltdat$lts=Ltemp(ltdat$W.m2, ltdat$WS.ms, ltdat$ta, ltdat$rh,ltdat$width.filled,  10,  ltdat$Conductance/1000, ltdat$Abs)

#Plot all observed leaf temperatures vs. predicted leaf temperatrues
par(mfcol=c(2,1), mar=c(4,4,2,1), mgp=c(2.5, .5, 0))
plot(ltdat$Tlo.C, ltdat$lts, bty="l", xlim=c(20,50), ylim=c(20,53), 
     xlab=expression("Observed T"[L]~degree~"C"), cex=1, cex.axis=1,
     ylab=expression("Predicted T"[L]~degree~"C"))
abline(0,1, lty=2, lwd=2)
abline(lm(ltdat$lts~ltdat$Tlo.C),  lty=1, lwd=2)
LT.mod=lm(lts~Tlo.C, data=ltdat)
#Prediction data frame for creating 95% CI of linear model
new.LT <- data.frame(LTmin=seq.int(min(na.omit(ltdat$Tlo.C)-10), max(na.omit(ltdat$Tlo.C+10)), length.out=1115))
pred.LT <- predict(LT.mod, newdata=list(Tlo.C=new.LT$LTmin), interval = 'confidence')
polygon(c(rev(new.LT$LTmin), (new.LT$LTmin)), c(((rev(pred.LT[,3]))), (pred.LT[ ,2])), col =alpha("black", 0.25), border = NA)
summary(lm(ltdat$lts~ltdat$Tlo.C))
text(21, 50, "a)")
text(40,23, expression("y= 1.05x + -2.97"), pos=4)
text(40,20, expression("r"^2~"= 0.57; p<0.01"), pos=4)

str(ltdat)
#Calculate maximum leaf temperature (Tmo) and maximum predicted leaf temperature (Predicted TL)
maxltemps=aggregate(ltdat[,c(11,19)], list("Latin"=ltdat$genspe, "leaf"=ltdat$leaf),  FUN=function(x)(quantile(na.omit(x), 0.975)))
meanmaxltemps=aggregate(maxltemps[,c(3,4)], list("Latin"=maxltemps$Latin),  FUN=function(x)(mean(na.omit(x))))

#Plot extreme leaf temperatures
plot(meanmaxltemps$Tlo.C, meanmaxltemps$lts, bty="l",cex=1, cex.axis=1,
     xlim=c(20,50), ylim=c(20,53), xlab=expression("T"[MO]~degree~"C"),
     ylab=expression("Predicted T"[MO]~degree~"C"))
cor.test(meanmaxltemps$lts, meanmaxltemps$Tlo.C)
maxLT.mod=lm(lts~Tlo.C, data=meanmaxltemps)
newmax.LT <- data.frame(LTmin=seq.int(min(na.omit(meanmaxltemps$lts)-300), max(na.omit(meanmaxltemps$Tlo.C+10)), length.out=1115))
predmax.LT <- predict(maxLT.mod, newdata=list(Tlo.C=newmax.LT$LTmin), interval = 'confidence')
polygon(c(rev(newmax.LT$LTmin), (newmax.LT$LTmin)), c(((rev(predmax.LT[,3]))), (predmax.LT[ ,2])), col =alpha("black", 0.25), border = NA)
abline(lm(meanmaxltemps$lts~meanmaxltemps$Tlo.C), lwd=2)
summary(lm(meanmaxltemps$lts~meanmaxltemps$Tlo.C))
abline(0, 1, lty=2)
text(21, 50, "b)")
text(40,23, expression("y= 0.75x + 10.01"), pos=4)
text(40,20, expression("r"^2~"= 0.76; p<0.01"), pos=4)


#----------------------------------------------------------------------------------------------------
# Heat tolerances and extreme leaf temperatures 
#----------------------------------------------------------------------------------------------------
htnlt=merge(maxltemps, HTs, by=c("Latin")) # maximum leaf temperatures (all leaves)
sphtnlt=merge(meanmaxltemps, HTs, by=c("Latin")) #mean maximum leaf temperatures 

par(mfcol=c(2,1), mar=c(4,4,2,1), mgp=c(2.5, .5, 0))
plot(sphtnlt$Tlo.C, sphtnlt$ctmax,  xlim=c(30, 53), ylim=c(35,53), bty="l", col=alpha("#575727", 0.75), pch=19, 
     xlab=expression('T'[MO] ~degree~"C"), ylab=expression("Photosynthetic heat tolerance" ~degree~"C"), cex=1, cex.axis=1)
points(htnlt$Tlo.C, htnlt$ctmax, col=alpha("#575727", 0.5), cex=1)
summary(lm(sphtnlt$ctmax~sphtnlt$Tlo.C))
cor.test(sphtnlt$Tlo.C,sphtnlt$ctmax)
abline(lm(sphtnlt$ctmax~sphtnlt$Tlo.C), lwd=2, col=alpha("#575727", 0.75))
#arrows(sphtnlt$mobs.tmax-sphtnlt$mobs.tmax.sd, sphtnlt$ctmax, sphtnlt$mobs.tmax+sphtnlt$mobs.tmax.sd, sphtnlt$ctmax, length=0.05, angle=90, code=3, col=alpha("#575727", 0.25))
#arrows(sphtnlt$mobs.tmax, sphtnlt$ctmax.lci, sphtnlt$mobs.tmax, sphtnlt$ctmax.uci, length=0.05, angle=90, code=3, col=alpha("#575727", 0.25))

points(sphtnlt$Tlo.C, sphtnlt$tcrit, col=alpha("#C4C41C",0.75), pch=17, cex=1)
points(htnlt$Tlo.C, htnlt$tcrit, col=alpha("#C4C41C", 0.75),pch=2, cex=1)
summary(lm(sphtnlt$tcrit~sphtnlt$Tlo.C))
abline(lm(sphtnlt$tcrit~sphtnlt$Tlo.C), lty=2, lwd=2, col=alpha("#C4C41C", 0.75), cex=1)
cor.test(sphtnlt$Tlo.C,sphtnlt$tcrit)
#arrows(datq$mobs.tmax-datq$mobs.tmax.sd, datq$tcrit, datq$mobs.tmax+datq$mobs.tmax.sd, datq$tcrit, length=0.05, angle=90, code=3, col=alpha("#C4C41C", 0.25))
#arrow(datq$mobs.tmax, datq$tcrit.lci, datq$mobs.tmax, datq$tcrit.uci, length=0.05, angle=90, code=3, col=alpha("#C4C41C", 0.25))
abline(0,1, lwd=2, lty=3,col="gray")
text(31, 53, "a)")
text(47.5, 36, expression(T[crit]~"="), pos=4);points(51, 36, col=alpha("#C4C41C", 0.75), pch=2, cex=1)
text(47.5, 38, expression(T[50]~"="), pos=4);points(51, 38, col=alpha("#575727", 0.75) ,pch=1, cex=1)


#Plot predicted in situ leaf temperatures vs. heat tolerances

#Retrieve microclimate data along with thermoregulatory traits
#Microclimate data was obtined using the NicheMapR R package 
lt.dist.dat=read.csv("NicheMapR.data.csv")
str(lt.dist.dat)

#Estimate in situ leaf temperature (Tmis)
lt.dist.dat$plt.gopmax=Ltemp(lt.dist.dat$SOLR, 0.1, lt.dist.dat$S.TA, lt.dist.dat$S.RH,  lt.dist.dat$width.filled, lt.dist.dat$Alt, lt.dist.dat$gopmax/1000, lt.dist.dat$Abs)

#Remove empty rows
lt.dist.dat=lt.dist.dat[-which(is.na(lt.dist.dat$S.TA)),]
lt.dist.dat2=lt.dist.dat[-which(is.na(lt.dist.dat$plt.gopmax)),]
unique(lt.dist.dat2$name)
str(lt.dist.dat2)

#Calculate the max leaf temperature for each leaf.
distltdat=aggregate(lt.dist.dat$plt.gopmax, list("genspe"=lt.dist.dat$name, "leaf"=lt.dist.dat$leaf), FUN=function(x)(quantile(na.omit(x), 0.975)))
colnames(distltdat)[3]="Tmis" #adjust names of column to  "maximum in situ" leaf temperature (TMIS)
unique(distltdat$genspe)

#Adjust names in the heat tolerance data to merge into 'Tmis' data
HTs[which(HTs$Latin=="Galphimia gracilis"),]$Latin="Galphimia speciosa"

#Merge Tmis data and heat tolerance data
Fdistltdat.ind.leaves=merge(distltdat, HTs, by.x=c("genspe"), by.y=c("Latin"))
#Remove data for spp with few than 15 occurrence records and the Baurhinia
Fdistltdat.ind.leaves=Fdistltdat.ind.leaves[-which((Fdistltdat.ind.leaves$genspe %in%  c("Eugenia coronata", "Gardenia taitensis", "Bauhinia divaricata"))==T),]

#Calculate the max leaf temperature per species.
mdistlts=aggregate(Fdistltdat.ind.leaves$Tmis, list("Latin"=Fdistltdat.ind.leaves$genspe), FUN=function(x)(mean(na.omit(x))))
colnames(mdistlts)[2]="Tmis" #adjust names of column to  "maximum in situ" leaf temperature (TMIS)
mdistlt=merge(mdistlts, HTs, all=T)

#
plot(jitter(mdistlt$Tmis), mdistlt$ctmax,  xlim=c(30, 52), ylim=c(35,53), bty="l", col=alpha("#575727", 0.75), pch=19, 
     xlab=expression('T'[MIS] ~degree~"C"), ylab=expression("Photoynthetic heat tolerance" ~degree~"C"))
points(Fdistltdat.ind.leaves$Tmis, Fdistltdat.ind.leaves$ctmax, col=alpha("#575727", 0.5))
summary(lm(mdistlt$ctmax~mdistlt$Tmis))
cor.test(mdistlt$Tmis, mdistlt$ctmax)
abline(lm(mdistlt$ctmax~mdistlt$Tmis), lwd=2, col=alpha("#575727", 0.75))

points(mdistlt$Tmis, mdistlt$tcrit,  col=alpha("#C4C41C",0.75), pch=17)
points(jitter(Fdistltdat.ind.leaves$Tmis), Fdistltdat.ind.leaves$tcrit, pch=2, col=alpha("#C4C41C", 0.75))
summary(lm( mdistlt$tcrit~mdistlt$Tmis))
abline(lm(mdistlt$tcrit~mdistlt$Tmis), lty=2, lwd=2, col=alpha("#C4C41C", 0.75))
cor.test(mdistlt$tcrit, mdistlt$Tmis)
abline(0,1, lwd=2, lty=3,col="gray")
text(31, 53, "b)")
text(47.5, 36, expression(T[crit]~"="), pos=4);points(51, 36, col=alpha("#C4C41C", 0.75), pch=2, cex=1)
text(47.5, 38, expression(T[50]~"="), pos=4);points(51, 38, col=alpha("#575727", 0.75) ,pch=1, cex=1)


#----------------------------------------------------------------------------------------------------
# Sensitivity analysis for Leaf temp model
#----------------------------------------------------------------------------------------------------
library(pse)

#Remove Altitude term and set to 10m for elevation in Miami & Fairchild Tropical BotanicGarden
Ltemp2=function( W.m2, ws, Ta, RH, WE, g, abs){
  u=ws #windspeedm/s
  A=10 #altitude above sea-level (m)
  d=(WE/100)*.72 #leaf width (m)
  #d=w* #characteristic dimension of leaf(m)
  cp=29.35 #specific heat of air @constant pressure(J*mol^-1*C^-1) -this should change with air temp, and could be better
  y=6.66*10^-4 #psychrometric constant is constants, but does change with temp, so the value could be better
  sb=(5.67*10^-8) #Stefan-Boltzman constant
  es=(0.611)*exp((17.502*Ta)/(Ta+240.97)) #saturation vapor pressure Kpa (eq. 3.8, Jones p. 41)
  ea=es*(RH/100) #vapor pressure of air (eq. 3.11, Jones p. 41)
  D= es-ea #vapor presure defecit (eq. 3.12, Jones p. 41)
  DELT=(17.502*240.97*es)/((240.97+Ta)^2) #slope of the saturation vapor pressure function (eq 3.9, Jones p. 41)
  pa= 101.3*exp(-A/8200) #atmospheric pressure (eq. 3.7, Jones p. 41)
  s=DELT/pa #slope of saturation mole fraction (eq. 3.10, Jones p. 41)
  gHa=1.4*0.135*sqrt(u/d) #boundary layer conductance for heat
  gr= (4*0.98*(5.67*10^-8)*((273.15+Ta)^3))/cp #radiative conductance can't find an equation for specific heat capacity of air in diferent temps or air pressures
  gva=1.4*0.147*sqrt(u/d) #boundary layer conductance for vapor
  gHr=gHa+gr #sum of boundary layer and radiative conductance (mol*m^-2*s^-1)
  gv=(0.5*g*gva)/(g+gva) #conductance to water vapor (eq. 14.2, Jones)
  ym= y*gHr/gv#modified psychrometric constant
  Rld.sky.lawn=0.98*(5.67*10^-8) * (((273.3+Ta)^4)+((273.3-20)^4)) ##(Jones pg. 27)
  #deg2rad <- function(deg) {(deg * pi) / (180)}
  #view=cos(deg2rad(angle))
  Rsw=W.m2*abs*(1+.26) #(Jones pg. 27)
  Rn=Rsw+Rld.sky.lawn-(2*(0.98*((5.67*10^-8)*((273.3+Ta)^4))))
  #
  TL=Ta+(ym/(s+ym))*((Rn/(gHr*cp))-(D/(pa*ym)))
  return(TL)
}

#Create a data frame for the sensitivity analysis based on oberved environmental variable
sens.dat=data.frame(
  Wm2=runif(1115, min(na.omit(ltdat$W.m2)), max(na.omit(ltdat$W.m2))),#1
  Wms=runif(1115, min(na.omit(ltdat$WS.ms)), max(na.omit(ltdat$WS.ms))),#2
  Ta=runif(1115, min(na.omit(ltdat$ta)), max(na.omit(ltdat$ta))),#3
  RH=runif(1115, min(na.omit(ltdat$rh)), max(na.omit(ltdat$rh))), #4
  We=runif(1115, min(na.omit(ltdat$width.filled)), max(na.omit(ltdat$width.filled))),#5
  gs=runif(1115, min(na.omit(ltdat$gopmax/1000)), max(na.omit(ltdat$gopmax/1000))),#6
  abs=runif(1115, min(na.omit(ltdat$Abs)), max(na.omit(ltdat$Abs))))#7
modelRun=function (sens.dat){return(mapply(Ltemp2, sens.dat[,1],sens.dat[,2],sens.dat[,3],sens.dat[,4], sens.dat[,5], sens.dat[,6], sens.dat[,7]))}
factors=c("Wm2", "WS", "Ta", "RH",  "We", "gs", "abs" )
q=c("qunif","qunif", "qunif", "qunif", "qunif", "qunif", "qunif")
q.arg=list(
  list(min(na.omit(ltdat$W.m2)), max(na.omit(ltdat$W.m2))), 
  list(min(na.omit(ltdat$WS.ms)), max(na.omit(ltdat$WS.ms))), 
  list(min(na.omit(ltdat$ta)), max(na.omit(ltdat$ta))),
  list(min(na.omit(ltdat$rh)), max(na.omit(ltdat$rh))), 
  list(min(na.omit(ltdat$width.filled)), max(na.omit(ltdat$width.filled))), 
  list(min(na.omit(ltdat$gopmax/1000)), max(na.omit(ltdat$gopmax/1000))), 
  list(min(na.omit(ltdat$Abs)), max(na.omit(ltdat$Abs))))

#Run latin hyperube sampling procedure: Note this can take a while, so there is a
#practice run
set.seed(238)

#Practice
myLHS=LHS(modelRun, factors, 12, q=q, q.arg=q.arg, nboot=100, opts="COR")
#More iterations 
myLHS=LHS(modelRun, factors, 1000, q=q, q.arg=q.arg, nboot=100000, opts="COR")

#Extract latin hyper cube results
get.data(myLHSpractice)
myLHS$prcc[[1]]
envrvals=c(myLHS$prcc[[1]]$PRCC[,1][1],
           myLHS$prcc[[1]]$PRCC[,1][2],
           myLHS$prcc[[1]]$PRCC[,1][3],
           myLHS$prcc[[1]]$PRCC[,1][4])
leafrvals=c(myLHS$prcc[[1]]$PRCC[,1][5], 
            myLHS$prcc[[1]]$PRCC[,1][6],
            myLHS$prcc[[1]]$PRCC[,1][7])

envrvals.minci=c(myLHS$prcc[[1]]$PRCC[,4][1],
                 myLHS$prcc[[1]]$PRCC[,4][2],
                 myLHS$prcc[[1]]$PRCC[,4][3],
                 myLHS$prcc[[1]]$PRCC[,4][4])
leafrvals.minci=c(myLHS$prcc[[1]]$PRCC[,4][5], 
                  myLHS$prcc[[1]]$PRCC[,4][6],
                  myLHS$prcc[[1]]$PRCC[,4][7])

envrvals.maxci=c(myLHS$prcc[[1]]$PRCC[,5][1],
                 myLHS$prcc[[1]]$PRCC[,5][2],
                 myLHS$prcc[[1]]$PRCC[,5][3],
                 myLHS$prcc[[1]]$PRCC[,5][4])
leafrvals.maxci=c(myLHS$prcc[[1]]$PRCC[,5][5], 
                  myLHS$prcc[[1]]$PRCC[,5][6],
                  myLHS$prcc[[1]]$PRCC[,5][7])

#Create data frame of sensitivity results
rvals=data.frame("factors"=factors, rvals=c(envrvals,leafrvals), 
                 minci=c(envrvals.minci,leafrvals.minci), 
                 maxci=c(envrvals.maxci,leafrvals.maxci))

#Plot Results of sensitivity LHS
par(mfcol=c(1,1), mai=c(1,1,1,1),mar=c(3,0,0,1),oma=c(1,0,0,6))
barplot(rvals$rvals, xlim=c(-1.25,1.5), axes=F, ylab="",
        col=c("burlywood4", "burlywood4","burlywood4","burlywood4","chartreuse4", "chartreuse4", "chartreuse4"), horiz=T, las=1,  xlab=NULL)
axis(4,outer=F, line=-2, at=seq(.625, 8.25, by = 1.25),cex=1.25, c(expression("Solar irradiance"), expression('Windspeed'), expression("Air temperature"), expression("Relative humidity"), expression("Effective leaf width"), expression("Stomatal conductance"), expression("Absorptivity")), las=1, lty=0, mgp=c(2,.5,1))
mtext("Effect on Leaf Temperture", side=1, at=0, line=3, cex=1)
axis(1, at=seq(-1, 1, by=0.2), cex.lab=1)
abline(v=0, col="black", lwd=1.3, lty=2)
