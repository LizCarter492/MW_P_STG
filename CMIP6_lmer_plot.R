library(corrplot)
library(multipanelfigure)
library(biwavelet)
library(lubridate)

library(raster)
library(ncdf4)
library(biwavelet)
library(lubridate)
setwd("H:/Research/MW_ClimateChange/data/climate_data")
ta_ras<-brick("air.mon.mean.nc", varname = "air")
ta_nc<-nc_open("air.mon.mean.nc")
time_nc<-ncvar_get(ta_nc, varid="time")
time<-as.POSIXct("1800-01-01 00:00")+as.difftime(time_nc,units="hours")
time<-as.Date(time, origin=as.Date("1800-01-01  00:00:0.0"))
nc_close(ta_nc)

PML_STGA_N.ext <- extent(c(360-135, 360-55, 70, 89))
PML_STGA_S.ext <- extent(c(360-135, 360-55, 30, 50))
WA_LSTCA_W.ext <- extent(c( 360-90, 360-80, 30, 40))
WA_LSTCA_E.ext <- extent(c( 360-75, 360-65, 30, 40))
EA_LSTCA_W.ext <- extent(c(360-15, 360-0, 7.5, 32))
EA_LSTCA_E.ext <- extent(c(360-30, 360-15, 7.5, 32))
EP_LSTCA_E.ext<- extent(360-110,360-105,25,32)
EP_LSTCA_W.ext<- extent(360-120,360-115,25,32)


PML_STGA_N<-as.numeric(extract(ta_ras, PML_STGA_N.ext, fun=mean, na.rm=T))
PML_STGA_S<-as.numeric(extract(ta_ras, PML_STGA_S.ext, fun=mean, na.rm=T))

WA_LSTCA_W<-as.numeric(extract(ta_ras, WA_LSTCA_W.ext, fun=mean, na.rm=T))
WA_LSTCA_E<-as.numeric(extract(ta_ras, WA_LSTCA_E.ext, fun=mean, na.rm=T))

EA_LSTCA_W<-as.numeric(extract(ta_ras, EA_LSTCA_W.ext, fun=mean, na.rm=T))
EA_LSTCA_E<-as.numeric(extract(ta_ras, EA_LSTCA_E.ext, fun=mean, na.rm=T))

EP_LSTCA_W<-as.numeric(extract(ta_ras, EP_LSTCA_W.ext, fun=mean, na.rm=T))
EP_LSTCA_E<-as.numeric(extract(ta_ras, EP_LSTCA_E.ext, fun=mean, na.rm=T))

PML_STG=PML_STGA_S-PML_STGA_N
WA_LSTC=WA_LSTCA_W-WA_LSTCA_E
EA_LSTC=EA_LSTCA_E-EA_LSTCA_W
EP_LSTC=EP_LSTCA_W-EP_LSTCA_E


time<-as.character(time)
CI<-data.frame(time,
               PML_STGA_N, 
               PML_STGA_S,
               PML_STG,
               
               WA_LSTCA_W,
               WA_LSTCA_E,
               WA_LSTC,
               
               EA_LSTCA_W,
               EA_LSTCA_E,
               EA_LSTC,
               
               EP_LSTCA_W,
               EP_LSTCA_E,
               EP_LSTC
               
)



P.nc<-nc_open("D:/DELL_PHD/E/GPCC precip/precip.mon.combined.total.v7.nc")
date<-ncvar_get(P.nc, "time")
date<-as.Date(date, origin="1800-01-01")
P<-ncvar_get(P.nc, "precip")
nc_close(P.nc)

gp_x<-rbind(c(360-105,360-90), c(38,48))
gp_x<-extent(gp_x)

DATE<-seq.Date(from=as.Date("1948-01-01"), to=as.Date("2017-12-31"), by="months")

CI<-CI[which(DATE%in%date),]
P<-brick("D:/DELL_PHD/E/GPCC precip/precip.mon.combined.total.v7.nc", varname="precip")
CI$P<-extract(P, gp_x,fun=mean)[date %in% DATE]

plot_sub<-CI
plot_sub[,-1]<-apply(plot_sub[,-1],2,function(x)(x-mean(x))/sd(x))
plot_sub<-plot_sub[625:744,]
par(mar=c(7.1, 4.1, 4.1, 2.1))
plot(-1*plot_sub$EA_LSTC, type="l", lwd=2, 
     col="gray",xlab="", xaxt="n", ylim=c(-2.5,3.5),
     ylab="Standardized variable")
points(plot_sub$WA_LSTC, type="l", lwd=2, col="purple")
points(-1*plot_sub$PML_STG, type="l", lwd=2, col="green")
points(plot_sub$P, type="l", lwd=3, col="blue")
axis(1, at=1:120,labels=plot_sub$time, las=2)
legend("topleft", 
       legend = c("MW Precip", "EA_LSTC", "WA_LSTC", "PML_STG"), 
       col = c("blue","gray", "purple", "green"), 
       lty = 1, 
       text.col = "black", 
       lwd=2,
       horiz = T , 
       inset = c(0.01, 0.01))


#read in lake Ontario lake levels
ll<-read.csv("H:/Research/MW_ClimateChange/data/OnatarioLakeLevels/miHuron1918.csv")
CI$ll<-rep(NA)
for(i in unique(year(CI$time))){
  for(j in unique(month(CI$time))){
    CI$ll[which(year(CI$time)==i & month(CI$time)==j)]<-ll[which(ll$year==i),(j+1)]
  }
}

CI$PML_EA<-CI$PML_STG*CI$EA_LSTC
CI$PML_WA<-CI$PML_STG*CI$WA_LSTC
CI$EA_WA<-CI$WA_LSTC*CI$EA_LSTC
m1<-lm(ll~EA_LSTC + WA_LSTC + PML_STG + PML_EA+PML_WA+EA_WA, data=CI)




###########################################################################################
####wavelet coherence######################################################################
###########################################################################################
t2 = cbind(CI$time, CI$EA_LSTC)
t1 = cbind(CI$time, CI$P)
# Specify the number of iterations. The more, the better (>1000).  For the
# purpose of this tutorial, we just set it = 10
nrands = 100

wtc.AB = wtc(t1, t2, nrands = nrands)

par(mfrow=c(1,3), oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12, ylab = "Scale", xlab = "Period", 
     plot.cb = TRUE, main = "P v EA_LSTC")


t2 = cbind(CI$time, CI$WA_LSTC)
# Specify the number of iterations. The more, the better (>1000).  For the
# purpose of this tutorial, we just set it = 10
nrands = 100

wtc.AB = wtc(t1, t2, nrands = nrands)

plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12,  xlab = "", 
     plot.cb = TRUE, main = "P v WA_LSTC")



t2 = cbind(CI$time, CI$PML_STG)
# Specify the number of iterations. The more, the better (>1000).  For the
# purpose of this tutorial, we just set it = 10
nrands = 100

wtc.AB = wtc(t1, t2, nrands = nrands)

plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12,  xlab = "", 
     plot.cb = TRUE, main = "P v PML_STG")

CI$time<-as.Date(CI$time, format=("%Y-%m-%d"))
main_z<-CI
for(i in 1:12){
  main_z[month(main_z$time) == i, -1]<-apply(main_z[month(main_z$time) == i, -1], 2, function(x) (x-mean(x))/sd(x))
}

t2 = cbind(main_z$time, main_z$EA_LSTC)
t1 = cbind(main_z$time, main_z$P)
# Specify the number of iterations. The more, the better (>1000).  For the
# purpose of this tutorial, we just set it = 10
nrands = 100

wtc.AB = wtc(t1, t2, nrands = nrands)

par(mfrow=c(1,3), oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12, ylab = "Scale", xlab = "Period", 
     plot.cb = TRUE, main = "P v EA_LSTC")


t2 = cbind(main_z$time, main_z$WA_LSTC)
# Specify the number of iterations. The more, the better (>1000).  For the
# purpose of this tutorial, we just set it = 10
nrands = 100

wtc.AB = wtc(t1, t2, nrands = nrands)

plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12,  xlab = "Period", 
     plot.cb = TRUE, main = "P v WA_LSTC")



t2 = cbind(main_z$time, main_z$PML_STG)
# Specify the number of iterations. The more, the better (>1000).  For the
# purpose of this tutorial, we just set it = 10
nrands = 100

wtc.AB = wtc(t1, t2, nrands = nrands)

plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12,  xlab = "Period", 
     plot.cb = TRUE, main = "P v PML_STG")


#####################################################################
##future climate
setwd("H:/Research/MW_ClimateChange/data/CMIP6_ssp370_ens68_tas_mon")
CI_p<-array(NA, dim=c(length(2030:3000), 4, length(list.files())))
tas<-nc_open(list.files()[3])
time<-ncvar_get(tas, varid="time")
time<-as.Date(time, origin=as.Date("1850-01-01", format="%Y-%m-%d"))
nc_close(tas)
remove(tas)
for(i in 3:length(list.files())){
  tas<-nc_open(list.files()[i])
  ta<-ncvar_get(tas, varid="tas")
  lat<-ncvar_get(tas, varid="lat")
  lon<-ncvar_get(tas, varid="lon")
  nc_close(tas)
  ta<-ta[,,2030:3000]
  ta_ras<-t(brick(ta[,seq(144,1),], ymn=min(lon), ymx=max(lon),xmn=min(lat), xmx=max(lat)))
  PML_STGA_N<-as.numeric(extract(ta_ras, PML_STGA_N.ext, fun=mean, na.rm=T))
  PML_STGA_S<-as.numeric(extract(ta_ras, PML_STGA_S.ext, fun=mean, na.rm=T))
  
  WA_LSTCA_W<-as.numeric(extract(ta_ras, WA_LSTCA_W.ext, fun=mean, na.rm=T))
  WA_LSTCA_E<-as.numeric(extract(ta_ras, WA_LSTCA_E.ext, fun=mean, na.rm=T))
  
  EA_LSTCA_W<-as.numeric(extract(ta_ras, EA_LSTCA_W.ext, fun=mean, na.rm=T))
  EA_LSTCA_E<-as.numeric(extract(ta_ras, EA_LSTCA_E.ext, fun=mean, na.rm=T))
  
  EP_LSTCA_W<-as.numeric(extract(ta_ras, EP_LSTCA_W.ext, fun=mean, na.rm=T))
  EP_LSTCA_E<-as.numeric(extract(ta_ras, EP_LSTCA_E.ext, fun=mean, na.rm=T))
  #remove(ta, ta_ras)
  CI_p[,1,i]=PML_STGA_S-PML_STGA_N
  CI_p[,2,i]=WA_LSTCA_W-WA_LSTCA_E
  CI_p[,3,i]=EA_LSTCA_E-EA_LSTCA_W
  CI_p[,4,i]=EP_LSTCA_W-EP_LSTCA_E
  remove(ta, ta_ras, PML_STGA_N, PML_STGA_S, EA_LSTCA_E,EA_LSTCA_W, WA_LSTCA_E,WA_LSTCA_W)
}
save(CI_p, file="PML_STG_WA_LSTC_EA_LSTC_EP_LSTC.RData")

#the model
library(lme4)
library(lmerTest)
CI$month<-as.factor(month(CI$time, label=T, abbr=T))
rs<-lmer(P~PML_STG + WA_LSTC + EA_LSTC+
           PML_STG:WA_LSTC+PML_STG:EA_LSTC+
           EA_LSTC:WA_LSTC+
           (1|month),
         data=CI)


pred<-array(NA, dim(CI_p)[c(1,3)])
for(i in 1:dim(CI_p)[3]){
  newx<-as.data.frame(CI_p[,,i])
  names(newx)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
  newx$month<-as.factor(month(time[2030:3000], label=T, abbr=T))
  pred[,i]<-predict(rs, newdata=newx)
}

CI_p_mean<-apply(CI_p[,,], c(1,2), mean, na.rm=T)
newx<-as.data.frame(CI_p_mean)
names(newx)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
newx$month<-as.factor(month(time[2030:3000], label=T, abbr=T))

#climatological mean end of 20th century
CI<-CI[,-16]
z20<-CI[1:12,]
for(i in 1:12){
  z20[i,-1]<-apply(CI[which(month(CI$time)==i & 
                              year(CI$time)>=1970 & 
                              year(CI$time)<2000),-1],2,mean)
}
z20<-as.data.frame(z20)
names(z20)<-names(CI)


#climatological mean of end of 21st century
z21<-CI_p_mean[1:12,]
for(i in 1:12){
  z21[i,]<-apply(CI_p_mean[which(month(time[2030:3000])==i & year(time[2030:3000])>=2070 & year(time[2030:3000])<2100),],2,mean)
} 



P<-predict(rs, newdata=newx)
z21<-as.data.frame(z21)
names(z21)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
z21$P<-rep(NA)
for(i in 1:12){
  z21$P[i]<-mean(P[which(month(time[2030:3000])==i &
                           year(time[2030:3000])>=2070 & 
                           year(time[2030:3000])<2100)])
} 

par(mfrow=c(4,1), mar=c(3, 4, 0.5, 0.5))

plot(z21$EA_LSTC, type="l", col="red", lwd=3, xaxt="n", xlab="",
     ylab="EA_LSTC (oC)", ylim=c(-10,6))
points(z20$EA_LSTC, type="l", lwd=3)
for(i in 2:dim(CI_p)[3]){
  df<-as.data.frame(CI_p[,,i])
  names(df)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
  df_m<-df[1:!2,]
  for(j in 1:12){
    df_m[j,]<-apply(df[which(month(time[2030:3000])==j & year(time[2030:3000])>=2070 & year(time[2030:3000])<2100),],2,mean,na.rm=T)
  }
  points(df_m$EA_LSTC, type="l", lwd=1, col="pink")
  
  
}
points(z20$EA_LSTC, type="l", lwd=3)
points(z21$EA_LSTC, type="l", lwd=3, col="red")
axis(1, at=1:12,labels=month(1:12,label=T, abbr=T),  las=2)
legend("top", 
       legend = c("20th observed", "21st CMIP6 mean", "21st members"), 
       col = c("black","red", "pink"), 
       lty = 1, 
       lwd=c(3,3,1),
       text.col = "black", 
       horiz = T , 
       inset = c(0.01, 0.01))



plot(z21$WA_LSTC, type="l", col="red", lwd=3, xaxt="n",ylim=c(-15,15), xlab="",ylab="WA_LSTC (oC)")
points(z20$WA_LSTC, type="l", lwd=3)
for(i in 1:dim(CI_p)[3]){
  df<-as.data.frame(CI_p[,,i])
  names(df)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
  df_m<-df[1:12,]
  for(j in 1:12){
    df_m[j,]<-apply(df[which(month(time[2030:3000])==j & year(time[2030:3000])>=2070 & year(time[2030:3000])<2100),],2,mean,na.rm=T)
  }
  points(df_m$WA_LSTC, type="l", lwd=1, col="pink")
  
  
}
points(z20$WA_LSTC, type="l", lwd=3)
points(z21$WA_LSTC, type="l", lwd=3, col="red")
axis(1, at=1:12,labels=month(1:12,label=T, abbr=T),  las=2)
legend("top", 
       legend = c("20th observed", "21st CMIP6 mean", "21st members"), 
       col = c("black","red", "pink"), 
       lty = 1, 
       lwd=c(3,3,1),
       text.col = "black",
       horiz = T , 
       inset = c(0.01, 0.01))



plot(z20$PML_STG, type="l", lwd=3, xaxt="n", xlab="",
     ylab="PML_STG (oC)", ylim=c(15,45))
points(z21$PML_STG, type="l", col="red", lwd=3)
for(i in 2:dim(CI_p)[3]){
  df<-as.data.frame(CI_p[,,i])
  names(df)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
  df_m<-df[1:!2,]
  for(j in 1:12){
    df_m[j,]<-apply(df[which(month(time[2030:3000])==j & year(time[2030:3000])>=2070 & year(time[2030:3000])<2100),],2,mean,na.rm=T)
  }
  points(df_m$PML_STG, type="l", lwd=1, col="pink")
  
  
}
points(z20$PML_STG, type="l", lwd=3)
points(z21$PML_STG, type="l", lwd=3, col="red")
axis(1, at=1:12,labels=month(1:12,label=T, abbr=T),  las=2)
legend("top", 
       legend = c("20th observed", "21st CMIP6 mean", "21st members"), 
       col = c("black","red", "pink"), 
       lty = 1, 
       lwd=c(3,3,1),
       text.col = "black",
       horiz = T , 
       inset = c(0.01, 0.01))



plot(z21$P, type="l", col="blue", lwd=3, 
     xaxt="n", xlab="",ylab="P (mm)", ylim=c(0,160))
points(z20$P, type="l", lwd=3)
for(i in 1:dim(CI_p)[3]){
  df<-as.data.frame(CI_p[,,i])
  names(df)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
  df$month<-as.factor(month(time[2030:3000], label=T, abbr=T))
  P<-predict(rs,df)
  P_m<-rep(NA,12)
  for(j in 1:12){
    P_m[j]<-mean(P[which(month(time[2030:3000])==j & year(time[2030:3000])>=2070 & year(time[2030:3000])<2100)])
  } 
  points(P_m, type="l", lwd=1, col="lightblue")
}
points(z21$P, type="l", col="blue", lwd=3)
points(z20$P, type="l", lwd=3)
axis(1, at=1:12,labels=month(1:12,label=T, abbr=T),  las=2)
legend("top", 
       legend = c("20th observed", "21st CMIP6 mean", "21st members"), 
       col = c("black","blue", "blue"), 
       lty = 1, 
       lwd=c(3,3,1),
       text.col = "black",
       horiz = T , 
       inset = c(0.01, 0.01))


r<-cor(main_z[,c("P", "EA_LSTC", "WA_LSTC", "PML_STG")])
corrplot(cor)


####################################################################################################
#future projections################################################################################
m1<-lmer(P~EA_LSTC+WA_LSTC+PML_STG+
           PML_STG:EA_LSTC + PML_STG:WA_LSTC+
           EA_LSTC:WA_LSTC+
           (1|month), data=dat)

pred<-array(NA, dim(CI_p)[c(1,3)])
for(i in 1:dim(CI_p)[3]){
  newx<-as.data.frame(CI_p[,,i])
  names(newx)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
  newx$month<-as.factor(month(time[2030:3000]))
  pred[,i]<-predict(m1, newdata=newx)
}

CI_p_mean<-apply(CI_p[,,], c(1,2), mean, na.rm=T)
newx<-as.data.frame(CI_p_mean)
names(newx)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
newx$month<-as.factor(month(time[2030:3000]))


z20<-CI[1:12,]
for(i in 1:12){
  z20[i,-1]<-apply(CI[which(month(CI$time)==i & 
                              year(CI$time)>=1970 & 
                              year(CI$time)<2000),-1],2,mean)
}
z20<-as.data.frame(z20)
names(z20)<-names(CI)

z21<-CI_p_mean[1:12,]
for(i in 1:12){
  z21[i,]<-apply(CI_p_mean[which(month(time[2030:3000])==i & year(time[2030:3000])>=2070 & year(time[2030:3000])<2100),],2,mean)
} 

P<-predict(elm, newdata=newx)
z21<-as.data.frame(z21)
names(z21)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
z21$P<-rep(NA)
for(i in 1:12){
  z21$P[i]<-mean(P[which(month(time[2030:3000])==i &
                           year(time[2030:3000])>=2070 & 
                           year(time[2030:3000])<2100)])
} 

plot(z21$P, type="l", col="red", lwd=3, xaxt="n", xlab="",ylab="P (mm)", ylim=c(10,140))
points(z20$P, type="l", lwd=3)
for(i in 1:dim(CI_p)[3]){
  df<-as.data.frame(CI_p[,,i])
  names(df)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
  df$month<-as.factor(month(time[2030:3000]))
  P<-predict(m3,df)
  for(j in 1:12){
    P_m[j]<-mean(P[which(month(time[2030:3000])==j & year(time[2030:3000])>=2070 & year(time[2030:3000])<2100)])
  } 
  points(P_m, type="l", lwd=1, col="pink")
}
points(z21$P, type="l", col="red", lwd=3)
points(z20$P, type="l", lwd=3)
axis(1, at=1:12,labels=month(1:12,label=T, abbr=T),  las=2)





plot(z21$EA_LSTC, type="l", col="red", lwd=3, xaxt="n", xlab="",ylab="EA_LSTC (oC)")
points(z20$EA_LSTC, type="l", lwd=3)
for(i in 2:dim(CI_p)[3]){
  df<-as.data.frame(CI_p[,,i])
  names(df)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
  df_m<-df[1:!2,]
  for(j in 1:12){
    df_m[j,]<-apply(df[which(month(time[2030:3000])==j & year(time[2030:3000])>=2070 & year(time[2030:3000])<2100),],2,mean,na.rm=T)
  }
  points(df_m$EA_LSTC, type="l", lwd=1, col="pink")
  
  
}
points(z20$EA_LSTC, type="l", lwd=3)
points(z21$EA_LSTC, type="l", lwd=3, col="red")

axis(1, at=1:12,labels=month(1:12,label=T, abbr=T),  las=2)

plot(z21$WA_LSTC, type="l", col="red", lwd=3, xaxt="n", xlab="",ylab="WA_LSTC (oC)")
points(z20$WA_LSTC, type="l", lwd=3)
for(i in 2:dim(CI_p)[3]){
  df<-as.data.frame(CI_p[,,i])
  names(df)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
  df_m<-df[1:!2,]
  for(j in 1:12){
    df_m[j,]<-apply(df[which(month(time[2030:3000])==j & year(time[2030:3000])>=2070 & year(time[2030:3000])<2100),],2,mean,na.rm=T)
  }
  points(df_m$WA_LSTC, type="l", lwd=1, col="pink")
  
  
}
points(z20$WA_LSTC, type="l", lwd=3)
points(z21$WA_LSTC, type="l", lwd=3, col="red")
axis(1, at=1:12,labels=month(1:12,label=T, abbr=T),  las=2)

plot(z20$PML_STG, type="l", lwd=3, xaxt="n", xlab="",ylab="PML_STG (oC)", ylim=c(15,35))
points(z21$PML_STG, type="l", col="red", lwd=3)
for(i in 2:dim(CI_p)[3]){
  df<-as.data.frame(CI_p[,,i])
  names(df)<-c("PML_STG", "WA_LSTC", "EA_LSTC", "EP_LSTC")
  df_m<-df[1:!2,]
  for(j in 1:12){
    df_m[j,]<-apply(df[which(month(time[2030:3000])==j & year(time[2030:3000])>=2070 & year(time[2030:3000])<2100),],2,mean,na.rm=T)
  }
  points(df_m$PML_STG, type="l", lwd=1, col="pink")
  
  
}
points(z20$PML_STG, type="l", lwd=3)
points(z21$PML_STG, type="l", lwd=3, col="red")
axis(1, at=1:12,labels=month(1:12,label=T, abbr=T),  las=2)


library(caret)

# Set training control
CI$month<-as.factor(month(CI$time))
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats = 5,
                              search = "random",
                              verboseIter = TRUE)

# Train the model
elm <- train(P~month + PML_STG + WA_LSTC + EA_LSTC,
                             #month:PML_STG+
                             #month:WA_LSTC+
                             #month:EA_LSTC, 
               data=CI,
                           method = "glmnet",
                           preProcess = c("center", "scale"),
                           tuneLength = 25,
                           trControl = train_control)
plot(coef(elm$finalModel, elm$finalModel$lambdaOpt))

# Check multiple R-squared
y_hat_elnet <- predict(elm, CI)
(rsq_enet <- cor(CI$P, y_hat_elnet)^2)


library("ggplot2"); theme_set(theme_bw())
## squash panels together
zero_margin <- theme(panel.spacing=grid::unit(0,"lines")) 
library("lattice")
library("dotwhisker")  ## coefficient plots


library(ggplot2)
CI$month<-as.factor(month(CI$time,label=T, abbr=T))
ggplot(CI, aes(x = PML_STG, y= P, color = month)) +
  geom_point() +
  geom_smooth(method='lm',formula=y~x)+
  facet_grid(. ~ (month)) 
  
ggplot(CI, aes(x = WA_LSTC, y= P, color = month)) +
  geom_point() +
  geom_smooth(method='lm',formula=y~x)+
  facet_grid(. ~ (month)) 

ggplot(CI, aes(x = EA_LSTC, y= P, color = month)) +
  geom_point() +
  geom_smooth(method='lm',formula=y~x)+
  facet_grid(. ~ (month)) 

  
names(CI)
CI$month<-as.factor(month(CI$time))
stdCI<-CI
stdCI[,2:14]<-apply(CI[,2:14],2,function(x)(x-mean(x))/sd(x))

ggplot(stdCI, aes(x = PML_STG, y= P, color = month)) +
  geom_point() +
  geom_smooth(method='lm',formula=y~x)+
  facet_grid(. ~ (month)) 

ggplot(stdCI, aes(x = WA_LSTC, y= P, color = month)) +
  geom_point() +
  geom_smooth(method='lm',formula=y~x)+
  facet_grid(. ~ (month)) 

ggplot(stdCI, aes(x = EA_LSTC, y= P, color = month)) +
  geom_point() +
  geom_smooth(method='lm',formula=y~x)+
  facet_grid(. ~ (month)) 

null<-lmer(P~(1|month), data=CI)
rs<-lmer(P~PML_STG + WA_LSTC + EA_LSTC+
           (1|month),
         data=CI)


varex<-data.frame(explained=c((809.0-15.3),(363.2-324.5),15.3, 324.5), category=c("monthly","residual", "monthly", "residual"), grp=c("explained", "explained", "unexplained", "unexplained")) 

library(ggplot2)
library(viridis)
library(hrbrthemes)


# Small multiple
ggplot(varex, aes(fill=grp, y=explained, x=category)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle("% variance explained by linear regression on surface temperature gradient") +
  theme_ipsum() +
  xlab("")

library(tidyverse) #for all data wrangling
library(cowplot) #for manuscript ready figures
library(lme4) #for lmer & glmer models
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
library(sjstats) #use for r2 functions

data(efc)
theme_set(theme_sjplot())