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

#read in lake Ontario lake levels
ll<-read.csv("H:/Research/MW_ClimateChange/data/OnatarioLakeLevels/ontario1918.csv")
CI$ll<-rep(NA)
for(i in unique(year(CI$time))){
        for(j in unique(month(CI$time))){
        CI$ll[which(year(CI$time)==i & month(CI$time)==j)]<-ll[which(ll$year==i),(j+1)]
        }
}

# CI$PML_EA<-CI$PML_STG*CI$EA_LSTC
# CI$PML_WA<-CI$PML_STG*CI$WA_LSTC
# CI$EA_WA<-CI$WA_LSTC*CI$EA_LSTC
# m1<-lm(ll~EA_LSTC + WA_LSTC + PML_STG + PML_EA+PML_WA+EA_WA, data=CI)

###### get coordinates definiting position of NASH, ERT, and NEPL####
library(lubridate)
library(ncdf4)
library(fields)
library(maps)
library(raster)
# setwd("D:/Users/ekc76/Box Sync/AFRI_Aexam/R scripts")
# source("filled.contour3.R")
# source("filled.legend.R")
library(sp)
library(rgeos)


NASH_e<-extent(260, 320, 15, 60)
ERT_e<-extent(360-114,360-94,32,42)
NEPL_e<-extent(260,300, 35, 60)

#update for August by downloading geopotential height monthly means file from here: https://www.esrl.noaa.gov/psd/data/gridded/data.gpcc.html
pre.nc <- nc_open("hgt.mon.mean.nc" )
time <- ncvar_get(pre.nc, varid="time")
time<-as.POSIXct("1800-01-01 00:00")+as.difftime(time,units="hours")
time<-as.Date(time, origin=as.Date("1800-01-01  00:00:0.0"))
mth<-month(time)
yr<-year(time)
nc_close(pre.nc)

hgt<-brick("hgt.mon.mean.nc" , varname="hgt", level=3)
u<-brick("uwnd.mon.mean.nc" , varname="uwnd", level=3)
v<-brick("vwnd.mon.mean.nc" , varname="vwnd", level=3)

#NASH_TS
hgt.r<-crop(hgt, NASH_e)
r.u<-crop(hgt, NASH_e)

int<-as.data.frame(1948:2019)
names(int)<-"year"
int$int<-NA
int$lat<-NA
int$lon<-NA




for(y in 1948:2017){
        hgt.AM<- mean(r.hgt[[which((mth==4 |mth==5) &yr %in% y)]])
        u.AM<-mean(r.u[[which((mth==4|mth==5)& yr %in% y)]])
        int$int[int$year==y]<-max(as.vector(hgt.AM))
        if(max(as.vector(hgt.AM))>1540){
                
                c.hgt<-rasterToContour(hgt.AM,levels=1540)
                c.u<-rasterToContour(u.AM,levels=0)
                int.pts <- gIntersection(c.hgt, c.u, byid = TRUE)
                if(is.null(int.pts)){
                        int$lat[which(int$year==y)] <- NA
                        int$lon[which(int$year==y)] <- NA
                }else{
                        int.coords <- int.pts@coords
                        if(dim(int.coords)[1]>1){
                                int$lat[which(int$year==y)] <- NA
                                int$lon[which(int$year==y)] <- NA
                        }else{
                                int$lon[which(int$year==y)] <- int.coords[1,1]
                                int$lat[which(int$year==y)] <- int.coords[1,2]
                        }
                }
                
        } else {
                int$lat[which(int$year==y)] <- NA
                int$lon[which(int$year==y)] <- NA
        }
}

y<-int$year[is.na(int$lat)]


write.table(int, "C:/Users/ekc76/Box Sync/Ont_flood_17/climatology/AM_NASH_1540_ridge_intensity.txt")




for(y in 1948:2017){
        hgt.AM<- mean(r.hgt[[which((mth==6 |mth==7) &yr %in% y)]])
        u.AM<-mean(r.u[[which((mth==6|mth==7)& yr %in% y)]])
        int$int[int$year==y]<-max(as.vector(hgt.AM))
        if(max(as.vector(hgt.AM))>1550){
                
                c.hgt<-rasterToContour(hgt.AM,levels=1550)
                c.u<-rasterToContour(u.AM,levels=0)
                int.pts <- gIntersection(c.hgt, c.u, byid = TRUE)
                if(is.null(int.pts)){
                        int$lat[which(int$year==y)] <- NA
                        int$lon[which(int$year==y)] <- NA
                }else{
                        int.coords <- int.pts@coords
                        if(dim(int.coords)[1]>1){
                                int$lat[which(int$year==y)] <- NA
                                int$lon[which(int$year==y)] <- NA
                        }else{
                                int$lon[which(int$year==y)] <- int.coords[1,1]
                                int$lat[which(int$year==y)] <- int.coords[1,2]
                        }
                }
                
        } else {
                int$lat[which(int$year==y)] <- NA
                int$lon[which(int$year==y)] <- NA
        }
}

y<-int$year[is.na(int$lat)]


write.table(int, "C:/Users/ekc76/Box Sync/Ont_flood_17/climatology/JJ_NASH_1550_ridge_intensity.txt")


for(y in 1948:2017){
        hgt.AM<- mean(r.hgt[[which((mth==8 |mth==9) &yr %in% y)]])
        u.AM<-mean(r.u[[which((mth==8|mth==9)& yr %in% y)]])
        int$int[int$year==y]<-max(as.vector(hgt.AM))
        if(max(as.vector(hgt.AM))>1560){
                
                c.hgt<-rasterToContour(hgt.AM,levels=1560)
                c.u<-rasterToContour(u.AM,levels=0)
                int.pts <- gIntersection(c.hgt, c.u, byid = TRUE)
                if(is.null(int.pts)){
                        int$lat[which(int$year==y)] <- NA
                        int$lon[which(int$year==y)] <- NA
                }else{
                        int.coords <- int.pts@coords
                        if(dim(int.coords)[1]>1){
                                int$lat[which(int$year==y)] <- NA
                                int$lon[which(int$year==y)] <- NA
                        }else{
                                int$lon[which(int$year==y)] <- int.coords[1,1]
                                int$lat[which(int$year==y)] <- int.coords[1,2]
                        }
                }
                
        } else {
                int$lat[which(int$year==y)] <- NA
                int$lon[which(int$year==y)] <- NA
        }
}

y<-int$year[is.na(int$lat)]



t2 = cbind(main_z$time, main_z$PML_STG)
t1 = cbind(main_z$time, main_z$P)
# Specify the number of iterations. The more, the better (>1000).  For the
# purpose of this tutorial, we just set it = 10
nrands = 100

wtc.AB = wtc(t1, t2, nrands = nrands)

par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12, ylab = "Scale", xlab = "Period", 
     plot.cb = TRUE, main = "Wavelet Coherence: P v PML_STG")


t2 = cbind(CI$time, CI$EA_LSTC)
t1 = cbind(CI$time, CI$ll)
# Specify the number of iterations. The more, the better (>1000).  For the
# purpose of this tutorial, we just set it = 10
nrands = 100

wtc.AB = wtc(t1, t2, nrands = nrands)

par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12, ylab = "Scale", xlab = "Period", 
     plot.cb = TRUE, main = "Wavelet Coherence: P v EA_LSTC")



t2 = cbind(CI$time, CI$PML_STG)
t1 = cbind(CI$time, CI$ll)
# Specify the number of iterations. The more, the better (>1000).  For the
# purpose of this tutorial, we just set it = 10
nrands = 100

wtc.AB = wtc(t1, t2, nrands = nrands, mother="dog")

par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12, ylab = "Scale", xlab = "Period", 
     plot.cb = TRUE, main = "Wavelet Coherence: P v WA_LSTC")


t2 = cbind(CI$time, CI$EP_LSTC)
t1 = cbind(CI$time, CI$ll)
# Specify the number of iterations. The more, the better (>1000).  For the
# purpose of this tutorial, we just set it = 10
nrands = 100

wtc.AB = wtc(t1, t2, nrands = nrands)

par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12, ylab = "Scale", xlab = "Period", 
     plot.cb = TRUE, main = "Wavelet Coherence: P v EP_LSTC")

cc<-ccf(main_z$ll,main_z$PML_STG, lag.max=120)
cc<-ccf(main_z$ll,main_z$EA_LSTC, lag.max=120)
cc<-ccf(main_z$ll,main_z$WA_LSTC, lag.max=240)

cc<-ccf(main_z$P,main_z$PML_STG, lag.max=12)
cc<-ccf(main_z$P,main_z$EA_LSTC, lag.max=12)
cc<-ccf(main_z$P,main_z$WA_LSTC, lag.max=12)

CI_z<-as.data.frame(apply(CI[,-1], 2, function(x)(x-mean(x))/sd(x)))
CI<-CI[,-19]
main_z<-CI

for(i in 1:12){
        main_z[month(CI$time)==i,-1]<-apply(CI[month(CI$time)==i,-1],2, function(x)(x-mean(x))/sd(x))
}

CI$month<-as.factor(month(CI$time))

m1<-lm(ll~EA_LSTC+WA_LSTC+PML_STG+PML_WA+PML_EA+EA_WA, data=main_z)
m2<-lmer(P~EA_LSTC+WA_LSTC+PML_STG+PML_WA+PML_EA+EA_WA+(1|month), data=dat)
m3<-lm(ll~month+EA_LSTC+WA_LSTC+PML_STG+month:PML_STG + month:WA_LSTC + month:EA_LSTC, data=CI)
library(glmnet)
dat<-main_z
dat$month<-month(dat$time)
jan<-dat[dat$month==1,]
feb<-dat[dat$month==2,]
mar<-dat[dat$month==3,]
apr<-dat[dat$month==4,]
may<-dat[dat$month==5,]
jun<-dat[dat$month==6,]
jul<-dat[dat$month==7,]
aug<-dat[dat$month==8,]
sep<-dat[dat$month==9,]
oct<-dat[dat$month==10,]
nov<-dat[dat$month==11,]
dec<-dat[dat$month==12,]
for(i in c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")){
        d<-get(i)
        m<-lm(P~EA_LSTC+WA_LSTC+PML_STG, data=d)
        assign(paste(i, "_m", sep=""), summary(m))
        print(summary(m))
}

plot(CI_z$P, type="l", ylim=c(-3,3))
points(CI_z$PML_STG, type="l", col="green")
points(CI_z$EA_LSTC, type="l", col="red")
points(CI_z$WA_LSTC, type="l", col="purple")

dat<-main_z
x<-as.matrix(dat[,c("PML_STG", "EA_LSTC", "WA_LSTC", "PML_EA", "PML_WA", "EA_WA")])
y<-as.matrix(dat$P)
f1<-glmnet(x,y)
plot(y, type="l")
points(predict(f1, newx=x,s=0), type="l", col="red")

monthly_means<-CI[1:12,]
for(i in 1:12){
        monthly_means[i,-1]<-apply(CI[which(month(CI$time)==i),-1],2,mean)
}


points(predict(m1), type="l", col="green")
plot(CI$P, type="l")
plot(resid(m1)~CI$WA_LSTC)

