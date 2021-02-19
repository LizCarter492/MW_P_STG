
library(lubridate)
library(ncdf4)
library(fields)
library(maps)
library(maptools)
library(raster)
library(data.table)
library(RColorBrewer)
library(viridis)
library(colorRamps)
setwd("D:/DELL_PHD/Box Sync/AFRI_Aexam/R scripts")
source("filled.contour3.R")
source("filled.legend.R")

in_path<-("D:/DELL_PHD/E/NASH_subseasonal_prediction/")
out_path<-("H:/MW_seasonal_warm_season_precip_forecast/")

######################################################
#Variable names:######################################
#gu: geostrophic u wind component
#gv: geostrophic v wind component
#agu: ageostrophic u wind component
#agv: ageostrophic v wind component
#p: precipitation
#mdf: moisture flux divergence
#ivt: integrated vapor transport
#q: specific humidity
#o: omega (surface lifting index)
#z: geopotential height
#h: sensible heat flux
#
#ltm prefix indicates long-term monthly mean
#######################################################
#######################################################


#(u,v); (gu, gv); (agu, agv): load total, gesotrophic, and ageostrophic wind####
gUV<-nc_open(paste(in_path,"geos_uv_1000_600mB.nc", sep=""))
gu<-ncvar_get(gUV, varid="u_g")
gv<-ncvar_get(gUV, varid="v_g")
nc_close(gUV)

unc<-nc_open(paste(in_path,"uwnd.mon.mean.nc", sep=""))
u<-ncvar_get(unc, varid="uwnd")
u<-u[,,1:5,]
nc_close(unc)

vnc<-nc_open(paste(in_path,"vwnd.mon.mean.nc", sep=""))
v<-ncvar_get(vnc, varid="vwnd")
v<-v[,,1:5,]
nc_close(vnc)

agu<-u-gu
agu<-apply(agu[,,1:4,1:840], c(1,2,4), mean)
agv<-v-gv
agv<-apply(agv[,,1:4,1:840], c(1,2,4), mean)

gu<-apply(gu[,,1:4,1:840], c(1,2,4), mean)
gv<-apply(gv[,,1:4,1:840], c(1,2,4), mean)

u<-apply(u[,,1:4,1:840], c(1,2,4), mean)
v<-apply(v[,,1:4,1:840], c(1,2,4), mean)

ltmagu<-ltmagv<-ltmgu<-ltmgv<-ltmu<-ltmv<-
  array(NA, dim=c(dim(agu)[1:2],12))
mthi<-rep(1:12, length(1948:2017))
yeari<-rep(1948:2017, each=12)


for(j in 1:12){
  ltmagu[,,j]<-apply(agu[,,which(mthi==j)],c(1,2), mean)
  ltmagv[,,j]<-apply(agv[,,which(mthi==j)],c(1,2), mean)
  ltmgu[,,j]<-apply(gu[,,which(mthi==j)],c(1,2), mean)
  ltmgv[,,j]<-apply(gv[,,which(mthi==j)],c(1,2), mean)
  ltmu[,,j]<-apply(u[,,which(mthi==j)],c(1,2), mean)
  ltmv[,,j]<-apply(v[,,which(mthi==j)],c(1,2), mean)
}


#(P): load GPCC precipitation####
setwd(in_path)
P.nc<-nc_open("precip.comb.v2018to2016-v6monitorafter.total.nc")
plat<-ncvar_get(P.nc, "lat")
plon<-ncvar_get(P.nc, "lon")
ptime<-ncvar_get(P.nc, "time")
ptime<-as.POSIXct("1800-01-01 00:00")+as.difftime(ptime,units="days")
ptime<-as.Date(ptime, origin=as.Date("1800-01-01  00:00:0.0"))
p<-ncvar_get(P.nc, "precip")
p<-p[,,which(year(ptime)>=1948 & year(ptime)<2018)]

nc_close(P.nc)

#(z): load geopotential height#####
setwd(in_path)
Z.nc<-nc_open("hgt.mon.mean.nc")
lat<-ncvar_get(Z.nc, "lat")
lon<-ncvar_get(Z.nc, "lon")
lev<-ncvar_get(Z.nc, "level")
z<-ncvar_get(Z.nc, "hgt")
time<-ncvar_get(Z.nc, "time")
time<-as.POSIXct("1800-01-01 00:00")+as.difftime(time,units="hours")
time<-as.Date(time, origin=as.Date("1800-01-01  00:00:0.0"))
z<-z[,,,year(time)<2018]
nc_close(Z.nc)

ltmz<-array(NA, dim=c(dim(z)[1:3],12))
mthi<-rep(1:12, length(1948:2017))
yeari<-rep(1948:2017, each=12)


for(j in 1:12){
  ltmz[,,,j]<-apply(z[,,,which(mthi==j)],c(1,2,3), mean)
}

#(MFD): load moisture flux divergence####
MFD<-array(NA, dim=c(144, 73, 12*length(1948:2018)))
mthi<-rep(1:12, length(1948:2017))
yeari<-rep(1948:2017, each=12)
setwd(paste(in_path,"NCAR_NCEP_Daily_moisture_budget/", sep=""))
for(g in (1948:2018)){
  mfd.nc <- nc_open(paste("MFD_daily_", g, ".nc", sep=""))
  time<-ncvar_get(mfd.nc, varid="time")
  time<-as.POSIXct("1800-01-01 00:00")+as.difftime(time,units="hours")
  time<-as.Date(time, origin=as.Date("1800-01-01  00:00:0.0"))
  mfd<-ncvar_get(mfd.nc, "Div_UQVQ")
  nc_close(mfd.nc)
  
  mth<-month(time)
  
  for(j in 1:12){
    MFD[,,which(mthi==j & yeari==g)]<-apply(mfd[,,mth==j], c(1,2), mean)
  }
  remove(mfd)
}


ltmMFD<-array(NA, dim=c(dim(MFD)[1:2],12))
for(j in 1:12){
  ltmMFD[,,j]<-apply(MFD[,,which(mthi==j)],c(1,2), mean)
}


#(UQ, VQ, IVT): load integrated vapor transport####
UQ<-VQ<-IVT<-array(NA, dim=c(144, 73, 12*length(1948:2017)))
mthi<-rep(1:12, length(1948:2017))
yeari<-rep(1948:2017, each=12)
setwd(paste(in_path,"NCAR_NCEP_Daily_moisture_budget/", sep=""))
for(g in (1948:2017)){
  VINT.nc <- nc_open(paste("VINT_daily_", g, ".nc", sep=""))
  time<-ncvar_get(VINT.nc, varid="time")
  time<-as.POSIXct("1800-01-01 00:00")+as.difftime(time,units="hours")
  time<-as.Date(time, origin=as.Date("1800-01-01  00:00:0.0"))
  uq<-ncvar_get(VINT.nc, "UQ")
  vq<-ncvar_get(VINT.nc, "VQ")
  nc_close(VINT.nc)
  mth<-month(time)
  
  for(j in 1:12){
    UQ[,,which(mthi==j & yeari==g)]<-apply(uq[,,mth==j], c(1,2), mean)
    VQ[,,which(mthi==j & yeari==g)]<-apply(vq[,,mth==j], c(1,2), mean)
    
  }
  remove(uq); remove(vq); remove(time); remove(mth)
}

IVT<-sqrt(UQ^2 + VQ^2)

ltmUQ<-ltmVQ<-ltmIVT<-array(NA, dim=c(dim(MFD)[1:2],12))
for(j in 1:12){
  ltmUQ[,,j]<-apply(UQ[,,which(mthi==j)],c(1,2), mean)
  ltmVQ[,,j]<-apply(VQ[,,which(mthi==j)],c(1,2), mean)
  ltmIVT[,,j]<-apply(IVT[,,which(mthi==j)],c(1,2), mean)
}


#(q): load specific humidity####
setwd(in_path)
Q.nc<-nc_open("shum.mon.mean.nc")
lat<-ncvar_get(Q.nc, "lat")
lon<-ncvar_get(Q.nc, "lon")
q<-ncvar_get(Q.nc, "shum")
time<-ncvar_get(Q.nc, "time")
time<-as.POSIXct("1800-01-01 00:00")+as.difftime(time,units="hours")
time<-as.Date(time, origin=as.Date("1800-01-01  00:00:0.0"))
q<-apply(q[,,1:4,year(time)<2018], c(1,2,4), mean)


ltmq<-array(NA, dim=c(dim(q)[1:2],12))
for(j in 1:12){
  ltmq[,,j]<-apply(q[,,which(mthi==j)],c(1,2), mean)
}



#(o): load omega####
setwd(in_path)
O.nc<-nc_open("omega.mon.mean.nc")
lat<-ncvar_get(O.nc, "lat")
lon<-ncvar_get(O.nc, "lon")
o<-ncvar_get(O.nc, "omega")
time<-ncvar_get(O.nc, "time")
time<-as.POSIXct("1800-01-01 00:00")+as.difftime(time,units="hours")
time<-as.Date(time, origin=as.Date("1800-01-01  00:00:0.0"))
o<-apply(o[,,1:4,year(time)<2018], c(1,2,4), mean)


ltmo<-array(NA, dim=c(dim(o)[1:2],12))
for(j in 1:12){
  ltmo[,,j]<-apply(o[,,which(mthi==j)],c(1,2), mean)
}

#
#####################################################################################
#set map graphical parameters####

#map dimensions
lat.min <- 10#min(lat)
lat.max <- 75#max(lat)
lon.min <- -132#min(lon-360)
lon.max <- max(lon-360)
lat.min.ar <- lat.min+1
lat.max.ar <- lat.max-1
lon.min.ar <- lon.min+1
lon.max.ar <- lon.max-1

lat_final <- rev(lat)
lon_final <- lon - 360
keep1 <- which(lon_final>=lon.min & lon_final<=lon.max)
keep2 <- which(lat_final>=lat.min & lat_final<=lat.max)

lon_all <- rep(lon_final,length(lat_final))
lat_all <- sort(rep(lat_final,length(lon_final)))

month_name<-format(ISOdate(2004,1:12,1),"%B")
mth<-rep(1:12, length(1948:2017))
yr<-rep(1948:2017, each=12)

#image size 
hgt<-500
wdt<-1600


p_c_df<-rbind(c(1),c(2),c(3), c(4), c(5), c(6), c(7), c(8), c(9), c(10), c(11), c(12))
lat.min <- 10#min(lat)
lat.max <- 75#max(lat)
lon.min <- -132#min(lon-360)
lon.max <- max(lon-360)
lat.min.ar <- lat.min+1
lat.max.ar <- lat.max-1
lon.min.ar <- lon.min+1
lon.max.ar <- lon.max-1

lat_final <- rev(lat)
lon_final <- lon - 360
keep1 <- which(lon_final>=lon.min & lon_final<=lon.max)
keep2 <- which(lat_final>=lat.min & lat_final<=lat.max)

lon_all <- rep(lon_final,length(lat_final))
lat_all <- sort(rep(lat_final,length(lon_final)))

month_name<-format(ISOdate(2004,1:12,1),"%B")
mth<-rep(1:12, length(1948:2017))
yr<-rep(1948:2017, each=12)

#image size 
hgt<-1500 #4 #
wdt<-3200 #12.8



for(m in 1:12){
  m_c<-p_c_df[m,]
  #mfd
  ltmUQ_sub<-apply(ltmUQ[,,m_c],c(1,2),mean)
  ltmUQ_sub<-t(ltmUQ_sub)
  ltmUQ_sub<-ltmUQ_sub[seq(dim(ltmUQ_sub)[1],1),]
  
  ltmVQ_sub<-apply(ltmVQ[,,m_c],c(1,2),mean)
  ltmVQ_sub<-t(ltmVQ_sub)
  ltmVQ_sub<-ltmVQ_sub[seq(dim(ltmVQ_sub)[1],1),]
  
  ltmmfd_sub<-apply(ltmMFD[,,m_c], c(1,2), mean)*30
  ltmmfd_sub<-t(ltmmfd_sub)
  ltmmfd_sub<-ltmmfd_sub[seq(dim(ltmmfd_sub)[1],1),]
  #geo
  
  ltmz850_sub<-apply(ltmz[,,3,m_c],c(1,2),mean)
  ltmz850_sub<-t(ltmz850_sub)
  ltmz850_sub<-ltmz850_sub[seq(dim(ltmz850_sub)[1],1),]
  
  gu_sub<-apply(gu[,,which(mth %in% (m_c))],c(1,2), mean)
  gu_sub<-t(gu_sub)
  gu_sub<-gu_sub[seq(dim(gu_sub)[1],1),]
  
  gv_sub<-apply(gv[,,which(mth %in% (m_c))],c(1,2), mean)
  gv_sub<-t(gv_sub)
  gv_sub<-gv_sub[seq(dim(gv_sub)[1],1),]
  
  
  #ageo     
  ltmo_sub<-apply(o[,,which (mth %in% m_c)], c(1,2), mean)
  ltmo_sub<-t(ltmo_sub)
  ltmo_sub<-ltmo_sub[seq(dim(ltmo_sub)[1],1),]
  
  agu_sub<-apply(agu[,,which(mth %in% (m_c))],c(1,2), mean)
  agu_sub<-t(agu_sub)
  agu_sub<-agu_sub[seq(dim(agu_sub)[1],1),]
  
  agv_sub<-apply(agv[,,which(mth %in% (m_c))],c(1,2), mean)
  agv_sub<-t(agv_sub)
  agv_sub<-agv_sub[seq(dim(agv_sub)[1],1),]
  
  png(paste("H:/Research/MW_ClimateChange/images/climatological maps/", month_name[m_c[1]],"_", month_name[m_c[2]],"climatologyREVISED.png",sep=""),
      res=300,
      #pointsize=15,
      type="cairo",
      width=wdt, height=hgt, 
  )

  par(mfrow=c(1,2)) 
  
  
  #geo
  
  uua<-as.vector(t(gu_sub))
  vva<-as.vector(t(gv_sub))
  
  uu<-uua
  vv<-vva

  speed <- sqrt(uu*uu+ vv*vv)
  speeda <- sqrt(uua*uua + vva*vva)
  uv <-which(speed>speeda & (1:length(uu))%in%seq(2,length(uu),by=4))
  uvall<-which(speed<=speeda & (1:length(uu))%in%seq(2,length(uu),by=4))
 
  mylevel<-c(1200, seq(1400, 1600, by=5), 1700)
  mycol<-topo.colors(length(mylevel)+1)
  
  mnm<-c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")[m]
  filled.contour3(lon_final,lat_final,t(ltmz850_sub),
                  xlim=c(lon.min,lon.max),
                  ylim=c(lat.min,lat.max),
                  level=mylevel,
                  col=mycol,
                  #frame.plot=FALSE,
                  main=title(mnm,cex.main=3),
                  
                  plot.axes={
                    axis(1, labels=FALSE, tick=FALSE);
                    axis(2, labels=FALSE, tick=FALSE);
                    if(m %in% 4:9){
                    arrow.plot(a1=lon_all[uv],a2=lat_all[uv],u=uua[uv],v=vva[uv],arrow.ex=0.5,length=0.06,xpd=FALSE,
                               xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3, col="gray40");
                    arrow.plot(a1=lon_all[uvall],a2=lat_all[uvall],u=uua[uvall],v=vva[uvall],arrow.ex=0.5,length=0.06,xpd=FALSE, 
                               xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3, col="gray40");
                    }else{
                    arrow.plot(a1=lon_all[uv],a2=lat_all[uv],u=uua[uv],v=vva[uv],arrow.ex=0.5,length=0.06,xpd=FALSE,
                               xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3, col="gray40");
                    arrow.plot(a1=lon_all[uvall],a2=lat_all[uvall],u=uua[uvall],v=vva[uvall],arrow.ex=0.5,length=0.06,xpd=FALSE, 
                               xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3, col="gray40");
                    }  
                    if(m_c[1]<5){ 
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1480,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1420,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1540,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub),lwd=1,lty=1, levels=c(1440,1500,1520,1540, 1560),
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max),add=TRUE);
                      
                    }else if(m_c[1]<6) {
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1480,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1420,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1540,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1500,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1440,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub),lwd=1,lty=1, levels=c(1440,1500,1520,1540, 1560),
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max),add=TRUE);
                      
                      
                    }else{
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1480,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1420,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1540,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1510,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1460,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, col="black",add=TRUE);
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1560,
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=4,lty=1, add=TRUE);
                      contour(x=lon_final, y=lat_final, z=t(ltmz850_sub),lwd=1,lty=1, levels=c(1440,1500,1520,1540),
                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max),  add=TRUE);
                      
                      
                    }
                    
                    map("world",add=T, lwd=3);
                  }
                  
  )
  
  # #ageo
  # uua<-as.vector(t(agu_sub))
  # vva<-as.vector(t(agv_sub))
  # uu<-uua
  # uua[uu>quantile(uu, 0.99, na.rm=T)]<-NA#quantile(uu, 0.95, na.rm=T)
  # #uua[uu<quantile(uu, 0.01, na.rm=T)]<-NA#quantile(uu, 0.05, na.rm=T)
  # vv<-vva
  # vva[vva>quantile(vv, 0.99, na.rm=T)]<-NA#quantile(vva, 0.95, na.rm=T)
  # #vva[vv<quantile(vv, 0.05, na.rm=T)]<-NA#quantile(vv, 0.05, na.rm=T)
  # #uua<-ifelse(uua>=0, sqrt(abs(uua)), -1*sqrt(abs(uua)))
  # #vva<-ifelse(vva>=0, sqrt(abs(vva)), -1*sqrt(abs(vva)))
  # speed <- sqrt(uua*uua+ vva*vva)
  # speeda <- sqrt(uua*uua + vva*vva)
  # uv <-which(speed>speeda & (1:length(uu))%in%seq(1,length(uu),by=2))
  # uvall<-which(speed<=speeda & (1:length(uu))%in%seq(1,length(uu),by=2))
  # 
  # mylevel<-c(-0.2, seq(-0.02, 0.02, by=0.001), 0.5)
  # mycol<-cm.colors(length(mylevel)+1)
  # 
  # filled.contour3(lon_final,lat_final,t(ltmo_sub),
  #                 xlim=c(lon.min,lon.max),
  #                 ylim=c(lat.min,lat.max),
  #                 level=mylevel,
  #                 col=mycol,
  #                 main=title(mnm,cex.main=3),
  #                 plot.axes={
  #                   axis(1, labels=FALSE, tick=FALSE);
  #                   axis(2, labels=FALSE, tick=FALSE);
  #                   arrow.plot(a1=lon_all[uv],a2=lat_all[uv],u=uua[uv],v=vva[uv],arrow.ex=0.5,length=0.06,xpd=FALSE,
  #                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3, col="gray40");
  #                   arrow.plot(a1=lon_all[uvall],a2=lat_all[uvall],u=uua[uvall],v=vva[uvall],arrow.ex=0.5,length=0.06,xpd=FALSE, 
  #                              xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3, col="gray40");
  #                   contour(x=lon_final, y=lat_final, z=t(ltmo_sub), levels=c(0.008),
  #                           xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=2,lty=1, col="black",add=TRUE);
  #                   contour(x=lon_final, y=lat_final, z=t(ltmo_sub), levels= c(-0.008),lty=2,
  #                           xlim=c(lon.min,lon.max),ylim=c(lat.min+5,lat.max), lwd=2, col="black",add=TRUE);
  #                   
  #                   map("world",add=T, lwd=3);
  #                 }
  #                 
  # )
  
  #mfd
  library(colorRamps)
  mylevel<-c(-7, seq(-3, 3, by=0.1), 10)
  mycol<-rev(matlab.like(length(mylevel)+1))
  uua<-as.vector(t(ltmUQ_sub))
  vva<-as.vector(t(ltmVQ_sub))
  #uua<-ifelse(uua>=0, sqrt(abs(uua)), -1*sqrt(abs(uua)))
  #vva<-ifelse(vva>=0, sqrt(abs(vva)), -1*sqrt(abs(vva)))
  
  uvall<-which((1:length(uua))%in%seq(1,length(uua),by=3))
  
  filled.contour3(lon_final,lat_final,t(ltmmfd_sub),
                  xlim=c(lon.min,lon.max),
                  ylim=c(lat.min,lat.max),
                  level=mylevel,
                  col=mycol,
                  #frame.plot=FALSE,
                  key.title=title(main=paste(j, k, sep=" "),cex=2),
                  
                  plot.axes={
                    axis(1, labels=FALSE, tick=FALSE);
                    axis(2, labels=FALSE, tick=FALSE);
                    arrow.plot(a1=lon_all[uvall],a2=lat_all[uvall],u=uua[uvall],v=vva[uvall],arrow.ex=0.4,length=0.06,xpd=FALSE,
                               xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3, col="grey40");
                    contour(x=lon_final, y=lat_final, z=t(ltmmfd_sub), levels= c(-1.5),lty=2,
                            xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3, col="black",add=TRUE);
                    
                    contour(x=lon_final, y=lat_final, z=t(ltmmfd_sub), levels= c(1.5),lty=1,
                            xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3, col="black",add=TRUE);
                    
                    
                    
                    
                    map("world",add=T, lwd=3);
                  }
                  
  )
  
  
  dev.off()
  
  
}


