
library(lubridate)
library(ncdf4)
library(fields)
library(maps)
library(raster)
setwd("C:/Users/ekc76/Box Sync/AFRI_Aexam/R scripts")
source("filled.contour3.R")
source("filled.legend.R")

#Variable names:
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

#ltm prefix indicates long-term monthly mean

#(u,v); (gu, gv); (agu, agv): load total, gesotrophic, and ageostrophic wind####
gUV<-nc_open("E:/NASH_subseasonal_prediction/geos_uv_1000_600mB.nc")
gu<-ncvar_get(gUV, varid="u_g")
gv<-ncvar_get(gUV, varid="v_g")
nc_close(gUV)

unc<-nc_open("E:/NASH_subseasonal_prediction/uwnd.mon.mean.nc")
u<-ncvar_get(unc, varid="uwnd")
u<-u[,,1:5,]
nc_close(unc)

vnc<-nc_open("E:/NASH_subseasonal_prediction/vwnd.mon.mean.nc")
v<-ncvar_get(vnc, varid="vwnd")
v<-v[,,1:5,]
nc_close(vnc)

agu<-u-gu
agu<-apply(agu[,,1:3,1:840], c(1,2,4), mean)
agv<-v-gv
agv<-apply(agv[,,1:3,1:840], c(1,2,4), mean)

gu<-apply(gu[,,1:3,1:840], c(1,2,4), mean)
gv<-apply(gv[,,1:3,1:840], c(1,2,4), mean)

u<-apply(u[,,1:3,1:840], c(1,2,4), mean)
v<-apply(v[,,1:3,1:840], c(1,2,4), mean)

ltmagu<-ltmagv<-ltmgu<-ltmgv<-ltmu<-ltmv<-array(NA, dim=c(dim(gu)[1:2],12))
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
setwd("E:/NASH_subseasonal_prediction/")
P.nc<-nc_open("precip.comb.v2018to2016-v6monitorafter.total.nc")
plat<-ncvar_get(P.nc, "lat")
plon<-ncvar_get(P.nc, "lon")
ptime<-ncvar_get(P.nc, "time")
ptime<-as.POSIXct("1800-01-01 00:00")+as.difftime(ptime,units="days")
ptime<-as.Date(ptime, origin=as.Date("1800-01-01  00:00:0.0"))
p<-ncvar_get(P.nc, "precip")
p<-p[,,year(ptime)<2018]
nc_close(P.nc)

#(z): load geopotential height#####
setwd("E:/NASH_subseasonal_prediction/")
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
setwd("E:/NASH_subseasonal_prediction/NCAR_NCEP_Daily_moisture_budget/")
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
UQ<-VQ<-IVT<-array(NA, dim=c(144, 73, 12*length(1948:2018)))
mthi<-rep(1:12, length(1948:2017))
yeari<-rep(1948:2017, each=12)
setwd("E:/NASH_subseasonal_prediction/NCAR_NCEP_Daily_moisture_budget/")
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
setwd("E:/NASH_subseasonal_prediction/")
Q.nc<-nc_open("shum.mon.mean.nc")
lat<-ncvar_get(Q.nc, "lat")
lon<-ncvar_get(Q.nc, "lon")
q<-ncvar_get(Q.nc, "shum")
time<-ncvar_get(Q.nc, "time")
time<-as.POSIXct("1800-01-01 00:00")+as.difftime(time,units="hours")
time<-as.Date(time, origin=as.Date("1800-01-01  00:00:0.0"))
q<-apply(q[,,1:5,year(time)<2018], c(1,2,4), mean)


ltmq<-array(NA, dim=c(dim(q)[1:2],12))
for(j in 1:12){
  ltmq[,,j]<-apply(q[,,which(mthi==j)],c(1,2), mean)
}



#(o): load omega####
#setwd("E:/NASH_subseasonal_prediction/")
setwd("G:/DELL_PHD/E/NASH_subseasonal_prediction/")

O.nc<-nc_open("omega.mon.mean.nc")
lat<-ncvar_get(O.nc, "lat")
lon<-ncvar_get(O.nc, "lon")
o<-ncvar_get(O.nc, "omega")
time<-ncvar_get(O.nc, "time")
time<-as.POSIXct("1800-01-01 00:00")+as.difftime(time,units="hours")
time<-as.Date(time, origin=as.Date("1800-01-01  00:00:0.0"))
o<-apply(o[,,1:3,year(time)<2018], c(1,2,4), mean)

ltmo<-array(NA, dim=c(dim(o)[1:2],12))
for(j in 1:12){
  ltmo[,,j]<-apply(o[,,which(mthi==j)],c(1,2), mean)
}

#(h): load sensible heat flux#####


#load climate variables and convert to factors####
library(lubridate)
CI<-read.table("C:/Users/ekc76/Box Sync/MW_seasonal_warm_season_precip_forecast/data/standardized_climate_indices.txt", header=T)
time<-seq.Date(as.Date("01-01-1948", format="%m-%d-%Y"), as.Date("12-31-2017", format="%m-%d-%Y"), by="month")
mth<-month(time)
yr<-year(time)
ss<-1948:2017
P<-CI$P[mth==7]
Pcat<-as.factor(ifelse(P<=(-0.3),1,ifelse(P>=0.3,3,2)))



SMI<-data.frame(matrix(rep(NA), nrow= length(ss), ncol=2))
names(SMI)<-c("SASMI", "NASMI")

for(i in 1:length(ss)){
  SMI[i,]<-apply(CI[which(mth %in% 6:8 & yr %in% ss[i]),c("SASMI", "NASMI")], 2, mean, na.rm=T)
}


CI_cat<-CI
for(j in 1:dim(CI)[2]){
  CI_cat[,j]<-ifelse(CI[,j]<=(-0.75),1,ifelse(CI[,j]>=0.75,3,2))
}


SMI_f<-SMI
for(j in 1:2){
  SMI_f[,j]<-ifelse(SMI[,j]<=(-0.66),1,ifelse(SMI[,j]>=0.66,3,2))
}
SMI_f <- as.data.frame(apply(SMI_f, 2, as.factor))

#set map graphical parameters

#map dimensions
lat.min <- 0#min(lat)
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


#lag between climate indices and composite
for(lagv in c(0,1,2)){

#absolute values maps####
for (i in 3:7){
  dat<-CI_cat[mth==i,-1]
  yr<-1948:2017
  for(j in names(dat)){
    #Map 1: JJA PJS, NASH anomaly (dark contour) mean position (gray contour)
    # MFD (color), geostrophic wind absolute (dark color when above average, gray when average/below average)
   
     pdf(paste("C:/Users/ekc76/Box Sync/MW_seasonal_warm_season_precip_forecast/images/climate_variables/climate_comp_maps/lag", lagv, "/ageostrophic/ageostrophic_abs_mfd_col_"
               ,j,"_month",i,"_lag",lagv,".pdf",sep=""),
        width=14, height=5)
    par(mfrow=c(1,3)) 
    nameslist<-(paste(c("low", "med", "high"), j, sep=" "))
    
    var<-dat[,j]
    L<-yr[var==1]
    M<-yr[var==2]
    H<-yr[var==3]
    for(k in c("L", "M", "H")){
      z300_sub<-apply(z[,,8,which(yr %in% get(k) & mth %in% (i+lagv))], c(1,2), mean)
      z300_sub<-t(z300_sub)
      z300_sub<-z300_sub[seq(dim(z300_sub)[1],1),]
      
      z850_sub<-apply(z[,,3,which(yr %in% get(k) & mth %in% (i+lagv))], c(1,2), mean)
      z850_sub<-t(z850_sub)
      z850_sub<-z850_sub[seq(dim(z850_sub)[1],1),]
     
      z925_sub<-apply(z[,,4,which(yr %in% get(k) & mth %in% (i+lagv))], c(1,2), mean)
      z925_sub<-t(z925_sub)
      z925_sub<-z925_sub[seq(dim(z925_sub)[1],1),]
      
      gu_sub<-apply(gu[,,which(yr %in% get(k) & mth %in% (i+lagv))],c(1,2), mean)
      gu_sub<-t(gu_sub)
      gu_sub<-gu_sub[seq(dim(gu_sub)[1],1),]
      
      gv_sub<-apply(gv[,,which(yr %in% get(k) & mth %in% (i+lagv))],c(1,2), mean)
      gv_sub<-t(gv_sub)
      gv_sub<-gv_sub[seq(dim(gv_sub)[1],1),]
      
      mfd_sub<-apply(MFD[,,which(yr %in% get(k) & mth %in% (i+lagv))],c(1,2), mean)
      mfd_sub<-t(mfd_sub)
      mfd_sub<-mfd_sub[seq(dim(mfd_sub)[1],1),]
      
      ltmz300_sub<-ltmz[,,8,i+lagv]
      ltmz300_sub<-t(ltmz300_sub)
      ltmz300_sub<-ltmz300_sub[seq(dim(ltmz300_sub)[1],1),]
      
      ltmz850_sub<-ltmz[,,3,i+lagv]
      ltmz850_sub<-t(ltmz850_sub)
      ltmz850_sub<-ltmz850_sub[seq(dim(ltmz850_sub)[1],1),]
      
      ltmz925_sub<-ltmz[,,4,i+lagv]
      ltmz925_sub<-t(ltmz925_sub)
      ltmz925_sub<-ltmz925_sub[seq(dim(ltmz925_sub)[1],1),]
      
      
      ltmgu_sub<-ltmgu[,,i+lagv]
      ltmgu_sub<-t(ltmgu_sub)
      ltmgu_sub<-ltmgu_sub[seq(dim(ltmgu_sub)[1],1),]
      
      ltmgv_sub<-ltmgv[,,i+lagv]
      ltmgv_sub<-t(ltmgv_sub)
      ltmgv_sub<-ltmgv_sub[seq(dim(ltmgv_sub)[1],1),]
      
      ltmmfd_sub<-ltmMFD[,,i+lagv]
      ltmmfd_sub<-t(ltmmfd_sub)
      ltmmfd_sub<-ltmmfd_sub[seq(dim(ltmmfd_sub)[1],1),]
      
      #gu_sub <- t(gu_sub)[keep1,keep2]
      #gv_sub <- t(gv_sub)[keep1,keep2]
      #ltmgu_sub <- t(ltmgu_sub)[keep1,keep2]
      #ltmgv_sub <- t(ltmgv_sub)[keep1,keep2]
      
      uu <- as.vector(t(gu_sub))
      vv <- as.vector(t(gv_sub))
      uua<-as.vector(t(ltmgu_sub))
      vva<-as.vector(t(ltmgv_sub))
      
      speed <- sqrt(uu*uu+ vv*vv)
      speeda <- sqrt(uua*uua + vva*vva)
      uv <-which(speed>speeda & (1:length(uu))%in%seq(1,length(uu),by=3))
      uvall<-which(speed<=speeda & (1:length(uu))%in%seq(1,length(uu),by=3))
      
      mylevel <- c(-5,seq(-0.4,0.4,length.out=18),5)
      mycol <- colorRampPalette(c("red","white","blue"))(length(mylevel)-1) 
      
     
      filled.contour3(lon_final,lat_final,t(mfd_sub),
                      xlim=c(lon.min,lon.max),
                      ylim=c(lat.min,lat.max),
                      level=mylevel,
                      col=mycol,
                      #frame.plot=FALSE,
                      key.title=title(main=paste(j, k, sep=" "),cex=2),
                      
                      plot.axes={               
                        axis(1, labels=FALSE, tick=FALSE);
                        axis(2, labels=FALSE, tick=FALSE);
                        arrow.plot(a1=lon_all[uv],a2=lat_all[uv],u=uu[uv],v=vv[uv],arrow.ex=0.17,length=0.2,xpd=FALSE,
                                   xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=2, col="gray40");
                        arrow.plot(a1=lon_all[uvall],a2=lat_all[uvall],u=uu[uvall],v=vv[uvall],arrow.ex=0.17,length=0.2,xpd=FALSE, 
                                   xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=2, col="gray75");
                       if(i<=5){ 
                        contour(x=lon_final, y=lat_final, z=t(z300_sub), levels=c(9060),col="green",
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3,lty=1, add=TRUE);
                        contour(x=lon_final, y=lat_final, z=t(ltmz300_sub), levels=9060,
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=0.8,lty=1, col="black",add=TRUE);
                        contour(x=lon_final, y=lat_final, z=t(z850_sub), levels=1540,col="green",
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3,lty=1, add=TRUE);
                        contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1540,
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=0.8,lty=1, col="black",add=TRUE);
                         contour(x=lon_final, y=lat_final, z=t(z850_sub),lwd=1,lty=1, levels=c(1480,1500,1520,1540),
                                 xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), col="green",add=TRUE);
                    
                         
                         
                         }else{
                         contour(x=lon_final, y=lat_final, z=t(z300_sub), levels=c(9100),col="green",
                                 xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3,lty=1, add=TRUE);
                         contour(x=lon_final, y=lat_final, z=t(ltmz300_sub), levels=9100,
                                 xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=0.8,lty=1, col="black",add=TRUE);
                         contour(x=lon_final, y=lat_final, z=t(z850_sub), levels=1560,col="green",
                                 xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3,lty=1, add=TRUE);
                         contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1560,
                                 xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=0.8,lty=1, col="green",add=TRUE);
                           contour(x=lon_final, y=lat_final, z=t(z850_sub),lwd=1,lty=1, levels=c(1480,1500,1520,1540),
                                   xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max),  col="green",add=TRUE);
                           
                           
                       }
                        map("world",add=T, lwd=1);
                        #map("state",add=T,lwd=1);
                        #map("lakes",add=T,lwd=1);
                      } 
                      
      )
      
      
      
   }
    dev.off()

    
    pdf(paste("C:/Users/ekc76/Box Sync/MW_seasonal_warm_season_precip_forecast/images/climate_variables/climate_comp_maps/lag", lagv, "/ageostrophic/ageostrophic_abs_z_anom_col_",j,"_month",i,"_lag",lagv,".pdf",sep=""),
        width=14, height=5)
    par(mfrow=c(1,3)) 
     
    filled.contour3(lon_final,lat_final,t(z850_anom),
                    xlim=c(lon.min,lon.max),
                    ylim=c(lat.min,lat.max),
                    level=mylevel,
                    col=mycol,
                    #frame.plot=FALSE,
                    key.title=title(main=paste(j, k, sep=" "),cex=2),
                    
                    plot.axes={               
                      axis(1, labels=FALSE, tick=FALSE);
                      axis(2, labels=FALSE, tick=FALSE);
                      arrow.plot(a1=lon_all[uv],a2=lat_all[uv],u=uu[uv],v=vv[uv],arrow.ex=0.17,length=0.2,xpd=FALSE,
                                 xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=2, col="gray40");
                      arrow.plot(a1=lon_all[uvall],a2=lat_all[uvall],u=uu[uvall],v=vv[uvall],arrow.ex=0.17,length=0.2,xpd=FALSE, 
                                 xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=2, col="gray75");
                      if(i<=5){ 
                        contour(x=lon_final, y=lat_final, z=t(z300_sub), levels=c(9060),col="green",
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3,lty=1, add=TRUE);
                        contour(x=lon_final, y=lat_final, z=t(ltmz300_sub), levels=9060,
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=0.8,lty=1, col="black",add=TRUE);
                        contour(x=lon_final, y=lat_final, z=t(z850_sub), levels=1540,col="green",
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3,lty=1, add=TRUE);
                        contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1540,
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=0.8,lty=1, col="black",add=TRUE);
                        contour(x=lon_final, y=lat_final, z=t(z850_sub),lwd=1,lty=1, levels=c(1480,1500,1520,1540),
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), col="green",add=TRUE);
                        
                        
                        
                      }else{
                        contour(x=lon_final, y=lat_final, z=t(z300_sub), levels=c(9100),col="green",
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3,lty=1, add=TRUE);
                        contour(x=lon_final, y=lat_final, z=t(ltmz300_sub), levels=9100,
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=0.8,lty=1, col="black",add=TRUE);
                        contour(x=lon_final, y=lat_final, z=t(z850_sub), levels=1560,col="green",
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=3,lty=1, add=TRUE);
                        contour(x=lon_final, y=lat_final, z=t(ltmz850_sub), levels=1560,
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max), lwd=0.8,lty=1, col="green",add=TRUE);
                        contour(x=lon_final, y=lat_final, z=t(z850_sub),lwd=1,lty=1, levels=c(1480,1500,1520,1540),
                                xlim=c(lon.min,lon.max),ylim=c(lat.min,lat.max),  col="green",add=TRUE);
                        
                        
                      }
                      map("world",add=T, lwd=1);
                      #map("state",add=T,lwd=1);
                      #map("lakes",add=T,lwd=1);
                    } 
                    
    )
    
    
    
  }
  dev.off()
  
    #Map 2: JJA omega (color) isoltherms (contour) agesotrophic wind absolute (dark color when above average)
    
    #Map 3: JJA MFD (color) q (contour) IVT (vector)
    
  }   
}
}


#anomalies maps
