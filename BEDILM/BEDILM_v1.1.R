

####################################
# DEPTH-INTEGRATED LAKE METABOLISM #
####################################
#
# Code for calculating daily rates of GPP, ER and NEP from high-frequency and depth-specific measurements of dissolved oxygen, temperature and light
#
# Code for Bayesian implementation written by Darren Giling (darren.giling@idiv.de)
#
# Physical fluxes and metabolic model following Staehr et al. 2012 L&O 57: 1317 and references therein
#
# Bayesian model and metabolic sub-models modified from Grace et al. 2015 L&O:M 13: 103 and Song et al. 2016 L&O:M doi: 10.1002/lom3.10112
#


### Please see the guide file accompanying this code ###


# packages to install
library(R2jags)
library(zoo)

#========================#
#   1. DEFINE INPUTS     # 
#========================#

### generally only need to update this section ###

# System
    dir_input <- "C:/example/BEDILM" # set 'dir_input' to folder with this code
    dir_input <- "C:/dg45koti/Dropbox/Code/Depth-integrated lake metabolism_v3/BEDILM" # set 'dir_input' to folder with this code
    
    data.file <- "Stechlin_example_2014_DOY_245-246.csv"
    
    data      <- read.csv(file.path(dir_input,"data",data.file)) # import from 'data' folder
    name      <- unique(data$lake) # can also provide in character

# Define the column numbers of required variables 
    wind.col     <- 5       # wind speed (m/s)
    pressure.col <- 6       # barometric pressure (mB)
    sal.col      <- 8       # salinity (ppt) - optional, provide a column of zeros if not available
    DO.cols      <- 9:17    # dissolved oxygen (mg/L or % sat - see below)
    TempC.cols   <- 18:26   # water temperature (deg C)
    PAR.cols     <- 27:35   # PAR (umol m^-2 s^-1)
    PAR0.col     <- 7       # PAR at surface (umol m^-2 s^-1)

# Define units of dissolved oxygen - set to 'FALSE' if measured in % saturation
    DO.in.mgL <- T

# Define whether to calculate Ds with measured or modelled DO.
    calc.Ds.with.modDO <- T

# Set profile measurement settings
    timestep        <- 1            # must be uniform (h)
    smooth.period   <- 4            # time period for smoothing (h)
    wind.height     <- 2            # height of wind sensor above lake surface (m)
    meas.depths     <- seq(1,17,2)  # measurement depths (m)
    layer.thickness <- rep(2, 9)    # thickness of each layer (m)

# Define physcial properties of the lake
    Zmax  <- 20  # maximum depth of profile (m): determines N2 for deepest layer
    Az    <- c(4054500,3729800,3469400,3320100,3176200,3022000,2863200,2716800,2565000) # surface area (m^2) of lake at depth z from hypsographic data
    p.gradient.threshold <- 0.08 # density gradient threshold for defining Zmix (kg m^-3 m^-1) (Sadro et al 2011 evaluated values between 0.005 and 0.10)

# Set the method to calculate the bottom of the metalimnion
    meta.symmetrical <- T

# Set number of iterations 
    n.temp.iter  <- 1000  # number of iterations for temperature curve model
    n.metab.iter <- 2000  # number of interations for metabolism model
    
    
#========================#
#        2. SETUP        # 
#========================#

# other variables
    num.layers <- length(meas.depths)
    sec.per.day<-86400


### Smoothing ###
# smoothing is applied across whole dataset, not within days, so this option requires continuous data (i.e. no missing days)
    
    smooth.rows <- (smooth.period/timestep)+1
    smooth <- function(x) (rollapply(x, smooth.rows, mean, na.rm=T, align="center"))
    
    # smoothing applied to DO, PAR and wind (following Obrador et al 2014):
    for (j in c(DO.cols, PAR0.col, PAR.cols, wind.col)) {
      data[((smooth.rows/2)+0.5):(nrow(data)-((smooth.rows/2)-0.5)),j] <- smooth(data[,j])
    }    

    data[,c(PAR.cols,PAR0.col)] <-  replace(data[,c(PAR.cols,PAR0.col)], data[,c(PAR.cols,PAR0.col)]<1,0) # set detection limit for PAR if wanted
    
### Define list of days (index k) ####
  
    k.list <- numeric()
    for (k in min(data$DOY):max(data$DOY)) { # check for full days with no NA values, add to k list  
          data.subset <- data[data$DOY==k,]
          check1<-if (nrow(data.subset)==24/timestep) T else F
          check2<-ifelse(is.na(min(data.subset[,c(5,6,7,DO.cols[1:length(meas.depths)],TempC.cols[1:length(meas.depths)]),])) , F , T) 
          check=all(check1,check2)
          ifelse(check==T,k.list<-c(k.list,k),k)  
    }


    
### Set up output table data frame  ###
    
    output.table<-NULL
    output.table<-data.frame(Lake=rep(name,length(k.list)*num.layers),DOY=rep(k.list, each=num.layers),layer.no=rep(1:num.layers, length(k.list)),depth=rep(meas.depths, length(k.list)),
                             layer.thickness=NA,layer.zone.mean=NA,TempC.mean=NA, TempC.sd=NA, PAR.max=NA, PAR.mean=NA, Zmix.mean=NA, metalim.max.mean=NA, thermocline.depth.mean=NA, mean.KD=NA, mean.Zeu=NA,  
                             Sc=NA, U10=NA, k600=NA, Ks=NA, Ds.sum=NA, N2=NA, Kv=NA, Dv.sum=NA, Dz.sum=NA, 
                             PPP=NA, R2=NA, pD=NA, DIC=NA, Rhat.test=NA, A.Rhat=NA, R.Rhat=NA, p.Rhat=NA, theta.Rhat=NA, A.auto.corr=NA, R.auto.corr=NA, p.auto.corr=NA, theta.auto.corr=NA,
                             p.mean=NA, p.sd=NA, A.mean=NA, A.sd=NA, R.mean=NA, R.sd=NA, theta.mean=NA, theta.sd=NA, ER.mean=NA, ER.sd=NA, GPP.mean=NA, GPP.sd=NA, NEP.mean=NA, NEP.sd=NA)

    
#========================#
#   3. LOOP OVER DAYS    # 
#========================#

  for (k in k.list)  {
      # k <- 246   # just for testing
      ptm <- proc.time()
      data.subset <- data[data$DOY==k,]  
      
      # define daily variables matrices
      DO.meas.matrix <- data.subset[,DO.cols]
      DO.meas.sat.matrix <- data.subset[,DO.cols]
      TempC.matrix <- data.subset[,TempC.cols]
      PAR.matrix <- data.subset[,PAR.cols]
      Sal.matrix <- data.subset[,sal.col]
      num.measurements <- nrow(data.subset)
      wind <- data.subset[, wind.col]  
      atmo.pressure <- data.subset[, pressure.col]
      
      is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
      data.subset[is.nan(data.subset)] <- NA
      
      # plotting
      #DO.plot<-data.matrix(DO.meas.matrix[1:nrow(DO.meas.matrix), rev(1:ncol(DO.meas.matrix))])  
      #x<-1:num.measurements
      #y<-rev(meas.depths)*(-1)
      #filled.contour(x, y, z=DO.plot, main=paste(c("DO meas","Day",k), sep=""), xlab="Timestep", ylab="Depth (m)") 
      
      #TempC.plot<-data.matrix(TempC.matrix[1:nrow(TempC.matrix), rev(1:ncol(TempC.matrix))])  
      #x<-1:num.measurements
      #y<-rev(meas.depths)*(-1)
      #filled.contour(x, y, z=TempC.plot, main=paste(c("TempC","Day",k), sep=""), xlab="Timestep", ylab="Depth (m)") 
      

      #-----------------------------------------------------------------------------------------------------
      #   3.1 Calculate mean daily light extension coefficient (KD; m^-1) from measured PARz if available
      #-----------------------------------------------------------------------------------------------------
      
      light.df <- data.frame(KD = rep(NA, num.measurements),
                             max.e = rep(NA, num.measurements),
                             e.r2 = rep(NA, num.measurements))
      e.measurements <- (which( is.na(data.subset[, PAR.cols[1]]) == FALSE))
      
      for (i in e.measurements) 
          {
          e <- as.numeric( t( data.subset[i, PAR.cols] ) )
          non0 <- which(e!=0)
          #plot(log(e[non0]+1)~meas.depths[non0])
          light.df$max.e[i] <- max(e, na.rm=T)
          
          if(length(non0)>1) 
            {
            light.model <- lm(log(e[non0]+1)~meas.depths[non0])
            light.df$KD[i] <- as.numeric(coefficients(light.model)[2])
            
            options(warn = -1)
            light.df$e.r2[i] <- cor(log(e[non0]+1), light.model$fitted.values)**2 
            options(warn = 0)
            }
          }
          
      # Mean of all profiles during day with log(PAR)~z r2>0.8 (following Obrador et al 2014)
      light.df.sub <- light.df[light.df$max.e>5,] 
      light.df.sub <- light.df.sub[light.df.sub$e.r2>0.8,]
      mean.KD <- mean(light.df.sub$KD, na.rm=T)
      mean.Zeu <- -4.6/mean.KD
          
    
      #---------------------------------------------------------
      # 3.2 Model temperature curve to define epi, meta and hypo
      #---------------------------------------------------------
     
      # temperature model: Rimmer et al. 2005
      
      temps <- NA; roc.t <- NA;roc.p <- NA; Zmix<-NA; p<-NA; threshold<-NA; thermocline<-NA; thermocline.depth<-NA; metalim.max<-NA
      roc.t.data <- NA; roc.p.data <- NA; p.data<- NA
      
      non.na.t <- which( is.na(data.subset[, TempC.cols[1]]) == FALSE)
      mod.temps.save <- matrix(NA, nrow=length(seq(0,Zmax,0.1)),ncol=length(non.na.t))
      
      for (i in non.na.t) {
          temps <- as.numeric(t( data.subset[i, TempC.cols] ))
          max.z <- length(temps)
          Tz <- as.numeric(temps)[1:max.z]
          Th <- min(as.numeric(temps), na.rm=T)
          Te <- max(as.numeric(temps), na.rm=T)
          z <- meas.depths[1:max.z]
          n.meas <- length(z)
          
          data.list <- list("z", "Th", "Te", "n.meas", "Tz")
          params=c("aa", "n")
          model.file <- 'temperature_curve_jags_code.txt'
            
          n.chains <- 3;  n.thin <- 3;  n.iter <- n.temp.iter*n.thin; n.burnin <- n.iter*0.5  # increase n.iter here
          
          temp.curve = NULL
          temp.curve <- do.call(jags.parallel,
                           list(data=data.list, inits=NULL, parameters.to.save=params, model.file = model.file,
                                n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                                n.thin = n.thin, n.cluster= n.chains, DIC = TRUE,
                                working.directory = dir_input, jags.seed = 123, digits=5))
          #traceplot(temp.curve)
          
          aa <- temp.curve$BUGSoutput$summary["aa","mean"]
          n <- temp.curve$BUGSoutput$summary["n","mean"]
          
          mod.depths <- seq(0,Zmax,0.1)
          mod.temps <- rep(NA,length(mod.depths))
             for (tt in 1:length(mod.temps)){
              mod.temps[tt] <- Th + (Te-Th) * (1/ (1+ ((aa*mod.depths[tt])^n))^(1-(1/n))) }
          
          mod.temps.save[,i]<-mod.temps
          
          p<-(1-((6.63e-06)*((mod.temps-4)^2)))*1000  # (kg m^-3)
          
          for (tt in 2:(length(mod.temps))) {       
            roc.p[tt] = (p[tt] - p[tt-1]) / (mod.depths[tt]-mod.depths[tt-1])
            roc.t[tt] = (mod.temps[tt] - mod.temps[tt-1]) / (mod.depths[tt]-mod.depths[tt-1])
          }
          
          thermocline<-which(roc.t==min(roc.t, na.rm=T))[1] 
          thermocline.depth[i]<-mod.depths[thermocline]
          
          threshold<-min(which(roc.p>p.gradient.threshold))
          if (abs(threshold)==Inf) Zmix[i] <- mod.depths[length(mod.depths)] else  Zmix[i] <- mod.depths[threshold]
          
          if (meta.symmetrical==T) {
              metalim.max[i]<- thermocline.depth[i] + (thermocline.depth[i] - Zmix[i])
              metalim.max[i] <- ifelse(metalim.max[i]>Zmax , Zmax, metalim.max[i])
          } else {
              threshold.max<-max(which(roc.p>p.gradient.threshold))
              if (abs(threshold.max)==Inf) metalim.max[i]<- mod.depths[length(mod.depths)]  else metalim.max[i]<- mod.depths[threshold.max]  
          }
        
      }
        
      mean.Zmix = mean(Zmix,na.rm=T)
      mean.meta.max = mean(metalim.max, na.rm=T)
      mean.thermocline = mean(thermocline.depth, na.rm=T)
      
      # surface and Zmax temps for N2 calcs - daily means
      temp.0 <- mod.temps.save[1,]
      temp.Zmax <- mod.temps.save[length(mod.temps),]
      
      ### write plot of measured and modelled temperature ###
        max.y <- Zmax+2
        x <- -meas.depths
        mean.mod.temps<-apply(mod.temps.save,1,FUN=mean)
        mean.meas.temps<-apply(TempC.matrix,2,FUN=mean)
        xlim=c(min(mean.meas.temps)-3,max(mean.meas.temps)+3)
        
        jpeg(file=file.path(dir_input,"/results/profile plots",paste(name,"_day_",k,".jpg", sep="")), width=500, height=500, pointsize=16, quality=200)
            plot(x ~ mean.meas.temps, las=1, ylim=c(-max.y,0), xlim=xlim,
                 ylab='Depth (m)', xlab= expression(paste("Water temp. (",degree,"C)")),yaxt='n', typ='n')
            axis(2, at=seq(0,-max.y,-4), labels=seq(0, max.y, 4), las=1, cex=0.8)
            title(paste0(name, " (day ", k,")"), line = 1, cex.main=1, font.main=1)
            rect(xlim[1], -mean.Zmix, xlim[2], -mean.meta.max, col='grey85', border=NA )
            segments(xlim[1],-mean.Zeu, xlim[2],  -mean.Zeu, lwd=2, lty=2)
            segments(xlim[1],-mean.thermocline, xlim[2],  -mean.thermocline, lwd=1)
            legend("topleft", legend=c("metalimnion", "mean thermocline", "mean Zeu", "mean meas. temp", "mean mod. temp"),
                   lwd=c(8,1,2,NA,2), col=c("grey","black", "black", "blue", "blue"), pt.lwd=2, lty=c(1,1,2,NA,1), pch=c(NA,NA,NA,21,NA), pt.cex=c(1,1,1,1.3,1), bty='n', cex=0.8)
            points(x ~ mean.meas.temps, pch=21, cex=1.3, col="blue")
            lines(-mod.depths~mean.mod.temps, lwd=2, col="blue")
        dev.off()
    
      # surface and Zmax temps for N2 calcs - daily means version
      #temp.0 <- mean.mod.temps[1]
      #temp.Zmax <- mean.mod.temps[length(mod.temps)]
      
      # assign depth zones (1 = epi; 2= meta; 3= hypo) - based on mean daily stratification
      zone<- matrix(NA, nrow = num.measurements, ncol= num.layers)
      
      for (i in 1:num.measurements) {
              for (j in 1:num.layers) {
                      ifelse(meas.depths[j]<Zmax, zone[i,j]<-1, 0)   
                      ifelse(meas.depths[j]>=Zmix[i], zone[i,j]<-2, 0)
                      ifelse(meas.depths[j]>=metalim.max[i], zone[i,j]<-3, 0)  
              }
      }
      mean.daily.zone <- round(apply(zone,2,FUN=mean))
      
      
      #--------------------------------------------------------------------
      # 3.3 Calculate DO mg/L from % sat or mg/L at saturation as required
      #--------------------------------------------------------------------
    
      # set up objects  
      kelvin<- matrix(NA, nrow = num.measurements, ncol= num.layers)
      DO.sat.pre<- matrix(NA, nrow = num.measurements, ncol= num.layers)
      DO.sat.matrix <- matrix(NA, nrow = num.measurements, ncol= num.layers)
      C <- matrix(NA, nrow = num.measurements, ncol= num.layers)
      deltaDO <- matrix(NA, nrow = num.measurements, ncol= num.layers)
      
      kelvin <- 273.15 + TempC.matrix
      
      # convert % sat to mg L-1 if required
      if(DO.in.mgL == F) { 
          inv.kelvin <- 1/kelvin
          p1 <- -862194900000*(inv.kelvin^4)+12438000000*(inv.kelvin^3)-66423080*(inv.kelvin^2)+157570.1*inv.kelvin-139.344
          p2 <- 2140.7*(inv.kelvin^2)-10.754*inv.kelvin+0.017674
          per.sat <- 0.01*exp(p1-Sal.matrix*p2)
          
          DO.meas.matrix <- DO.meas.sat.matrix * per.sat
          DO.sat.matrix <- 100 * per.sat
      } else {
          # calculate mg L-1 at saturation
          C <- (-173.4292 
                 + 249.6339 * (100 / kelvin) 
                 + 143.3483 * log(kelvin/ 100)  
                 -  21.8492 * (kelvin / 100) 
                 + Sal.matrix * ((-0.033096 + 0.014259 * kelvin/ 100) - 0.0017000 * (kelvin / 100)^2)  )
        
          DO.sat.pre <- exp(C) * 1.423
          DO.sat.matrix <- DO.sat.pre * ((atmo.pressure*0.0987 - 0.0112)/100) 
      }
      
  
      #--------------------------------------------------------
      # 3.4 Calculate rates of diffusive DO fluxes (Ds, Dv, Dz)
      #--------------------------------------------------------
      
          #---------------------------------------------
          # 3.4.1 Rate of exchange with atmosphere (Ds)
          #---------------------------------------------
          
          # epilimnion only (negative is fluxes into and positive is fluxes out of)
      
          # create objects and variables
          Sc <- matrix(NA, nrow = num.measurements, ncol= num.layers)
          U10 <- NA
          k600 <- NA
          Ks <- matrix(NA, nrow = num.measurements, ncol= num.layers)
          Ds.matrix <- matrix(NA, nrow = num.measurements, ncol= num.layers)
          alpha = 1.4125* (wind.height^(-0.15))
          
          # calculate Ds for each measurement time and layer
          for (i in 1:num.measurements) {
            for (j in 1:num.layers) {    
        
                  # schmidt coefficient
                  Sc[i,j] = -0.0476*(TempC.matrix[i,j]^3) + 3.7818*(TempC.matrix[i,j]^2) - 120.1*TempC.matrix[i,j] + 1800.6
      
                  # wind speed at 10m (m s-1)
                  U10[i] = wind[i] * alpha # Smith (1985)
      
                  # piston velocity (m h-1), following Cole and Caraco (1998)
                  # Other formulations have been tested with this model: see Equations 6a-d in Staehr et al. 2012 (and associated references)
                  k600[i] <- (2.07 + 0.215*(U10[i]^1.7))/100  
                  Ks[i,j] <- k600[i]*( (Sc[i,j]/600)^(-0.5)) # Jähne et al. (1987)
      
                  # Ds flux (mg L^-1 m^-3 timestep^-1)
                  Ds.matrix[i,j] <- (Ks[i,j]*(DO.meas.matrix[i,j] - DO.sat.matrix[i,j]))/Zmix[i] 
      
                  # only applied to epilimnion layers
                  Ds.matrix[i,j] <- ifelse(mean.daily.zone[j]==1, Ds.matrix[i,j], 0)
             }
          }
          
          Ds.matrix <- Ds.matrix*timestep  # (mg L^-1 m^-3 h^-1)
          
          # plotting
          #Ds.plot<-data.matrix(Ds.matrix[1:nrow(Ds.matrix), rev(1:ncol(Ds.matrix))])  
          #x<-1:num.measurements
          #y<-rev(meas.depths)*(-1)
          #filled.contour(x, y, z=Ds.plot, main=paste(c("Ds - flux to atmosphere","Day",k), sep=""), xlab="Timestep", ylab="Depth (m)")
          
          
          #-------------------------------------------------------------------
          # 3.4.2 Vertical transfer among layers due to eddy diffusivity (Dv)
          #-------------------------------------------------------------------
         
          # set up variables and empty matrices
          g <- 9.8 # m/s2
          p.matrix <- matrix(0, nrow = num.measurements, ncol= num.layers)
          N2 <- matrix(0, nrow = num.measurements, ncol= num.layers)
          Kv <- matrix(0, nrow = num.measurements, ncol= num.layers)
          Dv.matrix <- matrix(0, nrow = num.measurements, ncol= num.layers)
          p.Zmax <- NA
            
          # calculations for each timepoint 
          for (i in 1:num.measurements) {
                  
              # water density
              for (j in 1:num.layers) {     
                p.matrix[i,j]<-(1-((6.63e-06)*((TempC.matrix[i,j]-4)^2)))*1000  # water density (kg m^-3)
              }
              p.Zmax[i] <- (1-((6.63e-06)*((temp.Zmax[i]-4)^2)))*1000  # water density (kg m^-3)
                   
              # Brunt Väisälä bouyancy frequency (s^-2)  
              for (j in 1:(num.layers-1)) {
                N2[i,j] <-  - (  (g/p.matrix[i,j]) * ((p.matrix[i,j] - p.matrix[i,j+1])/(layer.thickness[j])) )
              }
              # N2 of deepest layer calculated with modelled density at Zmax
              N2[i,num.layers] <-  - ( (g/p.matrix[i,num.layers]) * ((p.matrix[i,num.layers] - p.Zmax[i])/(layer.thickness[num.layers])) )
              
              # vertical turbulent diffusivity
              for (j in 1:num.layers) {
                Kv[i,j] <- 2.941e-04 * ((Az[j] * 1e-06)^0.56)  * (N2[i,j]^(-0.43)) # vertical turbulent diffusivity; Hondzo and Stefan 1993
              }
           }
         
          max.finite <- function(X) max(X[is.finite(X)], na.rm=T)  # when there is no density gradient (or negative density gradient) between layers in surface waters - v1.1 option: replace with max
          Kv[,1:(num.layers-1)]<-sweep(Kv[,1:(num.layers-1)], MARGIN = 1, 
                STATS = apply(Kv[,1:(num.layers-1)], 1, max.finite),
                FUN =  function(x,s) ifelse(!is.finite(x), s, x)
          )
          Kv[ Kv=="NaN" | Kv=="Inf"] <- 0 # when there is negative density gradient between layers in last layer. This occurs when the modelled temp.Zmax is warmer than the deepest layer.
          Kv[Kv>1]<-1 # limit     
         
          # Dv rate (mg L^-1 m^-3 timestep^-1)
          for (i in 1:num.measurements) {  
            Dv.matrix[i,1] <- (Kv[i,1+1] * (DO.meas.matrix[i,1]-DO.meas.matrix[i,1+1]) / layer.thickness[1] ) * (Az[1] / (Az[1]*layer.thickness[1]))  # top layer
            Dv.matrix[i,num.layers] <- ((Kv[i,num.layers] * (DO.meas.matrix[i,num.layers]-DO.meas.matrix[i,num.layers-1]))  / layer.thickness[num.layers] ) * (Az[num.layers] / (Az[num.layers]*layer.thickness[num.layers])) # bottom layer
            
            for (j in 2:(num.layers-1)) {
                    Dv.matrix[i,j] <- (((Kv[i,j] * (DO.meas.matrix[i,j]-DO.meas.matrix[i,j-1])) + (Kv[i,j+1] * (DO.meas.matrix[i,j]-DO.meas.matrix[i,j+1]))) / layer.thickness[j] ) * (Az[j] / (Az[j]*layer.thickness[j]))  
            }
          }
          
          Dv.matrix <- Dv.matrix*timestep  # (mg L^-1 m^-3 h^-1)
          
          # plotting
          #N2.plot<-data.matrix(N2[1:nrow(N2), rev(1:ncol(N2))])  
          #x<-1:num.measurements
          #y<-rev(meas.depths)*(-1)
          #filled.contour(x, y, z=N2.plot, main=paste(c("N2 - buoyancy frequency","Day",k), sep=""), xlab="Timestep", ylab="Depth (m)")
        
          #Dv.plot<-data.matrix(Dv.matrix[1:nrow(Dv.matrix), rev(1:ncol(Dv.matrix))])  
          #x<-1:num.measurements
          #y<-rev(meas.depths)*(-1)
          #filled.contour(x, y, z=Dv.plot, main=paste(c("Dv - vertical transfer among layers","Day",k), sep=""), xlab="Timestep", ylab="Depth (m)") 
        
            
          #------------------------------------------------------------------
          # 3.4.3 Fluxes due to mixed-layer deepening (Dz)
          #------------------------------------------------------------------
          
          #  Following Bell et al. 2006,  Staehr et al. 2012
            
          Dz.matrix <- matrix(NA, nrow = num.measurements, ncol= num.layers)
          dZmix.dt <- 0; for (i in 2:num.measurements) { dZmix.dt[i]<- ((Zmix[i]-Zmix[i-1]) / timestep) } # Zmix deepening rate (m/s)
          
          # optional: set a maximum deepening rate to address issue of unrealisitic large Dz when there is short-term surface microstratification
          dZmix.dt[which(dZmix.dt>5)] <-0     # 
          dZmix.dt[which(dZmix.dt<(-5))] <-0
          
          # calculate Dz  (mg L^-1 m^-3 timestep-1)
          for (i in 2:num.measurements) 
              {    
              meta.plus <- which(meas.depths > (mean.Zmix-1) & meas.depths < (mean.meta.max+1)) # define layers (metalimnion and +/- 1 m above and below)
              
                      for (j in meta.plus)
                          {  
                          Dz.matrix[i,j] <- dZmix.dt[i] * ( (DO.meas.matrix[i-1,j]-DO.meas.matrix[i,j]) / layer.thickness[j] )
                          }
              }
          
          Dz.matrix[is.na(Dz.matrix)]<-0
          
          Dz.matrix <- Dz.matrix*timestep  # (mg L^-1 m^-3 h^-1)
          
          # plotting Dz
          #Dz.plot<-data.matrix(Dz.matrix[1:nrow(Dz.matrix), rev(1:ncol(Dz.matrix))])  
          #x<-1:num.measurements
          #y<-rev(meas.depths)*(-1)
          #filled.contour(x, y, z=Dz.plot, main=paste(c("Dz - mixed layer deepening flux","Day",k), sep=""), xlab="Timestep", ylab="Depth (m)") 
       
             
      #------------------------------------------------------------
      # 3.5 Model metabolism for each layer (index j) sequentially 
      #------------------------------------------------------------
    
      for (j in 1:num.layers) {
        
        # j <- 7  # just for testing
        
        ### select vectors for appropriate layer ###
          layer.mean <-  round(mean(zone[,j]),0)
          if(layer.mean==1) { epi.dummy <- 1 } else {epi.dummy=0}
          
          TempC <- TempC.matrix[,j]
          DO.sat <- DO.sat.matrix[,j]
          DO.meas <- DO.meas.matrix[,j]
          
          PAR <- PAR.matrix[,j] # to use directly measured or pre-calculated PARz
          
          Ds <- Ds.matrix[,j]  # passed to JAGS for DO measured version of Ds
          Ks.vect <- Ks[,j]    # passed to JAGS for DO modelled version of Ds
          Dv <- Dv.matrix[,j]
          Dz <- Dz.matrix[,j]
       
        ### package data for JAGS ###
          data.list <- list("num.measurements", "timestep", "DO.meas", "TempC", "PAR", "Ds", "Dv", "Dz", "Ks.vect", "DO.sat", "Zmix", "epi.dummy")
          params=c("A","R","p","theta","sd","ER","GPP","NEP","PPfit","DO.modelled")
          
        ### mcmc settings ###
          inits <- function() { list(sd=0.1) } # initial values (others are selected automatically)
          n.chains <- 3
          n.thin <- 100 
          n.iter <- n.metab.iter*n.thin          
          n.burnin <- n.iter*0.5
        
        ### model file ###
          model.file <- 'metab_model_jags_code.txt'
          # Contains the metablism model. The default uses mean temperature for ER sub-model (can be altered to 20 deg C, e.g.)
          # Prior distributions for GPP sub-model parameter p and temperature-dependence theta can be adjusted in the model file
          
        ### call JAGS ###
        metab = NULL
        metab <- do.call(jags.parallel,
                      list(data=data.list, inits=inits, parameters.to.save=params, model.file = model.file,
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                           n.thin = n.thin, n.cluster= n.chains, DIC = TRUE,
                           working.directory = dir_input, jags.seed = 123, digits=5))  
        
        ### inspect results ###
          # print(metab$BUGSoutput$summary, digits=3)                   
          # traceplot(metab)   # plot chains
        
        ### extract modelled DO ###
          DO.mod.means <- metab$BUGSoutput$summary[2:(num.measurements+1),"mean"]
          DO.mod.sd <- metab$BUGSoutput$summary[2:(num.measurements+1),"sd"]
        
        ### diagnostic multi-plot ###
          jpeg(file=file.path(dir_input,"/results/validation plots",paste(name,"_day",k,"_layer",j,".jpg", sep="")), width=1200, height=1200, pointsize=30, quality = 300)
            traceplot(metab, varname=c('A','p','R','theta'), ask=FALSE, mfrow=c(3,3), mar=c(2,2,0,8), new=FALSE)
            ymin=min(c(DO.meas,DO.mod.means))-0.1; ymax= max(c(DO.meas,DO.mod.means))+0.1
            plot(1:num.measurements,DO.mod.means, type="l", ylim=c(ymin,ymax),xlab="Timestep", lwd=2)
            points(1:num.measurements,DO.meas,pch=1,xlab="Timestep", col="grey60")  
            plot(1:num.measurements,TempC,pch=1,xlab="Timestep" , typ='p')
            plot(1:num.measurements,PAR,pch=1,xlab="Timestep" , typ='p')
          graphics.off()
        
        ### calculate chain convergence statistics ###
          srf<- metab$BUGSoutput$summary[,8] # select scale reduction factors (SRF) of each variable
          Rhat.test <- NULL
          Rhat.test <- ifelse(any(srf>1.1, na.rm=T)==TRUE,"Check convergence", "Convergence OK") # test if any are > 1.1
          metab.mcmc<-as.mcmc(metab)   # convert results to mcmc object
          ac.lag1 <- autocorr.diag(metab.mcmc, lags = 1) # call autocorrelation 
        
        ### add results to output.table ###
          row <- which(output.table$DOY==k & output.table$layer.no == j) 
          
          output.table$layer.thickness[row] <- layer.thickness[j]
          output.table$layer.zone.mean[row] <- layer.mean
          output.table$TempC.mean[row]<- mean(TempC)
          output.table$TempC.sd[row]<- sd(TempC)
          output.table$PAR.max[row]<- max(PAR) 
          output.table$PAR.mean[row]<- mean(PAR)
          output.table$Zmix.mean[row]<-mean.Zmix
          output.table$metalim.max.mean[row]<-mean.meta.max
          output.table$thermocline.depth.mean[row]<-mean.thermocline
          output.table$mean.KD[row]<-mean.KD
          output.table$mean.Zeu[row]<-mean.Zeu
          output.table$Sc[row] <- mean(Sc[,j])
          output.table$U10[row]<- mean(U10)
          output.table$k600[row]<- mean(k600)
          output.table$Ks[row]<- mean(Ks[,j])
          if (calc.Ds.with.modDO==T) {
            output.table$Ds.sum[row]<-sum(Ks.vect*(DO.mod.means - DO.sat)/Zmix ) *epi.dummy * timestep
          } else {
            output.table$Ds.sum[row]<- sum(Ds.matrix[,j]*timestep)
          }
          output.table$N2[row]<- mean(N2[,j])
          output.table$Kv[row]<- mean(Kv[,j])
          output.table$Dv.sum[row]<- sum(Dv.matrix[,j]*timestep)
          output.table$Dz.sum[row]<- sum(Dz.matrix[,j]*timestep)
          output.table$PPP[row] <- metab$BUGSoutput$summary["PPfit","mean"]
          output.table$R2[row] <-  cor(DO.mod.means,DO.meas)^2
          output.table$pD[row]<-metab$BUGSoutput$pD
          output.table$DIC[row]<- metab$BUGSoutput$DIC
          output.table$Rhat.test[row] <- Rhat.test
          output.table$A.Rhat[row]<- as.numeric(srf["A"])
          output.table$R.Rhat[row]<- as.numeric(srf["R"])
          output.table$p.Rhat[row]<- as.numeric(srf["p"])
          output.table$theta.Rhat[row]<- as.numeric(srf["theta"])
          output.table$A.auto.corr[row]<- ac.lag1[1,"A"]
          output.table$R.auto.corr[row]<- ac.lag1[1,"R"]
          output.table$p.auto.corr[row]<- ac.lag1[1,"p"]
          output.table$theta.auto.corr[row]<- ac.lag1[1,"theta"]
          output.table$p.mean[row] <- metab$BUGSoutput$mean$p
          output.table$p.sd[row]<- metab$BUGSoutput$sd$p
          output.table$A.mean[row]<- metab$BUGSoutput$mean$A
          output.table$A.sd[row]<- metab$BUGSoutput$sd$A
          output.table$R.mean[row]<- metab$BUGSoutput$mean$R
          output.table$R.sd[row]<- metab$BUGSoutput$sd$R
          output.table$theta.mean[row]<- metab$BUGSoutput$mean$theta
          output.table$theta.sd[row]<- metab$BUGSoutput$sd$theta
          output.table$ER.mean[row]<- metab$BUGSoutput$mean$ER
          output.table$ER.sd[row]<- metab$BUGSoutput$sd$ER
          output.table$GPP.mean[row]<- metab$BUGSoutput$mean$GPP
          output.table$GPP.sd[row]<- metab$BUGSoutput$sd$GPP
          output.table$NEP.mean[row]<- metab$BUGSoutput$mean$NEP
          output.table$NEP.sd[row]<- metab$BUGSoutput$sd$NEP
          
         ### overwrites output table to update for each layer ###
          write.csv(output.table, file=file.path(dir_input,"/results",paste(name, "_metab_model_results.csv", sep="")))  # output file
                    
     }  # end of layer loop
      
  } # end of day loop
    
  time<-proc.time() - ptm; time[3]/60; paste('mins')





