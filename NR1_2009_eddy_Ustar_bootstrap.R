
# This script:
# - reads a REddyProc input file with data for only 1 year and one site(the site must be included in "sites_df", currently only 7 sites)
# - gets if the time step from the data is hourly or half-hourly
# - estimates different Ustar thresholds based on Papale et al. 2006 (Ustar on original series, 5% bootstrap + median bootstrap + 95% bootstrap which gives 90% confidence interval)
# - does the NEE partition into GPP and Reco for all the different Ustar thresholds.
# - calculates aggregated GPP, NEE and Reco in gC/m2 for the whole year for all the different Ustar thresholds.
# Author: Francesc Montane
# Contact: fmontane@email.arizona.edu  


# IMPORTANT: Modify lines "setwd", "file_name", "Dir.s", and replace year 2009 with correct year.


setwd("D:/Sites_DOE/AmeriFlux/Niwot Ridge/L2_with_gaps/v0008/2009")


# remove all elements from workspace

rm(list=ls())


library(REddyProc)

file_name=paste("NR1_2009_2009_eddy.txt")
site=substr(file_name,1,3)
site_code<-paste("US-",site,sep="")
yearini=substr(file_name,5,8)
yearend=substr(file_name,5,8)

# create a data frame with site information (7 sites for now, include more sites if necessary)
sites_df<-data.frame(site=character(7),siteLat=numeric(7),siteLon=numeric(7),siteTime=numeric(7))
sites_df$site<-c("NR1","MMS","Ha1","Ho1","UMB","Vcm","Vcp")
sites_df$siteLat<-c(40.0329,39.3231,42.5378,45.2041,45.5598,35.8884,35.8624)
sites_df$siteLon<-c(-105.5464,-86.4131,-72.1715,-68.7402,-84.7138,-106.5321,-106.5974)
sites_df$siteTime<-c(-7,-5,-5,-5,-5,-7,-7)

mysite<-subset(sites_df[(sites_df$site==site),])

siteLat<-mysite$siteLat
siteLon<-mysite$siteLon
siteTime<-mysite$siteTime

Dir.s<-"D:/Sites_DOE/AmeriFlux/Niwot Ridge/L2_with_gaps/v0008/2009"
# Dir.s<-getwd()
EddyData.F <- fLoadTXTIntoDataframe(file_name, Dir.s)

# get time step from EddyData.F
# get if time step in theREddyProc input is hourly or half-hourly

time_step<-(EddyData.F$Hour[2]- EddyData.F$Hour[1])

# get the number of time steps per day

n_steps_day<-ifelse(time_step==1,24,48)

#+++ Add time stamp in POSIX time format
EddyDataWithPosix.F <- fConvertTimeToPosix(EddyData.F, 'YDH', Year.s='Year', Day.s='DoY', Hour.s='Hour')


#+++ Initialize new sEddyProc processing class
EddyProc.C <- sEddyProc$new(site_code, EddyDataWithPosix.F, c('NEE','Rg','Tair','VPD','Ustar'),DTS.n=n_steps_day)

# estimates Ustar on original series
Ustar <- EddyProc.C$sEstUstarThreshold()$UstarAggr
# estimates Ustar on original series + 5% bootstrap + median bootstrap + 95% bootstrap. This gives 90% confidence interval
# nSample is the number of repetitions in the bootstrap
# Papale et al. 2006 use nSample=100
# Papale et al. 2006 use 5% and 95% percentiles of the 100 bootstrapped thresholds estimates as 90% confidence interval boundaries.
res <- EddyProc.C$sEstUstarThresholdDistribution(nSample=100)


res_Ustar=res[1]
res_U05=res[2]
res_U50=res[3]
res_U95=res[4]

UstarThres.V.v<-c(res_Ustar,res_U05,res_U50,res_U95)
UstarThres.V.m<-matrix(UstarThres.V.v,nrow=1,byrow=TRUE)
UstarSuffix.V.s = c("Ustar", "U05", "U50", "U95")
EddyProc.C$sMDSGapFillAfterUStarDistr('NEE', UstarThres.m.n=UstarThres.V.m, UstarSuffix.V.s=UstarSuffix.V.s)
colnames(EddyProc.C$sExportResults())
resCols <- paste("NEE", UstarSuffix.V.s, "f", sep="_" )
cumNEE <- colSums( EddyProc.C$sExportResults()[,resCols])

#+++ Flux partitioning of one of the different gap filling setups, Note the Suffix.s
  EddyProc.C$sMDSGapFill('Tair', FillAll.b=FALSE)    # Gap-filled Tair needed for partitioning
  EddyProc.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='Ustar')
  EddyProc.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='U05')
  EddyProc.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='U50')
  EddyProc.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='U95')
  
  
  
  #+++ Export gap filled and partitioned data to standard data frame
  FilledEddyData.F <- EddyProc.C$sExportResults()
  #+++ Save results into (tab-delimited) text file in directory \out
  CombinedData.F <- cbind(EddyData.F, FilledEddyData.F)
  
  
  # time step in hours (either 1 hour or 0.5 hours)

h_time_step<-ifelse(time_step==1,1,0.5)

  
# change units from umol/m2/s to gC/m2 in each time step of the tower data (1 hour time step=3600 seconds) 
# conversion factors: 1000000 umol CO2 = 1mol CO2; 1mol C = 1 mol CO2; 1 mol C = 12 g C)
 
# calculate aggregated GPP
 
GPP_Ustar_f_gCm2<-data.frame(CombinedData.F$GPP_Ustar_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_Ustar_f_gCm2)<- c("GPP_Ustar_f_gCm2")
gpp1<-sum(GPP_Ustar_f_gCm2)  
 
 
GPP_U05_f_gCm2<-data.frame(CombinedData.F$GPP_U05_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_U05_f_gCm2)<- c("GPP_U05_f_gCm2")
gpp2<-sum(GPP_U05_f_gCm2) 

GPP_U50_f_gCm2<-data.frame(CombinedData.F$GPP_U50_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_U50_f_gCm2)<- c("GPP_U50_f_gCm2")
gpp3<-sum(GPP_U50_f_gCm2)

GPP_U95_f_gCm2<-data.frame(CombinedData.F$GPP_U95_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_U95_f_gCm2)<- c("GPP_U95_f_gCm2")
gpp4<-sum(GPP_U95_f_gCm2)

gpp_all<-as.data.frame(cbind(gpp1,gpp2,gpp3,gpp4)) 
gpp_all_time<-as.data.frame(cbind(GPP_Ustar_f_gCm2,GPP_U05_f_gCm2,GPP_U50_f_gCm2,GPP_U95_f_gCm2))

# get NEE in g C/m2

NEE_Ustar_f_gCm2<-data.frame(CombinedData.F$NEE_Ustar_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_Ustar_f_gCm2)<- c("NEE_Ustar_f_gCm2")
# calculate aggregated NEE for Up05 condition
NEE1<-sum(NEE_Ustar_f_gCm2)  
 
 
NEE_U05_f_gCm2<-data.frame(CombinedData.F$NEE_U05_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_U05_f_gCm2)<- c("NEE_U05_f_gCm2")
# calculate aggregated NEE for Up05 condition
NEE2<-sum(NEE_U05_f_gCm2) 

NEE_U50_f_gCm2<-data.frame(CombinedData.F$NEE_U50_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_U50_f_gCm2)<- c("NEE_U50_f_gCm2")
NEE3<-sum(NEE_U50_f_gCm2)

NEE_U95_f_gCm2<-data.frame(CombinedData.F$NEE_U95_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_U95_f_gCm2)<- c("NEE_U95_f_gCm2")
NEE4<-sum(NEE_U95_f_gCm2)

NEE_all<-as.data.frame(cbind(NEE1,NEE2,NEE3,NEE4)) 
NEE_all_time<-as.data.frame(cbind(NEE_Ustar_f_gCm2,NEE_U05_f_gCm2,NEE_U50_f_gCm2,NEE_U95_f_gCm2))
  
  
# get total respiration in gC/m2

Reco_Ustar_gCm2<-data.frame(CombinedData.F$Reco_Ustar * h_time_step* 3600 * 12 / 1000000)
names(Reco_Ustar_gCm2)<- c("Reco_Ustar_gCm2")
Reco1<-sum(Reco_Ustar_gCm2) 

Reco_U05_gCm2<-data.frame(CombinedData.F$Reco_U05 * h_time_step* 3600 * 12 / 1000000)
names(Reco_U05_gCm2)<- c("Reco_U05_gCm2")
Reco2<-sum(Reco_U05_gCm2) 

Reco_U50_gCm2<-data.frame(CombinedData.F$Reco_U50 * h_time_step* 3600 * 12 / 1000000)
names(Reco_U50_gCm2)<- c("Reco_U50_gCm2")
Reco3<-sum(Reco_U50_gCm2)

Reco_U95_gCm2<-data.frame(CombinedData.F$Reco_U95 * h_time_step* 3600 * 12 / 1000000)
names(Reco_U95_gCm2)<- c("Reco_U95_gCm2")
Reco4<-sum(Reco_U95_gCm2)

Reco_all<-as.data.frame(cbind(Reco1,Reco2,Reco3,Reco4)) 
Reco_all_time<-as.data.frame(cbind(Reco_Ustar_gCm2,Reco_U05_gCm2,Reco_U50_gCm2,Reco_U95_gCm2))

# get all data
# need to add NEE and Reco with correct units
  
alldataGPPunits<-cbind(CombinedData.F,gpp_all_time,NEE_all_time,Reco_all_time)
fWriteDataframeToFile(alldataGPPunits, 'uncertain_NR1_2009_Partition_Ustar_GPP_units.txt', 'out') 
fWriteDataframeToFile(gpp_all, 'uncertain_NR1_2009_GPP_4values_units.txt', 'out')  
fWriteDataframeToFile(NEE_all, 'uncertain_NR1_2009_NEE_4values_units.txt', 'out') 
fWriteDataframeToFile(Reco_all, 'uncertain_NR1_2009_Reco_4values_units.txt', 'out') 
fWriteDataframeToFile(UstarThres.V.m, 'uncertain_NR1_2009_Ustar_thresholds_4values.txt', 'out')
# write.table(UstarThres.V.m, './out/test_NR1_2009_Ustar_thresholds_10values.txt', col.names=TRUE,row.names=FALSE)
   
  
  
