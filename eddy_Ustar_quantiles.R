
# This script:
# - reads a REddyProc input file with data for only 1 year and one site(the site must be included in "sites_df", currently only 7 sites)
# - gets if the time step from the data is hourly or half-hourly
# - estimates different Ustar thresholds based on Ustar data quantiles on given probabilities,
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
EddySetups.C <- sEddyProc$new(site_code, EddyDataWithPosix.F, c('NEE','Rg','Tair','VPD','Ustar'),DTS.n=n_steps_day)

#+++ estimate different u* thresholds from the data
UstarThres_year.V.n <- quantile( subset(EddyDataWithPosix.F,Year==2009)$Ustar, probs=c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), na.rm=T)
# here we have only one year
UstarThres.V.m <- matrix( UstarThres_year.V.n, nrow=1, byrow=TRUE )
# suffixes to distinguish different Ustar setups
UstarSuffix.V.s <- c("Up05","Up10","Up15","Up20","Up25","Up30","Up35","Up40","Up45","Up50")    
# the next line takes some time to complete
EddySetups.C$sMDSGapFillAfterUStarDistr('NEE', UstarThres.m.n=UstarThres.V.m, UstarSuffix.V.s=UstarSuffix.V.s )
 colnames(EddySetups.C$sExportResults()) # Note the suffix in output columns
# inspect the mean across NEE estimates and uncertainty introduced by different uStar thresholds
resCols <- paste("NEE", UstarSuffix.V.s, "f", sep="_" )
cumNEE <- colSums( EddySetups.C$sExportResults()[,resCols])  

#+++ Flux partitioning of one of the different gap filling setups, Note the Suffix.s
  EddySetups.C$sMDSGapFill('Tair', FillAll.b=FALSE)    # Gap-filled Tair needed for partitioning
  EddySetups.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='Up05')
  EddySetups.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='Up10')
  EddySetups.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='Up15')
  EddySetups.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='Up20')
  EddySetups.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='Up25')
  EddySetups.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='Up30')
  EddySetups.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='Up35')
  EddySetups.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='Up40')
  EddySetups.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='Up45')
  EddySetups.C$sMRFluxPartition(Lat_deg.n=siteLat, Long_deg.n=siteLon, TimeZone_h.n=siteTime, Suffix.s='Up50')
  colnames(EddySetups.C$sExportResults())    # Note the suffix R and GPP output columns
  

  #+++ Export gap filled and partitioned data to standard data frame
  FilledEddyData.F <- EddySetups.C$sExportResults()
  #+++ Save results into (tab-delimited) text file in directory \out
  CombinedData.F <- cbind(EddyData.F, FilledEddyData.F)
  #fWriteDataframeToFile(CombinedData.F, 'NR1_1999_Partition_Ustar.txt', 'out')

  
  # time step in hours (either 1 hour or 0.5 hours)

h_time_step<-ifelse(time_step==1,1,0.5)

  
  # change units from umol/m2/s to gC/m2 in each time step of the tower data (1 hour time step=3600 seconds) 
# conversion factors: 1000000 umol CO2 = 1mol CO2; 1mol C = 1 mol CO2; 1 mol C = 12 g C)
 
GPP_Up05_f_gCm2<-data.frame(CombinedData.F$GPP_Up05_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_Up05_f_gCm2)<- c("GPP_Up05_f_gCm2")
# calculate aggregated GPP for Up05 condition
gpp1<-sum(GPP_Up05_f_gCm2) 

GPP_Up10_f_gCm2<-data.frame(CombinedData.F$GPP_Up10_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_Up10_f_gCm2)<- c("GPP_Up10_f_gCm2")
gpp2<-sum(GPP_Up10_f_gCm2)

GPP_Up15_f_gCm2<-data.frame(CombinedData.F$GPP_Up15_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_Up15_f_gCm2)<- c("GPP_Up15_f_gCm2")
gpp3<-sum(GPP_Up15_f_gCm2)

GPP_Up20_f_gCm2<-data.frame(CombinedData.F$GPP_Up20_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_Up20_f_gCm2)<- c("GPP_Up20_f_gCm2")
gpp4<-sum(GPP_Up20_f_gCm2)

GPP_Up25_f_gCm2<-data.frame(CombinedData.F$GPP_Up25_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_Up25_f_gCm2)<- c("GPP_Up25_f_gCm2")
gpp5<-sum(GPP_Up25_f_gCm2)

GPP_Up30_f_gCm2<-data.frame(CombinedData.F$GPP_Up30_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_Up30_f_gCm2)<- c("GPP_Up30_f_gCm2")
gpp6<-sum(GPP_Up30_f_gCm2)

GPP_Up35_f_gCm2<-data.frame(CombinedData.F$GPP_Up35_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_Up35_f_gCm2)<- c("GPP_Up35_f_gCm2")
gpp7<-sum(GPP_Up35_f_gCm2)

GPP_Up40_f_gCm2<-data.frame(CombinedData.F$GPP_Up40_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_Up40_f_gCm2)<- c("GPP_Up40_f_gCm2")
gpp8<-sum(GPP_Up40_f_gCm2)

GPP_Up45_f_gCm2<-data.frame(CombinedData.F$GPP_Up45_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_Up45_f_gCm2)<- c("GPP_Up45_f_gCm2")
gpp9<-sum(GPP_Up45_f_gCm2)

GPP_Up50_f_gCm2<-data.frame(CombinedData.F$GPP_Up50_f * h_time_step* 3600 * 12 / 1000000)
names(GPP_Up50_f_gCm2)<- c("GPP_Up50_f_gCm2")
gpp10<-sum(GPP_Up50_f_gCm2)

gpp_all<-as.data.frame(cbind(gpp1,gpp2,gpp3,gpp4,gpp5,gpp6,gpp7,gpp8,gpp9,gpp10)) 
gpp_all_time<-as.data.frame(cbind(GPP_Up05_f_gCm2,GPP_Up10_f_gCm2,GPP_Up15_f_gCm2,GPP_Up20_f_gCm2,GPP_Up25_f_gCm2,GPP_Up30_f_gCm2,GPP_Up35_f_gCm2,GPP_Up40_f_gCm2,GPP_Up45_f_gCm2,GPP_Up50_f_gCm2))

# get NEE in g C/m2

NEE_Up05_f_gCm2<-data.frame(CombinedData.F$NEE_Up05_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_Up05_f_gCm2)<- c("NEE_Up05_f_gCm2")
# calculate aggregated NEE for Up05 condition
NEE1<-sum(NEE_Up05_f_gCm2) 

NEE_Up10_f_gCm2<-data.frame(CombinedData.F$NEE_Up10_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_Up10_f_gCm2)<- c("NEE_Up10_f_gCm2")
NEE2<-sum(NEE_Up10_f_gCm2)

NEE_Up15_f_gCm2<-data.frame(CombinedData.F$NEE_Up15_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_Up15_f_gCm2)<- c("NEE_Up15_f_gCm2")
NEE3<-sum(NEE_Up15_f_gCm2)

NEE_Up20_f_gCm2<-data.frame(CombinedData.F$NEE_Up20_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_Up20_f_gCm2)<- c("NEE_Up20_f_gCm2")
NEE4<-sum(NEE_Up20_f_gCm2)

NEE_Up25_f_gCm2<-data.frame(CombinedData.F$NEE_Up25_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_Up25_f_gCm2)<- c("NEE_Up25_f_gCm2")
NEE5<-sum(NEE_Up25_f_gCm2)

NEE_Up30_f_gCm2<-data.frame(CombinedData.F$NEE_Up30_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_Up30_f_gCm2)<- c("NEE_Up30_f_gCm2")
NEE6<-sum(NEE_Up30_f_gCm2)

NEE_Up35_f_gCm2<-data.frame(CombinedData.F$NEE_Up35_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_Up35_f_gCm2)<- c("NEE_Up35_f_gCm2")
NEE7<-sum(NEE_Up35_f_gCm2)

NEE_Up40_f_gCm2<-data.frame(CombinedData.F$NEE_Up40_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_Up40_f_gCm2)<- c("NEE_Up40_f_gCm2")
NEE8<-sum(NEE_Up40_f_gCm2)

NEE_Up45_f_gCm2<-data.frame(CombinedData.F$NEE_Up45_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_Up45_f_gCm2)<- c("NEE_Up45_f_gCm2")
NEE9<-sum(NEE_Up45_f_gCm2)

NEE_Up50_f_gCm2<-data.frame(CombinedData.F$NEE_Up50_f * h_time_step* 3600 * 12 / 1000000)
names(NEE_Up50_f_gCm2)<- c("NEE_Up50_f_gCm2")
NEE10<-sum(NEE_Up50_f_gCm2)

NEE_all<-as.data.frame(cbind(NEE1,NEE2,NEE3,NEE4,NEE5,NEE6,NEE7,NEE8,NEE9,NEE10)) 
NEE_all_time<-as.data.frame(cbind(NEE_Up05_f_gCm2,NEE_Up10_f_gCm2,NEE_Up15_f_gCm2,NEE_Up20_f_gCm2,NEE_Up25_f_gCm2,NEE_Up30_f_gCm2,NEE_Up35_f_gCm2,NEE_Up40_f_gCm2,NEE_Up45_f_gCm2,NEE_Up50_f_gCm2)) 

# get total respiration in gC/m2

Reco_Up05_gCm2<-data.frame(CombinedData.F$Reco_Up05 * h_time_step* 3600 * 12 / 1000000)
names(Reco_Up05_gCm2)<- c("Reco_Up05_gCm2")
# calculate aggregated Reco for Up05 condition
Reco1<-sum(Reco_Up05_gCm2) 

Reco_Up10_gCm2<-data.frame(CombinedData.F$Reco_Up10 * h_time_step* 3600 * 12 / 1000000)
names(Reco_Up10_gCm2)<- c("Reco_Up10_gCm2")
Reco2<-sum(Reco_Up10_gCm2)

Reco_Up15_gCm2<-data.frame(CombinedData.F$Reco_Up15 * h_time_step* 3600 * 12 / 1000000)
names(Reco_Up15_gCm2)<- c("Reco_Up15_gCm2")
Reco3<-sum(Reco_Up15_gCm2)

Reco_Up20_gCm2<-data.frame(CombinedData.F$Reco_Up20 * h_time_step* 3600 * 12 / 1000000)
names(Reco_Up20_gCm2)<- c("Reco_Up20_gCm2")
Reco4<-sum(Reco_Up20_gCm2)

Reco_Up25_gCm2<-data.frame(CombinedData.F$Reco_Up25 * h_time_step* 3600 * 12 / 1000000)
names(Reco_Up25_gCm2)<- c("Reco_Up25_gCm2")
Reco5<-sum(Reco_Up25_gCm2)

Reco_Up30_gCm2<-data.frame(CombinedData.F$Reco_Up30 * h_time_step* 3600 * 12 / 1000000)
names(Reco_Up30_gCm2)<- c("Reco_Up30_gCm2")
Reco6<-sum(Reco_Up30_gCm2)

Reco_Up35_gCm2<-data.frame(CombinedData.F$Reco_Up35 * h_time_step* 3600 * 12 / 1000000)
names(Reco_Up35_gCm2)<- c("Reco_Up35_gCm2")
Reco7<-sum(Reco_Up35_gCm2)

Reco_Up40_gCm2<-data.frame(CombinedData.F$Reco_Up40 * h_time_step* 3600 * 12 / 1000000)
names(Reco_Up40_gCm2)<- c("Reco_Up40_gCm2")
Reco8<-sum(Reco_Up40_gCm2)

Reco_Up45_gCm2<-data.frame(CombinedData.F$Reco_Up45 * h_time_step* 3600 * 12 / 1000000)
names(Reco_Up45_gCm2)<- c("Reco_Up45_gCm2")
Reco9<-sum(Reco_Up45_gCm2)

Reco_Up50_gCm2<-data.frame(CombinedData.F$Reco_Up50 * h_time_step* 3600 * 12 / 1000000)
names(Reco_Up50_gCm2)<- c("Reco_Up50_gCm2")
Reco10<-sum(Reco_Up50_gCm2)

Reco_all<-as.data.frame(cbind(Reco1,Reco2,Reco3,Reco4,Reco5,Reco6,Reco7,Reco8,Reco9,Reco10)) 
Reco_all_time<-as.data.frame(cbind(Reco_Up05_gCm2,Reco_Up10_gCm2,Reco_Up15_gCm2,Reco_Up20_gCm2,Reco_Up25_gCm2,Reco_Up30_gCm2,Reco_Up35_gCm2,Reco_Up40_gCm2,Reco_Up45_gCm2,Reco_Up50_gCm2))

# get all data
# need to add NEE and Reco with correct units
  
alldataGPPunits<-cbind(CombinedData.F,gpp_all_time,NEE_all_time,Reco_all_time)
fWriteDataframeToFile(alldataGPPunits, 'NR1_2009_Partition_Ustar_GPP_units.txt', 'out') 
fWriteDataframeToFile(gpp_all, 'NR1_2009_GPP_10values_units.txt', 'out')  
fWriteDataframeToFile(NEE_all, 'NR1_2009_NEE_10values_units.txt', 'out') 
fWriteDataframeToFile(Reco_all, 'NR1_2009_Reco_10values_units.txt', 'out') 
fWriteDataframeToFile(UstarThres.V.m, 'NR1_2009_Ustar_thresholds_10values.txt', 'out')
# write.table(UstarThres.V.m, './out/test_NR1_2009_Ustar_thresholds_10values.txt', col.names=TRUE,row.names=FALSE)
 
