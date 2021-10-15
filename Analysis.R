##==============================================================================
##
## Script performs the analysis required for the paper:
## "Runoff Coefficients of High-flow Events in Undisturbed New England Basins"
## by Hosseini-Shakib et al.(2020)
##
## Author: Iman Hosseini-Shakib (ishakib@gmail.com)
##       
##==============================================================================
## Copyright 2021 Iman Hosseini-Shakib
## This file is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This file is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this file.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================

rm(list = ls())
graphics.off()

# Packages used
library(readr)
library(rstudioapi)
library(ggplot2)
library(Kendall)
library(xlsx)
library(RHydro)
library(highcharter)
library(dplyr)
library(htmlwidgets)
library(webshot)
library(beepr)

# Set working directory to source file location
wd=setwd(dirname(getActiveDocumentContext()$path)) 


# Assigning values to the following parameters:
nstations = 28 #number of USGS stream flow gauging stations to consider
Return = 2   #The assumed return period of the design flood which can be 2, 5 or 10
ndays = 5 #number of days to consider before (and after) which the flood event hydrograph tails are minimum of
baseflow_method = 3 #can be 0 for no baseflow separation, 1 for constant slope, 2 for constant discharge or "constant_slope" with slope=0 and 3 for "RLSWM" or concave method

#   1-1-Reading data
Sfolder=paste0(wd,"/Input_Data/snow_water-equivalent/") #directory of snow water equivalent timeseries
Tfolder=paste0(wd,"/Input_Data/temperature/") #directory of temperature timeseries
PRfolder = paste0(wd,"/Input_Data/precip_runoff/") #directory of precipitation & runoff timeseries
RPfolder = paste0(wd,"/Input_Data/") #directory of 2, 5 and 10 yr floods of stations
dir.create(paste0(wd,"/Output/")) #create a folder for outputs
OUTPUTfolder = paste0(wd,"/Output/") #directory of the output folder
 
RP = read.csv(paste0(RPfolder, "PeakFlowStatistics.csv")) #flood frequency of stations
USGS = as.character(paste0("0",RP$USGS.Gage.Station.Number))

for (f in 1:nstations) {
  print(f)
  if (Return == 2)
   design.flood = RP[f, 7] #calculated from daily data
  if (Return == 5)
   design.flood = RP[f, 8] #calculated from daily data
  if (Return == 10)
   design.flood = RP[f, 9] #calculated from daily data

    b <- read_delim(
    paste0(PRfolder, f, ".csv"),
    ",",
    escape_double = FALSE,
    trim_ws = TRUE,
    col_names = TRUE
    )
  b=b[-(1:273),] #to start from Oct 1, same as snow data
  SWE=read.csv(paste0(Sfolder,f,".csv"))
  b$S=SWE$SNDP[1:(nrow(b))]
  SWE$Date=as.Date(SWE$Date)
  TEMP=read.csv(paste0(Tfolder,f,".csv"))
  TEMP=TEMP[-(1:273),] #to start from Oct 1, same as snow data
  b$T=TEMP$TEMP[1:(nrow(b))]
  rm(SWE,TEMP)

    #   1-2-Separating year, month, day, runoff and precipitation columns
  
  date = as.Date(b$Date, "%m/%d/%Y")
  year = as.numeric(format(date, "%Y"))
  month = as.numeric(format(date, "%m"))
  day = as.numeric(format(date, "%d"))
  P = b[3]
  R = b[4]
  S = b[5]
  T = b[6]
  table <- data.frame(date, year, month, day, R, P,S,T)
  colnames(table) = c("date", "year", "month", "day", "RMCM", "PMCM","SMCM","T")
  rm(b,P, R,S,T,year, month, day, date)
  
  #   2- Finding trend in runoff coefficients of large events
  n = nrow(table)
  TotalR =matrix(0, n, 1)
  R = matrix(0, n, 1)
  P = matrix(0, n, 1)
  C = matrix(0, n, 1)
  BF = matrix(0, n, 1)
  Start=as.Date(matrix(0,n,1))
  End=as.Date(matrix(0,n,1))
  AvailableSWE = matrix(0, n, 1)
  AntecedentSWE = matrix(0, n, 1)
  SubsequentSWE = matrix(0, n, 1)
  AntecedentTEMP = matrix(0, n, 1)
  AntecedentBF = matrix(0, n, 1)
  conctime = matrix(0, n, 1)
  D = data.frame(table$date[1:n])
  z = table$RMCM
  y = table$PMCM
  x = table$date
  w = table$SMCM
  v = table$T
  path = paste0(OUTPUTfolder, f)
  dir.create(path)
  TailR=0
    for (i in 6:(n-6)) {
    while(i<=TailR) i=i+1
    max1 = max(z[(i-6):(i-1)])
    max2 = max(z[(i+1):(i+6)])
    max=max(max1,max2)
    if (z[i] >= max & z[i] >= design.flood) {
      j=i
      while(z[j]>0.5*max) j=j-1 
      while(z[j]> min(z[(j-ndays):(j)])) j=j-1
      TailL = j 
      j=i
      while(z[j]>0.5*max) j=j+1
      while(z[j]> min(z[(j):(j+ndays)])) j=j+1
      TailR = j 
      
      tc=time_of_conc(y[TailL:TailR],z[TailL:TailR]) #basin concentration time
      tl=max(1,ceiling(tc)) #basin lag time tl=ceiling(0.6*tc)
      
      hydrograph = z[(TailL-tl):(TailR+tl)]
      hyetograph = y[(TailL-tl):(TailR+tl)]
      xaxis_date = as.Date(x[(TailL-tl):(TailR+tl)],format=c("%Y,%m,%d"))
      snow = w[(TailL-tl):(TailR+tl)]
      temperature = v[(TailL-tl):(TailR+tl)]
      
      if (baseflow_method == 0)
        baseflow = rep(0, times=(TailR - TailL+2*tl + 1))
      if (baseflow_method == 1)
        baseflow = baseflow_sep(hydrograph,
                                method = "DFM",
                                parms = c(c = 0.99))
      if (baseflow_method == 2)
        baseflow = baseflow_sep(hydrograph,
                                method = "constant_slope",
                                parms = c(c_slope = 0))
      if (baseflow_method == 3)
        baseflow = baseflow_sep(hydrograph, method = "RLSWM")
      
      
      day1=as.Date(xaxis_date[1])
      dayn=as.Date(xaxis_date[length(xaxis_date)])
      AnteSWE.Date=xaxis_date[1]-1 #Antecedent snow water equivalent one day before the start of precipitation event
      SubseSWE.Date=xaxis_date[length(xaxis_date)-tl]+1 #Subsequent snow water equivalent one day after the end of precipitation event
      Stable=table[table$date==AnteSWE.Date,]
      Etable=table[table$date==SubseSWE.Date,]
      AnteSWE=Stable$SMCM
      SubseSWE=Etable$SMCM
      AnteTEMP.Date=xaxis_date[1]-1 #Antecedent Temperature one day before the start of event
      Ttable=table[table$date==AnteTEMP.Date,]
      AnteTEMP=Ttable$T
      AnteBF=Stable$RMCM
      Melt=max(0,(AnteSWE-SubseSWE))
      
      hydrograph[1:tl]=as.numeric(0)
      hyetograph[(length(hyetograph)-tl):(length(hyetograph))]=as.numeric(0)
      baseflow[1:tl]=as.numeric(0)
      
      
      
      chart <- highchart() %>%
              hc_yAxis_multiples(list(
         title = list(text = "rainfall (MCM)"),
         reversed = TRUE
       ),
       list(
         title = list(text = "flow (MCM)"),
         opposite = TRUE
       )) %>%
       hc_add_series(name = "Rainfall",
                     data = hyetograph ,
                      type = "column") %>% #%>% mutate(value = value) %>% .$value
       hc_add_series(
         name = "Runoff",
         data = hydrograph ,
         type = "spline",
         yAxis = 1
       ) %>% #%>% .$value
       hc_add_series(
         name = "Baseflow",
         data = baseflow ,
         type = "spline",
         yAxis = 1
       ) %>% #%>% .$value
       hc_xAxis(categories = xaxis_date, title = list(text = "date"))
      
      setwd(path)
      saveWidget(widget = chart,
                file = paste0("station ", f, " flood ", x[i], ".html"))
      webshot(
        url = paste0("station ", f, " flood ", x[i], ".html"),
        file = paste0("station ", f, " flood ", x[i], ".jpeg"),
        delay = 1
        )
      event_table = data.frame(xaxis_date, hyetograph, hydrograph, baseflow,snow,Melt,AnteSWE,SubseSWE,AnteTEMP,tl)
      colnames(event_table) = c(
        "Date",
        "Precipitation (MCM)",
        "Total Runoff (MCM)",
        "Baseflow (MCM)",
        "Snow Water Equivalent(MCM)",
        "Snowmelt (MCM)",
        "Antecedent Snow Water Equivalent (MCM)",
        "Subsequent Snow Water Equivalent (MCM)",
        "Antecedent Temperature (deg.C)",
        "Time of Concentration (day)"
      )
      write.xlsx2(event_table, paste0("station ", f," Event Table.xlsx"),
                  sheetName = paste0(x[i]),
                  append = T)
      Start[i]=day1
      End[i]=dayn
      AntecedentSWE[i]=mean(AnteSWE)
      SubsequentSWE[i]=mean(SubseSWE)
      AntecedentTEMP[i]=mean(AnteTEMP)
      AntecedentBF[i]=AnteBF
      conctime[i]=ceiling(mean(tl))
      P[i] = sum(hyetograph[1:(length(hyetograph)-tl)])
      BF[i] = sum(baseflow[tl:length(baseflow)])
      AvailableSWE[i]=Melt
      TotalR[i]=sum(hydrograph[tl:length(hydrograph)])
      R[i] = max(0,TotalR[i]-BF[i])
      C[i] = R[i] / (P[i]+AvailableSWE[i])
      
    }
  }
  large.events = data.frame(Start,D,End, TotalR,R, P, C,BF,AvailableSWE,AntecedentSWE,SubsequentSWE,AntecedentTEMP,AntecedentBF,conctime)
  large.events = subset(large.events,large.events$C >0&large.events$C<1)
  df = data.frame(large.events$table.date.1.n., large.events$C)
  x = large.events$table.date.1.n.
  y = large.events$C
  if (length(y) > 3) {
    MK = MannKendall(y)
    Trend = cbind(MK$tau, MK$sl, MK$S, MK$D, MK$varS)
    colnames(Trend) = c("Tau", "sl", "S", "D", "varS")
  } else {
    MK = 0
    Trend = 0
  }
  rm(df, i, max1,max2,max, n, R, P, C, D, z)
  large.events$month = as.numeric(format(large.events$table.date.1.n., "%m"))
  
  #   3- Exporting Files
  
  write.xlsx2(large.events,
              paste0(OUTPUTfolder, f, " Output.xlsx"),
              sheetName = "01-LargeEvents")

  rm(x, y)
}

#   5- Table of Trends & p-values
for (i in 1:nstations){
file1 = read.xlsx(paste0(OUTPUTfolder, i, " Output.xlsx"),sheetName = "01-LargeEvents")
x = file1$table.date.1.n.
y = file1$C
if (length(y) > 3) {
  MK = MannKendall(y)
  Trend = cbind(MK$tau, MK$sl, MK$S, MK$D, MK$varS)
  colnames(Trend) = c("Tau", "sl", "S", "D", "varS")
} else {
  MK = 0
  Trend = 0
}
write.xlsx2(Trend, file = paste0(OUTPUTfolder, "Slopes.xlsx"),sheetName = paste0(i),append = T)
}

slopes.table = as.data.frame(NULL)
for (i in 1:nstations) {
  file2 = read.xlsx(paste0(OUTPUTfolder, "Slopes.xlsx"), sheetName = paste0(i))
    if (ncol(file2) == 6) {
    slopes.table[i,1] = paste0(i)
    slopes.table[i, 2] = file2[1, 4]
    slopes.table[i, 3] = file2[1, 3]
    slopes.table[i, 4] = if (file2[1, 4] > 0)
      "+"
    else
      "-"
    slopes.table[i, 5] = if (file2[1, 3] <= 0.10)
      "Yes"
    else
      "No"
  } else {
    slopes.table[i,] = 0
  }
}
colnames(slopes.table) = c("Basin",
                           "Mann-Kendall S-value",
                           "Mann-Kendall p-value",
                           "Trend",
                           "Significance at 90% CL")

write.xlsx(slopes.table,
           file = paste0(OUTPUTfolder, "Slopes.xlsx"),
           sheetName = "Slopes", append = T)
setwd(wd)
rm(list = ls())
beep(sound = 2, expr = "All Done!")
