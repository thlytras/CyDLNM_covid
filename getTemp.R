library(FluDLNM)

cyStations <- NOAA_countryStations("CY", from=2010)
#cyStations <- subset(cyStations, format(end, "%Y")==2021)
cyStations <- subset(cyStations, station.name %in% c("PAFOS INTL", "AKROTIRI", "NICOSIA/ATHALASSA", "LARNACA"))

cyTemps <- NOAA_getGSOD(cyStations, 2010:2021)


save(cyStations, cyTemps, file="input/cyTemps.RData")

# with(subset(cyTemps, station.name=="AKROTIRI"), plot(date, temp, type="l", col="blue", lwd=0.3))
# with(subset(cyTemps, station.name=="LARNACA"), points(date, temp, type="l", col="orange", lwd=0.3))
# with(subset(cyTemps, station.name=="PAFOS INTL"), points(date, temp, type="l", col="green", lwd=0.3))
