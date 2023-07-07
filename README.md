# apollo_wxt_data_processing

## Overview
Description of data formats and processing of observational data recorded by APOLLO and WXT weather stations during field experiments FESST@HH 2020 and FESSTVaL 2021. Details on the instruments and experiments can be found in Kirsch et al. (2022) and Hohenegger et al. (2023).

## Files
### Code
- **apollo_level0_to_level1.py**: Converts APOLLO level 0 data to level 1 data
- **wxt_level0_to_level1.py**: Converts WXT level 0 data to level 1 data
- **fesstval_level1_to_level2.py**: Converts APOLLO and WXT level 1 data to level 2 data
- **fesstval_routines.py**: Routines for reading and processing APOLLO and WXT data

### Meta data
- **stations_fessthh.txt**: Meta data of APOLLO and WXT stations for FESST@HH 2020
- **stations_fesstval.txt**: Meta data of APOLLO and WXT stations for FESSTVaL 2021
- **APOLLO_serials.txt**: Serial numbers of APOLLO data loggers 
- **APOLLO_calibration.txt**: Calibration constants (bias) of APOLLO NTC thermometers

## APOLLO data formats
### Level 0
The raw measurement data as produced by instrument is called level 0 data. At the beginning of each measurement period (which usually corresponds to a restart of the logger) a new data file is created. This means, that a single file may either contain only a few seconds or up to a few weeks of measurement data. The files are in the ASCII format and named \<SSSSSS>\_\<YYMMDDHHMM>\.csv, according to the first six digits of the logger serial number and the current time in UTC (=CEST-2h). The data format looks like the following:
```
#START_DATE=20200610;NODE_ID=641064310410E30;Lat.: 53.37873;Long.: 10.19132;FW=v.114a9
#MMDD;HHMISS_LOG(UTC);MMDD_GPS;HHMISS_GPS(UTC);R_NTC(OHM);P_BME(Pa);STATUS(DEC(CGMWLT))
0610;185958;;;14241;101485;5
0610;185959;;;14242;101492;5
0610;190000;;;14235;101481;7
0610;190001;;;14261;101483;7
0610;190002;0610;190001;14261;101480;7
0610;190003;0610;190002;14266;101483;7
```
The first header line contains essential meta data of the data logger:
- START_DATE: Date string YYYYMMDD of start time
- NODE_ID: 13- or 14-digit alphanumeric serial number of the logger
- Lat./Long.: Geographical latitude and longitude of the logger at start time (may not be identical to the actual measurement position, if logger was activated prior to its final installation)
- FW: Firmware version of the logger software

The data format consists of the following columns:
- MMDD: Month and day of the current logger time in UTC (in firmware version 114a9 the day string may miss the leading zero for days <10)
- HHMISS_LOG: Hour, minute, and second of the current logger time in UTC
- MMDD_GPS: Month and day of current GPS time in UTC (only present during active GPS connection, usually for a few seconds at the beginning of each hour)
- HHMISS_GPS: Hour, minute, and second of current GPS time in UTC
- R_NTC: Resistance of NTC thermometer in Ohm
- P_BME: Air pressure of BME280 sensor in Pa
- STATUS: Operation status of the six logger functionalities <ins>C</ins>onsole, <ins>G</ins>PS, <ins>M</ins>easurement, <ins>W</ins>iFi, <ins>L</ins>oRa and <ins>T</ins>TN. The number is a decimal representation of a binary number corresponding to the string of status bits (1=On, 0=Off) in inverted (this results from a human communication error) order. Status 7 would translate into Console=On, GPS=On, Measurement=On, WiFi=Off, LoRa=Off, TTN=Off as its 6-digit binary representation is 000111 corresponding to the string TLWMGC.

### Level 1
Level 1 APOLLO data is measurement data in a standardized format, that is easily accessible for first preliminary analyses. The data has not yet passed any quality checks. The processing steps inlcude:
- converting resistance value of NTC thermometer into temperature
- correcting time drift of logger time using hourly GPS time stamps (level 0 data files without any valid GPS time stamp are ignored)
- creating regular 1-s time grid (i.e., filling up missing time steps with NaNs and removing duplicate time steps)
- Merging sub-daily level 0 files from same days and splitting multi-day files into daily data files      

This data is stored in ASCII files named APOLLO-\<AAA>\_\<YYYYMMDD>\.txt according to respective APOLLO instrument number and date. The data format looks like the following:
```
#DATE=20200610,APOLLO=001,SERIAL=641064310410E3,FW=v.114a9
#LAT=53.37869,LON=10.19198,ALT=3.0,H_T=3.0,H_P=2.0
#TIME(UTC);T_NTC(degC);P_BME(hPa);GPS;WIFI;LORA;TTN
185958;17.13;1014.92;0;0;0;0
185959;17.14;1014.81;1;0;0;0
190000;17.10;1014.83;1;0;0;0
190001;17.10;1014.80;1;0;0;0
190002;17.09;1014.83;1;0;0;0
190003;17.15;1014.82;1;0;0;0
```

The first two header lines contain essential meta data of the data logger and the measurement location:
- DATE: Date string YYYYMMDD 
- APOLLO: APOLLO instrument number (001-110)
- SERIAL: 13- or 14-digit alphanumeric serial number of the logger
- FW: Firmware version of the logger software
- LAT: Geographical latitude of the actual measurement location in °N
- LON: Geographical longitude of the actual measurement location in °E                        
- ALT: Altitude of the measurement location above sea level in m                        
- H_T: Approximate height of the NTC temperature sensor above ground in m
- H_P: Approximate height of the BME280 pressure sensor above ground in m
  
The data format consists of the following columns:
- TIME: Hour, minute, and second of the GPS-corrected and regularized logger time in UTC. Missing measurement time steps (e.g., during maintenance breaks) are filled up with empty data and (rare) duplicate time steps are removed. Files may begin after 00:00:00 UTC or end before 23:59:59 UTC, if the respective measurement data is not available or processed yet. 
- T_NTC: Air temperature of NTC thermometer in °C calculated from resistance *R* according to the simplified Steinhart-Hart equation $$T = \frac{1}{a_0 + a_1 \ln{R} + a_3 \ln^3{R}} - 273.15,$$ where $a_0 = 1.12881 \cdot 10^{-3}$, $a_1 = 2.34177 \cdot 10^{-4}$, and $a_3 = 8.74952 \cdot 10^{-8}$.
- P_BME: Air pressure of BME280 sensor in hPa						 
- GPS: Status of GPS (1=On, 0=Off)
- WIFI: Status of WiFi
- LORA: Status of LoRa
- TTN: Status of TTN

### Level 2
Level 2 APOLLO data has passed several quality and consistency checks and complies with the SAMD product standard (Jahnke-Bornemann, 2022). The performed quality checks are described in Kirsch et al. (2022). The data is stored in daily netCDF4 files named \<CCC>\_uhh_apollo_l2_\<VAR>\_v\<VERS>\_\<YYYYMMDD000000>\.nc
according to the campaign (fessthh=FESST@HH 2020, fval=FESSTVaL 2021), the measurement variable (ta=air temperature, 
pa=air pressure), the dataset version number, and the date. The data header looks like the following:
```
dimensions:
        time = UNLIMITED ; // (86400 currently)
        nv = 2 ;
        station = 82 ;
        character = UNLIMITED ; // (40 currently)
variables:
        int time(time) ;
                time:standard_name = "time" ;
                time:units = "seconds since 1970-01-01 00:00:00 UTC" ;
                time:bounds = "time_bnds" ;
                time:calendar = "standard" ;
        int time_bnds(time, nv) ;
                time_bnds:long_name = "start and end of time averaging intervals" ;
        char station_id(station, character) ;
                station_id:long_name = "station identifier code" ;
        char station_name(station, character) ;
                station_name:long_name = "station name" ;
        float lat(station) ;
                lat:standard_name = "latitude" ;
                lat:long_name = "latitude of instrument location" ;
                lat:units = "degrees_north" ;
        float lon(station) ;
                lon:standard_name = "longitude" ;
                lon:long_name = "longitude of instrument location" ;
                lon:units = "degrees_east" ;
        float zsl(station) ;
                zsl:standard_name = "altitude" ;
                zsl:long_name = "altitude of instrument location above mean sea level" ;
                zsl:comment = "altitude of ground level" ;
                zsl:units = "m" ;
        float zag(station) ;
                zag:standard_name = "height" ;
                zag:long_name = "height of sensor above ground" ;
                zag:units = "m" ;
        char lcz(station, character) ;
                lcz:long_name = "local climate zone of station" ;
                lcz:comments = "after Stewart and Oke (2012), doi:10.1175/BAMS-D-11-00019.1" ;
        float ta(time, station) ;
                ta:_FillValue = NaNf ;
                ta:least_significant_digit = 2 ;
                ta:standard_name = "air_temperature" ;
                ta:long_name = "air temperature" ;
                ta:units = "K" ;
                ta:comment = "Data of station 113PGa smoothed because of higher noise level than usual" ;

// global attributes:
                :_NCProperties = "version=2,netcdf=4.7.3,hdf5=1.10.4" ;
                :Title = "Air temperature observed by APOLLO station network" ;
                :Institution = "Meteorological Institute, University of Hamburg (UHH), Germany" ;
                :Contact_Person = "Prof. Dr. Felix Ament (felix.ament@uni-hamburg.de)" ;
                :Source = "APOLLO (Autonomous cold POoL LOgger), logger software version 114a9" ;
                :History = "Data processed with fesstval_level1_to_level2.py, version 1.0" ;
                :Conventions = "CF-1.7 where applicable" ;
                :Processing_date = "2021-03-18 18:11:00 UTC" ;
                :Author = "Bastian Kirsch (bastian.kirsch@uni-hamburg.de)" ;
                :Comments = "FESST@HH field experiment (June - August 2020), doi:10.25592/uhhfdm.8967" ;
                :Licence = "This data is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License (CC-BY-SA-4.0)." ;
```

The data format contains the following coordinates and variables:
- time: time in seconds since 1970-01-01 00:00:00 UTC (leap seconds are ignored)
- time_bnds: begin and end of measurement interval in seconds since 1970-01-01 00:00:00 UTC
- station_id: Identifier string of stations (007WSa-113PGa for FESST@HH 2020 and 001Sa-102Sa for FESSTVaL 2021)
- station_name: Descriptive name string of stations
- lat: Geographical latitude of stations in °N 
- lon: Geographical longitude of stations in °E
- zsl: Altitude of stations above sea level in m
- zag: Height of sensors above ground in m
- lcz: Local climate zones of stations after Stewart and Oke (2012) 
- ta, pa: Air temperature in K or air pressure in Pa

## WXT data formats
### Level 0
The raw measurement data as produced by instrument is called level 0 data. The data is stored in daily ASCII files named \<WWYYMMDD>\.wxt according to 2-digit WXT station number (01-21) and the date. The data format looks like the following:
```
DATE;TIME;TIMEZONE;VERSION;STATION;FIX;MSG;WXTDATA
16.06.2020;00:00:00;UTC;2.2b;01;1;9;0R9,Fx=1,La=5329.7214N,Lo=01005.1144E,Al=-1.9M
16.06.2020;00:00:02;UTC;2.2b;01;1;3;0r3,Rc=0.00M,Rd=0s,Ri=0.0M,Hc=0.0M,Hd=0s,Hi=0.0M,Rp=0.7M,Hp=0.0MC}J
16.06.2020;00:00:03;UTC;2.2b;01;1;4;0r4,Tr=13.0C,Ra=0.0M,Sl=0.000733V,Rt=1051.0R,Sr=-0.213999VFRr
16.06.2020;00:00:04;UTC;2.2b;01;1;1;0r1,Dn=052D,Dm=057D,Dx=064D,Sn=0.3M,Sm=0.4M,Sx=0.5MI[B
16.06.2020;00:00:07;UTC;2.2b;01;1;2;0r2,Ta=13.8C,Tp=13.9C,Ua=75.2P,Pa=1015.7HGZZ
16.06.2020;00:00:10;UTC;2.2b;01;1;9;0R9,Fx=1,La=5329.7214N,Lo=01005.1145E,Al=-2.1M
16.06.2020;00:00:12;UTC;2.2b;01;1;3;0r3,Rc=0.00M,Rd=0s,Ri=0.0M,Hc=0.0M,Hd=0s,Hi=0.0M,Rp=0.7M,Hp=0.0MC}J
16.06.2020;00:00:13;UTC;2.2b;01;1;4;0r4,Tr=13.1C,Ra=0.0M,Sl=0.001388V,Rt=1051.4R,Sr=-0.059149VNlZ
16.06.2020;00:00:13;UTC;2.2b;01;1;5;0r5,Th=13.4C,Vh=0.0V,Vs=13.0V,Vr=3.618V,Id=N4450549G@R
16.06.2020;00:00:14;UTC;2.2b;01;1;1;0r1,Dn=055D,Dm=070D,Dx=078D,Sn=0.2M,Sm=0.3M,Sx=0.4MNuZ
16.06.2020;00:00:17;UTC;2.2b;01;1;2;0r2,Ta=13.8C,Tp=13.9C,Ua=75.1P,Pa=1015.7HKV_
```

The data format consists of the following columns:
- DATE: GPS-synchronized date in the format DD.MM.YYYY
- TIME: GPS-synchronized time in the format HH:MM:SS
- TIMEZONE: Timezone of time stamp (always UTC)
- VERSION: Version number of data logger software
- STATION: 2-digit WXT station number (01-21)
- FIX: GPS fix, i.e., availability of valid GPS position data (0=no, 1=yes)
- MSG: Number of data message (1, 2, 3, 4, 5, or 9)
- WXTDATA: Data message containing WXT measurement data and GPS information

The meteorological measurement data of the WXT module and geographical position of the GPS receiver are consecutively sent to the data logger as message telegrams. There are six different types of messages (message 1-5 for WXT data and message 9 for GPS data), that are recorded at a fixed sampling interval of 10 s (15 s for message 5), but not at a fixed second. Therefore, the time stamp appears to be irregularly at first glance. The five WXT message strings always start with the message identifier, followed by the comma-separated measurement data and end with a 3-digit check sum. The following measurement parameters are available:

Message 1:
- Dn: Minimum wind direction in ° (referring to the 4-Hz raw sonic anemometer data within the sampling interval)
- Dm: Average wind direction in °
- Dx: Maximum wind direction in °
- Sn: Minimum wind speed in m/s
- Sm: Average wind speed in m/s
- Sx: Maximum wind speed in m/s

Message 2:
- Ta: Air temperature in °C
- Tp: Internal temperature in °C
- Ua: Relative humidity of air in %
- Pa: Air pressure in hPa

Message 3:
- Rc: Rainfall accumulation in mm
- Rd: Rainfall duration in s
- Ri: Rainfall intensity in mm/h
- Hc: Hail accumulation in hits/cm^2
- Hd: Hail duration in s
- Hi: Hail intensity in hits/(cm^2 h)
- Rp: Rainfall peak intensity in mm/h (identical to Ri in the used configuration)
- Hp: Hail peak intensity in hits/(cm^2 h) (identical to Hi in the used configuration)

Message 4: 
- Tr: Pt1000 air temperature in °C
- Ra: Tipping bucket rainfall accumulation in mm (unused)
- Sl: Voltage of ultrasonic level sensor in V (unused)
- Rt: Pt1000 resistance in Ohm
- Sr: Solar radiation in V (unused)

Message 5:
- Th: Heating temperature in °C
- Vh: Heating voltage in V
- Vs: Supply voltage in V
- Vr: Reference voltage in V
- Id: 8-digit alphanumeric serial number of WXT sensor

Message 9:
- Fx: GPS fix
- La: Geographical latitude in format ddmm.mmmmu °N, where dd, mm.mmmm and u denote integer degrees, decimal minutes, and the hemisphere, respectively. The value 5330.7500N would translate into 53°30.75'N=53°30'45''N=53.5125°N
- Lo: Geographical longitude in format dddmm.mmmm °E
- Al: Altitude above sea level in m

### Level 1
Similar to the APOLLO data format, level 1 data is measurement data in a standardized format, that is easily accessible, but has not yet passed any quality checks. The following processing steps are performed:
- extracting measurement data from message telegrams
- rounding time stamp of individual messages to next 10-s time step
- creating regular 10-s time grid (i.e., filling up missing time steps with NaNs and removing duplicate or overlapping time steps in case of corrupt GPS data)

Accordingly, the data is stored in ASCII files named WXT-\<WW>\_\<YYYYMMDD>\.txt, that look like the following:
```
#DATE=20200616,WXT=01,SERIAL=N4450549,FW=2.2b
#LAT=53.49537N,LON=10.08528E,ALT=2.0M,H_WXT=3.0M
#TIME(UTC);TT(degC);TT_PT1000(degC);PP(hPa);RH(%);FF(m/s);FB(m/s);DD(deg);RR(mm);HAIL_HITS(1/cm2);U_SUPPLY(V)
000000;;;;;;;;;;
000010;13.8;13.0;1015.7;75.2;0.4;0.5;57;0.0;0.0;
000020;13.8;13.1;1015.7;75.1;0.3;0.4;70;0.0;0.0;13.0
000030;13.9;13.0;1015.7;74.9;0.3;0.3;72;0.0;0.0;13.0
000040;13.9;13.1;1015.7;75.0;0.3;0.4;63;0.0;0.0;
000050;13.9;13.0;1015.7;75.3;0.4;0.5;59;0.0;0.0;13.0
000100;13.9;13.0;1015.7;75.2;0.4;0.5;65;0.0;0.0;13.0
000110;13.9;13.0;1015.7;75.2;0.3;0.4;76;0.0;0.0;
```

While the meta data in the first two header lines correspond to these in the APOLLO level 1 format, additional data columns for selected WXT measurement parameters are included:
- TT: WXT air temperature in °C (message parameter Ta)
- TT_PT1000: Pt1000 air temperature in °C (Tr)
- PP: Air pressure in hPa (Pa)
- RH: Relative humidity of air in % (Ua)
- FF: Wind speed in m/s (Sm)
- FB: Maximum wind speed (gust) in m/s (Sx)
- DD: Wind direction in ° (Dm)
- RR: Rainfall accumulation in mm (Rc)
- HAIL_HITS: Hail accumulation in hits/cm^2 (Hc)
- U_SUPPLY: Supply voltage in V (Vs)
  
### Level 2
The WXT level 2 data format is defined analogically to the APOLLO data. Due to internal quality checks implemented by the manufacturer of the sensor, only basic plausibility checks are performed on the data. The data is stored in daily netCDF4 files named \<CCC>\_uhh_wxt_l2_\<VAR>\_v\<VERS>\_\<YYYYMMDD000000>\.nc, where the following variables are available:
- ta: Pt1000 air temperature in K
- pa: Air temperature in Pa
- hur: Relative humidity of air in 1
- wspeed: Wind speed in m/s
- wspeed_max: Maximum wind speed (gust) in m/s
- wdir : Wind direction in °
- precip: Rainfall accumulation in mm
- hail: Hail accumulation in hits/m^2

## References
- Hohenegger, C. and the FESSTVaL team (2023): *FESSTVaL: the Field Experiment on Submesoscale Spatio-Temporal Variability in Lindenberg*, Bull. Am. Meteorol. Soc. (in review)
- Jahnke-Bornemann, A. (2022): *The SAMD product standard (Standardized Atmospheric Measurement Data) version 2.2*, http://doi.org/10.25592/UHHFDM.10416.
- Kirsch, B., Hohenegger, C., Klocke, D., Senke, R., Offermann, M., and Ament, F. (2022): *Sub-mesoscale observations of convective cold pools with a dense station network in Hamburg, Germany*, Earth Syst. Sci. Data, 14, 3531–3548, https://doi.org/10.5194/ESSD-14-3531-2022.
- Stewart, I. D., and T. R. Oke (2012): *Local climate zones for urban temperature studies*, Bull. Am. Meteorol. Soc., 93, 1879-1900, https://doi.org/10.1175/BAMS-D-11-00019.1.

## Data sets
- APOLLO and WXT level 2 data of FESST@HH 2020: https://doi.org/10.25592/UHHFDM.10172
- APOLLO and WXT level 2 data of FESSTVaL 2021: https://doi.org/10.25592/UHHFDM.10179

## Data policy
The data policy for FESSTVaL campaign data applies (https://doi.org/10.25592/UHHFDM.10181).

## Contact
Bastian Kirsch (bastian.kirsch@uni-hamburg.de)<br>
Felix Ament (felix.ament@uni-hamburg.de)<br>
Meteorologisches Institut, Universität Hamburg, Germany

Last updated: 5 July 2023
