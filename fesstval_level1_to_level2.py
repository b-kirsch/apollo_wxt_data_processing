# -*- coding: utf-8 -*-
"""
@author: Bastian Kirsch (bastian.kirsch@uni-hamburg.de)

Script to convert level1 APOLLO (Autonomous cold POoL LOgger) and 
WXT data to level2 daily .nc-files complying with the SAMD data standard

Data structure needed in datadir: ./<YYYY>/<MM>/<DD>/<APOLLO,WXT>-<SSS,SS>_<YYYYMMDD>.txt
Output data structure in outdir: ./<YYYY>/<MM>/<DD>/fessthh_uhh_<instr>_l2_<var>_v00_<YYYYMMDDhhmmss>.nc

Processing steps:
    - removing values outside plausible ranges
    - removing spikes in NTC temperature during active Wifi
    - removing implausible values compared to the network mean
    - removing periods at beginning and end of measurement periods
    - manually removing and correcting selected stations and periods
    - merging data from different instruments at same station
    - merging daily data from all stations
    - calibrating APOLLO temperature sensors
    - writing level2 data as SAMD-conform netCDF4 files
    
Dependences on non-standard software:
- fesstval_routines.py 

Required meta data files:
- stations_fessthh.txt
- APOLLO_serials.txt    
- APOLLO_calibration.txt
    
Last updated: 5 February 2021
"""

import numpy as np
import pandas as pd
import datetime as dt
from netCDF4 import Dataset
from netCDF4 import stringtochar
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import fesstval_routines as fst

t_run = dt.datetime.now()
#----------------------------------------------------------------------------
# General settings
maindir          = '.'
datadir_apo      = maindir+'APOLLO_data/level1/'
outdir_apo       = maindir+'APOLLO_data/level2/'
datadir_wxt      = maindir+'WXT_data/level1/'
datadir_wxt_perm = maindir+'Wettermast_data/2020_FESSTHH/'
outdir_wxt       = maindir+'WXT_data/level2/'
plotdir          = maindir+'Cold-Pools/Plots/APOLLO_Quicklooks/'
meta_file        = maindir+'FESSTVaL/FESSTHH/stations_fessthh.txt'
serial_file      = maindir+'FESSTVaL/APOLLO/apollo_serials.txt'
cali_file        = maindir+'FESSTVaL/APOLLO/apollo_calibration.txt'

start_date = dt.date(2020,6,1)
end_date   = dt.date(2020,8,31)

#type_proc = 'APOLLO'
type_proc = 'WXT'

write_data  = False

#----------------------------------------------------------------------------
# Processing settings

wxt_permanent = {'002MIw':'WXT_WMH',
                 '004MIw':'WXT_BBG'}
time_res      = {'APOLLO':1, 
                 'WXT':10}
outdirs      = {'APOLLO':outdir_apo, 
                 'WXT':outdir_wxt}

cali_apollo   = fst.apollo_calibration(califile=cali_file)


quicklook = False

if quicklook:
    fs  = 12 #fontsize of plot
    mpl.rcParams['font.size'] = fs
    mpl.rcParams['font.sans-serif'] = ['Tahoma']
    mpl.rcParams['axes.labelpad'] = 8
    mpl.rcParams['axes.autolimit_mode'] = 'round_numbers'
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['grid.color'] = 'black'
    mpl.rcParams['grid.linestyle'] = 'dotted'
    mpl.rcParams['grid.alpha'] = 0.25
    mpl.rcParams['xtick.labelsize'] = fs
    mpl.rcParams['ytick.labelsize'] = fs
    mpl.rcParams['legend.frameon'] = True
    mpl.rcParams['legend.handlelength'] = 1.5
    mpl.rcParams['legend.handletextpad'] = 0.5 


#----------------------------------------------------------------------------
# Data corrections 

# Removing values outside physically plausible values  
plausibility_check = True      # APOLLO and WXT
TT_min,TT_max      = 0,40      # (°C)   
PP_min,PP_max      = 950,1050  # (hPa)
RH_min,RH_max      = 0,100     # (%)
FF_min,FF_max      = 0,40      # (m/s)
DD_min,DD_max      = 0,360     # (°)
RR_min,RR_max      = 0,10      # (mm/10s)
HA_min             = 0         # (1/cm2)

# Removing start and end of station specific measurement period (may be disturbed by installation)
startend_removal   = True      # APOLLO and WXT
mins_remove        = 15        # Remove first and last xx minutes

# Removing spikes (TT: only for periods of activated Wifi)
spike_removal      = True      #
TT_diff_max        = 0.5       # (K) Max. allowed difference from running median
TT_n_run           = 15        # (s) Time window of running median filter
PP_diff_max        = 1         # (hPa) Max. allowed difference from running median
PP_n_run           = 30        # (s) Time window of running median filter

# Removing implausible values compared to network mean
space_consistency   = True       # APOLLO
TT_diff_mean        = 15         # (K) Max. allowed difference from network mean
apollo_erroneous_tt = [11,15,26] # Erroneous loggers
TT_diff_mean_err    = 2          # (K) Max. allowed difference from network mean
                                 # for erroneous instruments
t_windows_err       = 900        # (s) Time window removed around values detected
                                 # by TT_diff_mean_err                     


# Calibration of temperature sensors
calibration        = True


# Smoothing of time series for selected sensors
data_smoothing     = True      # Only APOLLO
t_smooth           = 5          # s
smooth_tt          = ['113PGa'] # Senors of single stations      


# Manually removing specific periods with erroneous measurements
periods_remove_tt   = {'093PGa':[dt.datetime(2020,8,15,11,0),dt.datetime(2020,8,15,18,0)],
                       '051OGa':[dt.datetime(2020,7,30,9,55),dt.datetime(2020,7,30,10,5)],
                      }

periods_remove_pp   = {'040UHa':[dt.datetime(2020,9,4,13,30),dt.datetime(2020,9,5,0,0)],
                      }


# periods influenced by warm air of local ventilation system
periods_remove_039  = {dt.date(2020,6,12): [dt.time(6,30),dt.time(9,0)],
                       dt.date(2020,6,15): [dt.time(5,0),dt.time(6,30)],
                       dt.date(2020,6,16): [dt.time(4,0),dt.time(5,0)],
                       dt.date(2020,6,18): [dt.time(4,0),dt.time(4,30)],
                       dt.date(2020,6,19): [dt.time(4,0),dt.time(18,30)],
                       dt.date(2020,6,22): [dt.time(4,0),dt.time(13,30)],
                       dt.date(2020,6,23): [dt.time(5,0),dt.time(16,0)],
                       dt.date(2020,6,25): [dt.time(4,0),dt.time(9,30)],
                       dt.date(2020,6,26): [dt.time(4,30),dt.time(6,30)],
                       dt.date(2020,6,29): [dt.time(4,0),dt.time(18,30)],
                       dt.date(2020,6,30): [dt.time(4,0),dt.time(18,30)],
                       dt.date(2020,7,1):  [dt.time(4,0),dt.time(6,30)],
                       dt.date(2020,7,3):  [dt.time(4,0),dt.time(8,30)],
                       dt.date(2020,7,4):  [dt.time(5,0),dt.time(15,30)],
                       dt.date(2020,7,6):  [dt.time(4,0),dt.time(16,30)],
                       dt.date(2020,7,7):  [dt.time(4,0),dt.time(18,30)],
                       dt.date(2020,7,8):  [dt.time(4,0),dt.time(18,30)],
                       dt.date(2020,7,9):  [dt.time(4,0),dt.time(18,30)],
                       dt.date(2020,7,10): [dt.time(7,0),dt.time(14,30)],
                       dt.date(2020,7,13): [dt.time(4,30),dt.time(5,30)],
                       dt.date(2020,7,14): [dt.time(4,0),dt.time(6,30)],
                       dt.date(2020,7,20): [dt.time(4,0),dt.time(6,0)],
                       dt.date(2020,7,21): [dt.time(4,0),dt.time(16,30)],
                       dt.date(2020,7,22): [dt.time(4,0),dt.time(15,0)],
                       dt.date(2020,7,23): [dt.time(4,0),dt.time(14,30)],
                       dt.date(2020,7,24): [dt.time(4,30),dt.time(13,0)],
                       dt.date(2020,7,27): [dt.time(4,0),dt.time(6,0)],
                       dt.date(2020,7,28): [dt.time(4,0),dt.time(16,0)],
                       dt.date(2020,7,29): [dt.time(4,0),dt.time(18,30)],
                       dt.date(2020,7,30): [dt.time(4,0),dt.time(13,30)],
                       dt.date(2020,8,7):  [dt.time(4,0),dt.time(5,0)],
                       dt.date(2020,8,10): [dt.time(5,30),dt.time(10,30)],
                       dt.date(2020,8,11): [dt.time(4,0),dt.time(13,30)],
                       dt.date(2020,8,12): [dt.time(4,0),dt.time(9,0)],
                       dt.date(2020,8,13): [dt.time(4,0),dt.time(12,0)],
                       dt.date(2020,8,14): [dt.time(6,30),dt.time(12,0)],
                       dt.date(2020,8,17): [dt.time(4,30),dt.time(8,30)],
                       dt.date(2020,8,18): [dt.time(5,30),dt.time(14,0)],
                       dt.date(2020,8,19): [dt.time(4,30),dt.time(9,0)],
                       dt.date(2020,8,20): [dt.time(4,30),dt.time(13,0)],
                       dt.date(2020,8,27): [dt.time(4,0),dt.time(14,0)],
                       dt.date(2020,8,28): [dt.time(4,0),dt.time(13,30)],
                       dt.date(2020,8,31): [dt.time(4,0),dt.time(9,30)],
                       }

#----------------------------------------------------------------------------
# Definition of level 2 data format

meas_type       = 'fessthh'  # Campaign
institution     = 'uhh'      # Responsible Institution
version_dataset = 0          # In case of new dataset version number difference 
                             # needs to explained in History section
version_proc    = '1.0'      # Version number of this python script                                           


def write_netcdf_file(varstr,dataframe_write,meta_write,type_instr=type_proc,
                      kkk=meas_type,sss=institution,vers_ds=version_dataset,
                      vers_proc=version_proc,t_res=time_res[type_proc],
                      odir=outdirs[type_proc],nchar=40,c_level=1):
    
    # Permission denied error = trying to create an already existing file
    # Invalid ID error = trying to access an unclosed file
    
    #c_level = Compression level (1-9)
    
    if dataframe_write.notnull().sum().sum() == 0:
        print('No '+varstr+' data found to be written!')
        return
    
    nstat    = dataframe_write.shape[1]
    vnn      = {'TT':'ta',
                'PP':'pa',
                'RH':'hur',
                'FF':'wspeed',
                'FB':'wspeed_max',
                'DD':'wdir',
                'RR':'precip',
                'HA':'hail',
                }
    
    varnames = {'TT':'Air temperature',
                'PP':'Air pressure',
                'RH':'Relative humidity of air',
                'FF':'Wind speed',
                'FB':'Wind speed maximum (gust)',
                'DD':'Wind direction',
                'RR':'Precipitation amount',
                'HA':'Hail amount',
                }
    
    if type_instr == 'APOLLO': 
        ht_str = 'HT_'+varstr 
        source = 'APOLLO (Autonomous cold POoL LOgger), '+ \
                 'logger software version 114a9'
        lsd_tt = 2
    if type_instr == 'WXT': 
        ht_str = 'HT_WXT'
        source = 'Vaisala Weather Transmitter WXT536 with external Pt1000 thermometer, '+ \
                 'logger software version 2.2b'
        lsd_tt = None         
                 
    t_ave = t_res #if varstr not in ['TT','PP','RH'] else 0     
    
    # File naming and creation
    dtime    = dataframe_write.index[0]
    writedir = odir+dtime.strftime('%Y/%m/%d/') 
    dtstr    = dtime.strftime('%Y%m%d%H%M%S')
    if os.path.isdir(writedir) == False: os.makedirs(writedir)
    filename = writedir+kkk+'_'+sss+'_'+type_instr.lower()+'_l2_'+vnn[varstr]+\
               '_v'+str(vers_ds).zfill(2)+'_'+dtstr+'.nc'
    if os.path.isfile(filename): os.remove(filename)

    ncfile = Dataset(filename,'w',format='NETCDF4')
    
    # Dimensions
    ncfile.createDimension('time',None)
    ncfile.createDimension('nv',2)
    ncfile.createDimension('station',nstat)
    ncfile.createDimension('character',None) # Only used for creation of char string
    
    #Dimension Variables
    utime              = fst.datetime_to_unixtime(dataframe_write.index)
    time               = ncfile.createVariable('time', 'i4', ('time',))
    time[:]            = utime
    time.standard_name = 'time'
    time.units         = 'seconds since 1970-01-01 00:00:00 UTC' 
    time.bounds        = 'time_bnds'
    time.calendar      = 'standard'
    
    nv                 = ncfile.createVariable('time_bnds', 'i4', ('time','nv',))
    nv[:,:]            = np.column_stack((utime-t_ave,utime))
    nv.long_name       = 'start and end of time averaging intervals'

    # Instrument meta-data        
    station_id              = ncfile.createVariable('station_id', 'S1', ('station','character',))
    station_id[:]           = stringtochar(meta_write['STATION'].to_numpy(dtype='S'+str(nchar)))
    station_id.long_name    = 'station identifier code'
    
    station_name            = ncfile.createVariable('station_name', 'S1', ('station','character',))
    station_name[:]         = stringtochar(meta_write['NAME'].to_numpy(dtype='S'+str(nchar)))
    station_name.long_name  = 'station name'

    lat                     = ncfile.createVariable('lat', 'f4', ('station',)) 
    lat[:]                  = meta_write['LAT'].values
    lat.standard_name       = 'latitude'
    lat.long_name           = 'latitude of instrument location'
    lat.units               = 'degrees_north'
    
    lon                     = ncfile.createVariable('lon', 'f4', ('station',)) 
    lon[:]                  = meta_write['LON'].values
    lon.standard_name       = 'longitude'
    lon.long_name           = 'longitude of instrument location'
    lon.units               = 'degrees_east'
    
    zsl                     = ncfile.createVariable('zsl', 'f4', ('station',)) 
    zsl[:]                  = meta_write['ALT'].values
    zsl.standard_name       = 'altitude'
    zsl.long_name           = 'altitude of instrument location above mean sea level'
    zsl.comment             = 'altitude of ground level'
    zsl.units               = 'm'
    
    zag                     = ncfile.createVariable('zag', 'f4', ('station',)) 
    zag[:]                  = meta_write[ht_str].values
    zag.standard_name       = 'height'
    zag.long_name           = 'height of sensor above ground'
    zag.units               = 'm'
     
    lcz                     = ncfile.createVariable('lcz', 'S1', ('station','character',)) 
    lcz[:]                  = stringtochar(meta_write['LCZ'].to_numpy(dtype='S'+str(nchar)))
    lcz.long_name           = 'local climate zone of station'
    lcz.comments            = 'after Stewart and Oke (2012), doi:10.1175/BAMS-D-11-00019.1'
    
    
    # Measured Variables
    #zlib=True,complevel=4,least_significant_digit=3
    if varstr == 'TT':
        ta                       = ncfile.createVariable(vnn[varstr],'f4',
                                                         ('time','station',),
                                                         fill_value='nan',
                                                         zlib=True,
                                                         complevel=c_level,
                                                         least_significant_digit=lsd_tt)
        ta[:,:]                  = fst.TC_to_TK(dataframe_write)
        ta.standard_name         = 'air_temperature'
        ta.long_name             = 'air temperature'
        ta.units                 = 'K'
        if type_instr == 'APOLLO':
            ta.comment           = 'Data of station 113PGa smoothed because of '+\
                                   'higher noise level than usual'
    
    if varstr == 'PP':
        pa                       = ncfile.createVariable(vnn[varstr],'f4',
                                                         ('time','station',),
                                                         fill_value='nan',
                                                         zlib=True,
                                                         complevel=c_level,
                                                         least_significant_digit=0)
        pa[:,:]                  = dataframe_write*100
        pa.standard_name         = 'air_pressure '
        pa.long_name             = 'air pressure'
        pa.units                 = 'Pa'
        pa.comment               = 'unnormalized air pressure'
        
    if varstr == 'RH':
        hur                      = ncfile.createVariable(vnn[varstr],'f4',
                                                         ('time','station',),
                                                         fill_value='nan',
                                                         zlib=True,
                                                         complevel=c_level)
        hur[:,:]                 = dataframe_write/100
        hur.standard_name        = 'relative_humidity'
        hur.long_name            = 'relative humidity of air'
        hur.units                = '1' 
        
    if varstr == 'FF':    
        wspeed                   = ncfile.createVariable(vnn[varstr],'f4',
                                                         ('time','station',),
                                                         fill_value='nan',
                                                         zlib=True,
                                                         complevel=c_level)
        wspeed[:,:]              = dataframe_write
        wspeed.standard_name     = 'wind_speed'
        wspeed.long_name         = 'wind speed'
        wspeed.comment           = 'average of 4-Hz data within measurement interval'
        wspeed.units             = 'm s-1'
        
    if varstr == 'FB':    
        wspeed_max               = ncfile.createVariable(vnn[varstr],'f4',
                                                         ('time','station',),
                                                         fill_value='nan',
                                                         zlib=True,
                                                         complevel=c_level)
        wspeed_max[:,:]          = dataframe_write
        wspeed_max.standard_name = 'wind_speed_of_gust'
        wspeed_max.long_name     = 'wind speed maximum (gust)'
        wspeed_max.comment       = 'maximum 3-s average of 4-Hz data within measurement interval'
        wspeed_max.units         = 'm s-1'    
        
    if varstr == 'DD':    
        wdir                     = ncfile.createVariable(vnn[varstr],'f4',
                                                         ('time','station',),
                                                         fill_value='nan',
                                                         zlib=True,
                                                         complevel=c_level)
        wdir[:,:]                = dataframe_write
        wdir.standard_name       = 'wind_from_direction'
        wdir.long_name           = 'wind direction'
        wdir.comment             = 'average of 4-Hz data within measurement interval'
        wdir.units               = 'degree' 
        
    if varstr == 'RR':    
        precip                   = ncfile.createVariable(vnn[varstr],'f4',
                                                         ('time','station',),
                                                         fill_value='nan',
                                                         zlib=True,
                                                         complevel=c_level)
        precip[:]                = dataframe_write # mm = kg m-2
        precip.standard_name     = 'precipitation_amount'
        precip.long_name         = 'amount of precipitation'
        precip.comment           = 'accumulated amount of precipitation within measurement interval'
        precip.units             = 'kg m-2'
        
    if varstr == 'HA':  
        hail                     = ncfile.createVariable(vnn[varstr],'f4',
                                                         ('time','station',),
                                                         fill_value='nan',
                                                         zlib=True,
                                                         complevel=c_level)
        hail[:]                  = dataframe_write*10**4 # cm-2 -> m-2
        hail.long_name           = 'amount of hail'
        hail.comment             = 'number of hail hits within measurement interval'
        hail.units               = 'm-2'
    
    # Global attributes
    ncfile.Title             = varnames[varstr]+' observed by '+type_instr+' station network'
    ncfile.Institution       = 'Meteorological Institute, University of Hamburg (UHH), Germany'
    ncfile.Contact_Person    = 'Prof. Dr. Felix Ament (felix.ament@uni-hamburg.de)'
    ncfile.Source            = source
    ncfile.History           = 'Data processed with fesstval_level1_to_level2.py,'+\
                               ' version '+vers_proc
    #ncfile.Dependencies      = 'external'
    ncfile.Conventions       = 'CF-1.7 where applicable'
    ncfile.Processing_date   = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')    
    ncfile.Author            = 'Bastian Kirsch (bastian.kirsch@uni-hamburg.de)'     
    ncfile.Comments          = 'FESST@HH field experiment (June - August 2020),'+\
                               ' doi:10.25592/uhhfdm.8967'
    ncfile.Licence           = 'This data is licensed under a '+\
                               'Creative Commons Attribution-ShareAlike 4.0 '+\
                               'International License (CC-BY-SA-4.0).'   
    ncfile.close()
    return


#----------------------------------------------------------------------------
print('Converting '+type_proc+' level 1 to level 2 data')
print('Start Date: '+start_date.strftime('%Y-%m-%d'))
print('End Date  : '+end_date.strftime('%Y-%m-%d'))

# Station meta data
meta_data = fst.fessthh_stations('all',include_ht=True,include_time=True,
                                 include_serial=True,include_lcz=True,
                                 metafile=meta_file,serialfile=serial_file)

ii_type        = meta_data['STATION'].str.endswith(type_proc[0].lower())  
ii_type_unique = ii_type & (meta_data['STATION'].duplicated() == False)

meta_data['NAME'] = meta_data['NAME'].str.replace('ä','ae')
meta_data['NAME'] = meta_data['NAME'].str.replace('ö','oe')
meta_data['NAME'] = meta_data['NAME'].str.replace('ü','ue')
meta_data['NAME'] = meta_data['NAME'].str.replace('Ö','Oe')
meta_data['NAME'] = meta_data['NAME'].str.replace('ß','ss')


days = pd.date_range(start=start_date,end=end_date,freq='d')
time_res = {'APOLLO':1,'WXT':10}

TT_removed   = pd.DataFrame(index=meta_data[ii_type_unique]['STATION'],columns=days)
PP_removed   = pd.DataFrame(index=meta_data[ii_type_unique]['STATION'],columns=days)
TTPT_removed = pd.DataFrame(index=meta_data[ii_type_unique]['STATION'],columns=days)
RH_removed   = pd.DataFrame(index=meta_data[ii_type_unique]['STATION'],columns=days)
FF_removed   = pd.DataFrame(index=meta_data[ii_type_unique]['STATION'],columns=days)
FB_removed   = pd.DataFrame(index=meta_data[ii_type_unique]['STATION'],columns=days)
DD_removed   = pd.DataFrame(index=meta_data[ii_type_unique]['STATION'],columns=days)
RR_removed   = pd.DataFrame(index=meta_data[ii_type_unique]['STATION'],columns=days)
HA_removed   = pd.DataFrame(index=meta_data[ii_type_unique]['STATION'],columns=days)

for iday,d in enumerate(days):
    print(' ')
    print(d.strftime('%Y-%m-%d'))
    
    #--------------------
    print('Reading level 1 data')
    
    time = pd.date_range(start=d,periods=86400/time_res[type_proc],
                         freq=str(time_res[type_proc])+'s')        
    
    TT_data   = pd.DataFrame(index=time)
    PP_data   = pd.DataFrame(index=time)
    
    if type_proc == 'APOLLO':
        WIFI_data = pd.DataFrame(index=time)
    
    if type_proc == 'WXT':
        TTPT_data = pd.DataFrame(index=time)
        RH_data   = pd.DataFrame(index=time)
        FF_data   = pd.DataFrame(index=time)
        FB_data   = pd.DataFrame(index=time)
        DD_data   = pd.DataFrame(index=time)
        RR_data   = pd.DataFrame(index=time)
        HA_data   = pd.DataFrame(index=time)
        
    stat_instr    = {} # Daily instrument numbers per station
    stat_startend = {} # Daily start and end index per station (only if
                       # measurment period starts or ends)
    
    for stat in meta_data[ii_type_unique]['STATION']:        
        # Merging data of different instruments at same station
        sdata_list = []
        stat_instr[stat] = []
        stat_startend[stat] = [np.nan,np.nan]
        
        for ii in meta_data[(meta_data['STATION'] == stat)].index:
            stat_start = meta_data.loc[ii]['START']
            stat_end   = meta_data.loc[ii]['END']
            #stat_name  = meta_data.loc[ii]['NAME']
            
            if (d.date() >= stat_start.date()) and (d.date() <= stat_end.date()): 
                s = meta_data.loc[ii][type_proc]
                stat_instr[stat].append(s)
                if d.date() == stat_start.date(): 
                    stat_startend[stat][0] = np.where(time == stat_start)[0][0]
                if d.date() == stat_end.date(): 
                    stat_startend[stat][1] = np.where(time == stat_end)[0][0]    
                
                if stat not in wxt_permanent.keys():
                    sdata = fst.read_fesstval_level1(type_proc,s,d.year,d.month,d.day,
                                                     datadir_apollo=datadir_apo,
                                                     datadir_wxt=datadir_wxt,
                                                     include_monitoring_data=True,
                                                     mute=True)
                else:
                    sdata = fst.read_fesstval_wettermast(wxt_permanent[stat],
                                                         d.year,d.month,d.day,
                                                         datadir=datadir_wxt_perm)

            else:
                sdata = pd.DataFrame()
            
            sdata_list.append(sdata)
        
        stat_data = pd.concat(sdata_list)
        empty = False if len(stat_data.columns) > 0 else True
        
        TT_data[stat] = stat_data['TT'].reindex(time) if not empty else np.nan 
        PP_data[stat] = stat_data['PP'].reindex(time) if not empty else np.nan
        
        if type_proc == 'APOLLO':
            WIFI_data[stat] = stat_data['WIFI'].reindex(time) if not empty else np.nan
        
        if type_proc == 'WXT':
            TTPT_data[stat] = stat_data['TT_PT1000'].reindex(time) if not empty else np.nan
            RH_data[stat]   = stat_data['RH'].reindex(time) if not empty else np.nan
            FF_data[stat]   = stat_data['FF'].reindex(time) if not empty else np.nan
            FB_data[stat]   = stat_data['FB'].reindex(time) if not empty else np.nan
            DD_data[stat]   = stat_data['DD'].reindex(time) if not empty else np.nan
            RR_data[stat]   = stat_data['RR'].reindex(time) if not empty else np.nan
            HA_data[stat]   = stat_data['HA'].reindex(time) if not empty else np.nan
    
    
    instr_day = np.array([stat_instr[s][-1] if len(stat_instr[s]) > 0 else -1 \
                         for s in meta_data[ii_type_unique]['STATION'].values])
    meta_data[type_proc+'_DAY'] = -1    
    pd.set_option('mode.chained_assignment',None)
    meta_data[type_proc+'_DAY'].loc[ii_type_unique] = instr_day    
    pd.set_option('mode.chained_assignment','warn')
    
    #--------------------        
    print('Correcting data')
    
    TT_nan_pre = TT_data.isnull().sum()
    PP_nan_pre = PP_data.isnull().sum()
    if type_proc == 'WXT':
        TTPT_nan_pre = TTPT_data.isnull().sum()
        RH_nan_pre   = RH_data.isnull().sum()
        FF_nan_pre   = FF_data.isnull().sum()
        FB_nan_pre   = FB_data.isnull().sum()
        DD_nan_pre   = DD_data.isnull().sum()
        RR_nan_pre   = RR_data.isnull().sum()
        HA_nan_pre   = HA_data.isnull().sum()
    
    if plausibility_check:
        TT_data.mask((TT_data < TT_min) | (TT_data > TT_max),inplace=True)
        PP_data.mask((PP_data < PP_min) | (PP_data > PP_max),inplace=True)
        if type_proc == 'WXT': 
            TTPT_data.mask((TTPT_data < TT_min) | (TTPT_data > TT_max),inplace=True)
            RH_data.mask((RH_data < RH_min) | (RH_data > RH_max),inplace=True)
            FF_data.mask((FF_data < FF_min) | (FF_data > FF_max),inplace=True)
            FB_data.mask((FB_data < FF_min) | (FB_data > FF_max),inplace=True)
            DD_data.mask((DD_data < DD_min) | (DD_data > DD_max),inplace=True)
            RR_data.mask((RR_data < RR_min) | (RR_data > RR_max),inplace=True)
            HA_data.mask((HA_data < HA_min),inplace=True)                       
            
    if startend_removal:
        i_mins_remove = int(mins_remove*(60/time_res[type_proc]))
        for col in meta_data[ii_type_unique]['STATION']:
            i_start,i_end = stat_startend[col]
            ii_start = np.zeros_like(time,dtype=bool)
            ii_end   = np.zeros_like(time,dtype=bool)
            if np.isfinite(i_start):
                ii_start = (time >= time[i_start]) & (time < time[i_start+i_mins_remove])
                # [i_start:i_start+i_mins_remove]
            if np.isfinite(i_end):
                ii_end = (time > time[i_end-i_mins_remove]) & (time <= time[i_end])
                #[i_end-i_mins_remove:i_end]
                
            ii_startend = ii_end | ii_start    
            
            TT_data[col].mask(ii_startend,inplace=True)
            PP_data[col].mask(ii_startend,inplace=True)
            if type_proc == 'WXT':
                TTPT_data[col].mask(ii_startend,inplace=True)
                RH_data[col].mask(ii_startend,inplace=True)
                FF_data[col].mask(ii_startend,inplace=True)
                FB_data[col].mask(ii_startend,inplace=True)
                DD_data[col].mask(ii_startend,inplace=True)
                RR_data[col].mask(ii_startend,inplace=True)
                HA_data[col].mask(ii_startend,inplace=True)

                
    if spike_removal & (type_proc == 'APOLLO'):        
        WIFI_sum = WIFI_data.sum()
        #ii_spike_sum = 0
        for col in WIFI_sum[WIFI_sum > 0].index:
            TT_run_median = TT_data[col].rolling(TT_n_run,center=True,min_periods=3).median()
            ii_spike = ((TT_run_median - TT_data[col]).abs() > TT_diff_max) & \
                       (WIFI_data[col] == 1)
            TT_data[col].mask(ii_spike,inplace=True) 
            #ii_spike_sum += ii_spike.sum()
        
        PP_spikes = (PP_data.diff().abs() > PP_diff_max).sum()
        for col in PP_spikes[PP_spikes > 0].index:
            PP_run_median = PP_data[col].rolling(PP_n_run,center=True,min_periods=6).median()
            ii_spike = ((PP_run_median - PP_data[col]).abs() > PP_diff_max)
            PP_data[col].mask(ii_spike,inplace=True) 
            
        #PP_data.mask(PP_data.diff().abs() > PP_diff_max,inplace=True)    
        
    if space_consistency & (type_proc == 'APOLLO'):
        TT_mean = TT_data.mean(axis=1)
        for col in meta_data[ii_type_unique]['STATION']:
            TT_data[col].mask((TT_data[col]-TT_mean).abs()>TT_diff_mean,inplace=True)  
            instr_err = set(apollo_erroneous_tt).intersection(stat_instr[col])
            if len(instr_err) > 0:
                mask = TT_data[col].isnull()
                ii = np.where((TT_data[col]-TT_mean).abs()>TT_diff_mean_err)[0]
                for i in ii: mask[i-t_windows_err:i+t_windows_err] = True
                TT_data[col].mask(mask,inplace=True)
                
                
    if calibration & (type_proc == 'APOLLO'):
        for col in meta_data[ii_type_unique]['STATION']:
            n_instr = len(stat_instr[col])
            # Loop over instruments per station (mostly 1)
            for i_apo,apo in enumerate(stat_instr[col]):
                cali = cali_apollo[apo]
                if np.isnan(cali):
                    print('*** No calibration constant found for APOLLO '+str(apo)+' ***')
                    continue   
                if n_instr == 1: 
                    TT_data[col] -= cali
                if (n_instr == 2) & (i_apo == 0):
                    TT_data[col].iloc[:stat_startend[col][1]+1] -= cali
                if (n_instr == 2) & (i_apo == 1):
                    TT_data[col].iloc[stat_startend[col][0]:] -= cali    
            
            
    if len(periods_remove_tt.keys()) > 0:
        for col in periods_remove_tt.keys():
            if col in meta_data[ii_type_unique]['STATION'].to_list():
                ii1,ii2 = periods_remove_tt[col]
                TT_data[col].mask((TT_data.index>=ii1)&((TT_data.index<=ii2)),inplace=True)
     
    if len(periods_remove_pp.keys()) > 0:
        for col in periods_remove_pp.keys():
            if col in meta_data[ii_type_unique]['STATION'].to_list():
                ii1,ii2 = periods_remove_pp[col]
                PP_data[col].mask((PP_data.index>=ii1)&((PP_data.index<=ii2)),inplace=True)
    
    # Single data at WXT Holm            
    if ('047OGw' in meta_data[ii_type_unique]['STATION'].to_list()) & (d == dt.date(2020,6,20)):
        col = '047OGw'
        TT_data[col].mask(TT_data[col].notnull(),inplace=True)
        PP_data[col].mask(PP_data[col].notnull(),inplace=True)
        TTPT_data[col].mask(TTPT_data[col].notnull(),inplace=True)
        RH_data[col].mask(RH_data[col].notnull(),inplace=True)
        FF_data[col].mask(FF_data[col].notnull(),inplace=True)
        FB_data[col].mask(FB_data[col].notnull(),inplace=True)
        DD_data[col].mask(DD_data[col].notnull(),inplace=True)
        RR_data[col].mask(RR_data[col].notnull(),inplace=True)
        HA_data[col].mask(HA_data[col].notnull(),inplace=True)    
                
    if '039UHa' in meta_data[ii_type_unique]['STATION'].to_list():
        for dd in periods_remove_039.keys():
            if dd == d:
                ii1,ii2 = periods_remove_039[dd]
                TT_data['039UHa'].mask((TT_data.index.time>=ii1)&\
                                       ((TT_data.index.time<=ii2)),inplace=True)                     
                    
        
    if data_smoothing:
        for col in smooth_tt:
            if col in meta_data[ii_type_unique]['STATION'].to_list():
                TT_data[col] = TT_data[col].rolling(2*t_smooth+1,center=True,
                                                    min_periods=t_smooth,axis=0).mean()
    #--------------------        
    TT_removed[d] = TT_data.isnull().sum() - TT_nan_pre   
    PP_removed[d] = PP_data.isnull().sum() - PP_nan_pre
    
    if type_proc == 'WXT':
        TTPT_removed[d] = TTPT_data.isnull().sum() - TTPT_nan_pre
        RH_removed[d]   = RH_data.isnull().sum() - RH_nan_pre
        FF_removed[d]   = FF_data.isnull().sum() - FF_nan_pre
        FB_removed[d]   = FB_data.isnull().sum() - FB_nan_pre
        DD_removed[d]   = DD_data.isnull().sum() - DD_nan_pre
        RR_removed[d]   = RR_data.isnull().sum() - RR_nan_pre
        HA_removed[d]   = HA_data.isnull().sum() - HA_nan_pre
               
        
            
    #--------------------
    if write_data:
        print('Writing level 2 data')
        if type_proc == 'APOLLO':
            write_netcdf_file('TT',TT_data,meta_data[ii_type_unique])
            write_netcdf_file('PP',PP_data,meta_data[ii_type_unique])
        if type_proc == 'WXT':
            write_netcdf_file('TT',TTPT_data,meta_data[ii_type_unique])
            write_netcdf_file('PP',PP_data,meta_data[ii_type_unique]) 
            write_netcdf_file('RH',RH_data,meta_data[ii_type_unique]) 
            write_netcdf_file('FF',FF_data,meta_data[ii_type_unique]) 
            write_netcdf_file('FB',FB_data,meta_data[ii_type_unique])
            write_netcdf_file('DD',DD_data,meta_data[ii_type_unique]) 
            write_netcdf_file('RR',RR_data,meta_data[ii_type_unique]) 
            write_netcdf_file('HA',HA_data,meta_data[ii_type_unique])
           
            
    if quicklook:
        
        fig,ax = plt.subplots(2,1,figsize=(9,8),dpi=150)
    
        ax[0].plot(TT_data,linewidth=0.7)
        ax[0].set_title(d.strftime('APOLLO Level 2, %d %b %Y'))
        ax[0].set_ylabel('Temperature (°C)')
        
        ax[1].plot(PP_data,linewidth=0.5)
        ax[1].set_ylabel('Pressure (hPa)')
        
        for x in [0,1]:
            ax[x].xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
                
        plotname = plotdir+type_proc+d.strftime('_level2_%Y%m%d')+'.png'
        fig.savefig(plotname,bbox_inches='tight')
        plt.close()    
        print('Plot done!') 


def print_statistics(removed_data,var):
    n_total = len(days) * (86400/time_res[type_proc]) * ii_type_unique.sum()
    n_removed = removed_data.sum().sum()
    frac_removed = (n_removed/n_total) * 100
    
    if len(days) > 1:
        n_stat_removed = ((removed_data > 0).sum(axis=1) > 0).sum()
    else:
        n_stat_removed = (removed_data > 0).sum().sum()
    
    print(var+': {:1.0f}'.format(n_removed)+' ({:4.2f} %) at '.format(frac_removed)+\
          '{:1.0f} station(s)'.format(n_stat_removed)) 
    return

print(' ')
print('*********')
print('Total number of single measurements removed:')
print_statistics(TT_removed,'TT')
print_statistics(PP_removed,'PP')           
    
if type_proc == 'WXT':
    print_statistics(TTPT_removed,'TTPT')
    print_statistics(RH_removed,'RH')
    print_statistics(FF_removed,'FF')
    print_statistics(FB_removed,'FB')
    print_statistics(DD_removed,'DD')
    print_statistics(RR_removed,'RR')
    print_statistics(HA_removed,'HA')
  
#----------------------------------------------------------------------------
print('*********')
fst.print_runtime(t_run)