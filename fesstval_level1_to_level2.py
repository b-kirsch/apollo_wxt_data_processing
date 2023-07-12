# -*- coding: utf-8 -*-
"""
@author: Bastian Kirsch (bastian.kirsch@uni-hamburg.de)

Script to convert level1 APOLLO (Autonomous cold POoL LOgger) and 
WXT data to level2 daily .nc-files complying with the SAMD data standard

Data structure needed in data_dir: ./<YYYY>/<MM>/<DD>/<APOLLO,WXT>-<SSS,SS>_<YYYYMMDD>.txt
Output data structure in out_dir: ./<YYYY>/<MM>/<DD>/fval_uhh_<instr>_l2_<var>_v00_<YYYYMMDDhhmmss>.nc

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
- stations_fesstval.txt
- APOLLO_serials.txt  
- APOLLO_calibration.txt  
    
Last updated: 21 April 2022
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
maindir     = '.'
datadir_apo = maindir+'APOLLO_data/level1/'
datadir_wxt = maindir+'WXT_data/level1/'
outdir_apo  = maindir+'APOLLO_data/level2/' 
outdir_wxt  = maindir+'WXT_data/level2/' 
plotdir     = maindir+'Cold-Pools/Plots/FESSTVaL2021_Quicklooks/'
logdir      = maindir+'FESSTVaL/'
meta_file   = maindir+'FESSTVaL/stations_fesstval.txt'
serial_file = maindir+'FESSTVaL/APOLLO/apollo_serials.txt'
cali_file   = maindir+'FESSTVaL/APOLLO/apollo_calibration.txt'

start_date = dt.date(2021,5,17)
end_date   = dt.date(2021,8,27)

#type_proc = 'APOLLO'
type_proc  = 'WXT'

write_data = False
quicklook  = False
log_stats  = False

#----------------------------------------------------------------------------
# Processing settings
time_res      = {'APOLLO':1, 
                 'WXT':10}
outdirs       = {'APOLLO':outdir_apo, 
                 'WXT':outdir_wxt}
cali_apollo   = fst.apollo_calibration(califile=cali_file)


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
apollo_erroneous_tt = [3,9]      # Erroneous loggers
TT_diff_mean_err    = 2          # (K) Max. allowed difference from network mean
                                 # for erroneous instruments
t_windows_err       = 900        # (s) Time window removed around values detected
                                 # by TT_diff_mean_err 
                                 
# Removing constant values (only PP)  
constant_removal    = True       # APOLLO 
n_unique_pp         = 10         # Minimum number of unique values per day
t_windows_const     = 10         # (s) Time window for calculating rolling var  
                         

# Calibration of temperature sensors
calibration        = True


# Smoothing of time series for selected sensors
data_smoothing     = False      # Only APOLLO
t_smooth           = 5          # s
smooth_tt          = [''] # Senors of single stations      


# Manually removing specific periods with erroneous measurements
periods_remove_tt   = {'028Ea':[dt.datetime(2021,7,10,16,30),dt.datetime(2021,7,11,13,50)],
                       '051Fa':[dt.datetime(2021,6,24,15,30),dt.datetime(2021,6,24,21,10)],
                      }

periods_remove_pp   = {'050Ea':[dt.datetime(2021,6,16,10,50),dt.datetime(2021,6,16,11,10)],
                       '028Ea':[dt.datetime(2021,7,10,16,30),dt.datetime(2021,7,11,13,50)],
                      }

#----------------------------------------------------------------------------
# Definition of level 2 data format

meas_type       = 'fval'     # Campaign
institution     = 'uhh'      # Responsible Institution
version_dataset = 0          # In case of new dataset version number difference 
                             # needs to explained in History section
version_proc    = '2.0'      # Version number of this python script
                             # 1.0 = fessthh, 2.0=fesstval                                           


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
                'HA':'Number of hail hits',
                }
    
    if type_instr == 'APOLLO': 
        ht_str = 'HT_'+varstr 
        source = 'APOLLO (Autonomous cold POoL LOgger), '+ \
                 'logger software version 115'
        lsd_tt = 2
    if type_instr == 'WXT': 
        ht_str = 'HT_WXT'
        source = 'Vaisala Weather Transmitter WXT536 with external Pt1000 thermometer, '+ \
                 'logger software version 2.2c'
        lsd_tt = None         
                 
    t_ave = t_res #if varstr not in ['TT','PP','RH'] else 0     
    
    # File naming and creation
    dtime    = dataframe_write.index[0]
    writedir = odir+dtime.strftime('%Y/%m/')  #'%Y/%m/%d/'
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
            ta.comment           = 'No calibration constants applied for '+\
                                   'stations 072Ga, 075Ha and 092Ha' 
    
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
        hail.long_name           = 'number of hail hits'
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
    ncfile.Comments          = 'FESSTVaL field experiment (May - August 2021),'+\
                               ' doi:10.25592/uhhfdm.9766'
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
meta_data = fst.fesstval_stations('all',metafile=meta_file,serialfile=serial_file,
                                  include_time=True,include_serial=True,
                                  include_lcz=True)

ii_type        = meta_data['STATION'].str.endswith(type_proc[0].lower())  
ii_type_unique = ii_type & (meta_data['STATION'].duplicated() == False)


if type_proc == 'APOLLO':
    apo_miss_cali = []
    for apo in meta_data[ii_type_unique]['APOLLO']:
        if not np.isfinite(cali_apollo[apo]): apo_miss_cali.append(apo)
    apo_miss_cali.sort()
    # APOLLO  81 = 072Ga
    # APOLLO  89 = 075Ha
    # APOLLO 100 = 092Ha 

meta_data['NAME'] = meta_data['NAME'].str.replace('ö','oe')
meta_data['NAME'] = meta_data['NAME'].str.replace('ü','ue')
meta_data['NAME'] = meta_data['NAME'].str.replace('ß','ss')

if type_proc == 'APOLLO':
    meta_data['HT_TT'] = 3 #(m)
    meta_data['HT_PP'] = 2 #(m)
    
if type_proc == 'WXT': 
    meta_data['HT_WXT'] = 3 #(m)

days = pd.date_range(start=start_date,end=end_date,freq='d')

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
    print('Reading level 1 '+type_proc+' data')
    
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
                
                sdata = fst.read_fesstval_level1(type_proc,s,d.year,d.month,d.day,
                                                 datadir_apollo=datadir_apo,
                                                 datadir_wxt=datadir_wxt,
                                                 include_monitoring_data=True,
                                                 mute=True)
                
                if d.date() == stat_start.date(): 
                    stat_startend[stat][0] = np.where(time == stat_start)[0][0]
                    # Needed for APOLLO 9 on 12 Aug 2021 (moved from 053Fa to 047Fa)
                    sdata = sdata.loc[stat_start:dt.datetime(d.year,d.month,d.day,23,59,59)]
                    
                if d.date() == stat_end.date(): 
                    stat_startend[stat][1] = np.where(time == stat_end)[0][0]  
                    # ...
                    sdata = sdata.loc[dt.datetime(d.year,d.month,d.day,0,0,0):stat_end]
                
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
                
    if constant_removal & (type_proc == 'APOLLO'):
        ii_unique = (PP_data.nunique() >= 1) & (PP_data.nunique() < n_unique_pp)
        for col in ii_unique[ii_unique].index: PP_data[col] = np.nan
        
        ii_const = PP_data.rolling(t_windows_const,min_periods=5).var().round(6) == 0
        PP_data.mask(ii_const,inplace=True)
                
    if calibration & (type_proc == 'APOLLO'):
        for col in meta_data[ii_type_unique]['STATION']:
            n_instr = len(stat_instr[col])
            # Loop over instruments per station (mostly 1)
            for i_apo,apo in enumerate(stat_instr[col]):
                cali = cali_apollo[apo]
                if np.isnan(cali):
                    #print('*** No calibration constant found for APOLLO '+str(apo)+' ***')
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
           
            
    if quicklook & (type_proc == 'APOLLO'):
        
        stat_id = TT_data.columns
        
        cmap = mpl.cm.tab10#nipy_spectral
        norm = mpl.colors.BoundaryNorm(np.arange(11),cmap.N,clip=True)
    
        
        fig,ax = plt.subplots(2,1,figsize=(7,6),dpi=150)
        
        for i in range(len(stat_id)):
            c = cmap(norm(int(stat_id[i][:2]))) 
            ax[0].plot(TT_data[stat_id[i]],linewidth=0.7,color=c,label=stat_id[i])
            ax[1].plot(PP_data[stat_id[i]],linewidth=0.7,color=c)
            
        ax[0].set_title(d.strftime('APOLLO Level 2 data, %d %B %Y'),loc='left',
                        fontsize=12)
        ax[0].set_ylabel('Temperature (°C)',fontsize=12)
        ax[0].legend(fontsize=8,ncol=3,bbox_to_anchor=(1.05, 1), loc='upper left')
        ax[1].set_ylabel('Pressure (hPa)',fontsize=12)
        ax[1].set_xlabel('Time (UTC)',fontsize=12)
        
        
        for x in [0,1]:
            ax[x].xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
            ax[x].set_xlim([TT_data.index[0],TT_data.index[-1]])
        
        plotname = plotdir+d.strftime('apollo_level2_%Y%m%d.png')
        fig.savefig(plotname,bbox_inches='tight')
        plt.close()    
        print('Plot done!') 
        
    if quicklook & (type_proc == 'WXT'):
        
        stat_id = TT_data.columns
        
        cmap = mpl.cm.tab20
        norm = mpl.colors.BoundaryNorm(np.arange(19),cmap.N,clip=True)
        
        lw = 0.8
    
        fig,ax = plt.subplots(3,2,figsize=(12,8),dpi=150)
        
        for i in range(len(stat_id)):
            c = cmap(norm(i)) 
            ax[0,0].plot(TT_data[stat_id[i]],linewidth=lw,color=c)
            ax[1,0].plot(PP_data[stat_id[i]],linewidth=lw,color=c)
            ax[2,0].plot(RH_data[stat_id[i]],linewidth=lw,color=c)
            ax[0,1].plot(FF_data[stat_id[i]],linewidth=lw,color=c,label=stat_id[i])
            ax[1,1].plot(FB_data[stat_id[i]],linewidth=lw,color=c)
            ax[2,1].plot(RR_data[stat_id[i]],linewidth=lw,color=c)
            
        ax[0,0].set_title(d.strftime('WXT Level 2 data, %d %B %Y'),loc='left',
                          fontsize=12)
        ax[0,0].set_ylabel('Temperature (°C)',fontsize=12)
        ax[1,0].set_ylabel('Pressure (hPa)',fontsize=12)
        ax[2,0].set_ylabel('Relative Humidity (%)',fontsize=12)
        ax[0,1].set_ylabel('Wind Speed (m/s)',fontsize=12)
        ax[1,1].set_ylabel('Wind Gust (m/s)',fontsize=12)
        ax[2,1].set_ylabel('Rain Amount (mm)',fontsize=12)
        
        ax[0,1].legend(fontsize=8,ncol=2,bbox_to_anchor=(1.05, 1), loc='upper left')
        
        for i in [0,1]: ax[2,i].set_xlabel('Time (UTC)',fontsize=12)
        for i in [0,1,2]: ax[i,1].set_ylim([0,ax[i,1].get_ylim()[1]])
        
        for x in [0,1]:
            for y in [0,1,2]:
                ax[y,x].xaxis.set_major_formatter(mpl.dates.DateFormatter('%H'))
                ax[y,x].set_xlim([TT_data.index[0],TT_data.index[-1]])
        
        plt.tight_layout()
        plotname = plotdir+d.strftime('wxt_level2_%Y%m%d.png')
        fig.savefig(plotname,bbox_inches='tight')
        plt.close()    
        print('Plot done')     


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
if type_proc == 'APOLLO':
    print_statistics(TT_removed,'TT')
    print_statistics(PP_removed,'PP')           
    
if type_proc == 'WXT':
    print_statistics(TTPT_removed,'TTPT')
    print_statistics(PP_removed,'PP')
    print_statistics(RH_removed,'RH')
    print_statistics(FF_removed,'FF')
    print_statistics(FB_removed,'FB')
    print_statistics(DD_removed,'DD')
    print_statistics(RR_removed,'RR')
    print_statistics(HA_removed,'HA')
    
if log_stats:
    print(' ')
    print('Logging statistics of removed measurements to file')   
    
    
    def log_statistics(removed_data,var,type_data=type_proc,log=logdir,base=maindir):

        filename = log+'data_removed_fesstval_level2_'+type_data.lower()+'_'+var.lower()+'.txt'
        
        header1 = 'FESSTVaL '+type_data+' level 2 removed single '+var+' measurements,'+\
                  ' data base = '+base+dt.datetime.now().strftime(' (%Y-%m-%d %H:%M LT)')
        header2 = 'DATE'
        for s in removed_data.index: header2 = header2+';'+s 
        
        wfile = open(filename,'w')
        wfile.write(header1+'\n')
        wfile.write(header2+'\n')
        wfile.close()
        
        removed_data.transpose().to_csv(filename,mode='a',index=True,
                                        header=False,sep=';',na_rep='',
                                        date_format='%Y%m%d')
        return
    
    if type_proc == 'APOLLO':
        log_statistics(TT_removed,'TT')
        log_statistics(PP_removed,'PP')           
        
    if type_proc == 'WXT':
        log_statistics(TTPT_removed,'TTPT')
        log_statistics(PP_removed,'PP')
        log_statistics(RH_removed,'RH')
        log_statistics(FF_removed,'FF')
        log_statistics(FB_removed,'FB')
        log_statistics(DD_removed,'DD')
        log_statistics(RR_removed,'RR')
        log_statistics(HA_removed,'HA')
    
   
#----------------------------------------------------------------------------
print(' ')
print('*** Finshed! ***')
fst.print_runtime(t_run)