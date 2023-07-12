# -*- coding: utf-8 -*-
"""
@author: Bastian Kirsch (bastian.kirsch@uni-hamburg.de)

Script to convert raw (=level0) APOLLO (Autonomous cold POoL LOgger) data 
(Software Version 114a9, FESST@HH 2020) to level1 daily .txt-files

Data structure needed in datadir: ./<SSS>_<first 4 characters of serial number>/<YYMMDDHHMM>.csv
Output data structure in outdir: ./<YYYY>/<MM>/<DD>/APOLLO-<SSS>_<YYYYMMDD>.txt

Dependences on non-standard software:
- fesstval_routines.py 

Required meta data files:
- stations_fessthh.txt
- APOLLO_serials.txt

Last updated: 16 October 2020
"""

import numpy as np
import pandas as pd
import datetime as dt
import os
import glob
import fesstval_routines as fst

t_run = dt.datetime.now()

#----------------------------------------------------------------------------
# Paths and basic settings
maindir     = '.'
datadir     = maindir+'APOLLO_data/level0/'
outdir      = maindir+'APOLLO_data/level1/'
meta_file   = maindir+'FESSTVaL/FESSTHH/stations_fessthh.txt'
serial_file = maindir+'FESSTVaL/APOLLO/apollo_serials.txt'

overwrite_data  = False                  # Delete all existing files prior to writing !!!

apollo_start = 1   # 1
apollo_end   = 110 # 110

fessthh_start = dt.datetime(2020,5,28,0,0,0)
fessthh_end   = dt.datetime(2020,9,4,23,59,59)

# Setting station coordinates outside of FESST@HH period 
only_fessthh = True
lat_set,lon_set,alt_set = 53.56809,9.97489, 16 #Geomatikum 18th floor
ht_tt_set,ht_pp_set     = 72,72

start_time_tolerance = 1  # (h) Max. allowed time difference between start time
                          # of file and official start time of measurements
                          # for files to be processed (to allow for loggers
                          # started before final installation)

# Time drift correction
t_drift_tolerance  = 86400  # (s) Max. time drift between logger time and GPS time considered as valid

#----------------------------------------------------------------------------
print('Processing APOLLO Level 0 data')
print('APOLLO Start:',apollo_start)
print('APOLLO End:',apollo_end)
if overwrite_data:
    print('Overwriting is enabled')
else:
    print('Overwriting is disabled')

filename_len = len('ssssss-yymmddhhmm.csv')
header3      = '#TIME(UTC);T_NTC(degC);P_BME(hPa);GPS;WIFI;LORA;TTN'

def wdata_tp(tp):
    return '{:2.2f}'.format(tp) if pd.notnull(tp) else '' 

def wdata_s(s):
    try:
        return '{:d}'.format(int(s)) 
    except ValueError:
        return ''

protected_files = [] # Needed to overwrite existing files only once

for apo in range(apollo_start,apollo_end+1):
    serial  = fst.apollo_serial(apo,serialfile=serial_file)
    station_info = fst.fessthh_stations(apo,find='APOLLO',metafile=meta_file,
                                        serialfile=serial_file,
                                        include_ht=True,include_time=True)
    
    if not station_info.empty:
        station_name  = station_info['NAME'].values[-1]
        station_lat   = station_info['LAT'].values[-1]
        station_lon   = station_info['LON'].values[-1]
        station_alt   = station_info['ALT'].values[-1]
        station_ht_tt = station_info['HT_TT'].values[-1]
        station_ht_pp = station_info['HT_PP'].values[-1]
        station_start = pd.to_datetime(station_info['START'].values)[0]
        station_end   = pd.to_datetime(station_info['END'].values)[0]
        if only_fessthh & (station_end.year == 2099): station_end = fessthh_end
    else:
        continue

    print('')
    print('APOLLO '+str(apo)+' ('+serial[:6]+') '+station_name)
    
    if len(station_info.index) > 1:
        print('*** Multiple station information found for this device ***')
        continue
    
    apollo_dir = datadir
        
    if os.path.isdir(apollo_dir) == False:
        print('Directory '+apollo_dir+' does not exist!')
        continue
    
    file_list_all = os.listdir(apollo_dir) 
    
    if len(file_list_all) == 0:
        print('No data files to process available')
        continue 
    
    # Rename files missing the leading serial-string in its name (produced by
    # logger software versions <114a9)
    for f in file_list_all:
        if len(f) != filename_len: fst.rename_apollo_file(apollo_dir+f)
    # Save new filelist in case of renamed files
    file_list_all = os.listdir(apollo_dir)  
    
    # Check which files are already processed in case existing data should
    # not be overwritten    
    processed_data = []
    for f in glob.iglob(outdir+'**/APOLLO-'+str(apo).zfill(3)+'*.txt',recursive=True):
        processed_data.append(f.replace('\\','/'))
    processed_data.sort()    
        
    if (len(processed_data) > 0) & (overwrite_data == False):    
        last_processed_datetime = fst.last_datetime_of_file(processed_data[-1])
    else:
        last_processed_datetime = station_start   
    
    # Remove all pre-existing level1 files prior to processing    
    if overwrite_data: 
        for p in processed_data: os.remove(p)      
        
    file_list = [] 
    
    for f in file_list_all:
        # Sort out invalid files or filenames
        if os.path.getsize(apollo_dir+f) == 0: continue
        if len(f) != filename_len: continue
        # Select file dates to be processed
        if f.startswith(serial[:6]) == False: continue
        ftime = pd.to_datetime(f.split('-')[-1][:10],format='%y%m%d%H%M')
        if (overwrite_data == False) & (ftime < last_processed_datetime): continue
        if only_fessthh & (ftime < station_start-dt.timedelta(hours=start_time_tolerance)): continue
        if only_fessthh & (ftime > station_end): continue
        if ftime > dt.datetime.now(): continue
        file_list.append(f)
                
    file_list.sort()
    nfiles = len(file_list)
    if nfiles == 0: print('No valid files to process found!')
    
    for nf,afile in enumerate(file_list):
        filename = apollo_dir+afile
        filesize = os.path.getsize(filename)
        if filesize > 1e5:
            sizestr = ' ({:3.1f}'.format(filesize/(1024**2))+' MB)'
        else :
            sizestr = ' ({:3.1f}'.format(filesize/(1024**1))+' KB)'
        print('File '+str(nf+1)+' of '+str(nfiles)+': '+afile+sizestr)
        
        if filename.endswith('0000.csv'):
            print('*** File probably has a wrong time stamp '+\
                 '(starts at 00:00 UTC) ***')
            continue
        
        # Reading data
        data = fst.read_apollo_level0(filename)
        ntime = len(data.index)
        if data.empty: 
            print('*** Data file '+afile+' is empty ***')
            continue
        
        serial_data = data.serial
        status_data = data.status
        lat_data    = data.lat
        lon_data    = data.lon
        fw_data     = data.firmware
        
        if serial[:6] != serial_data:
            print('*** Serial in file does not match nominal serial of APOLLO ***')
            continue
        
        if (data['DTLOG'][0].date() < dt.date(2020,1,1)) or \
            (data['DTLOG'][ntime-1].date() > dt.datetime.now().date()):
            print('*** Data time starts before 2020 or ends in the future ***')
            continue
            
        dt_corr = fst.correct_time_drift(pd.to_datetime(data['DTLOG'].values),
                                         pd.to_datetime(data['DTGPS'].values),
                                         data['GPS'].values,
                                         t_drift_tolerance)
        if len(dt_corr) == 0: continue
        
        data['DTCORR'] = dt_corr
        ii_dupli = data['DTCORR'].duplicated()
        n_dupli = ii_dupli.sum()
        data = data[~ii_dupli]
        
        if n_dupli >= 60:
            print('*** '+str(n_dupli)+' duplicate timesteps removed ***')
        
        dt_grid = pd.date_range(start=pd.to_datetime(data['DTCORR'].values[0]),
                                end=pd.to_datetime(data['DTCORR'].values[-1]),
                                freq='s')
        data.set_index('DTCORR',inplace=True)
        data = data.reindex(dt_grid)
        
        # Drop data before start or after end of meausurements
        ii_invalid = (data.index < station_start) | (data.index > station_end)
        if ii_invalid.sum() > 0: data.drop(data.index[ii_invalid],inplace=True)
        
        if len(data.index) == 0: 
            print('*** No valid data to process found ***')
            continue
        
        days = pd.date_range(start=data.index[0].date(),end=data.index[-1].date())
        
        # Write the data into separate files for each day
        for d in days:
            print('Writing data for '+d.strftime('%Y-%m-%d'))
            
            writedir = outdir+d.strftime('%Y/%m/%d/') 
            if os.path.isdir(writedir) == False: os.makedirs(writedir)
            
            writefile = writedir+'APOLLO-'+str(apo).zfill(3)+'_'+d.strftime('%Y%m%d')+'.txt'
            if os.path.isfile(writefile) == False: protected_files.append(writefile)
            
            if (os.path.isfile(writefile) == True) & (overwrite_data == True) & \
               (writefile not in protected_files):
                protected_files.append(writefile)
                os.remove(writefile)
            
            if os.path.isfile(writefile) == False:
                # Writing file header if file does not already exist
                if (d.date() < station_start.date()) or (d.date() > station_end.date()):
                    lat_write,lon_write,alt_write = lat_set,lon_set,alt_set
                    ht_tt_write,ht_pp_write       = ht_tt_set,ht_pp_set
                    print('*** WARNING: Default station coordinates are used ***')
                else:
                    lat_write,lon_write,alt_write = station_lat,station_lon,station_alt
                    ht_tt_write,ht_pp_write       = station_ht_tt,station_ht_pp
                    
                header1 = '#DATE='+d.strftime('%Y%m%d')+',APOLLO='+str(apo).zfill(3)+ \
                          ',SERIAL='+serial+',FW='+str(fw_data)  
                header2 = '#LAT={:.5f}N,'.format(lat_write)+'LON={:.5f}E,'.format(lon_write)+\
                          'ALT={:.1f}M,'.format(alt_write)+'H_T={:.1f}M,'.format(ht_tt_write)+\
                          'H_P={:.1f}M'.format(ht_pp_write)

                wfile = open(writefile,'w')
                wfile.write(header1+'\n')
                wfile.write(header2+'\n')
                wfile.write(header3+'\n')
                wfile.close()
                
            else:
                # Create filling data if continue writing into exisiting file 
                last_timestep_old = fst.last_datetime_of_file(writefile)
                
                if last_timestep_old == dt.datetime(2000,1,1,0,0,0): 
                    print('*** Output file exists but is empty! ***')
                    continue
                
                if last_timestep_old.strftime('%H%M%S') == '235959': 
                    print('*** Output file is already filled! ***')
                    continue

                dt_grid_fill = pd.date_range(start=last_timestep_old,end=dt_grid[-1],\
                                             freq='s')[1:]
                data = data.reindex(dt_grid_fill)
                nfill = len(dt_grid_fill) - len(dt_grid)
                if nfill < 0:
                    # Due to time drift and recalibration of internal clock after restart, 
                    # first time step of new file may be before last timestep of previous 
                    # file. Cut new data where previous data ended
                    # Remark: Time drift should be usually smaller than 1 min
                    print(str(np.abs(nfill))+' timestep(s) removed due to overlapping '+\
                          'time stamps of consecutive data sets!')
                    print('Last old time step:',last_timestep_old)
                    print('First new time step:',dt_grid[0])
                    
            iday = (data.index.date == d.date())
            
            wdata = np.column_stack([data.index[iday].strftime('%H%M%S'), 
                                     data['TT_NTC'][iday].apply(wdata_tp),
                                     data['PP_BME'][iday].apply(wdata_tp),
                                     data['GPS'][iday].apply(wdata_s),
                                     data['WIFI'][iday].apply(wdata_s),
                                     data['LORA'][iday].apply(wdata_s),
                                     data['TTN'][iday].apply(wdata_s),
                                     ]) 
            
            wfile = open(writefile,'ab') 
            np.savetxt(wfile,wdata,delimiter=';',fmt=['%s']*wdata.shape[1])
            wfile.close()   
            
#----------------------------------------------------------------------------
print('*********')
fst.print_runtime(t_run)