# -*- coding: utf-8 -*-
"""
@author: Bastian Kirsch (bastian.kirsch@uni-hamburg.de)

Script to convert raw (=level0) APOLLO (Autonomous cold POoL LOgger) data 
(Software Version 115, FESSTVaL 2021) to level1 daily .txt-files

Data structure needed in datadir: ./<SSS>_<first 4 characters of serial number>/<YYMMDDHHMM>.csv
Output data structure in outdir: ./<YYYY>/<MM>/<DD>/APOLLO-<SSS>_<YYYYMMDD>.txt

Dependences on non-standard software:
- fesstval_routines.py 

Required meta data files:
- stations_fesstval.txt
- APOLLO_serials.txt

last updated: 28 October 2021
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
logdir      = maindir+'APOLLO_data/log/'
meta_file   = maindir+'FESSTVaL/stations_fesstval.txt'
serial_file = maindir+'FESSTVaL/APOLLO/apollo_serials.txt'

apollo_start = 1   # 1
apollo_end   = 110 # 110

start_date   = dt.date(2021,4,28) #dt.date(2021,4,28)
end_date     = dt.date(2021,9,1)  #dt.date(2021,9,1)

overwrite_data = True # Delete all existing level 1 files prior to writing
log_to_file    = True
only_fesstval  = True

min_filesize   = 20000 # 20 KB = ~10 min


# Files without valid GPS time stamp but trusted logger time stamp
# (=no significant drift of logger time in previous/following files)
files_no_gps = ['04E437-2108240812.csv',  # APOLLO 44 
                '58F237-2107011327.csv', # APOLLO 51
                '007443-2107051356.csv', # APOLLO 64
                '007343-2106251536.csv', # APOLLO 69
                ]

n_header_remove = 10 # Number of lines removed from file header due to
                     # unreliable logger time stamp

'''
Manually corrected files (corrupt lines deleted):

- APOLLO 5   087343-2107220814.csv 
- APOLLO 6   58E437-2107261515.csv 
- APOLLO 12  4CD103-2104281736.csv 
- APOLLO 13  40E237-2104280945.csv
- APOLLO 28  507343-2105070914.csv 
- APOLLO 31  C4E137-2104291730.csv 
- APOLLO 32  E47143-2104301043.csv 
- APOLLO 44  04E437-2104292312.csv 
- APOLLO 45  44E437-2104300854.csv 
- APOLLO 48  104E53-2105190937.csv 
- APOLLO 53  114E13-2104291420.csv 
- APOLLO 54  C41094-2105031439.csv 
- APOLLO 55  60E337-2105031541.csv 
- APOLLO 57  D8EC37-2104300741.csv 
- APOLLO 66  1C9938-2104291703.csv 
- APOLLO 69  007343-2104281033.csv 
- APOLLO 76  08E437-2105031300.csv 
- APOLLO 101 487943-2104300754.csv 

'''
#----------------------------------------------------------------------------
print('Processing APOLLO Level 0 data')
print('APOLLO Start:',apollo_start)
print('APOLLO End  :',apollo_end)
print('Date Start  :',start_date)
print('Date End    :',end_date)
if overwrite_data:
    print('Overwriting is enabled')
else:
    print('Overwriting is disabled')
    
if log_to_file:
    print('Logging to file is enabled')
    if not os.path.isdir(logdir): os.mkdir(logdir)
    log_file = dt.datetime.now().strftime(logdir+'log_apollo_l0_l1_%Y%m%d_%H%M.txt')
    
    def print_to_file(message,log=log_file):
        with open(log,'a') as f:
            print(message,file=f)
        return 
else:
    print('Logging to file is disabled') 
          

filename_len = len('ssssss-yymmddhhmm.csv')
header3      = '#TIME(UTC);T_NTC(degC);P_BME(hPa);GPS;WIFI;LORA;TTN'

ht_tt = 3 # m above ground
ht_pp = 2 # m above ground

t_drift_tol = 300 # (s) Max. time drift between logger time and GPS time considered as valid

protected_files = [] # Needed to overwrite existing files only once

stations = fst.fesstval_stations('all',metafile=meta_file,serialfile=serial_file,
                                  include_time=True,include_serial=True)
ii_apo   = stations['STATION'].str.endswith('a')  


for apo in range(apollo_start,apollo_end+1):
    iapo = stations[ii_apo]['APOLLO'] == apo
    station_info = stations[ii_apo][iapo]
    
    nstat_apo = len(station_info.index)
    istat_apo = 0
    
    if station_info.empty: continue
    
    serial        = station_info['SERIAL'].iloc[istat_apo]
    station_name  = station_info['NAME'].iloc[istat_apo]
    station_start = station_info['START']
    station_end   = station_info['END']

    if nstat_apo == 2:
        station_name = station_info['NAME'].iloc[0]+' / '+\
                       station_info['NAME'].iloc[1]
            
    print('')
    print('APOLLO '+str(apo)+' ('+serial[:6]+') '+station_name)
    
    apo_dir = datadir+str(apo).zfill(3)+'_'+serial[:6]+'/'
        
    if os.path.isdir(apo_dir) == False:
        print('Directory '+apo_dir+' does not exist!')
        continue
    
    file_list_apo = os.listdir(apo_dir) 
    
    if len(file_list_apo) == 0:
        print('No data files to process available')
        continue 
    
    # Check which files are already processed in case existing data should
    # not be overwritten    
    processed_data = []
    for f in glob.iglob(outdir+'**/APOLLO-'+str(apo).zfill(3)+'*.txt',recursive=True):
        processed_data.append(f.replace('\\','/'))
    processed_data.sort()    
        
    if (len(processed_data) > 0) & (not overwrite_data):    
        last_processed_datetime = fst.last_datetime_of_file(processed_data[-1])
    else:
        last_processed_datetime = station_start.min()
    
    # Remove all pre-existing level1 files prior to processing    
    if overwrite_data: 
        for p in processed_data: os.remove(p) 
        
    # Remove small level0 files (and 'original' files)
    size_list_apo = [os.path.getsize(apo_dir+f) for f in file_list_apo]
    remove_list = []
    for i in range(len(file_list_apo)):
        if (size_list_apo[i] < min_filesize) or ('ori' in file_list_apo[i]): 
            remove_list.append(file_list_apo[i])
    for r in remove_list: file_list_apo.remove(r)
        
    nfiles = len(file_list_apo)
    if nfiles > 0:
        print(str(nfiles)+' files to process found')
    else: 
        print('No files to process found!')
    
    # Sort files after time stamp    
    ii_sort_time = np.argsort([f.split('-')[-1] for f in file_list_apo])
    file_list_apo_sorted = [file_list_apo[i] for i in ii_sort_time]          
        
    # Loop over files
    for nf,afile in enumerate(file_list_apo_sorted):
        
        filename   = apo_dir+afile
        filesize   = os.path.getsize(filename)
        fileserial = afile.split('-')[0]
        
        size_str  = ' ({:3.1f}'.format(filesize/(1024**2))+' MB)'
        print_str = 'APOLLO '+str(apo)+' '+afile+size_str+' '
            
        print(str(nf+1)+'/'+str(nfiles)+': '+print_str)
        
        # Check conditions to sort out invalid files
        err_msg = ''
        
        try:
            filetime = pd.to_datetime(afile.split('-')[1].split('.')[0],
                                      format='%y%m%d%H%M')
        except ValueError:
            err_msg = 'Invalid time string detected in file name'
        
        if (len(err_msg) == 0) & (os.path.getsize(apo_dir+afile) == 0):
            err_msg = 'File is empty'
            
        if (len(err_msg) == 0) & (len(afile) not in [filename_len,filename_len-1]):
            err_msg = 'Invalid file name (length)'
            
        if (len(err_msg) == 0) & (fileserial != serial[:6]):
            if (fileserial.upper() != serial[:6]) & (fileserial.upper() != serial[1:6]):
                err_msg = 'File name does not start with correct serial'
            
        if (len(err_msg) == 0) & filename.endswith('0000.csv'):
            err_msg = 'File probably has a wrong time stamp (starts at 0 UTC)'
        
        if only_fesstval:    
            if (len(err_msg) == 0) & (filetime < station_start.min()-dt.timedelta(hours=1)):
                err_msg = 'File starts before installation of instrument'
                                      
            if (len(err_msg) == 0) & (filetime > station_end.max()):
                err_msg = 'File starts after de-installation of instrument' 
             
        if len(err_msg) > 0:
            print(err_msg)
            if log_to_file: print_to_file(print_str+err_msg)
            continue    
        

        if (not overwrite_data) & (filetime < last_processed_datetime): 
            print('Skipping file (overwriting is disabled)')
            continue
        
        if (filetime.date() < start_date) or (filetime.date() > end_date): 
            print('Skipping file (outside of processing period)')
            continue
 
        # Reading data
        data,read_msg = fst.read_apollo_level0(filename,return_msg=True)
        if log_to_file & (len(read_msg) > 0): print_to_file(print_str+read_msg)
        if data.empty: continue
        
        ntime       = len(data.index)
        serial_data = data.serial
        status_data = data.status
        lat_data    = data.lat
        lon_data    = data.lon
        fw_data     = data.firmware
        
        if serial[:6] != serial_data:
            msg = 'Serial in file header does not match nominal serial of APOLLO'
            if log_to_file: print_to_file(print_str+msg)
            print(msg)
            continue
        
        # Remove first lines of file due to unreliable logger time stamp
        data = data.iloc[n_header_remove:]
        
        # Correct logger time stamp for drift
        dt_log_invalid = data['DTLOG'].isnull()
        data.drop(data.index[dt_log_invalid],inplace=True)
        
        dt_corr,dt_corr_msg = fst.correct_time_drift(data['DTLOG'],data['DTGPS'],
                                                     tdrift_tol=t_drift_tol,
                                                     return_msg=True)
        
        if log_to_file & (len(dt_corr_msg) > 0): print_to_file(print_str+dt_corr_msg)
        if (len(dt_corr) == 0) & (afile not in files_no_gps): continue
        
        if afile not in files_no_gps:
            data['DTCORR'] = dt_corr
        else:
            data['DTCORR'] = data['DTLOG']
            msg = 'Logger time stamp applyed without GPS correction'
            if log_to_file: print_to_file(print_str+msg)
            print(msg)
            
        ii_dupli = data['DTCORR'].duplicated()
        n_dupli = ii_dupli.sum()
        data.drop(data.index[ii_dupli],inplace=True)
        
        if n_dupli >= 100:
            msg = str(n_dupli)+' duplicate timesteps removed'
            if log_to_file: print_to_file(print_str+msg)
            print(msg)
        
        fvi = data['DTCORR'].first_valid_index()
        lvi = data['DTCORR'].last_valid_index()
        dt_grid = pd.date_range(start=pd.to_datetime(data['DTCORR'].loc[fvi]),
                                end=pd.to_datetime(data['DTCORR'].loc[lvi]),
                                freq='s')
        data.set_index('DTCORR',inplace=True)
        data = data.reindex(dt_grid)
        
        istat_apo = 0
        if nstat_apo == 2:
            ii_stat = (station_info['START'].dt.date <= dt_grid[0].date()) & \
                      (station_info['END'].dt.date > dt_grid[0].date())
            istat_apo = np.where(ii_stat)[0][0]
            
        station_lat = station_info['LAT'].iloc[istat_apo]    
        station_lon = station_info['LON'].iloc[istat_apo]
        station_alt = station_info['ALT'].iloc[istat_apo]
        
        if only_fesstval:
            # Drop data before start or after end of meausurements
            ii_invalid = (data.index < station_start.iloc[istat_apo]) | \
                         (data.index > station_end.iloc[istat_apo])
            if ii_invalid.sum() > 0: 
                data.drop(data.index[ii_invalid],inplace=True)
                
            # Compare coordinates in header with nominal coordinates
            if (np.abs(lat_data-station_lat) > 0.005) or (np.abs(lon_data-station_lon) > 0.01):
                msg = 'Coordinates in file header deviate from station coordinates'    
                if log_to_file: print_to_file(print_str+msg)
                print(msg)
                print('Nominal Coordinates:',str(station_lat),str(station_lon))
                print('File Coordinates   :',str(lat_data),str(lon_data))
        
        if len(data.index) == 0: 
            msg = 'No valid data to process found after correcting time stamp'
            if log_to_file: print_to_file(print_str+msg)
            print(msg)
            continue
        
        # Convert Integer data to nullable data type to ensure correct writing format
        int_list = ['GPS','WIFI','LORA','TTN']
        data[int_list] = data[int_list].astype('Int8')
        
        days = pd.date_range(start=data.index[0].date(),end=data.index[-1].date())
        
        # Write the data into separate files for each day
        for d in days:
            print(d.strftime('%Y-%m-%d'))
            
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
                lat_write = station_lat if only_fesstval else lat_data
                lon_write = station_lon if only_fesstval else lon_data
                alt_write = station_alt if only_fesstval else 0
                    
                header1 = '#DATE='+d.strftime('%Y%m%d')+',APOLLO='+str(apo).zfill(3)+ \
                          ',SERIAL='+serial+',FW='+str(fw_data)  
                header2 = '#LAT={:.5f}N,'.format(lat_write)+ \
                          'LON={:.5f}E,'.format(lon_write)+ \
                          'ALT={:.1f}M,'.format(alt_write)+ \
                          'H_T={:.1f}M,'.format(ht_tt)+ \
                          'H_P={:.1f}M'.format(ht_pp)

                wfile = open(writefile,'w')
                wfile.write(header1+'\n')
                wfile.write(header2+'\n')
                wfile.write(header3+'\n')
                wfile.close()
                
            else:
                # Create filling data if continue writing into exisiting file 
                last_timestep_old = fst.last_datetime_of_file(writefile)
                
                if last_timestep_old == dt.datetime(2000,1,1,0,0,0):
                    msg = d.strftime('Level 1 file for %Y-%m-%d exists but is empty')
                    if log_to_file: print_to_file(print_str+msg)
                    print(msg)
                    continue
                
                if last_timestep_old.strftime('%H%M%S') == '235959': 
                    msg = d.strftime('Level 1 file for %Y-%m-%d is already filled')
                    if log_to_file: print_to_file(print_str+msg)
                    print(msg)
                    continue

                dt_grid_fill = pd.date_range(start=last_timestep_old,
                                             end=dt_grid[-1],
                                             freq='s')[1:]
                data = data.reindex(dt_grid_fill)
                nfill = len(dt_grid_fill) - len(dt_grid)
                if nfill < 0:
                    # Due to time drift and recalibration of internal clock after restart, 
                    # first time step of new file may be before last timestep of previous 
                    # file. Cut new data where previous data ended
                    # Remark: Time drift should be usually smaller than 1 min
                    msg = str(np.abs(nfill))+' timesteps removed due to overlapping '+\
                          'time stamps of consecutive data sets!'
                    if log_to_file: print_to_file(print_str+msg)
                    print(msg)
                    print('Last old time step:',last_timestep_old)
                    print('First new time step:',dt_grid[0])  
        
            iday = (data.index.date == d.date())
            wvar = ['TT_NTC','PP_BME']+int_list
            
            data[wvar].loc[iday].to_csv(writefile,mode='a',
                                        index=True,header=False,
                                        sep=';',na_rep='',
                                        date_format='%H%M%S',
                                        float_format='%1.2f')
          
#----------------------------------------------------------------------------
print(' ')
print('*** Finshed! ***')
runtime_str = fst.print_runtime(t_run,return_str=True)
if log_to_file: print_to_file('*** '+runtime_str+' ***')