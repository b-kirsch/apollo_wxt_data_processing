# -*- coding: utf-8 -*-
"""
@author: Bastian Kirsch (bastian.kirsch@uni-hamburg.de)

Script to read FESSTVaL WXT level0 data, extract meteorological variables 
and GPS data from WXT messages and write them into level1 daily .txt-files

Input data:  <maindir>/WXT_data/level0/WXT_<SS>/<YYYYMM>/<SSYYMMDD>.wxt
Output data: <maindir>/WXT_data/level1/<YYYY>/<MM>/<DD>/WXT-<SS>_<YYYYMMDD>.txt

Dependences on non-standard software:
- fesstval_routines.py 

Required meta data files:
- stations_fesstval.txt

last updated: 19 October 2021

"""

import fesstval_routines as fst
import numpy as np
import datetime as dt
import pandas as pd
import os
import glob

t_run = dt.datetime.now()
#----------------------------------------------------------------------------
# Paths and basic settings
maindir   = '.'
datadir   = maindir+'WXT_data/level0/'
outdir    = maindir+'WXT_data/level1/' 
meta_file = maindir+'FESSTVaL/stations_fesstval.txt'

start_date = dt.date(2021,4,29)
end_date   = dt.date(2021,9,1)

wxt_start  = 1  #1
wxt_end    = 21 #21

overwrite_data = True
log_to_file    = True
only_fesstval  = True


'''
Manually corrected data due to incorrect time stamps:
- WXT01 after 01.07.2021 02:18:25 (2 timesteps corrected)
- WXT02 after 11.08.2021 23:11:03 (61 timesteps corrected)
- WXT03 after 08.07.2021 13:20:12 (4 timesteps deleted)
- WXT13 after 18.05.2021 06:55:48 (4 timesteps deleted)
- WXT16 after 28.08.2021 11:00:44 (1 timestep deleted)
'''

#----------------------------------------------------------------------------
print('Processing WXT level 0 data')
print('WXT Start :',wxt_start)
print('WXT End   :',wxt_end)
print('Date Start:',start_date)
print('Date End  :',end_date)
if overwrite_data:
    print('Overwriting is enabled')
else:
    print('Overwriting is disabled')

if log_to_file:
    print('logging to file is enabled')
    if not os.path.isdir('log'): os.mkdir('log')
    log_file = dt.datetime.now().strftime('log/log_wxt_l0_l1_%Y%m%d_%H%M.txt')
    
    def print_to_file(message,log=log_file):
        with open(log,'a') as f:
            print(message,file=f)
        return 
else:
    print('logging to file is disabled')    
        

header3 = '#TIME(UTC);TT(degC);TT_PT1000(degC);PP(hPa);RH(%);'+\
          'FF(m/s);FB(m/s);DD(deg);RR(mm);HAIL_HITS(1/cm2);U_SUPPLY(V)'
          
ht_wxt = 3 # m above ground          
        
def wdata_dd(f):
    try:
        return '{:d}'.format(int(f)) 
    except ValueError:
        return ''      
    
stations  = fst.fesstval_stations('all',metafile=meta_file,include_time=True)
ii_wxt    = stations['STATION'].str.endswith('w')  

for wxt in range(wxt_start,wxt_end+1):
    iwxt = stations[ii_wxt]['WXT'] == wxt
    station_info = stations[ii_wxt][iwxt]
    nstat = len(station_info.index)
    
    if station_info.empty: continue
        
    station_start = station_info['START'].squeeze()
    station_end   = station_info['END'].squeeze()
    station_lat   = station_info['LAT'].squeeze()
    station_lon   = station_info['LON'].squeeze()
    station_alt   = station_info['ALT'].squeeze()

    print('')
    print('WXT '+str(wxt).zfill(2)+' '+station_info['NAME'].squeeze())
    
    wxtdir = datadir+'WXT_'+str(wxt).zfill(2)+'/'
    month_list = os.listdir(wxtdir) 
    
    if not os.path.isdir(wxtdir) or (len(month_list) == 0):
        print('No data available!')
        continue
    
    available_data = []
    for mm in month_list:
        file_list = os.listdir(wxtdir+mm+'/')
        for f in file_list: available_data.append(wxtdir+mm+'/'+f)
    available_dates = [f.split('/')[-2]+os.path.splitext(os.path.basename(f))[0][-2:] \
                       for f in available_data] 
        
    # Correct filename extension to lower case if upper case
    for f in available_data:
        if f.endswith('.WXT'): os.rename(f,f.replace('.WXT','.wxt'))
    
    # Check which files are already processed in case existing data should
    # not be overwritten    
    processed_data = []
    for f in glob.iglob(outdir+'**/WXT-'+str(wxt).zfill(2)+'*.txt',recursive=True):
        processed_data.append(f.replace('\\','/'))
    processed_data.sort()  
    
    if overwrite_data:
        start_date_read = start_date
    else:
        if len(processed_data) == 0:
            start_date_read = start_date
        else:
            last_processed_datetime = fst.last_datetime_of_file(processed_data[-1])
            if last_processed_datetime.strftime('%H%M%S') == '235950':
                start_date_read = last_processed_datetime.date()+dt.timedelta(days=1)
            else:
                start_date_read = last_processed_datetime.date()   
    
    days = pd.date_range(start=start_date_read,end=end_date,freq='d') 
    
    if len(days) == 0: print('No files to process found!')

    for d in days:
        if d.strftime('%Y%m%d') not in available_dates: 
            continue
  
        if (d.date() < station_start.date()) or (d.date() > station_end.date()):
            if only_fesstval: continue
        
        print(d.strftime('%Y-%m-%d'))
        data = fst.read_wxt_level0(wxt,d.year,d.month,d.day,datadir=datadir)
        if data.empty: 
            msg = 'WXT'+str(wxt).zfill(2)+d.strftime(' %Y-%m-%d: ')+\
                  'File could not be read or does not exist'
            if log_to_file: print_to_file(msg)
            continue
        
        # Remove potential negative jumps in timestamp
        data.set_index('DTIME',inplace=True)
        data_ut   = fst.datetime_to_unixtime(data.index)
        ut_diff   = np.diff(data_ut)
        nt_remove = 0
        
        while np.all(ut_diff >= 0) == False:
            nt_remove += 1
            ii_diff_neg = np.where(ut_diff < 0)[0]
            if nt_remove == 1:
                dt_jump = data.index[ii_diff_neg][0]
                if pd.isnull(dt_jump): dt_jump = data.index[ii_diff_neg-1][0]
                    
            data.drop(data.index[ii_diff_neg+1],axis=0,inplace=True) 
            data_ut = fst.datetime_to_unixtime(data.index)
            ut_diff = np.diff(data_ut)
            
        if nt_remove >= 5:
            msg = 'WXT'+str(wxt).zfill(2)+d.strftime(' %Y-%m-%d: ')+\
                  str(nt_remove)+' time steps removed due to negative jumps at '+\
                      dt_jump.strftime('%H:%M:%S UTC')
            print(msg)
            if log_to_file: print_to_file(msg)
            
        # Create regular time grid (shift data to next full 10s time step)
        dt_grid = pd.date_range(start=d,freq='10s',periods=8640)    
        data_wxt = pd.DataFrame(np.zeros(len(dt_grid)),index=dt_grid,columns=['x'])
        
        for col in data.columns:
            if col in ['FIX','MSG','WXTDATA']: continue
            ii_nan   = data[col].isnull()
            data_fin = data[col][~ii_nan]
            ii_dupli = data_fin.index.duplicated() # in case duplicate timestamps remain 
                                                   # after jump removal
            data_wxt[col] = (data_fin[~ii_dupli]).reindex(dt_grid,method='pad',limit=1)
        data_wxt.drop('x',1,inplace=True)
        
        lat = np.nanmedian(data_wxt['LAT']) if data_wxt['LAT'].notnull().sum() > 0 else np.nan
        lon = np.nanmedian(data_wxt['LON']) if data_wxt['LON'].notnull().sum() > 0 else np.nan
        alt = np.nanmedian(data_wxt['ALT']) if data_wxt['ALT'].notnull().sum() > 0 else np.nan
        fw  = data.firmware
        
        write_lat = station_lat if only_fesstval else lat
        write_lon = station_lon if only_fesstval else lon
        write_alt = station_alt if only_fesstval else alt
        
        if (np.abs(lat-station_lat) > 0.001) or (np.abs(lon-station_lon) > 0.002):
            msg = 'WXT'+str(wxt).zfill(2)+d.strftime(' %Y-%m-%d: ')+\
                  'Geographical coordinates in file deviate from nominal coordinates'
            if only_fesstval:      
                print(msg)
                print('Nominal Coordinates:',str(station_lat),str(station_lon))
                print('File Coordinates   :',str(lat),str(lon))
                if log_to_file: print_to_file(msg)
        
        if data_wxt['TT'].isnull().all():
            print('No valid TT data found!')
            continue
        
        serial_msg = data[data['MSG'] == '5']['WXTDATA'].iloc[-1] 
        iid        = serial_msg.find('Id=')
        serial     = serial_msg[iid+3:iid+11]
        
        # Convert Integer data to nullable data type to ensure correct writing format
        data_wxt['DD'] = data_wxt['DD'].astype('Int16')
        
        # Drop data before start or after end of meausurements
        if only_fesstval:
            ii_invalid = (data_wxt.index < station_start) | (data_wxt.index > station_end)
            if ii_invalid.sum() > 0: data_wxt.drop(data_wxt.index[ii_invalid],inplace=True)
        
        # Writing data to file        
        writedir = outdir+d.strftime('%Y/%m/%d/')
        if os.path.isdir(writedir) == False: os.makedirs(writedir)
        
        writefile = writedir+'WXT-'+str(wxt).zfill(2)+'_'+d.strftime('%Y%m%d')+'.txt'
        if os.path.isfile(writefile):
            if overwrite_data: 
                os.remove(writefile)
            else: 
                print('File already exisists and overwriting is disabled!')
                continue
            
        header1 = '#DATE='+d.strftime('%Y%m%d')+',WXT='+str(wxt).zfill(2)+\
                  ',SERIAL='+serial+',FW='+fw  
        header2 = '#LAT={:.5f}N,'.format(write_lat)+\
                  'LON={:.5f}E,'.format(write_lon)+\
                  'ALT={:.1f}M,'.format(write_alt)+\
                  'H_WXT={:.1f}M'.format(ht_wxt)
        
        wfile = open(writefile,'w')
        wfile.write(header1+'\n')
        wfile.write(header2+'\n')
        wfile.write(header3+'\n')
        wfile.close()
        
        wvar = ['TT','TT_PT1000','PP','RH','FF','FB','DD','RR','HA','VS']
        data_wxt[wvar].to_csv(writefile,mode='a',index=True,header=False,
                              sep=';',na_rep='',date_format='%H%M%S',
                              float_format='%1.1f')

#----------------------------------------------------------------------------
print(' ')
print('*** Finshed! ***')
runtime_str = fst.print_runtime(t_run,return_str=True)
if log_to_file: print_to_file('*** '+runtime_str+' ***')