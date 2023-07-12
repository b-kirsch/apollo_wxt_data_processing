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
- stations_fessthh.txt

Last updated: 8 July 2020

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
meta_file = maindir+'FESSTVaL/FESSTHH/stations_fessthh.txt'

start_date = dt.date(2020,4,29)
end_date   = dt.date(2020,9,4)

wxt_start  = 1  #1
wxt_end    = 21 #21

overwrite_data  = True

ht_wxt_set = np.nan

# Manually corrected data due to incorrect time stamps:
# - WXT02, 16.06.2020 after 11:59:21 UTC
# - WXT02, 29.06.2020 after 21:31:35 UTC
# - WXT08, 22.06.2020 after 17:08:55 UTC 
# - WXT12, 08.06.2020 after 06:34:11 UTC

# CR manually inserted:
# - WXT09, 06.08.2020 after 16:31:50 UTC

# Corrupt line remove:
# - WXT06, 16.08.2020 after 07:17:00 UTC


#----------------------------------------------------------------------------
print('Processing WXT level 0 data')
if overwrite_data:
    print('Overwriting is enabled')
else:
    print('Overwriting is disabled')

header3 = '#TIME(UTC);TT(degC);TT_PT1000(degC);PP(hPa);RH(%);'+\
          'FF(m/s);FB(m/s);DD(deg);RR(mm);HAIL_HITS(1/cm2);U_SUPPLY(V)'
        
def wdata_wxt(f):
    return str(f).replace('nan','')

def wdata_dd(f):
    try:
        return '{:d}'.format(int(f)) 
    except ValueError:
        return ''         

for wxt in range(wxt_start,wxt_end+1):
    wxt_dir = datadir + 'WXT_'+str(wxt).zfill(2)+'/'
    station_info = fst.fessthh_stations(wxt,find='WXT',include_ht=True,
                                        include_time=True,metafile=meta_file)
    nstat = len(station_info.index)
    
    if station_info.empty == False:
        station_name    = station_info['NAME'].values
        station_lat     = station_info['LAT'].values
        station_lon     = station_info['LON'].values
        station_alt     = station_info['ALT'].values
        station_ht_wxt  = station_info['HT_WXT'].values
        station_start   = pd.to_datetime(station_info['START'].values)
        station_end     = pd.to_datetime(station_info['END'].values)
        print_stat_info = ''
        for s in station_name: print_stat_info += s+' / '
        print_stat_info = print_stat_info[:-3]
    else:
        print_stat_info = ''

    print('')
    print('WXT '+str(wxt).zfill(2)+' '+print_stat_info)
    
    if nstat > 1:
        print('*** Multiple station information found for this device! ***')

    if os.path.isdir(wxt_dir) == False:
        print('No data available!')
        continue
    
    month_list = os.listdir(wxt_dir) 
    
    if len(month_list) == 0:
        print('No data available!')
        continue  
    
    available_data = []
    for mm in month_list:
        file_list = os.listdir(wxt_dir+mm+'/')
        for f in file_list: available_data.append(wxt_dir+mm+'/'+f)
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
    
    days_stat = [pd.date_range(start=station_start[s].date(),\
                               end=station_end[s].date(),freq='d') \
                 for s in range(nstat)]
    
    if len(days) == 0: print('No files to process found!')

    for d in days:
        if d.strftime('%Y%m%d') not in available_dates: continue
    
        istat = -1
        for s in range(nstat):
            if d.date() in days_stat[s]:
                istat = s
                break
            
        if (d.date() < station_start[istat].date()) or (d.date() > station_end[istat].date()):
            continue
        
        print('Writing data for '+d.strftime('%Y-%m-%d'))
        data = fst.read_wxt_level0(wxt,d.year,d.month,d.day)
        if data.empty: continue
        
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
            
        if nt_remove >= 2:
            print('*** WARNING: '+str(nt_remove)+' time steps removed due to '+\
                  'negative jumps at '+dt_jump.strftime('%Y-%m-%d %H:%M:%S UTC')+' ***')
            
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
        
        if data_wxt['TT'].isnull().all():
            print('No valid TT data found!')
            continue
        
        serial_msg = data[data['MSG'] == '5']['WXTDATA'].iloc[-1] 
        iid        = serial_msg.find('Id=')
        serial     = serial_msg[iid+3:iid+11]
        
        # Drop data before start or after end of meausurements
        ii_invalid = (data_wxt.index < station_start[istat]) | (data_wxt.index > station_end[istat])
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
            
        if istat == -1:
            lat_write,lon_write,alt_write,ht_wxt_write = lat,lon,alt,ht_wxt_set
            print('*** WARNING: Station coordinates from file are used ***')
        else:
            lat_write,lon_write    = station_lat[istat],station_lon[istat]
            alt_write,ht_wxt_write = station_alt[istat],station_ht_wxt[istat]
            
        header1 = '#DATE='+d.strftime('%Y%m%d')+',WXT='+str(wxt).zfill(2)+\
                  ',SERIAL='+serial+',FW='+fw  
        header2 = '#LAT={:.5f}N,'.format(lat_write)+'LON={:.5f}E,'.format(lon_write)+\
                  'ALT={:.1f}M,'.format(alt_write)+'H_WXT={:.1f}M'.format(ht_wxt_write)
        
        wfile = open(writefile,'w')
        wfile.write(header1+'\n')
        wfile.write(header2+'\n')
        wfile.write(header3+'\n')
        wfile.close()
        
        wdata = np.column_stack([data_wxt.index.strftime('%H%M%S'),
                                 data_wxt['TT'].fillna(''),
                                 data_wxt['TT_PT1000'].fillna(''),
                                 data_wxt['PP'].fillna(''),
                                 data_wxt['RH'].fillna(''),
                                 data_wxt['FF'].fillna(''),
                                 data_wxt['FB'].fillna(''),
                                 data_wxt['DD'].apply(wdata_dd),
                                 data_wxt['RR'].fillna(''),
                                 data_wxt['HA'].fillna(''),
                                 data_wxt['VS'].fillna(''),
                                ])
        
        fmt = ['%s']*wdata.shape[1]
        wfile = open(writefile,'ab')
        np.savetxt(wfile,wdata,delimiter=';',fmt=fmt)
        wfile.close()    

#----------------------------------------------------------------------------
print('*********')
fst.print_runtime(t_run)