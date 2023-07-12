# -*- coding: utf-8 -*-
"""
@author: Bastian Kirsch (bastian.kirsch@uni-hamburg.de)

Functions used for processing FESST@HH 2020 and FESSTVaL 2021 data

Last updated: 6 July 2023
"""

import numpy as np
import datetime as dt
import pandas as pd
import xarray as xr
import os
import ftplib
import glob

maindir = '.'

#----------------------------------------------------------------------------
def print_runtime(dt_start,return_str=False):
    t_run  = np.round((dt.datetime.now()-dt_start).total_seconds())
    dt_ref = dt.datetime(2000,1,1,0,0,0) + dt.timedelta(seconds=t_run)
    if t_run < 3600: 
        print_str = dt_ref.strftime('Run time: %M:%S min')
    if (t_run >= 3600) & (t_run < 86400): 
        print_str = dt_ref.strftime('Run time: %H:%M h')
    if t_run >= 86400: 
        print_str = dt_ref.strftime('Run time: %d day(s) %H:%M h')
    print(print_str)    
    if return_str:
        return print_str
    else:
        return    

def TK_to_TC(TTK):
    return TTK - 273.15

def datetime_to_unixtime(dtime):
    return np.array(dtime.astype(np.int64)/1e9,dtype=int)

def unixtime_to_datetime(utime):
    return pd.to_datetime(utime*1e9)  

def R_to_T(R):
   #Steinhart-Hart-Equation for NTC (Resistance to Temperature)
   alpha = 1.12881283*10**(-3)
   beta  = 2.34176548*10**(-4)
   gamma = 8.74951992*10**(-8)
   T = (alpha+(beta*np.log(R)+gamma*np.log(R)**3))**(-1)
   return TK_to_TC(T)

def correct_apollo_date(datestr):
    # Fill leading zero for day if missing (065=5 June)
    if len(datestr) == 4: 
        return datestr
    else:
        return datestr[:2]+'0'+datestr[2:] 

def correct_apollo_date_advanced(datestr):
    # Fill leading zero for day or month if missing (065=5 June, 71=1 July,620=20 June)
    if type(datestr) != str: return datestr
    if len(datestr) == 4: 
        return datestr
    if len(datestr) == 3:
        if datestr.startswith('0') == False: 
            return datestr.zfill(4)
        else:
            return datestr[:2]+'0'+datestr[2:] 
    if len(datestr) == 2:
        return '0'+datestr[0]+'0'+datestr[1]
    
def correct_apollo_time_advanced(timestr):
    if type(timestr) != str: return timestr
    return timestr.zfill(6)
    
def correct_time_drift(dt_log,dt_gps,tdrift_tol=300,return_msg=False):
    
    if dt_gps.shape != dt_log.shape:
        msg = 'DTLOG and DTGPS do not have the same shape'
        print(msg)
        if not return_msg:
            return pd.DataFrame()
        else:
            return pd.DataFrame(),msg
        
    # Convert logger time and GPS time to unixtime
    ut_log = pd.Series(datetime_to_unixtime(dt_log))
    ut_gps = pd.Series(datetime_to_unixtime(dt_gps))
    
    # Calculate time drift in seconds
    drift_sec = ut_log-ut_gps
    
    # Define valid time drifts by maximum (absoulte) tolerance
    valid_drift = (drift_sec.abs() <= tdrift_tol)
    
    # Calculate synchronity (equality of time steps)
    diff_sync = ut_gps[valid_drift].diff()/ut_log[valid_drift].diff()
    #diff_sync = ut_gps.diff()/ut_log.diff()
    
    # Define synchronized time steps
    valid_sync = (diff_sync == 1)
    
    # Determine valid GPS time steps
    valid_gps = (valid_drift & valid_sync)
    n_valid   = valid_gps.sum()
    n_tot     = dt_log.shape[0]
    
    if n_valid == 0:
        if n_tot >= 3600:
            msg = 'No valid GPS time for {:d} hour(s) of data'.format(n_tot//3600)
            print(msg)
        else:
            msg = ''
        if not return_msg:
            return pd.DataFrame()
        else:
            return pd.DataFrame(),msg
    
    msg = ''
    drift_max = drift_sec[valid_gps].max()
    if drift_max > 60:
        msg = 'Maximum time drift of {:d} seconds'.format(drift_max)+\
              ' (valid time steps = {:d})'.format(n_valid)
        print(msg) 
        
    # Create array of time drift-correcting seconds by padding (=forward-filling) 
    # valid time drifts to next valid time drift (with backward-filling at beginning)    
    drift_sec[~valid_gps] = np.nan   
    corr_sec = drift_sec.fillna(method='ffill')
    fvi = corr_sec.first_valid_index()
    corr_sec.iloc[:fvi] = corr_sec.iloc[fvi]
        
    # Correct logger time 
    dt_corr = unixtime_to_datetime(ut_log-corr_sec)
    
    if not return_msg:
        return dt_corr
    else:
        return dt_corr,msg 
    

def extract_val(sstr,varstr,varendstr,return_str=False):
    #if type(sstr) != str: print(sstr)
    i_start = sstr.find(varstr)
    if i_start == -1: return np.nan
    i_end = sstr[i_start+len(varstr):].find(varendstr)
    if varendstr == '': i_end = len(sstr)-1
    if i_end == -1: return np.nan
    valstr = sstr[i_start+len(varstr)+1:i_end+i_start+len(varstr)]
    if return_str: return valstr
    #if valstr.replace('.','').isdecimal(): return float(valstr)
    #else: return np.nan 
    try:
        return float(valstr)
    except ValueError:
        return np.nan 
    
def gps_to_dec(gps_val):
    degrees = gps_val//100.
    minutes = (gps_val%100.)/60.
    return degrees+minutes   

def apollo_status(status_float,status_len):
    try:
        # *** Invert status_str since status is calculated according to binary exponents, 
        # not positions in string ***
        return format(int(status_float),'0'+str(status_len)+'b')[::-1] 
    except ValueError:
        return ''
    

def read_apollo_level0(filename,status_str='CGMWLT',header_only=False,
                       return_msg=False):
    
    if os.path.isfile(filename) == False:
        print('File '+filename+' does not exist!')
        return pd.DataFrame()   
    
    data_columns_type = {'DATELOG':str,
                         'TIMELOG':str,
                         'DATEGPS':str,                     
                         'TIMEGPS':str,
                         'R_NTC'  :float,
                         'PP_BME' :float,
                         'STATUS' :float,
                         }
    
    dict_status = {'C':'CONSOLE',
                   'G':'GPS',
                   'M':'MEASURE',
                   'W':'WIFI',
                   'L':'LORA',
                   'T':'TTN',
                   }
    
    n_status   = len(status_str)
    
    # Read header
    try:
        header_str = pd.read_csv(filename,header=None,nrows=1)[0][0]
    except ValueError:
        msg = 'Header of file could not be correctly read'
        print(msg)
        if not return_msg:
            return pd.DataFrame()
        else:
            return pd.DataFrame(),msg     
        
    i_year     = header_str.find('START_DATE=') 
    year_str   = header_str[i_year+11:i_year+15]
    i_serial   = header_str.find('NODE_ID=')         
    serial_str = header_str[i_serial+8:i_serial+14]
    i_lat      = header_str.find('Lat.:')
    lat        = float(header_str[i_lat+6:i_lat+14]) if i_lat > 0 else np.nan
    i_lon      = header_str.find('Long.:')
    lon        = float(header_str[i_lon+7:i_lon+15]) if i_lon > 0 else np.nan
    i_fw       = header_str.find('FW=')
    fw_str     = header_str[i_fw+3:i_fw+11] if i_fw > 0 else 'v.114a7'
    
    if header_only:
        return {'YYSTR':year_str,'SERIAL':serial_str,'LAT':lat,'LON':lon,'FW':fw_str}
    
    # Read data
    try:
        data = pd.read_csv(filename,header=1,sep=';',names=list(data_columns_type.keys()),
                           dtype=data_columns_type)
    except ValueError:
        msg= 'File could not be correctly read'
        print(msg)
        if not return_msg:
            return pd.DataFrame()
        else:
            return pd.DataFrame(),msg       
    
    if data.empty: 
        msg= 'File is empty'
        print(msg)
        if not return_msg:
            return data
        else:
            return data,msg  
    
    data.status   = status_str
    data.serial   = serial_str
    data.lat      = lat
    data.lon      = lon
    data.firmware = fw_str
    
    # Correction of date and time strings if needed 
    if data['DATELOG'][0][0] == '0':
        data['DATELOG'] = data['DATELOG'].apply(correct_apollo_date)
    else:
        print('File probably misses leading zeros, advanced correction applied')
        data['DATELOG'] = data['DATELOG'].apply(correct_apollo_date_advanced) 
        data['TIMELOG'] = data['TIMELOG'].apply(correct_apollo_time_advanced)
        data['DATEGPS'] = data['DATEGPS'].apply(correct_apollo_date_advanced) 
        data['TIMEGPS'] = data['TIMEGPS'].apply(correct_apollo_time_advanced)
        
    if data['DATELOG'][0] == '0101':
        print('File probably has a wrong time stamp (starts on 1 Jan)')    
    
    try:        
        data['DTLOG']   = pd.to_datetime(year_str+data['DATELOG']+data['TIMELOG'],\
                                         format='%Y%m%d%H%M%S',errors='coerce')
        data['DTGPS']   = pd.to_datetime((year_str+data['DATEGPS']+data['TIMEGPS']),\
                                         format='%Y%m%d%H%M%S',errors='coerce').\
                                         fillna(dt.datetime(2000,1,1,0,0,0))   
    except ValueError:
        msg = 'Timestamp of file could not be correctly read'
        print(msg)
        if not return_msg:
            return pd.DataFrame()
        else:
            return pd.DataFrame(),msg    
                                     
        
    data['PP_BME'] = data['PP_BME']/100.
    data['TT_NTC'] = R_to_T(data['R_NTC'])
    
    data['STATUS_BIN'] = data['STATUS'].apply(lambda x:apollo_status(x,n_status))
    for i,s in enumerate(status_str): 
        data[dict_status[s]] = data['STATUS_BIN'].str[i].astype(float)   
    
    if not return_msg:
        return data
    else:
        return data,''


def read_wxt_level0(wxt,yy,mm,dd,datadir=maindir+'WXT_data/level0/'):
    filedate = dt.date(yy,mm,dd)
    filedir = datadir+'WXT_'+str(wxt).zfill(2)+filedate.strftime('/%Y%m/')
    
    filename = filedir+str(wxt).zfill(2)+filedate.strftime('%y%m%d.wxt')  
    
    if os.path.isfile(filename) == False:
        # Trying old file naming convention
        filename = filedir+filedate.strftime('%Y%m%d.wxt') 
        
    if os.path.isfile(filename) == False:
        print('File '+filename+' does not exist!')
        return pd.DataFrame({})  
    
    data_columns_type = {'DATE'    :str,
                         'TIME'    :str,
                         'TIMEZONE':str,                     
                         'VERSION' :str,
                         'STATION' :str,
                         'FIX'     :str,
                         'MSG'     :str,
                         'WXTDATA' :str,
                         }
    
    try:
        data = pd.read_csv(filename,header=0,sep=';',names=list(data_columns_type.keys()),
                            dtype=data_columns_type,encoding='latin1')
    except ValueError:
        print('File '+filename+' could not be correctly read!')
        return pd.DataFrame()
    

    data['FIX']=data['FIX'].fillna(-1).astype(int)
    
    data.timezone = data['TIMEZONE'][0]
    data.station  = int(data['STATION'][0])
    data.firmware = data['VERSION'][0]
    
    data.drop('TIMEZONE',1,inplace=True)
    data.drop('STATION',1,inplace=True)
    data.drop('VERSION',1,inplace=True)
    
    data['DTIME'] = pd.to_datetime(data['DATE']+data['TIME'],format='%d.%m.%Y%H:%M:%S')
    data.drop('DATE',1,inplace=True)
    data.drop('TIME',1,inplace=True)
    
    if wxt != data.station:
        print('*** WARNING: Input station number ('+str(wxt)+') does not match '+\
              'file station number ('+str(data.station)+') ***')
            
    data['WXTDATA'].mask(data['WXTDATA'].isnull(),other='',inplace=True)        
    
    # MSG 1        
    data['DD'] = data['WXTDATA'].apply(lambda x:extract_val(x,'Dm','D')) # DD (°)
    data['FF'] = data['WXTDATA'].apply(lambda x:extract_val(x,'Sm','M')) # FF (m/s)
    data['FB'] = data['WXTDATA'].apply(lambda x:extract_val(x,'Sx','M')) # Gusts (m/s)
    
    # MSG 2
    data['TT'] = data['WXTDATA'].apply(lambda x:extract_val(x,'Ta','C')) # TT_WXT (°C)
    data['RH'] = data['WXTDATA'].apply(lambda x:extract_val(x,'Ua','P')) # RH (%)
    data['PP'] = data['WXTDATA'].apply(lambda x:extract_val(x,'Pa','H')) # PP (hPa)
    
    # MSG 3
    data['RR'] = data['WXTDATA'].apply(lambda x:extract_val(x,'Rc','M')) # RR_accu (mm)
    data['HA'] = data['WXTDATA'].apply(lambda x:extract_val(x,'Hc','M')) # Hail Hits (1/cm2)
    
    # MSG 4
    data['TT_PT1000'] = data['WXTDATA'].apply(lambda x:extract_val(x,'Tr','C')) # TT_Pt1000 (°C)
    
    # MSG 5
    data['VS'] = data['WXTDATA'].apply(lambda x:extract_val(x,'Vs','V')) # Supply Voltage (V)
    
    # MSG 9 (GPS)
    data['LAT'] = data['WXTDATA'].apply(lambda x:gps_to_dec(extract_val(x,'La','N')))
    data['LON'] = data['WXTDATA'].apply(lambda x:gps_to_dec(extract_val(x,'Lo','E')))
    data['ALT'] = data['WXTDATA'].apply(lambda x:extract_val(x,'Al','M'))
    
    data['LAT'].mask(data['FIX']==0,inplace=True)  
    data['LON'].mask(data['FIX']==0,inplace=True)
    data['ALT'].mask(data['FIX']==0,inplace=True)   
    
    # data.lat = data['LAT'].median()
    # data.lon = data['LON'].median()
    # data.alt = data['ALT'].median()
    
    # data.drop('LAT',1,inplace=True)
    # data.drop('LON',1,inplace=True)
    # data.drop('ALT',1,inplace=True)
            
    return data


def last_datetime_of_file(filename):
    # Only for APOLLO/WXT level 1 data
    filedata = pd.read_csv(filename,sep=';',header=2,dtype='str')['#TIME(UTC)']
    if filedata.empty: return dt.datetime(2000,1,1,0,0,0)
    last_str = filedata.values[-1]
    date_str = os.path.splitext(os.path.basename(filename))[0][-8:]  
    return pd.to_datetime(date_str+last_str) 


def rename_apollo_file(filename):
    # Rename APOLLO level0 data files produced by logger software versions <114a9
    # to match the naming convention 'ssssss-yymmddhhmm.csv'
    fname_old = os.path.basename(filename)
    if len(fname_old) != len('yymmddhhmm.csv'): return
    
    filedir = os.path.dirname(filename)+'/'
    
    header_info = read_apollo_level0(filename,header_only=True)
    serial = header_info['SERIAL']
    fname_new = serial[:6]+'-'+fname_old
    
    os.rename(filename,filedir+fname_new)
    
    print('File '+fname_old+' renamed to '+fname_new)
    return


def read_fesstval_level1(stat_type,stat,yy,mm,dd,
                         include_monitoring_data=False,
                         return_filename=False,mute=False,
                         datadir_apollo=maindir+'APOLLO_data/level1/',
                         datadir_wxt=maindir+'WXT_data/level1/',
                         ):
    # Read FESST@HH level1 APOLLO and WXT data
    stat_type = stat_type.upper()
    if stat_type == 'A': stat_type = 'APOLLO'
    if stat_type == 'W': stat_type = 'WXT'
    if stat_type not in ['APOLLO','WXT']: stat_type = 'APOLLO'
    n_header = 3
    
    if stat_type == 'APOLLO':
        datadir = datadir_apollo
        z_fill = 3
        data_columns_type = {'DTIME':str,
                             'TT'   :float,
                             'PP'   :float,                     
                             'GPS'  :float,
                             'WIFI' :float,
                             'LORA' :float,
                             'TTN'  :float,
                             }
        
    if stat_type == 'WXT':
        datadir = datadir_wxt
        z_fill = 2
        data_columns_type = {'DTIME'    :str,
                             'TT'       :float,
                             'TT_PT1000':float,
                             'PP'       :float, 
                             'RH'       :float,
                             'FF'       :float,
                             'FB'       :float,
                             'DD'       :float,
                             'RR'       :float,
                             'HA'       :float,
                             'VS'       :float,
                             }
    
    filedate = dt.date(yy,mm,dd)
    filename = datadir+filedate.strftime('%Y/%m/%d/')+\
               stat_type+'-'+str(stat).zfill(z_fill)+'_'+filedate.strftime('%Y%m%d')+'.txt'
              
    if return_filename: return filename          
              
    if os.path.isfile(filename) == False:
        if mute == False: print('File '+filename+' does not exist!')
        return pd.DataFrame({})  
    
    try:
        data = pd.read_csv(filename,header=n_header-1,sep=';', \
                           names=list(data_columns_type.keys()), \
                           dtype=data_columns_type)
    except ValueError:
        if mute == False: print('File '+filename+' could not be correctly read!')
        return pd.DataFrame({}) 
    
    data['DTIME'] = pd.to_datetime(filedate.strftime('%Y%m%d')+data['DTIME'], \
                                   format='%Y%m%d%H%M%S')
    data.set_index('DTIME',inplace=True)
    
    if include_monitoring_data == False:
        for col in ['GPS','WIFI','LORA','TTN','VS']:
            if col in data.columns: data.drop(col,1,inplace=True)

    header  = pd.read_csv(filename,sep=';',header=None,nrows=2)
    header1 = header[0][0]
    header2 = header[0][1]
    
    fw_str = extract_val(header1,'FW','',return_str=True)
    lat    = extract_val(header2,'LAT',',')
    lon    = extract_val(header2,'LON',',')
    alt    = extract_val(header2,'ALT',',')
    
    data.firmware = fw_str
    data.lat      = lat
    data.lon      = lon
    data.alt      = alt
    
    if stat_type == 'APOLLO':
        serial_str = extract_val(header1,'SERIAL',',',return_str=True)
        ht_tt      = extract_val(header2,'H_T',',')
        ht_pp      = extract_val(header2,'H_P','')
        
        data.serial = serial_str
        data.ht_tt  = ht_tt
        data.ht_pp  = ht_pp
        
    if stat_type == 'WXT': 
        ht_wxt = extract_val(header2,'H_WXT','')
        data.ht_wxt  = ht_wxt
            
    return data

def read_fesstval_level2(stat_type,var_str,yy,mm,dd,dataset_version=0,
                         mute=False,return_ds=False,return_meta=False,
                         datadir_apollo=maindir+'APOLLO_data/level2/',
                         datadir_wxt=maindir+'WXT_data/level2/'):
    
    # Read FESST@HH level2 APOLLO and WXT data
    stat_type = stat_type.upper()
    if stat_type == 'A': stat_type = 'APOLLO'
    if stat_type == 'W': stat_type = 'WXT'
    if stat_type not in ['APOLLO','WXT']: stat_type = 'APOLLO'
    campaign = 'fval' 
    if yy == 2019: campaign = 'prefval'
    if yy == 2020: campaign = 'fessthh'
    
    var_dict = {'TT':'ta',
                'PP':'pa',
                'RH':'hur',
                'FF':'wspeed',
                'FB':'wspeed_max',
                'DD':'wdir',
                'RR':'precip',
                'HA':'hail',
                }
    
    if stat_type == 'APOLLO':
        datadir  = datadir_apollo
        var_list = ['TT','PP']

    if stat_type == 'WXT':
        datadir  = datadir_wxt
        var_list = list(var_dict.keys())
        
    if var_str not in var_list:
        print('Variable '+var_str+' is not available! Pick one of '+str(var_list)) 
        return pd.DataFrame()
        
    filedate = dt.date(yy,mm,dd)  # ('%Y/%m/%d/')
    filename = datadir+filedate.strftime('%Y/%m/')+campaign+'_uhh_'+stat_type.lower()+\
               '_l2_'+var_dict[var_str]+'_v'+str(dataset_version).zfill(2)+'_'+\
               filedate.strftime('%Y%m%d000000.nc')
                    
    if os.path.isfile(filename) == False:
        if mute == False: print('File '+filename+' does not exist!')
        return pd.DataFrame()  
    
    try:
        ds = xr.open_dataset(filename)
    except ValueError:
        if mute == False: print('File '+filename+' could not be correctly read!')
        return pd.DataFrame() 
    
    if return_ds: return ds
        
    val = ds[var_dict[var_str]].values
    station_ids = ds['station_id'].values.astype(str).astype(object)
    
    if var_str == 'TT': val  = TK_to_TC(val)
    if var_str == 'PP': val *= 0.01
    if var_str == 'RH': val *= 100
    
    data = pd.DataFrame(val,columns=station_ids,index=ds['time'].values)

    if not return_meta:
        return data#.round(2)
    else:
        meta = pd.DataFrame(index=station_ids)
        meta['LAT']   = ds['lat'].values
        meta['LON']   = ds['lon'].values
        meta['ALT']   = ds['zsl'].values
        meta['HT']    = ds['zag'].values
        meta['NAME']  = ds['station_name'].values.astype(str).astype(object)
        #meta['INSTR'] = ds['instrument_id'].values
        meta['LCZ']   = ds['lcz'].values.astype(str).astype(object)
        
        return data,meta #.round(2)
  
  
def read_fesstval_wettermast(stat_str,yy,mm,dd,
                             datadir=maindir+'Wettermast_data/2020_FESSTHH/'):
    
    stat_list = ['WMH','NIK','WXT_WMH','WXT_BBG']
    
    if stat_str not in stat_list:
        print('Station '+stat_str+' not available! Pick one of '+stat_list)
        return pd.DataFrame({}) 
    
    filedate = dt.date(yy,mm,dd)
    week = filedate.isocalendar()[1]
    n_header = 7
    
    if stat_str == 'WMH':
        column_names = {'TT002' :'TT002',
                        'TT010' :'TT010',
                        'TT050' :'TT050',
                        'TT110' :'TT110',
                        'TT175' :'TT175',
                        'TT250' :'TT250',
                        'TT280' :'TT280',
                        'DT002' :'TD002',
                        'DT010' :'TD010',
                        'DT050' :'TD050',
                        'DT110' :'TD110',
                        'DT175' :'TD175',
                        'DT250' :'TD250',
                        'DT280' :'TD280',
                        'P002'  :'PP',
                        'FF010' :'FF010',
                        'FF050' :'FF050',
                        'FF110' :'FF110',
                        'FF175' :'FF175',
                        'FF250' :'FF250',
                        'FF280' :'FF280',
                        'FB010' :'FB010', #strongest gust
                        'FB050' :'FB050',
                        'FB110' :'FB110',
                        'FB175' :'FB175',
                        'FB250' :'FB250',
                        'FB280' :'FB280',
                        'DD010' :'DD010',
                        'DD050' :'DD050',
                        'DD110' :'DD110',
                        'DD175' :'DD175',
                        'DD250' :'DD250',
                        'DD280' :'DD280',
                        'RR'    :'RR',
                        'NHU'   :'CBH', # cloud base height (m)
                        'NHO'   :'CTH', # cloud top height (m)
                        'G'     :'GS',  # Globalstrahlung/kurzwellige Einstrahlung (W/m2)
                        'L'     :'LS',  # Langwellige Einstrahlung (W/m2)
                        'R'     :'SO',  # Kurzwellige Ausstrahlung (W/m2)
                        'E'     :'LO',  # Langwellige Ausstrahlung (W/m2)
                        'U_BHF' : 'SF',  # Sensibler Wärmestrom (ohne Querwindkorrektur; W/m2)
                        'U_LF'  : 'LF',  # Latenter Wärmestrom (W/m2)
                        'WU_BW' : 'WW',  # Vertikale Windkomponente aus USAT (m/s)
                        }
    
    
    if stat_str == 'WXT_WMH':
        column_names = {'WXT_T'  :'TT',
                        'WXT_TR' :'TT_PT1000',
                        'WXT_P'  :'PP',
                        'WXT_RH' :'RH',
                        'WXT_VEL':'FF',
                        'WXT_MXV':'FB',
                        'WXT_DIR':'DD',
                        'WXT_R'  :'RR',
                        'WXT_H'  :'HA',
                        'WXT_VS' :'VS',
                        } 
        
    if stat_str == 'WXT_BBG':
        column_names = {'WXT7_T'  :'TT',
                        'WXT7_TR' :'TT_PT1000',
                        'WXT7_P'  :'PP',
                        'WXT7_RH' :'RH',
                        'WXT7_VEL':'FF',
                        'WXT7_MXV':'FB',
                        'WXT7_DIR':'DD',
                        'WXT7_R'  :'RR',
                        'WXT7_H'  :'HA',
                        'WXT7_VS' :'VS',
                        }    
    
    if stat_str == 'NIK':
        column_names = {'XNIK_TK096' :'TT096',
                        'XNIK_TK111' :'TT111',
                        'XNIK_TK123' :'TT123',
                        }
        
    fname = stat_str
    if stat_str == 'WXT_WMH': fname = 'WXT'
    if stat_str == 'WMH'    : fname = 'MASTER'
    
    filename = datadir+'KW'+str(week)+'_'+fname+'.txt'
    
    if os.path.isfile(filename) == False:
        print('File '+filename+' does not exist!')
        return pd.DataFrame({})  
    
    try:
        data = pd.read_csv(filename,header=n_header-1,sep=';',na_values=99999)
    except ValueError:
        print('File '+filename+' could not be correctly read!')
        return pd.DataFrame({})     
    
    time_fmt = '%H:%M'
    timeres_sec = 60
    if 'WXT' in stat_str:
        time_fmt = '%H:%M:%S'
        data['TIME'] = data['TIME'].apply(lambda x:correct_wm_time(x))
        timeres_sec = 10
    
    dtime_to_utc = dt.timedelta(hours=1)
    data['DTIME'] = pd.to_datetime(data['$Names=DATE']+data['TIME'], \
                                   format='%d.%m.%Y'+time_fmt) - dtime_to_utc
    data.set_index('DTIME',inplace=True)
    
    datatime = pd.date_range(start=filedate, freq=str(timeres_sec)+'s',
                             periods=int(86400/timeres_sec))
    
    data = data.reindex(datatime)

    data.drop('$Names=DATE',1,inplace=True)
    data.drop('TIME',1,inplace=True)
    for c in data.columns:
        if c not in column_names.keys():
            data.drop(c,1,inplace=True)
            
    data.rename(columns=column_names,inplace=True)        
        
    lat,lon,alt = np.nan,np.nan,np.nan
    
    if stat_str == 'WMH': lat,lon,alt = 53.51990,10.10516,1.0
    if stat_str == 'NIK': lat,lon,alt = 53.54750, 9.99065,6.0
    if stat_str in ['WXT_WMH','WXT_BBG']: 
        if stat_str == 'WXT_WMH': sstr = 'Wettermast'
        if stat_str == 'WXT_BBG': sstr = 'Baursberg'
        stat_info = fessthh_stations(sstr,include_ht=True)
        lat       = stat_info['LAT'].to_numpy(float)[0]
        lon       = stat_info['LON'].to_numpy(float)[0]
        alt       = stat_info['ALT'].to_numpy(float)[0]
        ht_wxt    = stat_info['HT_WXT'].to_numpy(float)[0]
        
        data.ht_wxt  = ht_wxt
  
    data.lat      = lat
    data.lon      = lon
    data.alt      = alt
    
    return data        


def remove_implausible_values(data,varstr):
    
    var_thresholds = {'TT':[0, 40],     # (°C)
                      'PP':[950,1050],  # (hPa)
                      'RH':[0,100],     # (%)
                      'FF':[0,40],      # (m/s)
                      'DD':[0,360],     # (°)
                      'RR':[0,10],      # (mm/10s)
                      }
    
    if varstr not in list(var_thresholds.keys()):
        print('Variable '+varstr+' not available for removal of implausible values!')
        return data
    
    vmin,vmax = var_thresholds[varstr]
    
    with np.errstate(invalid='ignore'):
        ii_implaus = (data < vmin) | (data > vmax)
        
    n_implaus = ii_implaus.sum() if len(ii_implaus.shape) == 1 else ii_implaus.sum().sum()
        
    if n_implaus > 0:
        print(str(n_implaus)+' implausible '+varstr+' values removed!')
    
    return data.mask(ii_implaus)
    
def apollo_serial(apollo,serialfile=maindir+\
                  'FESSTVaL/APOLLO/apollo_serials.txt'):
    data = pd.read_csv(serialfile,sep=';').set_index('#APOLLO')
    if apollo in data.index:
        return data['Serial'][apollo]
    else:
        print('APOLLO '+str(apollo)+' not available!')
        return ''
    
def serial_apollo(serial,serialfile=maindir+\
                  'FESSTVaL/APOLLO/apollo_serials.txt'):
    data = pd.read_csv(serialfile,sep=';').set_index('#APOLLO')
    if type(serial) != str: serial = str(serial)
    ii_find = data['Serial'].str.startswith(serial).values
    return np.array(data.index[ii_find])  

def apollo_calibration(califile=maindir+'FESSTVaL/APOLLO/apollo_calibration.txt'):
    data = pd.read_csv(califile,sep=';',comment='#')
    cali = dict(zip(data['APOLLO'],data['TT_BIAS(K)']))
    return cali

def fessthh_stations(nn,find='NAME',metafile=maindir+'FESSTVaL/FESSTHH/stations_fessthh.txt',
                     serialfile=maindir+'FESSTVaL/APOLLO/apollo_serials.txt',
                     include_ht=False,include_time=False,include_serial=False,
                     include_lcz=False):

    if find.upper() == 'A': find = 'APOLLO'
    if find.upper() == 'W': find = 'WXT'
    if find.upper() == 'N': find = 'NAME'
    if find.upper() == 'S': find = 'SERIAL'
    if find.upper() == 'L': find = 'LCZ'
    
    if find == 'SERIAL': include_serial = True
    if find == 'LCZ': include_lcz = True
    
    if find not in ['STATION','NAME','APOLLO','WXT','SERIAL','LCZ']: find = 'NAME'
    if (find == 'NAME' or find == 'SERIAL' or find == 'LCZ') and (type(nn) != str):
        try:
            nn = str(nn)
        except ValueError:     
            print(find+' search item has to be a string!')
            return pd.DataFrame()
    if (find in ['STATION','APOLLO','WXT']) and (type(nn) != int):
        try:
            nn = int(nn)
        except ValueError:    
            print(find+' search item has to be an integer!')
            return pd.DataFrame()
    
    names = ['STATION','NAME','LAT','LON','ALT','APOLLO','WXT',
             'HT_TT','HT_PP','HT_WXT','START','END','LCZ']
    dtypes = [str,str,float,float,float,str,str,float,float,float,str,str,str]
    data = pd.read_csv(metafile,sep=';',header=0,names=names,dtype=dict(zip(names,dtypes)))

    data['APOLLO'] = data['APOLLO'].apply(lambda x: int(x) if type(x) == str else 0)
    data['WXT']    = data['WXT'].apply(lambda x: int(x) if type(x) == str else 0)
    data['START']  = data['START'].apply(lambda x: pd.to_datetime(x,format='%Y%m%d %H%M'))
    data['END']    = data['END'].apply(lambda x: pd.to_datetime(x,format='%Y%m%d %H%M'))
    data['END'].fillna(dt.datetime(2099,12,31,23,59,59),inplace=True)
    data['SERIAL'] = data['APOLLO'].apply(lambda x: apollo_serial(x,serialfile=serialfile) \
                                          if x > 0 else '')
    data['LCZ']    = data['LCZ'].apply(lambda x: x if x != np.nan else '')
        
    if find == 'STATION': 
        ii_find = data[find].str.startswith(str(nn).zfill(3)).values 
    if find == 'NAME': 
        if nn == '':
            # all measurement locations without duplicates due to instrument replacements
            # (last status of network)
            ii_find = (data['STATION'].str[:3].to_frame().duplicated(keep='last') == False)
        if nn == 'l2':
            # all measurement location without duplicates due to instrument replacements
            # that have existed once (level 2 meta data)
            ii_find = (data['STATION'].duplicated() == False)
        elif nn == 'all':
            # all measurement locations with duplicates due to instrument replacements
            ii_find = data['STATION'].notnull()
        else:    
            ii_find = data[find].str.contains(nn,case=False).values
    if find in ['APOLLO','WXT']:
        ii_find = (data[find] == nn).values
    if find == 'SERIAL':
        ii_find = data[find].str.startswith(nn).values
    if find == 'LCZ':
        ii_find = (data[find] == nn).values
    
    if include_ht == False:
        data.drop('HT_TT',1,inplace=True)
        data.drop('HT_PP',1,inplace=True)
        data.drop('HT_WXT',1,inplace=True)
    if include_time == False:
        data.drop('START',1,inplace=True)
        data.drop('END',1,inplace=True) 
    if include_serial == False:
        data.drop('SERIAL',1,inplace=True)    
    if include_lcz == False:
        data.drop('LCZ',1,inplace=True)    
    if find in ['APOLLO','SERIAL']:    
        data.drop('WXT',1,inplace=True) 
    if find == 'WXT':    
        data.drop('APOLLO',1,inplace=True)     

    return data[ii_find]


def fesstval_stations(nn,find='STATION',metafile=maindir+'FESSTVaL/stations_fesstval.txt',
                      serialfile=maindir+'FESSTVaL/APOLLO/apollo_serials.txt',
                      include_time=False,include_serial=False,include_lcz=False):
    
    
    if find.upper() == 'A': find = 'APOLLO'
    if find.upper() == 'W': find = 'WXT'
    if find.upper() == 'S': find = 'SERIAL'
    if find.upper() == 'L': find = 'LCZ'
    if find.upper() == 'N': find = 'NAME'
    
    if find == 'SERIAL': include_serial = True
    if find == 'LCZ': include_lcz = True
    
    if find not in ['STATION','NAME','APOLLO','WXT','SERIAL','LCZ']: find = 'STATION'
    
    names_dtypes = {'STATION':str,
                    'NAME'   :str,
                    'LAT'    :float,
                    'LON'    :float,
                    'ALT'    :float,
                    'APOLLO' :str,
                    'WXT'    :str,
                    'START'  :str,
                    'END'    :str,
                    'LCZ'    :str,
                    'MAINT'  :str,
                   }
    
    names = list(names_dtypes.keys())

    data = pd.read_csv(metafile,sep=';',header=0,names=names,dtype=names_dtypes)
    
    
    data['APOLLO'] = data['APOLLO'].apply(lambda x: int(x) if type(x) == str else 0)
    data['WXT']    = data['WXT'].apply(lambda x: int(x) if type(x) == str else 0)
    data['START']  = data['START'].apply(lambda x: pd.to_datetime(x,format='%Y%m%d %H%M'))
    data['END']    = data['END'].apply(lambda x: pd.to_datetime(x,format='%Y%m%d %H%M'))
    data['SERIAL'] = data['APOLLO'].apply(lambda x: apollo_serial(x,serialfile=serialfile) \
                                      if x > 0 else '')   
    
    if find == 'STATION': 
        #ii_find = data[find].str.startswith(str(nn).zfill(3)).values 
        ii_find = data[find].str.contains(nn,case=False).values
        
    if find == 'NAME':
         ii_find = data[find].str.contains(nn,case=False).values
         
    if find in ['APOLLO','WXT']:
        ii_find = (data[find] == nn).values 
        
    if find == 'SERIAL':
        ii_find = data[find].str.startswith(nn).values
        
    if nn == '': 
        ii_find = (data['STATION'].to_frame().duplicated(keep='last') == False)
        
    if nn == 'all':
        ii_find = data['STATION'].notnull() 
        
    if not include_time:
        data.drop('START',1,inplace=True)
        data.drop('END',1,inplace=True)  
    if include_serial == False:
        data.drop('SERIAL',1,inplace=True)    
    if include_lcz == False:
        data.drop('LCZ',1,inplace=True)    
    if find in ['APOLLO','SERIAL']:    
        data.drop('WXT',1,inplace=True) 
    if find == 'WXT':    
        data.drop('APOLLO',1,inplace=True) 
        
    data.drop('MAINT',1,inplace=True)  
    
    return data[ii_find]


def wxt_extent_filedir(wxt_file):
    if len(wxt_file) != 12:
        print('WXT filename'+wxt_file+' has wrong length!')
        return ''
    statstr,datestr = wxt_file[:2],wxt_file[2:8]
    filedate = pd.to_datetime(datestr,format='%y%m%d')
    
    wxt_file_ext = 'WXT_'+statstr+filedate.strftime('/%Y%m/')+wxt_file
    return wxt_file_ext