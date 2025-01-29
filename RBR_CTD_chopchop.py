# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:24:01 2024

@author: OllieTwigge
"""
#%%
# Match up each profile to CTD count - maybe using the time [filename, start, end]
# add attributes from arnos table to the netcdf 
# including Ship station number 
# change the filename of each Profile to the one
# CTD: <CruiseID>_CTD_RBR[unit-type]_<SN[X][X][X][X][X][X]>_<YYYYMMDD>_<HHMM>.[rsk][.txt] |
# Cruise ID is HYD2401
#
# # add new parameter that says up or downcast see below
# To make a new 
# xdata['gg'] = ('time',np.ones(len(xdata.depth)))

import xarray as xr
import pandas as pd
import numpy as np
from glob import glob
from datetime import datetime, timedelta
from pathlib import Path
import os
#%% just some functions

def extract_RBR_times_ID(fullpath): 
    # extracts the start and end time of each nc file that has come from ruskin
    # also now extracts the instrument type (platformID)
    pstart = []
    pend = []
    pid = []
    for i in fullpath:
        xdata = xr.open_dataset(i)   
        epstart = pd.to_datetime(xdata.attrs['time coverage start'])
        epend = pd.to_datetime(xdata.attrs['time coverage end'])
        epid = xdata.attrs['platform id']
        pstart.append(epstart)
        pend.append(epend)
        pid.append(epid)
        xdata.close
        
    return pd.DataFrame({'fullpath':fullpath, 'start':pstart, 'end':pend, 'pid': pid})

def new_file_names(fdat, cid): 
    # takes the info from inside each file and makes new file names for the inkfish standard
    # the fdat needs to have profile count/CTD number in it 
    nfms = []    
    
    for i in range(len(fdat)):
        #print(i)
        if np.isnan(fdat['PROFILE_COUNT'][i]):# == np.nan:
            nfms.append(os.path.basename(fdat['fullpath'][i]))
        else:
            fpath = fdat['fullpath'][i]
            ofname = Path(fpath).stem #original filename
            ofext = Path(fpath).suffixes[0]
            fsn = ofname.split('_')[0] # filename serial number
            
            if not any(i in fdat.columns for i in ['PROFILE_COUNT', 'CTD_Number']):
            # set(['PROFILE_COUNT', 'CTD_Number']) & set(list(fdat.columns)):
                raise Exception('No CTD numbering provided or column name needs adjusting')
            
            profile_num = int(fdat['PROFILE_COUNT'][i])
            ctd_num = f"{profile_num:03d}"
            
            ipid = fdat['pid'][i].split(' ')[0]#.replace(' ', '_') #instrument platform ID
            isn = fdat['pid'][i].split(' ')[1] # instrument serial number
            
            if fsn != isn:
                raise Exception(f'Mismatching of Serial Numbers between filename and filedata at {ofname}')
            
            ist = fdat['start'][i].strftime('%Y%m%d_%H%M') # instrument start time
            #cid = 'HYD2401' #cruise ID
            # CTD: <CruiseID>_CTD_RBR[unit-type]_<SN[X][X][X][X][X][X]>_<YYYYMMDD>_<HHMM>.[rsk][.txt]
            
            nfm = f'{cid}_CTD{ctd_num}_{ipid}_SN{isn}_{ist}{ofext}' #new_file_name
            nfms.append(nfm)
    print(nfms)
    return(nfms)

def trigger_index_morethan(df_name, col_name, 
                       trigger_float = 0.0, buffer = 0, window_size = 10):
    df = pd.DataFrame(df_name[col_name])
    df['is_trig'] = df > trigger_float
    df['is_streak'] = df['is_trig'].rolling(window=window_size).sum()
    df_index = df[buffer:].index[df[buffer:]['is_streak'] == window_size][0]
    print(f'Trig: {df_index} at {df_name.loc[df_index,"pressure"]}')
    return df_index

def trigger_index_lessthan(df_name, col_name, 
                       trigger_float = 0.0, buffer = 0, window_size = 10):
    df = pd.DataFrame(df_name[col_name])
    df['is_trig'] = df < trigger_float
    df['is_streak'] = df['is_trig'].rolling(window=window_size).sum()
    df_index = df[buffer:].index[df[buffer:]['is_streak'] == window_size][0]
    print(f'Trig: {df_index} at {df_name.loc[df_index,"pressure"]}')
    return df_index
#%%
metdat = pd.read_excel('RBR-Metadata.xlsx', header= 1, parse_dates=True)
metdat['RAW FILE'] = metdat['RAW FILE'].ffill()
metdat['DATE_TIME (UTC)'] = pd.to_datetime(metdat['DATE_TIME (UTC)'])

#%%
mainpath = r"C:\temp\Ollie-RBR"
testpath = mainpath+r"\*.nc"

fullpath= list(glob(testpath, recursive=False))
#%%
fdat = extract_RBR_times_ID(fullpath)
fdat = fdat.merge(metdat[['fullpath','PROFILE_COUNT']], how = 'left', on = 'fullpath')
#fdat = fdat[fdat['fullpath'].isin(metdat['fullpath'])].reset_index()
fdat['newfn'] = new_file_names(fdat, cid = 'HYD2401')
fdat['newfp'] = fdat['newfn'].apply(lambda x: os.path.join(mainpath, x))
metdat = metdat.merge(fdat[['fullpath','newfn','newfp']], how = 'left', on = 'fullpath')

#%% This is the renaming
for i in range(len(fdat)):
    os.rename(fdat['fullpath'][i], fdat['newfp'][i])
    
#%% this is the big loop
for i in range(len(metdat)):
    if metdat['newfp'][i] == np.nan:
        continue
    with xr.open_dataset(metdat['newfp'][i]) as xdata:
        print('\nProcessing:', os.path.basename(metdat['newfp'][i]))
        
        yd = xdata.pressure.to_dataframe()
        yd['pchange'] = yd['pressure'].diff()
        yd['pchange'] = yd['pchange'].replace(0, np.nan)
        yd['tchange'] = yd.index.diff().total_seconds()
        yd['descent_rate'] = yd['pchange']/yd['tchange']
        
        bot_index = trigger_index_lessthan(yd, 'descent_rate', 
                               trigger_float = -0.4, buffer = 100, window_size = 10)
        max_index = yd['pressure'].idxmax()
        print(f'Max: {max_index} at {yd.loc[max_index,"pressure"]}')
        bot_index = max(bot_index, max_index)
        
        yd['downcast'] = False
        yd.loc[yd.index[0]:bot_index,'downcast'] = True
        bot_pres = yd.loc[bot_index,'pressure']
        print('Bottom Pressure =', bot_pres)
        
        xdata.attrs['Bottom Pressure'] = bot_pres
        xdata['descent_rate'] = ('time', yd['descent_rate'])
        xdata['downcast'] = ('time', yd['downcast'])
        
        row = metdat[metdat['newfp'] == metdat['newfp'][i]] #change this to fullpath

        xdata.attrs['EXPEDITION'] = 'HYD2401'
        xdata.attrs['CTD_NUMBER'] = int(row['PROFILE_COUNT'].item())
        xdata.attrs['OBSERVATION_TYPE'] = row['OBSERVATION_TYPE'].item()
        xdata.attrs['OBSERVATION_TIMESTAMP'] = str(row['DATE_TIME (UTC)'].item()) 
        xdata.attrs['NOTES'] = row['NOTES'].item()
        xdata.attrs['OPERATOR'] = row['RECORDED_BY'].item()
        xdata.attrs['LATITUDE'] = row['y'].item()
        xdata.attrs['LONGITUDE'] = row['x'].item()
        xdata.attrs['REGION'] = row['Region'].item()
        xdata.attrs['VESSEL'] = row['Vessel'].item()
        xdata.attrs['METHOD'] = row['Method'].item()
        xdata.attrs['original raw filename'] = os.path.basename(row['fullpath'].item())
        xdata.attrs['new expedition filename'] = os.path.basename(row['newfp'].item())
        
        xdata.to_netcdf(os.path.join(mainpath,r"with_meta",os.path.basename(metdat['newfp'][i])))
        
#%%
fno = fullpath[56]#os.path.join(r"C:\temp\Ollie-RBR\testing.nc")
#%%
xdata = xr.open_dataset(r"C:\temp\Ollie-RBR\with_meta\HYD2401_CTD003_RBRconcerto3_SN214481_20240911_2049.nc")
#%%
print(list(xdata.keys()))
print(xdata.attrs)
#%%
xdata.close()
#%%
xdata.to_netcdf(r"C:\temp\Ollie-RBR\testing_new.nc")


#%%
yd = xdata.pressure.to_dataframe()
yd['pchange'] = yd['pressure'].diff()
yd['pchange'] = yd['pchange'].replace(0, np.nan)
yd['tchange'] = yd.index.diff().total_seconds()
yd['descent_rate'] = yd['pchange']/yd['tchange']

#%%
bot_index = trigger_index_lessthan(yd, 'descent_rate', 
                       trigger_float = -0.4, buffer = 100, window_size = 10)
max_index = yd['pressure'].idxmax()
print(f'Max: {max_index} at {yd.loc[max_index,"pressure"]}')
bot_index = max(bot_index, max_index)

#%%
yd['downcast'] = False
yd.loc[yd.index[0]:bot_index,'downcast'] = True
bot_pres = yd.loc[bot_index,'pressure']
print('Bottom Pressure =', bot_pres)
#%%
xdata.attrs['Bottom Pressure'] = bot_pres
xdata['descent_rate'] = ('time', yd['descent_rate'])
xdata['downcast'] = ('time', yd['downcast'])

#%% NEED TO CHANGE THIS TO A DICTIONARY
row = metdat[metdat['newfp'] == fno] #change this to fullpath
xdata.attrs['EXPEDITION'] = 'HYD2401'
xdata.attrs['CTD_NUMBER'] = int(row['PROFILE_COUNT'].item())
xdata.attrs['OBSERVATION_TYPE'] = row['OBSERVATION_TYPE'].item()
xdata.attrs['OBSERVATION_TIMESTAMP'] = str(row['DATE_TIME (UTC)'].item()) 
xdata.attrs['NOTES'] = row['NOTES'].item()
xdata.attrs['OPERATOR'] = row['RECORDED_BY'].item()
xdata.attrs['LATITUDE'] = row['y'].item()
xdata.attrs['LONGITUDE'] = row['x'].item()
xdata.attrs['REGION'] = row['Region'].item()
xdata.attrs['VESSEL'] = row['Vessel'].item()
xdata.attrs['METHOD'] = row['Method'].item()
xdata.attrs['original raw filename'] = os.path.basename(row['fullpath'].item())
xdata.attrs['new expedition filename'] = os.path.basename(row['newfp'].item())
#%%
xdatafin = xdata.copy()
#%%
xdatafin.to_netcdf('C:\\temp\\Ollie-RBR\\HYD2401_CTD059_RBRconcerto3_SN214481_20240921_0327.nc')

#%%
# Check which range each datetime falls into
results = []
ftak = []
for dt in metdat['Profile_start']: #metdat['DATE_TIME (UTC)']:#
    for fullpath, start, end in zip(fdat['fullpath'],fdat['start'], fdat['end']):
        if (start - timedelta(seconds=0)) <= (dt) <= end:
            results.append((dt, (fullpath, start, end)))
            ftak.append(fullpath)
            break
    else:
        results.append((dt, None))  # No matching range
        ftak.append(None)
# Output the results
for dt, match in results:
    if match:
        print(f"{dt} falls in range {match}")
    else:
        print(f"{dt} does not fall in any range")
#%%
metdat['fullpath'] = ftak
fixing = metdat[['PROFILE_COUNT','DATE_TIME (UTC)','fullpath']]
#%% Notepad
metdat.loc[metdat.query('63 <= PROFILE_COUNT <= 89').index,
           'DATE_TIME (UTC)'] = metdat.loc[metdat.query('63 <= PROFILE_COUNT <= 89').index,
                                           'DATE_TIME (UTC)'].apply(lambda x:
                                                                    x + timedelta(seconds=60))
                                                                    
descent_trigger = int(-0.4) #in m/s
surface_buffer = 300
bot_index = yd[surface_buffer:].index[yd[surface_buffer:]['descent_rate'] < - descent_trigger][0]
# maybe -0.4 as a cuttoff                                                                    

window_size = 10
yd['is_negative'] = yd['descent_rate'] < descent_trigger
yd['negative_streak'] = yd['is_negative'].rolling(window=window_size).sum()
bot_index = yd[surface_buffer:].index[yd[surface_buffer:]['negative_streak'] == window_size][0]






