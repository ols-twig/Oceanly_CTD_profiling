# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 09:22:05 2024

@author: OllieTwigge
"""
# renames each nc file taken from ruskin ctd to the standard required for the NAS
# at this point it needs to have a metatable also provided to assign the CTD number 
# TODO make it detect the last CTD number from the NAS and then just increment
# TODO make into a complete script that keeps the OriginalRaw then makes new folder with the new names 


import xarray as xr
import pandas as pd
import numpy as np
import os 
from pathlib import Path

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
            
            if any(i in fdat.columns for i in ['PROFILE_COUNT', 'CTD_Number']):
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

#%% inputting the meta file to assign CTD number
metdat = pd.read_excel('RBR-Metadata.xlsx', header= 1, parse_dates=True)
metdat['RAW FILE'] = metdat['RAW FILE'].ffill()
metdat['DATE_TIME (UTC)'] = pd.to_datetime(metdat['DATE_TIME (UTC)'])

#%%
mainpath = r"C:\temp\Ollie-RBR"
testpath = mainpath+r"\*.nc"

fullpath= list(glob(testpath, recursive=False))
fdat = extract_RBR_times_ID(fullpath)
fdat = fdat.merge(metdat[['fullpath','PROFILE_COUNT']], how = 'left', on = 'fullpath')
#fdat = fdat[fdat['fullpath'].isin(metdat['fullpath'])].reset_index()
fdat['newfn'] = new_file_names(fdat, cid = 'HYD2401')
fdat['newfp'] = fdat['newfn'].apply(lambda x: os.path.join(mainpath, x))

#%% This is the renaming
for i in range(len(fdat)):
    os.rename(fdat['fullpath'][i], fdat['newfp'][i])
    