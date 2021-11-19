# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19, 2021
Hydrograph Prep (v1
@author: tclarkin (USBR 2021)

...

"""
### Import Libraries ###
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from functions import nwis_import,analyze_voldur,init_hydrograph_plot

### Begin User Input ###
# Set Working Directory
os.chdir("C://Users//tclarkin//Documents//Projects//Anderson_Ranch_Dam//hydrographs")
site = "ard"
site_source = "13186000"

# Hydrograph Settings
duration = 120  # int, days
count = 20       # int, number of hydrographs to extract
sel_wy = None # list or None

### Begin Script ###
# Check for directory
if not os.path.isdir("hydro_prep"):
    os.mkdir("hydro_prep")

# Import daily data
site_daily = nwis_import(site_source,"dv")

# Remove data prior to 1989 (instantaneous data)
site_daily = site_daily.loc[site_daily.index.year>=1989,:]

# ID annual maxima for duration
site_evs = analyze_voldur(site_daily,duration)

# Select top hydrographs corresponding to user specified count
site_evs_sort = site_evs["avg_flow"].sort_values(ascending=False)
top_evs = site_evs.loc[site_evs_sort[0:count].index]

# If sel_wy specified, add to list
if sel_wy!=None:
    # If sel_wy is longer than count, return message
    if len(sel_wy)<count:
        print(f'{len(sel_wy)} events specified, while only {count} events are expected')
    sel_evs = pd.DataFrame(columns=site_evs.columns)
    for wy in sel_wy:
        if wy in top_evs.index:
            print(f"{wy} already in top events")
            continue
        else:
            print(f"{wy} not in top events; adding")
            sel_evs.loc[wy,:] = site_evs.loc[wy,:]
    if len(sel_evs)>0:
        top_evs = top_evs.drop(index=top_evs.tail(len(sel_evs)).index)
        for wy in sel_evs.index:
            top_evs.loc[wy,:] = sel_evs.loc[wy,:]

# Retrieve instantaneous data
hydro = pd.DataFrame(index=range(0,duration*24))
init_hydrograph_plot(duration*24)
for wy in top_evs.index:
    print(wy)
    start = top_evs.loc[wy,"start"].strftime("%Y-%m-%d")
    end = top_evs.loc[wy,'end'].strftime('%Y-%m-%d')
    site_inst = nwis_import(site_source,"iv",start=start,end=end)

    # Convert to hourly
    site_hour = site_inst["flow"].resample("H").mean()
    site_hour = site_hour.interpolate()
    site_hour = site_hour.dropna()
    site_hour = site_hour.reset_index(drop=True)

    if len(site_hour) < len(hydro):
        print("Adding timesteps...")
        n = 1
        while len(site_hour) < len(hydro):
            site_hour = site_hour.append(site_hour.tail(1))
            n += 1
        print(f"{n} timesteps added")
    if len(site_hour) > len(hydro):
        print(f"Clipping {len(site_hour)-len(hydro)} timesteps")
        site_hour = site_hour[0:duration*24]

    hydro[wy] = site_hour.values.astype("int")
    plt.plot(hydro[wy],label=wy)

# Save figure and files
plt.legend()
plt.savefig(f"hydro_prep/{site}_input_hydrographs.jpg", bbox_inches='tight', dpi=300)
hydro.to_csv(f"hydro_prep/{site}_input_hydrographs.csv")
top_evs.to_csv(f"hydro_prep/{site}_input_hydrographs_summary.csv")
